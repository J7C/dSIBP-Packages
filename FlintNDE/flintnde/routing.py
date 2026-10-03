"""统一普通点/奇点调度与命名直线路径规划。

端点必须显式命名，并可用字面量 ``"inf"`` 表示无穷远。无穷远在局部变量
``sinv=1/s`` 下处理，系数矩阵为 ``-sinv**(-2) A(1/sinv)``。调度器在普通点调用
Cauchy--DFT；奇点依次尝试 exact power-log、Fuchsian reduction 和认证的指数广义幂级数。
路径计划报告两端及路径上的全部奇点；内部奇点由相邻普通匹配点和局部基跨点元数据表示。
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Any

from flint import acb, acb_mat, arb

from .core import arb_ball_from_json, arb_ball_to_json, require_exact_keys
from .exact_gaussian import gaussian_rational
from .local_solutions import LocalReductionError, LocalSolutionBasis, build_local_solution_basis
from .singularities import (
    RationalMatrixSystem,
    SingularityInventory,
    SingularityRecord,
    analyze_singularities,
)
from .transport import build_straight_path

def _normalize_message_language(value: str) -> str:
    """规范化路径规划提示语言；只接受公开合同中的 EN/CN。"""

    if not isinstance(value, str):
        raise TypeError("message_language must be a string")
    language = value
    if language not in {"EN", "CN"}:
        raise ValueError('message_language must be "EN" or "CN"')
    return language


def _normalize_singularity_mode(value: str) -> str:
    """校验公开奇点模式是否属于允许集合。"""

    if not isinstance(value, str):
        raise TypeError("singularity_mode must be a string")
    if value not in {"avoid", "singularity_jump"}:
        raise ValueError(
            'singularity_mode must be exactly "avoid" or "singularity_jump"'
        )
    return value

class AdaptivePathSingularityError(ValueError):
    """缺省避奇点模式发现内部奇点时返回结构化路段报告。"""

    def __init__(self, message: str, report: dict[str, Any]) -> None:
        super().__init__(message)
        self.report = report
        self.singular_path_pairs = tuple(
            (item["start"], item["target"]) for item in report["segments"]
        )



@dataclass(frozen=True)
class NamedPoint:
    """保存既有奇点路由使用的诊断名称及坐标。"""

    name: str
    value: Any

    def __post_init__(self) -> None:
        """拒绝空名称和含混的无穷远拼写。"""

        if not self.name.strip():
            raise ValueError("point name must not be empty")
        if isinstance(self.value, str) and self.value.lower() == "inf" and self.value != "inf":
            raise ValueError('infinity must be written exactly as "inf"')

    @property
    def is_infinity(self) -> bool:
        """返回端点是否表示无穷远。"""

        return isinstance(self.value, str) and self.value == "inf"

    def finite_acb(self) -> acb:
        """把有限端点转为 Acb，同时保留已有 Acb ball。"""

        if self.is_infinity:
            raise ValueError(f"{self.name} is infinity, not a finite point")
        if isinstance(self.value, acb):
            return acb(self.value)
        if isinstance(self.value, complex):
            return acb(str(self.value.real), str(self.value.imag))
        try:
            return gaussian_rational(self.value).to_acb()
        except (TypeError, ValueError):
            return acb(str(self.value))

    def exact_gaussian_value(self) -> str | None:
        """当端点不是浮点输入时，返回可移植的 Q(i) 精确值。"""

        if self.is_infinity or isinstance(self.value, (float, complex, acb, arb)):
            return None
        try:
            return str(gaussian_rational(self.value))
        except (TypeError, ValueError):
            return None

@dataclass(frozen=True)
class PointClassification:
    """保存一个命名点的分类和收敛半径规则。"""

    name: str
    value: acb | None
    value_label: str
    kind: str
    singularity_identifier: str | None
    pole_order: int
    convergence_radius: arb | None

    def to_json(self) -> dict[str, Any]:
        """序列化点分类和半径，不引入 binary64 精度损失。"""

        return {
            "name": self.name,
            "value": self.value_label,
            "classification": self.kind,
            "singularity_identifier": self.singularity_identifier,
            "pole_order": self.pole_order,
            "convergence_radius": (
                None if self.convergence_radius is None else self.convergence_radius.str(40)
            ),
            "radius_rule": (
                "nearest other finite singularity"
                if self.kind != "ordinary"
                else "nearest finite singularity"
            ),
        }


def _nearest_distance(point: acb, records: list[SingularityRecord]) -> arb | None:
    """返回最近有限记录的距离；局部圆无界时返回 None。"""

    distances = [abs(point - record.location) for record in records if record.location is not None]
    if not distances:
        return None
    return min(distances, key=lambda value: float(value.mid()))


def classify_point(
    endpoint: NamedPoint,
    inventory: SingularityInventory,
) -> PointClassification:
    """使用已验证清单分类命名有限点或无穷远。"""

    if endpoint.is_infinity:
        record = inventory.infinity
        return PointClassification(
            endpoint.name,
            None,
            "inf",
            record.kind,
            record.identifier if record.kind != "ordinary" else None,
            record.pole_order,
            None,
        )
    point = endpoint.finite_acb()
    matches = [
        record
        for record in inventory.finite
        if record.location is not None and abs(point - record.location).contains(0)
    ]
    if len(matches) > 1:
        raise ValueError(f"{endpoint.name}: point overlaps more than one isolated singularity")
    if not matches:
        radius = _nearest_distance(point, list(inventory.finite))
        return PointClassification(
            endpoint.name,
            point,
            point.str(40),
            "ordinary",
            None,
            0,
            radius,
        )
    record = matches[0]
    others = [candidate for candidate in inventory.finite if candidate is not record]
    radius = _nearest_distance(point, others)
    return PointClassification(
        endpoint.name,
        point,
        record.location_exact or point.str(40),
        record.kind,
        record.identifier,
        record.pole_order,
        radius,
    )


@dataclass(frozen=True)
class LocalExpansion:
    """统一局部方法调度器的返回结果。"""

    classification: PointClassification
    working_variable: str
    method: str
    convergence_radius: arb | None
    matrix_coefficients: tuple[acb_mat, ...] | None
    power_log_basis: LocalSolutionBasis | None
    manifest: dict[str, Any]


def _working_system_and_point(
    system: RationalMatrixSystem, endpoint: NamedPoint
) -> tuple[RationalMatrixSystem, NamedPoint, bool]:
    """把无穷远端点映射为反演变量的零点。"""

    if endpoint.is_infinity:
        return system.inverted("sinv"), NamedPoint(endpoint.name, 0), True
    return system, endpoint, False


def prepare_local_expansion(
    system: RationalMatrixSystem,
    endpoint: NamedPoint,
    *,
    order: int,
    radius_fraction: float = 0.60,
    radius: acb | None = None,
    sample_count: int | None = None,
    output_layout: Any | None = None,
) -> LocalExpansion:
    """按 exact 结构选择普通点、Fuchsian 化或指数广义幂级数引擎。

    高阶 pole 会先做 exact shearing；仍保留二阶 pole 时还可接受单重不同主导根的
    start-only 形式渐近 sector。普通点若没有有限奇点限制局部圆，则必须显式提供
    ``radius``；形式渐近基的返回收敛半径为 ``None``。
    """

    if order <= 0:
        raise ValueError("local expansion order must be positive")
    if not 0 < radius_fraction < 1:
        raise ValueError("radius_fraction must lie between zero and one")
    working_system, working_endpoint, inverted = _working_system_and_point(system, endpoint)
    inventory = analyze_singularities(
        working_system,
        output_layout=output_layout,
        filename=("singularity_inventory_sinv.json" if inverted else "singularity_inventory.json"),
    )
    classification = classify_point(working_endpoint, inventory)
    if classification.kind == "ordinary":
        sampling_radius = radius
        if sampling_radius is None:
            if classification.convergence_radius is None:
                raise ValueError(
                    f"{endpoint.name}: no finite singularity bounds the Cauchy disk; provide radius"
                )
            sampling_radius = acb(
                classification.convergence_radius * arb(str(radius_fraction))
            )
        coefficients = working_system.to_analytic_system(inventory).taylor_matrix_coefficients(
            working_endpoint.finite_acb(),
            order,
            radius=sampling_radius,
            sample_count=sample_count,
        )
        manifest = {
            "schema": "flintnde_local_dispatch_v1",
            "point": classification.to_json(),
            "working_variable": working_system.variable_name,
            "infinity_transformation": inverted,
            "method": "ordinary_cauchy_dft",
            "order": order,
            "sample_count": sample_count or max(32, 2 * order),
            "sampling_radius": sampling_radius.str(40),
        }
        if output_layout is not None:
            output_layout.write_json(
                "transport", f"{endpoint.name}_ordinary_expansion.json", manifest
            )
        return LocalExpansion(
            classification,
            working_system.variable_name,
            "ordinary_cauchy_dft",
            classification.convergence_radius,
            tuple(coefficients),
            None,
            manifest,
        )

    exact_point = working_endpoint.exact_gaussian_value() or classification.value_label
    try:
        gaussian_rational(exact_point)
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"{endpoint.name}: automatic exact Frobenius extraction currently requires a Q(i) singularity"
        ) from error
    basis = build_local_solution_basis(working_system, exact_point, order)
    manifest = {
        "schema": "flintnde_local_dispatch_v1",
        "point": classification.to_json(),
        "working_variable": working_system.variable_name,
        "infinity_transformation": inverted,
        "method": basis.method,
        "order": order,
        "frobenius": basis.manifest.get("frobenius"),
        "local_basis": basis.manifest,
    }
    if output_layout is not None:
        output_layout.write_json("frobenius", f"{endpoint.name}_frobenius.json", manifest)
    return LocalExpansion(
        classification,
        working_system.variable_name,
        basis.method,
        (
            classification.convergence_radius
            if basis.continuation_ready
            else None
        ),
        None,
        basis,
        manifest,
    )


@dataclass(frozen=True)
class PathPlan:
    """保存命名路径点、所选局部方法和跨奇点要求。"""

    name: str
    original_variable: str
    working_variable: str
    infinity_transformation: bool
    max_step_over_radius: float
    start: PointClassification
    target: PointClassification
    internal_singularities: tuple[PointClassification, ...]
    points: tuple[dict[str, Any], ...]
    segments: tuple[dict[str, Any], ...]
    planning_action: str
    singularity_mode: str
    message_language: str
    continuation_ready: bool
    messages: tuple[str, ...]

    def to_json(self) -> dict[str, Any]:
        """序列化完整路径决策和生成点。"""

        return {
            "schema": "flintnde_adaptive_path_plan_v1",
            "path_name": self.name,
            "original_variable": self.original_variable,
            "working_variable": self.working_variable,
            "infinity_transformation": self.infinity_transformation,
            "max_step_over_convergence_radius": self.max_step_over_radius,
            "start": self.start.to_json(),
            "target": self.target.to_json(),
            "internal_singularities": [item.to_json() for item in self.internal_singularities],
            "points": list(self.points),
            "segments": list(self.segments),
            "planning_action": self.planning_action,
            "singularity_mode": self.singularity_mode,
            "message_language": self.message_language,
            "continuation_ready": self.continuation_ready,
            "messages": list(self.messages),
        }


class AdaptivePath(list[acb]):
    """可直接输运且携带奇点跨越、半径和逐步比值的路径变量。

    列表元素只包含可实际求值的普通匹配点；内部奇点作为相邻两点之间的 transition
    元数据保存，避免把未必有限的 ``Y(singularity)`` 冒充普通 snapshot。
    """

    def __init__(
        self,
        points: list[acb],
        *,
        working_system: RationalMatrixSystem,
        transitions: list[dict[str, Any]],
        step_reports: list[dict[str, Any]],
        internal_singularities: tuple[PointClassification, ...],
        start_classification: PointClassification,
        target_classification: PointClassification,
        working_variable: str,
        infinity_transformation: bool,
        save_requests: tuple[dict[str, Any], ...] = (),
        formal_asymptotic_match_estimate: dict[str, Any] | None = None,
        planning_action: str = "raw_points_automatic_plan",
        singularity_mode: str = "avoid",
        message_language: str = "EN",
    ) -> None:
        """保存直接执行点和与其严格对齐的 transition/step 诊断。"""

        super().__init__(acb(point) for point in points)
        if len(transitions) != max(0, len(points) - 1):
            raise ValueError("adaptive path transition count does not match point count")
        self._working_system = working_system
        self._transitions = tuple(transitions)
        self.step_reports = tuple(step_reports)
        self.internal_singularities = internal_singularities
        self.start_classification = start_classification
        self.target_classification = target_classification
        self.working_variable = working_variable
        self.infinity_transformation = infinity_transformation
        self.save_requests = save_requests
        self.formal_asymptotic_match_estimate = formal_asymptotic_match_estimate
        self.planning_action = planning_action
        self.singularity_mode = singularity_mode
        self.message_language = message_language
        finite_ratios = [
            report["step_over_convergence_radius"]
            for report in step_reports
            if report["step_over_convergence_radius"] is not None
        ]
        self.max_step_over_convergence_radius = max(finite_ratios, default=None)


    def to_json(self, digits: int = 80) -> dict[str, Any]:
        """序列化可执行通用路径，供另一进程恢复后直接输运。"""

        return adaptive_path_to_json(self, digits=digits)


def _adaptive_acb_record(value: acb, digits: int) -> dict[str, Any]:
    """保存通用路径节点的展示中点及可严格恢复的 Arb 球。"""

    return {
        "real": value.real.mid().str(digits, radius=False, more=True),
        "imag": value.imag.mid().str(digits, radius=False, more=True),
        "real_ball": arb_ball_to_json(value.real, digits),
        "imag_ball": arb_ball_to_json(value.imag, digits),
    }


def _adaptive_acb_from_record(record: Any, field_name: str) -> acb:
    """从当前计划 JSON 的严格 Arb 球字段恢复 Acb。"""

    require_exact_keys(
        record,
        {"real", "imag", "real_ball", "imag_ball"},
        field_name,
    )
    return acb(
        arb_ball_from_json(record["real_ball"], f"{field_name}.real_ball"),
        arb_ball_from_json(record["imag_ball"], f"{field_name}.imag_ball"),
    )

def _adaptive_classification_record(
    classification: PointClassification, digits: int
) -> dict[str, Any]:
    """保存执行所需的完整点分类，不依赖展示字符串反解析。"""

    return {
        "name": classification.name,
        "value": (
            None
            if classification.value is None
            else _adaptive_acb_record(classification.value, digits)
        ),
        "value_label": classification.value_label,
        "kind": classification.kind,
        "singularity_identifier": classification.singularity_identifier,
        "pole_order": classification.pole_order,
        "convergence_radius": (
            None
            if classification.convergence_radius is None
            else classification.convergence_radius.mid().str(
                digits, radius=False, more=True
            )
        ),
    }


def _adaptive_classification_from_json(
    record: dict[str, Any],
    inventory: SingularityInventory,
    field_name: str,
) -> PointClassification:
    """按当前系统重新分类序列化坐标，并核对计划未被用于另一系统。"""

    require_exact_keys(
        record,
        {
            "name",
            "value",
            "value_label",
            "kind",
            "singularity_identifier",
            "pole_order",
            "convergence_radius",
        },
        field_name,
    )
    if record["value"] is None:
        raise ValueError(f"{field_name} lacks a finite working coordinate")
    if not isinstance(record["name"], str) or not isinstance(record["kind"], str):
        raise ValueError(f"{field_name} name and kind must be strings")
    identifier = record["singularity_identifier"]
    if identifier is not None and not isinstance(identifier, str):
        raise ValueError(f"{field_name} singularity_identifier must be a string or null")
    pole_order = record["pole_order"]
    if isinstance(pole_order, bool) or not isinstance(pole_order, int):
        raise ValueError(f"{field_name} pole_order must be an integer")
    point = _adaptive_acb_from_record(record["value"], f"{field_name}.value")
    actual = classify_point(NamedPoint(record["name"], point), inventory)
    expected = (record["kind"], identifier, pole_order)
    observed = (
        actual.kind,
        actual.singularity_identifier,
        actual.pole_order,
    )
    if observed != expected:
        raise ValueError(
            f"{field_name} classification does not match the current rational system"
        )
    return actual

def adaptive_path_to_json(path: AdaptivePath, *, digits: int = 80) -> dict[str, Any]:
    """序列化 AdaptivePath 的节点与局部基 bridge。

    该记录面向两阶段接口：规划进程产生它，执行进程只恢复节点、核对奇点清单并
    输运，不重新选择路径。
    """

    if digits < 20:
        raise ValueError("adaptive-path serialization needs at least 20 decimal digits")
    transitions: list[dict[str, Any]] = []
    for transition in path._transitions:
        method = transition.get("method")
        if method == "ordinary_taylor":
            transitions.append({"method": method})
        elif method == "regular_singular_bridge":
            transitions.append(
                {
                    "method": method,
                    "singularity": _adaptive_classification_record(
                        transition["singularity"], digits
                    ),
                }
            )
        else:
            raise ValueError(f"unsupported adaptive transition method: {method}")
    return {
        "schema": "flintnde_adaptive_path_serialized_v1",
        "planning_precision_digits": digits,
        "points": [_adaptive_acb_record(point, digits) for point in path],
        "transitions": transitions,
        "step_reports": list(path.step_reports),
        "internal_singularities": [
            _adaptive_classification_record(item, digits)
            for item in path.internal_singularities
        ],
        "start_classification": _adaptive_classification_record(
            path.start_classification, digits
        ),
        "target_classification": _adaptive_classification_record(
            path.target_classification, digits
        ),
        "working_variable": path.working_variable,
        "infinity_transformation": path.infinity_transformation,
        "save_requests": list(path.save_requests),
        "formal_asymptotic_match_estimate": path.formal_asymptotic_match_estimate,
        "planning_action": path.planning_action,
        "singularity_mode": path.singularity_mode,
        "message_language": path.message_language,
    }


def adaptive_path_from_json(
    system: RationalMatrixSystem,
    record: dict[str, Any],
) -> AdaptivePath:
    """恢复并校验当前通用有理矩阵计划；本函数不调用路径规划器。"""

    if not isinstance(system, RationalMatrixSystem):
        raise TypeError("adaptive path restoration requires a RationalMatrixSystem")
    require_exact_keys(
        record,
        {
            "schema",
            "planning_precision_digits",
            "points",
            "transitions",
            "step_reports",
            "internal_singularities",
            "start_classification",
            "target_classification",
            "working_variable",
            "infinity_transformation",
            "save_requests",
            "formal_asymptotic_match_estimate",
            "planning_action",
            "singularity_mode",
            "message_language",
        },
        "adaptive path",
    )
    if record["schema"] != "flintnde_adaptive_path_serialized_v1":
        raise ValueError("unsupported adaptive-path schema")
    planning_digits = record["planning_precision_digits"]
    if isinstance(planning_digits, bool) or not isinstance(planning_digits, int):
        raise ValueError("adaptive-path planning_precision_digits must be an integer")
    if planning_digits < 20:
        raise ValueError("adaptive-path planning_precision_digits must be at least 20")
    inverted = record["infinity_transformation"]
    if not isinstance(inverted, bool):
        raise ValueError("adaptive-path infinity_transformation must be Boolean")
    working = system.inverted("sinv") if inverted else system
    inventory = analyze_singularities(working)
    point_records = record["points"]
    if not isinstance(point_records, list):
        raise ValueError("adaptive-path points must be a list")
    points = [
        _adaptive_acb_from_record(item, f"points[{index}]")
        for index, item in enumerate(point_records)
    ]
    if len(points) < 2:
        raise ValueError("a serialized adaptive path needs at least two points")
    transition_records = record["transitions"]
    if not isinstance(transition_records, list):
        raise ValueError("adaptive-path transitions must be a list")
    if len(transition_records) != len(points) - 1:
        raise ValueError("serialized adaptive transition count does not match its points")
    transitions: list[dict[str, Any]] = []
    for index, item in enumerate(transition_records):
        if not isinstance(item, dict) or "method" not in item:
            raise ValueError(f"transitions[{index}] must contain method")
        method = item["method"]
        if method == "ordinary_taylor":
            require_exact_keys(item, {"method"}, f"transitions[{index}]")
            transitions.append({"method": method})
        elif method == "regular_singular_bridge":
            require_exact_keys(
                item,
                {"method", "singularity"},
                f"transitions[{index}]",
            )
            singularity = _adaptive_classification_from_json(
                item["singularity"], inventory, f"transitions[{index}].singularity"
            )
            transitions.append({"method": method, "singularity": singularity})
        else:
            raise ValueError(f"unsupported adaptive transition method: {method}")
    start_classification = _adaptive_classification_from_json(
        record["start_classification"], inventory, "start_classification"
    )
    target_classification = _adaptive_classification_from_json(
        record["target_classification"], inventory, "target_classification"
    )
    internal_records = record["internal_singularities"]
    if not isinstance(internal_records, list):
        raise ValueError("adaptive-path internal_singularities must be a list")
    internal = tuple(
        _adaptive_classification_from_json(
            item, inventory, f"internal_singularities[{index}]"
        )
        for index, item in enumerate(internal_records)
    )
    step_reports = record["step_reports"]
    if not isinstance(step_reports, list):
        raise ValueError("serialized adaptive step reports must be a list")
    save_requests = record["save_requests"]
    if not isinstance(save_requests, list):
        raise ValueError("serialized adaptive save requests must be a list")
    working_variable = record["working_variable"]
    planning_action = record["planning_action"]
    if not isinstance(working_variable, str) or not isinstance(planning_action, str):
        raise ValueError("adaptive-path working_variable and planning_action must be strings")
    formal_estimate = record["formal_asymptotic_match_estimate"]
    if formal_estimate is not None and not isinstance(formal_estimate, dict):
        raise ValueError("adaptive-path formal_asymptotic_match_estimate must be an object or null")
    language = _normalize_message_language(record["message_language"])
    try:
        singularity_mode = _normalize_singularity_mode(record["singularity_mode"])
    except (TypeError, ValueError) as error:
        raise ValueError("serialized adaptive path has an invalid singularity mode") from error
    return AdaptivePath(
        points,
        working_system=working,
        transitions=transitions,
        step_reports=step_reports,
        internal_singularities=internal,
        start_classification=start_classification,
        target_classification=target_classification,
        working_variable=working_variable,
        infinity_transformation=inverted,
        save_requests=tuple(save_requests),
        formal_asymptotic_match_estimate=formal_estimate,
        planning_action=planning_action,
        singularity_mode=singularity_mode,
        message_language=language,
    )

def _map_path_endpoints_at_infinity(
    system: RationalMatrixSystem, start: NamedPoint, target: NamedPoint
) -> tuple[RationalMatrixSystem, NamedPoint, NamedPoint, bool]:
    """把含一个无穷远端点的路径映射到倒数变量。"""

    if start.is_infinity and target.is_infinity:
        raise ValueError("a path cannot have infinity as both endpoints")
    if not start.is_infinity and not target.is_infinity:
        return system, start, target, False
    working = system.inverted("sinv")

    def mapped(endpoint: NamedPoint) -> NamedPoint:
        if endpoint.is_infinity:
            return NamedPoint(endpoint.name, 0)
        value = endpoint.finite_acb()
        if abs(value).contains(0):
            raise ValueError("a finite zero endpoint maps to infinity in sinv; split the route explicitly")
        return NamedPoint(endpoint.name, 1 / value)

    return working, mapped(start), mapped(target), True


def _on_segment_parameter(
    point: acb, start: acb, target: acb, tolerance: float
) -> arb | None:
    """用 Arb 球几何返回与路径内部相交的奇点参数。"""

    direction = target - start
    denominator = (
        direction.real * direction.real + direction.imag * direction.imag
    )
    if denominator.contains(0):
        return None
    offset = point - start
    parameter = (
        offset.real * direction.real + offset.imag * direction.imag
    ) / denominator
    parameter_midpoint = parameter.mid()
    tolerance_ball = arb(str(tolerance))
    if not tolerance_ball < parameter_midpoint < arb(1) - tolerance_ball:
        return None
    projected = start + direction * acb(parameter)
    point_radius = max(point.real.rad(), point.imag.rad())
    separation = abs(point - projected)
    if separation.contains(0) or separation < tolerance_ball + point_radius:
        return parameter_midpoint
    return None

def _point_record(name: str, point: acb, classification: PointClassification) -> dict[str, Any]:
    """构造一个可移植的生成路径点记录。"""

    return {
        "name": name,
        "value": point.str(40),
        "classification": classification.kind,
        "singularity_identifier": classification.singularity_identifier,
        "convergence_radius": (
            None
            if classification.convergence_radius is None
            else classification.convergence_radius.str(40)
        ),
    }


def _map_detour_point(endpoint: NamedPoint, inverted: bool) -> NamedPoint:
    """把用户绕行点映射到当前工作变量；无穷远绕行点和反演后的零点均不允许。"""

    if endpoint.is_infinity:
        raise ValueError("detour points must be finite")
    if not inverted:
        return endpoint
    value = endpoint.finite_acb()
    if abs(value).contains(0):
        raise ValueError("a zero detour point maps to infinity in sinv")
    return NamedPoint(endpoint.name, 1 / value)


def _internal_singularities_on_leg(
    inventory: SingularityInventory,
    start: acb,
    target: acb,
    path_tolerance: float,
) -> list[PointClassification]:
    """按路径顺序返回一条有限线段内部的全部奇点。"""

    found: list[tuple[arb, PointClassification]] = []
    for record in inventory.finite:
        if record.location is None:
            continue
        parameter = _on_segment_parameter(record.location, start, target, path_tolerance)
        if parameter is None:
            continue
        endpoint = NamedPoint(record.identifier, record.location)
        found.append((parameter, classify_point(endpoint, inventory)))
    found.sort(key=lambda item: item[0])
    return [classification for _, classification in found]


def _ordinary_match_point(
    point: acb,
    classification: PointClassification,
    direction: acb,
    interval_length: arb,
    max_step_over_radius: float,
    preferred_distance: arb | None = None,
) -> acb:
    """奇点端点向区间内部移动到局部收敛圆中的普通匹配点。"""

    if classification.kind == "ordinary":
        return acb(point)
    distance = interval_length / 3
    if preferred_distance is not None:
        distance = min(
            distance,
            preferred_distance,
            key=lambda value: float(value.mid()),
        )
    if classification.convergence_radius is not None:
        distance = min(
            distance,
            classification.convergence_radius * arb(str(max_step_over_radius)),
            key=lambda value: float(value.mid()),
        )
    return point + direction * distance


def _route_point(value: Any, role: str) -> tuple[NamedPoint, bool]:
    """把裸坐标或 ``(coordinate,"save")`` 转为内部端点；名称不属于 save 合同。"""

    tagged = isinstance(value, (list, tuple)) and len(value) == 2 and value[1] == "save"
    coordinate = value[0] if tagged else value
    if isinstance(coordinate, NamedPoint):
        return coordinate, tagged
    return NamedPoint(role, coordinate), tagged


def _preflight_local_basis(
    system: RationalMatrixSystem,
    inventory: SingularityInventory,
    classification: PointClassification,
    order: int,
) -> tuple[SingularityRecord, LocalSolutionBasis]:
    """对一个非普通点执行正式建路所需的同一局部基能力检查。

    输入必须是 inventory 中的奇点分类。当前 exact 局部 dispatcher 只接受 Q(i)
    中心；位置或局部谱超出能力时统一抛出 ``LocalReductionError``，调用方决定
    是把失败记录进计划，还是立即停止正式建路。
    """

    record = inventory.find_finite(classification.singularity_identifier)
    if record is None or record.location_exact is None:
        raise LocalReductionError(
            f"{classification.name}: {classification.kind} local basis requires "
            "an exact Q(i) singular location"
        )
    try:
        basis = build_local_solution_basis(system, record.location_exact, order)
    except (ValueError, NotImplementedError, LocalReductionError) as error:
        raise LocalReductionError(
            f"{classification.name}: unsupported {classification.kind} local basis: {error}"
        ) from error
    return record, basis


def build_adaptive_path(
    system: RationalMatrixSystem,
    start: Any,
    target: Any,
    *,
    detour_points: tuple[Any, ...] = (),
    path_name: str | None = None,
    max_step_over_radius: float = 0.45,
    path_tolerance: float = 1.0e-24,
    formal_asymptotic_order: int | None = None,
    formal_minimum_order_factor: float = 3.0,
    singularity_mode: str = "avoid",
    message_language: str = "EN",
    output_layout: Any | None = None,
) -> AdaptivePath:
    """生成可直接作为 ``transport_path`` 参数的自描述 Acb 路径。

    调用本函数就表示从原始端点自动规划；已有路径应直接交给输运函数。缺省
    ``singularity_mode="avoid"`` 拒绝穿过内部奇点，用户可给出明确绕行点；只有显式
    选择 ``"singularity_jump"`` 才建立局部基 bridge，且其多值分支必须由用户确认。
    步长限制始终使用程序内部发现的最近奇点和收敛半径。
    """

    language = _normalize_message_language(message_language)
    resolved_singularity_mode = _normalize_singularity_mode(singularity_mode)
    start, start_save = _route_point(start, "start")
    target, target_save = _route_point(target, "target")
    parsed_detours = tuple(
        _route_point(point, f"detour_{index + 1:03d}")
        for index, point in enumerate(detour_points)
    )
    detour_points = tuple(point for point, _save in parsed_detours)
    save_flags = (start_save, *(save for _point, save in parsed_detours), target_save)
    if not 0 < max_step_over_radius < 1:
        raise ValueError("max_step_over_radius must lie between zero and one")
    if path_tolerance <= 0:
        raise ValueError("path_tolerance must be positive")
    if formal_asymptotic_order is not None and formal_asymptotic_order <= 0:
        raise ValueError("formal_asymptotic_order must be positive")
    if formal_minimum_order_factor <= 1:
        raise ValueError("formal_minimum_order_factor must exceed one")
    resolved_name = path_name or f"{start.name}_to_{target.name}"
    working, working_start, working_target, inverted = _map_path_endpoints_at_infinity(
        system, start, target
    )
    mapped_detours = tuple(
        _map_detour_point(endpoint, inverted) for endpoint in detour_points
    )
    original_inventory = analyze_singularities(system, output_layout=output_layout)
    working_inventory = (
        analyze_singularities(
            working,
            output_layout=output_layout,
            filename="singularity_inventory_sinv.json",
        )
        if inverted
        else original_inventory
    )
    route = (working_start, *mapped_detours, working_target)
    route_classifications = tuple(
        classify_point(endpoint, working_inventory) for endpoint in route
    )
    bad_middle_points = [
        item.name
        for item, save_requested in zip(
            route_classifications[1:-1], save_flags[1:-1]
        )
        if item.kind != "ordinary" and not save_requested
    ]
    if bad_middle_points:
        raise ValueError(
            "singular intermediate points must carry the save tag; unsupported points: "
            + ", ".join(bad_middle_points)
        )

    internal_by_leg: list[tuple[int, list[PointClassification]]] = []
    checkpoints: list[tuple[NamedPoint, PointClassification]] = []
    for leg_index, (left, right) in enumerate(zip(route[:-1], route[1:]), start=1):
        left_value = left.finite_acb()
        right_value = right.finite_acb()
        if abs(right_value - left_value).contains(0):
            raise ValueError(f"path leg {leg_index} has identical endpoints")
        internal = _internal_singularities_on_leg(
            working_inventory,
            left_value,
            right_value,
            path_tolerance,
        )
        if internal:
            internal_by_leg.append((leg_index, internal))
        if not checkpoints:
            checkpoints.append((left, route_classifications[leg_index - 1]))
        checkpoints.extend(
            (NamedPoint(item.name, item.value), item) for item in internal
        )
        checkpoints.append((right, route_classifications[leg_index]))

    internal_classifications = tuple(
        item for _, records in internal_by_leg for item in records
    )
    bridge_classifications = tuple(
        classification
        for _point, classification in checkpoints[1:-1]
        if classification.kind != "ordinary"
    )
    if internal_by_leg and resolved_singularity_mode == "avoid":
        singular_segments = [
            {
                "leg_index": leg_index,
                "start": route[leg_index - 1].finite_acb().str(40),
                "target": route[leg_index].finite_acb().str(40),
                "singularities": [
                    {
                        "identifier": item.singularity_identifier or item.name,
                        "coordinate": item.value_label,
                        "classification": item.kind,
                        "pole_order": item.pole_order,
                    }
                    for item in records
                ],
            }
            for leg_index, records in internal_by_leg
        ]
        segment_summary = "; ".join(
            f"leg {item['leg_index']} ({item['start']} -> {item['target']}): "
            + ", ".join(
                record["identifier"] for record in item["singularities"]
            )
            for item in singular_segments
        )
        message = (
            "缺省避开奇点模式拒绝生成可执行路径。经过奇点的路段为："
            f"{segment_summary}。请通过 detour_points 给出明确绕行点；也可显式设置 "
            "singularity_mode='singularity_jump' 使用奇点折跃。奇点折跃选择的多值分支等价于某一"
            "绕行路径，必须由用户确认。"
            if language == "CN"
            else
            "Default avoid-singularity mode refuses to build an executable path. "
            f"Singular legs: {segment_summary}. Supply explicit detour_points, or set "
            "singularity_mode='singularity_jump' explicitly. A singularity jump selects a multivalued branch "
            "equivalent to a detour path; the user must confirm it."
        )
        raise AdaptivePathSingularityError(
            message,
            {
                "schema": "flintnde_singular_path_report_v1",
                "singularity_mode": "avoid",
                "message_language": language,
                "segments": singular_segments,
            },
        )
    formal_start_distance: arb | None = None
    formal_start_estimate: dict[str, Any] | None = None
    # 所有非普通点先走同一个局部 dispatcher；形式渐近基只允许从起点向外初始化。
    for role, classification in (
        ("start", route_classifications[0]),
        *(("internal", item) for item in bridge_classifications),
        ("target", route_classifications[-1]),
    ):
        if classification.kind == "ordinary":
            continue
        _record, local_basis = _preflight_local_basis(
            working, working_inventory, classification, 4
        )
        if not local_basis.continuation_ready and role != "start":
            raise LocalReductionError(
                f"{classification.name}: {local_basis.method} is start-only; "
                "internal/target continuation requires Stokes connection data"
            )
        if local_basis.method == "formal_exponential_asymptotic" and role == "start":
            if formal_asymptotic_order is not None:
                formal_start_distance, formal_start_estimate = (
                    local_basis.suggest_match_distance(
                        formal_asymptotic_order,
                        formal_minimum_order_factor,
                    )
                )
    internal_identifiers = [
        item.singularity_identifier or item.name
        for _, records in internal_by_leg
        for item in records
    ]
    internal_method_label = (
        "Frobenius bridges"
        if all(item.kind == "regular_singular" for item in bridge_classifications)
        else "automatic high-pole local-basis attempts"
    )

    if language == "CN":
        messages = [
            "输入为原始端点，已自动规划单变量路径；已有计划应直接交给执行函数。",
            (
                "当前显式使用奇点折跃模式；所选多值分支等价于某一绕行路径，"
                "必须由用户确认。"
                if resolved_singularity_mode == "singularity_jump"
                else "当前使用避开奇点模式（缺省）。"
            ),
            f"起点 '{start.name}' 分类为 {route_classifications[0].kind}",
            f"终点 '{target.name}' 分类为 {route_classifications[-1].kind}",
            f"步长不超过收敛半径的 {max_step_over_radius} 倍",
        ]
    else:
        messages = [
            "Raw endpoints were supplied, so the one-variable path was planned automatically; "
            "pass an existing plan directly to the execution function.",
            (
                "Singularity mode: singularity jump (explicit). The selected multivalued branch is "
                "equivalent to a detour path and must be confirmed by the user."
                if resolved_singularity_mode == "singularity_jump"
                else "Singularity mode: avoid (default)."
            ),
            f"start '{start.name}' is {route_classifications[0].kind}",
            f"target '{target.name}' is {route_classifications[-1].kind}",
            f"step length is limited by {max_step_over_radius} times the convergence radius",
        ]
    if inverted:
        messages.append("path points use sinv=1/s")
    if route_classifications[0].kind == "non_fuchsian_input_basis":
        messages.append(
            "the high-pole start requires the five-order asymptotic convergence report"
        )
        if formal_start_estimate is not None:
            messages.append(
                "formal start match distance uses the nearest exponential-root gap and places "
                "the requested order at one third of the estimated least-term degree"
            )
    if mapped_detours:
        messages.append(
            ("用户指定的绕行点：" if language == "CN" else "user-selected detours: ")
            + ", ".join(point.name for point in detour_points)
        )
    elif internal_by_leg:
        messages.append(
            (
                "内部奇点 "
                + ", ".join(internal_identifiers)
                + " 使用显式奇点折跃所允许的局部基桥接（"
                + internal_method_label
                + "）"
                if language == "CN"
                else "internal singularities "
                + ", ".join(internal_identifiers)
                + f" use local bridges ({internal_method_label})"
            )
        )
    else:
        messages.append(
            "直线路径不含内部奇点"
            if language == "CN"
            else "the straight path contains no internal singularity"
        )

    analytic_system = working.to_analytic_system(working_inventory)
    generated: list[acb] = []
    generated_classifications: list[PointClassification] = []
    transitions: list[dict[str, Any]] = []
    endpoint_matches: list[dict[str, Any]] = []

    def append_unique(
        point: acb,
        classification: PointClassification,
        transition: dict[str, Any] | None = None,
    ) -> None:
        """追加普通执行点，并让 transition 与新增线段一一对应。"""

        if generated and abs(generated[-1] - point).contains(0):
            return
        if generated:
            transitions.append(transition or {"method": "ordinary_taylor"})
        generated.append(acb(point))
        generated_classifications.append(classification)

    previous_right_classification: PointClassification | None = None
    for interval_index, ((left, left_classification), (right, right_classification)) in enumerate(
        zip(checkpoints[:-1], checkpoints[1:])
    ):
        left_value = left.finite_acb()
        right_value = right.finite_acb()
        delta = right_value - left_value
        interval_length = abs(delta)
        unit = delta / interval_length
        left_work = _ordinary_match_point(
            left_value,
            left_classification,
            unit,
            interval_length,
            max_step_over_radius,
            (
                formal_start_distance
                if interval_index == 0
                and left_classification is route_classifications[0]
                else None
            ),
        )
        right_work = _ordinary_match_point(
            right_value,
            right_classification,
            -unit,
            interval_length,
            max_step_over_radius,
        )
        for endpoint, classification, match in (
            (left, left_classification, left_work),
            (right, right_classification, right_work),
        ):
            if classification.kind != "ordinary":
                endpoint_matches.append(
                    {
                        "endpoint_name": endpoint.name,
                        "classification": classification.kind,
                        "singular_point": classification.value_label,
                        "ordinary_match_point": match.str(40),
                        "convergence_radius": (
                            None
                            if classification.convergence_radius is None
                            else classification.convergence_radius.str(40)
                        ),
                        "match_distance_over_convergence_radius": (
                            None
                            if classification.convergence_radius is None
                            else float((abs(match - endpoint.finite_acb()) / classification.convergence_radius).mid())
                        ),
                        "formal_asymptotic_match_estimate": (
                            formal_start_estimate
                            if interval_index == 0 and endpoint.name == route[0].name
                            else None
                        ),
                    }
                )
        if abs(right_work - left_work).contains(0):
            raise ValueError(
                f"path interval {interval_index + 1} has no ordinary region after matching"
            )
        left_work_classification = classify_point(
            NamedPoint(f"{resolved_name}_match_{interval_index + 1:03d}_left", left_work),
            working_inventory,
        )
        bridge_transition = None
        if generated and previous_right_classification is not None:
            bridge_transition = {
                "method": "regular_singular_bridge",
                "singularity": previous_right_classification,
            }
        append_unique(left_work, left_work_classification, bridge_transition)
        if working_inventory.finite:
            leg_points = build_straight_path(
                analytic_system,
                left_work,
                right_work,
                step_fraction=max_step_over_radius,
            )
        else:
            leg_points = [acb(left_work), acb(right_work)]
        for point in leg_points[1:]:
            local_classification = classify_point(
                NamedPoint(f"{resolved_name}_ordinary", point), working_inventory
            )
            append_unique(point, local_classification, {"method": "ordinary_taylor"})
        previous_right_classification = (
            right_classification if right_classification.kind != "ordinary" else None
        )

    if len(generated) < 2:
        raise ValueError("adaptive path requires at least two ordinary path points")
    path_records: list[dict[str, Any]] = []
    segment_records: list[dict[str, Any]] = []
    for index, (point, classification) in enumerate(
        zip(generated, generated_classifications), start=1
    ):
        path_records.append(
            {
                "point_index": index,
                "value": point.str(40),
                "convergence_radius": (
                    None
                    if classification.convergence_radius is None
                    else classification.convergence_radius.str(40)
                ),
            }
        )
        if index == len(generated):
            continue
        transition = transitions[index - 1]
        target_point = generated[index]
        if transition["method"] == "regular_singular_bridge":
            singularity = transition["singularity"]
            center = singularity.value
            if center is None:
                raise ValueError("finite Frobenius bridge is missing its center")
            for side, step_start, step_target in (
                ("match_to_singularity", point, center),
                ("singularity_to_match", center, target_point),
            ):
                step_length = abs(step_target - step_start)
                radius = singularity.convergence_radius
                ratio = None if radius is None else float((step_length / radius).mid())
                segment_records.append(
                    {
                        "step_index": len(segment_records) + 1,
                        "execution_segment_index": index,
                        "method": "local_singular_basis",
                        "side": side,
                        "singularity_identifier": singularity.singularity_identifier,
                        "step_length": step_length.str(40),
                        "controlling_convergence_radius": (
                            None if radius is None else radius.str(40)
                        ),
                        "radius_owner": "singularity",
                        "step_over_convergence_radius": ratio,
                    }
                )
        else:
            step_length = abs(target_point - point)
            radius = classification.convergence_radius
            ratio = None if radius is None else float((step_length / radius).mid())
            segment_records.append(
                {
                    "step_index": len(segment_records) + 1,
                    "execution_segment_index": index,
                    "method": "ordinary_taylor",
                    "side": None,
                    "singularity_identifier": None,
                    "step_length": step_length.str(40),
                    "controlling_convergence_radius": (
                        None if radius is None else radius.str(40)
                    ),
                    "radius_owner": "start_ordinary_point",
                    "step_over_convergence_radius": ratio,
                }
            )
    finite_ratios = [
        record["step_over_convergence_radius"]
        for record in segment_records
        if record["step_over_convergence_radius"] is not None
    ]
    actual_maximum_ratio = max(finite_ratios, default=None)
    if actual_maximum_ratio is not None and actual_maximum_ratio > max_step_over_radius + 1.0e-12:
        raise ArithmeticError("generated path violates max_step_over_radius")
    ratio_message = (
        "step/R ratios in path order: "
        + repr([record["step_over_convergence_radius"] for record in segment_records])
        + f"; maximum step/R: {actual_maximum_ratio}"
    )
    messages.append(ratio_message)
    save_requests: list[dict[str, Any]] = []
    original_route = (start, *detour_points, target)
    for route_index, (original_point, working_point, classification, save_requested) in enumerate(
        zip(original_route, route, route_classifications, save_flags)
    ):
        if not save_requested:
            continue
        request: dict[str, Any] = {
            "coordinate": "inf" if original_point.is_infinity else original_point.finite_acb().str(40),
            "working_coordinate": classification.value_label,
            "classification": classification.kind,
            "singularity_identifier": classification.singularity_identifier,
            "role": (
                "start" if route_index == 0 else
                "target" if route_index == len(original_route) - 1 else
                "detour"
            ),
            "execution_index": None,
        }
        if classification.kind == "ordinary":
            working_value = working_point.finite_acb()
            request["execution_index"] = next(
                (
                    index
                    for index, generated_point in enumerate(generated)
                    if abs(generated_point - working_value).contains(0)
                ),
                None,
            )
            if request["execution_index"] is None:
                raise ArithmeticError(
                    "saved ordinary route coordinate is absent from the execution path"
                )
        save_requests.append(request)
    warnings.warn("; ".join(messages), UserWarning, stacklevel=2)
    manifest = {
        "schema": "flintnde_executable_path_v1",
        "planning_action": "raw_points_automatic_plan",
        "singularity_mode": resolved_singularity_mode,
        "message_language": language,
        "path_name": resolved_name,
        "original_variable": system.variable_name,
        "working_variable": working.variable_name,
        "infinity_transformation": inverted,
        "max_step_over_convergence_radius": max_step_over_radius,
        "start": route_classifications[0].to_json(),
        "target": route_classifications[-1].to_json(),
        "detour_points": [item.to_json() for item in route_classifications[1:-1]],
        "internal_singularities": [item.to_json() for item in internal_classifications],
        "endpoint_matches": endpoint_matches,
        "path_points": path_records,
        "steps": segment_records,
        "step_over_convergence_radius_list": [
            record["step_over_convergence_radius"] for record in segment_records
        ],
        "actual_max_step_over_convergence_radius": actual_maximum_ratio,
        "formal_asymptotic_match_estimate": formal_start_estimate,
        "save_requests": save_requests,
        "direct_transport_path_type": "list[acb]",
        "messages": messages,
    }
    if output_layout is not None:
        output_layout.write_json("transport", f"{resolved_name}_path.json", manifest)
    return AdaptivePath(
        generated,
        working_system=working,
        transitions=transitions,
        step_reports=segment_records,
        internal_singularities=bridge_classifications,
        start_classification=route_classifications[0],
        target_classification=route_classifications[-1],
        working_variable=working.variable_name,
        infinity_transformation=inverted,
        save_requests=tuple(save_requests),
        formal_asymptotic_match_estimate=formal_start_estimate,
        planning_action="raw_points_automatic_plan",
        singularity_mode=resolved_singularity_mode,
        message_language=language,
    )


def build_adaptive_path_plan(
    system: RationalMatrixSystem,
    start: NamedPoint,
    target: NamedPoint,
    *,
    path_name: str | None = None,
    max_step_over_radius: float = 0.45,
    path_tolerance: float = 1.0e-24,
    local_analysis_order: int = 4,
    singularity_mode: str = "avoid",
    message_language: str = "EN",
    output_layout: Any | None = None,
) -> PathPlan:
    """生成命名直线路径报告，并明确其是否可继续执行。

    缺省 ``singularity_mode="avoid"`` 在直线含内部奇点时保留逐点诊断，但把
    ``continuation_ready`` 置为 False；显式 ``"singularity_jump"`` 才预检 power-log、
    shearing 或指数 sector 局部基，规划器从不在 pole 上计算 Taylor 级数，也不把
    未通过门禁的高阶 pole 自动认证为真正不规则奇点。
    """

    language = _normalize_message_language(message_language)
    resolved_singularity_mode = _normalize_singularity_mode(singularity_mode)
    if not 0 < max_step_over_radius < 1:
        raise ValueError("max_step_over_radius must lie between zero and one")
    if path_tolerance <= 0:
        raise ValueError("path_tolerance must be positive")
    if local_analysis_order <= 0:
        raise ValueError("local_analysis_order must be positive")
    resolved_name = path_name or f"{start.name}_to_{target.name}"
    working, working_start, working_target, inverted = _map_path_endpoints_at_infinity(
        system, start, target
    )
    original_inventory = analyze_singularities(system, output_layout=output_layout)
    working_inventory = (
        analyze_singularities(
            working,
            output_layout=output_layout,
            filename="singularity_inventory_sinv.json",
        )
        if inverted
        else original_inventory
    )
    start_classification = classify_point(working_start, working_inventory)
    target_classification = classify_point(working_target, working_inventory)
    start_value = working_start.finite_acb()
    target_value = working_target.finite_acb()
    if abs(target_value - start_value).contains(0):
        raise ValueError("start and target must be distinct points")

    internal_with_parameter: list[tuple[arb, PointClassification]] = []
    for record in working_inventory.finite:
        if record.location is None:
            continue
        parameter = _on_segment_parameter(
            record.location, start_value, target_value, path_tolerance
        )
        if parameter is None:
            continue
        endpoint = NamedPoint(record.identifier, record.location)
        internal_with_parameter.append((parameter, classify_point(endpoint, working_inventory)))
    internal_with_parameter.sort(key=lambda item: item[0])
    internal = tuple(item[1] for item in internal_with_parameter)

    internal_identifiers = [
        item.singularity_identifier or item.name for item in internal
    ]
    if language == "CN":
        messages = [
            "输入为原始端点，已生成单变量路径计划；已有计划应直接交给执行函数。",
            (
                "当前显式使用奇点折跃模式；所选多值分支等价于某一绕行路径，"
                "必须由用户确认。"
                if resolved_singularity_mode == "singularity_jump"
                else "当前使用避开奇点模式（缺省）。"
            ),
            f"起点 '{start.name}' 分类为 {start_classification.kind}",
            f"终点 '{target.name}' 分类为 {target_classification.kind}",
        ]
    else:
        messages = [
            "Raw endpoints were supplied, so a one-variable path plan was generated; "
            "pass an existing plan directly to the execution function.",
            (
                "Singularity mode: singularity jump (explicit). The selected multivalued branch is "
                "equivalent to a detour path and must be confirmed by the user."
                if resolved_singularity_mode == "singularity_jump"
                else "Singularity mode: avoid (default)."
            ),
            f"start '{start.name}' is {start_classification.kind}",
            f"target '{target.name}' is {target_classification.kind}",
        ]
    if inverted:
        messages.append(
            "无穷远使用 sinv=1/s 与 d/dsinv=-sinv^(-2) d/ds"
            if language == "CN"
            else "infinity handled with sinv=1/s and d/dsinv=-sinv^(-2) d/ds"
        )
    local_methods: dict[str, str] = {}
    unresolved_local_points: dict[str, str] = {}
    classified_roles = (
        (("start", start_classification),)
        + tuple(("internal", item) for item in internal)
        + (("target", target_classification),)
    )
    for role, classification in classified_roles:
        if classification.kind == "ordinary":
            continue
        if role == "internal" and resolved_singularity_mode == "avoid":
            continue
        identifier = classification.singularity_identifier or classification.name
        try:
            _record, local_basis = _preflight_local_basis(
                working, working_inventory, classification, local_analysis_order
            )
        except LocalReductionError as error:
            unresolved_local_points[identifier] = str(error)
        else:
            local_methods[identifier] = local_basis.method
            if not local_basis.continuation_ready and role != "start":
                unresolved_local_points[identifier] = (
                    f"{local_basis.method} is start-only; Stokes connection data are unavailable"
                )

    if internal:
        messages.append(
            (
                "直线路径上的内部奇点：" + ", ".join(internal_identifiers)
                if language == "CN"
                else "internal singularities on the straight path: "
                + ", ".join(internal_identifiers)
            )
        )
        if resolved_singularity_mode == "avoid":
            messages.append(
                (
                    "该计划不可继续执行；请给出明确绕行路径，或显式选择 singularity_jump 奇点折跃模式。"
                    if language == "CN"
                    else "This plan is not continuation-ready; supply an explicit detour path "
                    "or select singularity_jump mode explicitly."
                )
            )
        elif not unresolved_local_points:
            messages.append(
                "所有内部奇点均有认证的自动局部基"
                if language == "CN"
                else "all internal singularities have an automatic certified local basis"
            )
    else:
        messages.append(
            "直线路径不含内部奇点"
            if language == "CN"
            else "no internal singularity lies on the straight path"
        )
    if unresolved_local_points:
        messages.append(
            (
                "不支持的局部奇点："
                if language == "CN"
                else "unsupported local singularities: "
            )
            + "; ".join(
                f"{identifier}: {reason}"
                for identifier, reason in unresolved_local_points.items()
            )
        )
    warnings.warn("; ".join(messages), UserWarning, stacklevel=2)

    checkpoints: list[tuple[float, NamedPoint, PointClassification]] = [
        (0.0, working_start, start_classification)
    ]
    for parameter, classification in internal_with_parameter:
        checkpoints.append(
            (parameter, NamedPoint(classification.name, classification.value), classification)
        )
    checkpoints.append((1.0, working_target, target_classification))
    direction = target_value - start_value
    total_distance = abs(direction)
    unit = direction / total_distance
    analytic_system = working.to_analytic_system(working_inventory)
    generated: list[tuple[str, acb, PointClassification]] = []
    ordinary_counter = 0

    def append_point(name: str, point: acb, classification: PointClassification) -> None:
        if generated and abs(generated[-1][1] - point).contains(0):
            return
        generated.append((name, acb(point), classification))

    for interval_index, (left, right) in enumerate(zip(checkpoints[:-1], checkpoints[1:])):
        left_parameter, left_named, left_classification = left
        right_parameter, right_named, right_classification = right
        left_point = start_value + direction * left_parameter
        right_point = start_value + direction * right_parameter
        interval_length = abs(right_point - left_point)
        if interval_length.contains(0):
            continue
        append_point(left_named.name, left_point, left_classification)

        left_work = acb(left_point)
        if left_classification.kind != "ordinary":
            candidate = interval_length / 3
            if left_classification.convergence_radius is not None:
                candidate = min(
                    candidate,
                    left_classification.convergence_radius * arb(str(max_step_over_radius)),
                    key=lambda value: float(value.mid()),
                )
            left_work = left_point + unit * candidate

        right_work = acb(right_point)
        if right_classification.kind != "ordinary":
            candidate = interval_length / 3
            if right_classification.convergence_radius is not None:
                candidate = min(
                    candidate,
                    right_classification.convergence_radius * arb(str(max_step_over_radius)),
                    key=lambda value: float(value.mid()),
                )
            right_work = right_point - unit * candidate

        if not abs(left_work - left_point).contains(0):
            ordinary_counter += 1
            local = classify_point(NamedPoint(f"{resolved_name}_ordinary_{ordinary_counter:03d}", left_work), working_inventory)
            append_point(local.name, left_work, local)
        if not abs(right_work - left_work).contains(0):
            ordinary_path = build_straight_path(
                analytic_system,
                left_work,
                right_work,
                step_fraction=max_step_over_radius,
            )
            for point in ordinary_path[1:]:
                ordinary_counter += 1
                local = classify_point(
                    NamedPoint(f"{resolved_name}_ordinary_{ordinary_counter:03d}", point),
                    working_inventory,
                )
                append_point(local.name, point, local)
        if not abs(right_work - right_point).contains(0):
            append_point(right_named.name, right_point, right_classification)
        elif interval_index == len(checkpoints) - 2:
            append_point(right_named.name, right_point, right_classification)

    points = tuple(_point_record(name, point, classification) for name, point, classification in generated)
    segments: list[dict[str, Any]] = []
    for index, (left, right) in enumerate(zip(generated[:-1], generated[1:]), start=1):
        if left[2].kind != "ordinary":
            method = (
                "regular_singular_power_log_launch"
                if left[2].kind == "regular_singular"
                else "automatic_high_pole_local_basis"
            )
        elif right[2].kind != "ordinary":
            method = (
                "regular_singular_power_log_match"
                if right[2].kind == "regular_singular"
                else "automatic_high_pole_local_basis"
            )
        else:
            method = "ordinary_taylor"
        controlling_radius = (
            left[2].convergence_radius
            if left[2].kind != "ordinary"
            else right[2].convergence_radius
            if right[2].kind != "ordinary"
            else left[2].convergence_radius
        )
        step_length = abs(right[1] - left[1])
        segments.append(
            {
                "segment_index": index,
                "start_name": left[0],
                "target_name": right[0],
                "method": method,
                "step_length": step_length.str(40),
                "controlling_convergence_radius": (
                    None if controlling_radius is None else controlling_radius.str(40)
                ),
                "step_over_convergence_radius": (
                    None
                    if controlling_radius is None
                    else float((step_length / controlling_radius).mid())
                ),
            }
        )
    continuation_ready = (
        not unresolved_local_points
        and not (internal and resolved_singularity_mode == "avoid")
    )
    if local_methods:
        messages.append(
            "certified local methods: "
            + ", ".join(
                f"{identifier}={method}" for identifier, method in local_methods.items()
            )
        )
    plan = PathPlan(
        resolved_name,
        system.variable_name,
        working.variable_name,
        inverted,
        max_step_over_radius,
        start_classification,
        target_classification,
        internal,
        points,
        tuple(segments),
        "raw_points_automatic_plan",
        resolved_singularity_mode,
        language,
        continuation_ready,
        tuple(messages),
    )
    if output_layout is not None:
        output_layout.write_json("transport", f"{resolved_name}_path.json", plan.to_json())
    return plan
