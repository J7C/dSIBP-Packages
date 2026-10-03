"""MadStree 到 FlintNDE 的数值适配器。

MadStree 只提交连续复仿射单变量段、exact dlog 拉回、边界和 master 顺序。
FlintNDE 在本进程内完成边界输运、各段自动规划或顺序直输、dense 多点求值和
奇点领头阶解析。正规化级数控制协议接收 MadStree 的符号最低阶证书，生成自适应
exact ep 网格，并对 MadStree 回传的 Acb 终点球作 Laurent 重构和独立验证。
"""

from __future__ import annotations

from multiprocessing import get_context
from decimal import Decimal
from fractions import Fraction
import json
import math
import os
import re
import sys
import time
import traceback
import warnings
from pathlib import Path
from typing import Any


EVALUATE_SCHEMA = "madstree_flintnde_evaluate_v1"
EVALUATE_BATCH_SCHEMA = "madstree_flintnde_ep_batch_v1"
SERIES_CONTROL_SCHEMA = "madstree_flintnde_ep_series_control_v1"
DEFAULT_PARALLEL_TASK_COUNT = 12
REQUEST_KEYS = {
    "schema", "backendPackagePath", "masterDigest", "dimension", "segments",
    "pathPlanning", "singularityMode", "boundary",
    "workingPrecisionDigits", "primaryOrder", "referenceOrder",
    "targetRelativeError", "certificationMode", "messageLanguage",
    "columnVectorConvention", "dlogStatus",
}
MADSTREE_CITATION_NOTICE = {
    "intro": (
        "To support MadStree, please acknowledge our work when you use MadStree "
        "in your research. You can do so by citing the following papers:"
    ),
    "entries": [
        {
            "authors": "Jiaqi Chen and Bo Feng",
            "title": (
                "Towards Systematic Evaluation of de Sitter Correlators via "
                "Generalized Integration-By-Parts Relations"
            ),
            "arXiv": "2401.00129",
            "url": "https://arxiv.org/abs/2401.00129",
        },
        {
            "authors": "Jiaqi Chen, Bo Feng and Yi-Xiao Tao",
            "title": (
                "Multivariate hypergeometric solutions of cosmological (dS) "
                "correlators by d log-form differential equations"
            ),
            "arXiv": "2411.03088",
            "url": "https://arxiv.org/abs/2411.03088",
        },
    ],
    "packagePaper": "MadStree package paper, arXiv identifier pending.",
}

def _citation_text() -> str:
    """生成与 Wolfram 入口同序的人工 CLI 引用文本。"""

    lines = [MADSTREE_CITATION_NOTICE["intro"]]
    for index, entry in enumerate(MADSTREE_CITATION_NOTICE["entries"], start=1):
        lines.append(
            f'{index}. {entry["authors"]}, "{entry["title"]}", '
            f'arXiv:{entry["arXiv"]}, {entry["url"]}.'
        )
    lines.append(f'3. {MADSTREE_CITATION_NOTICE["packagePaper"]}')
    return "\n".join(lines)


def _with_citation_metadata(result: dict[str, Any]) -> dict[str, Any]:
    """为每个文件协议响应附加稳定引用元数据。"""

    enriched = dict(result)
    enriched["metadata"] = {"citationNotice": MADSTREE_CITATION_NOTICE}
    return enriched


def _exact_keys(value: Any, expected: set[str], label: str) -> dict[str, Any]:
    """要求对象字段与当前唯一 schema 完全一致。"""

    if not isinstance(value, dict):
        raise TypeError(f"{label} must be a JSON object")
    actual = set(value)
    if actual != expected:
        raise ValueError(
            f"{label} fields do not match the current schema; "
            f"missing={sorted(expected - actual)}, unexpected={sorted(actual - expected)}"
        )
    return value


def _integer(value: Any, label: str, minimum: int = 1) -> int:
    """读取严格整数，拒绝 bool 冒充整数。"""

    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise TypeError(f"{label} must be an integer >= {minimum}")
    return value


def _string(value: Any, label: str) -> str:
    """读取非空字符串。"""

    if not isinstance(value, str) or not value:
        raise TypeError(f"{label} must be a nonempty string")
    return value


def _complex_record(value: Any, label: str) -> None:
    """验证 Q(i) 或十进制复数记录。"""

    item = _exact_keys(value, {"real", "imag"}, label)
    _string(item["real"], f"{label}.real")
    _string(item["imag"], f"{label}.imag")


def _matrix(value: Any, dimension: int, label: str, validator: Any) -> None:
    """验证指定元素合同的方阵。"""

    if not isinstance(value, list) or len(value) != dimension:
        raise ValueError(f"{label} must be a {dimension} by {dimension} matrix")
    for row_index, row in enumerate(value):
        if not isinstance(row, list) or len(row) != dimension:
            raise ValueError(f"{label}[{row_index}] must have length {dimension}")
        for column_index, entry in enumerate(row):
            validator(entry, f"{label}[{row_index}][{column_index}]")


def _rational_function_record(value: Any, label: str) -> None:
    """验证一个 Q(i)(s) 有理函数系数记录。"""

    item = _exact_keys(value, {"numerator", "denominator"}, label)
    for key in ("numerator", "denominator"):
        coefficients = item[key]
        if not isinstance(coefficients, list) or not coefficients:
            raise ValueError(f"{label}.{key} must be a nonempty list")
        for index, coefficient in enumerate(coefficients):
            _complex_record(coefficient, f"{label}.{key}[{index}]")


def _letters(value: Any, dimension: int, label: str) -> None:
    """验证仿射 dlog letters。"""

    if not isinstance(value, list):
        raise TypeError(f"{label} must be a list")
    for index, record in enumerate(value):
        item = _exact_keys(record, {"alpha", "beta", "residue"}, f"{label}[{index}]")
        _complex_record(item["alpha"], f"{label}[{index}].alpha")
        _complex_record(item["beta"], f"{label}[{index}].beta")
        _matrix(item["residue"], dimension, f"{label}[{index}].residue", _complex_record)


def _validate_segment(value: Any, dimension: int, label: str) -> None:
    """验证一个连续复仿射段。"""

    item = _exact_keys(
        value, {"start", "points", "letters", "fromUserIndex", "userIndices"}, label
    )
    _string(item["start"], f"{label}.start")
    _integer(item["fromUserIndex"], f"{label}.fromUserIndex", 0)
    if not isinstance(item["points"], list) or not item["points"]:
        raise ValueError(f"{label}.points must be a nonempty list")
    if not isinstance(item["userIndices"], list) or len(item["userIndices"]) != len(item["points"]):
        raise ValueError(f"{label}.userIndices must match points")
    for index, point in enumerate(item["points"]):
        _complex_record(point, f"{label}.points[{index}]")
        _integer(item["userIndices"][index], f"{label}.userIndices[{index}]")
    if len(set(item["userIndices"])) != len(item["userIndices"]):
        raise ValueError(f"{label}.userIndices must be distinct")
    _letters(item["letters"], dimension, f"{label}.letters")


def _validate_boundary(value: Any, dimension: int) -> None:
    """验证有限边界或正则奇点边界。"""

    if not isinstance(value, dict):
        raise TypeError("boundary must be a JSON object")
    kind = value.get("kind")
    if kind == "finite":
        item = _exact_keys(value, {"kind", "values"}, "boundary")
        if not isinstance(item["values"], list) or len(item["values"]) != dimension:
            raise ValueError("finite boundary vector dimension mismatch")
        for index, record in enumerate(item["values"]):
            _complex_record(record, f"boundary.values[{index}]")
        return
    item = _exact_keys(
        value,
        {"kind", "variable", "matrix", "branches", "weights", "start", "target"},
        "boundary",
    )
    if kind != "regular_singular":
        raise ValueError('boundary.kind must be "finite" or "regular_singular"')
    _string(item["variable"], "boundary.variable")
    _string(item["start"], "boundary.start")
    _string(item["target"], "boundary.target")
    _matrix(item["matrix"], dimension, "boundary.matrix", _rational_function_record)
    if not isinstance(item["branches"], list) or not item["branches"]:
        raise ValueError("boundary.branches must be a nonempty list")
    if not isinstance(item["weights"], list) or len(item["weights"]) != len(item["branches"]):
        raise ValueError("boundary.weights must match branches")
    for index, branch in enumerate(item["branches"]):
        record = _exact_keys(branch, {"a", "b", "C"}, f"boundary.branches[{index}]")
        _complex_record(record["a"], f"boundary.branches[{index}].a")
        _integer(record["b"], f"boundary.branches[{index}].b", 0)
        if not isinstance(record["C"], list) or len(record["C"]) != dimension:
            raise ValueError("boundary branch vector dimension mismatch")
        for column, coefficient in enumerate(record["C"]):
            _complex_record(coefficient, f"boundary.branches[{index}].C[{column}]")
        _complex_record(item["weights"][index], f"boundary.weights[{index}]")


def _validate_request(value: Any) -> dict[str, Any]:
    """验证唯一请求格式；旧 schema 一律拒绝。"""

    data = _exact_keys(value, REQUEST_KEYS, "request")
    if data["schema"] != EVALUATE_SCHEMA:
        raise ValueError("unsupported MadStree-FlintNDE schema")
    _string(data["backendPackagePath"], "backendPackagePath")
    _string(data["masterDigest"], "masterDigest")
    dimension = _integer(data["dimension"], "dimension")
    if not isinstance(data["pathPlanning"], bool):
        raise TypeError("pathPlanning must be a JSON boolean")
    if data["singularityMode"] not in {"avoid", "singularity_jump"}:
        raise ValueError('singularityMode must be "avoid" or "singularity_jump"')
    if not data["pathPlanning"] and data["singularityMode"] != "avoid":
        raise ValueError("direct user-point transport supports avoid mode only")
    if data["messageLanguage"] not in {"EN", "CN"}:
        raise ValueError('messageLanguage must be "EN" or "CN"')
    _integer(data["workingPrecisionDigits"], "workingPrecisionDigits", 20)
    primary = _integer(data["primaryOrder"], "primaryOrder")
    reference = _integer(data["referenceOrder"], "referenceOrder")
    if reference <= primary:
        raise ValueError("referenceOrder must exceed primaryOrder")
    _string(data["targetRelativeError"], "targetRelativeError")
    if data["certificationMode"] not in {"embedded", "certified"}:
        raise ValueError('certificationMode must be "embedded" or "certified"')
    if data["columnVectorConvention"] != "Y'=A(s)Y":
        raise ValueError("columnVectorConvention mismatch")
    if data["dlogStatus"] != "certifiedByFormulaChecks":
        raise ValueError("uncertified dlog input")
    if not isinstance(data["segments"], list) or not data["segments"]:
        raise ValueError("segments must be a nonempty list")
    for index, segment in enumerate(data["segments"]):
        _validate_segment(segment, dimension, f"segments[{index}]")
    _validate_boundary(data["boundary"], dimension)
    return data


def _acb_record(value: Any, digits: int) -> dict[str, Any]:
    """同时写出用户可读中点和用于正规化认证的完整 Arb 球。"""

    from flintnde.core import arb_ball_to_json

    return {
        "real": value.real.mid().str(digits, radius=False, more=True),
        "imag": value.imag.mid().str(digits, radius=False, more=True),
        "real_ball": arb_ball_to_json(value.real, digits),
        "imag_ball": arb_ball_to_json(value.imag, digits),
    }


def _point_record(value: Any, digits: int) -> dict[str, Any]:
    """序列化实际输运节点坐标。"""

    return _acb_record(value, digits)


def _resolve_backend_path(data: dict[str, Any]) -> Path:
    """解析 FlintNDE 包并把它加入本进程导入路径。"""

    backend_path = Path(data["backendPackagePath"]).resolve()
    package_file = backend_path / "flintnde" / "__init__.py"
    if not package_file.is_file():
        raise FileNotFoundError(f"FlintNDE package not found: {package_file}")
    if str(backend_path) not in sys.path:
        sys.path.insert(0, str(backend_path))
    return backend_path


def _letters_to_system(records: list[dict[str, Any]], dimension: int, name: str) -> Any:
    """把 MadStree 仿射 dlog records 转成 FlintNDE 极点--留数系统。"""

    from flint import acb_mat
    from flintnde import PartialFractionSystem, gaussian_rational

    merged: list[tuple[Any, Any]] = []
    for record in records:
        beta = gaussian_rational(record["beta"])
        if beta.is_zero:
            continue
        alpha = gaussian_rational(record["alpha"])
        if alpha.is_zero:
            raise ValueError("dlog letter vanishes at the segment anchor")
        from flint import acb

        residue = acb_mat(
            [[acb(item["real"], item["imag"]) for item in row] for row in record["residue"]]
        )
        pole = -alpha / beta
        for index, (old_pole, old_residue) in enumerate(merged):
            if old_pole == pole:
                merged[index] = (old_pole, old_residue + residue)
                break
        else:
            merged.append((pole, residue))
    poles, residues = zip(*merged) if merged else ((), ())
    return PartialFractionSystem(
        constant=acb_mat(dimension, dimension),
        residues=tuple(residues),
        poles=tuple(pole.to_acb() for pole in poles),
        name=name,
    )


def _rational_system(boundary: dict[str, Any], dimension: int, name: str) -> Any:
    """重建正则奇点边界的一般 Q(i)(t) 系统。"""

    from flintnde import RationalMatrixSystem, rational_function

    matrix = tuple(
        tuple(rational_function(entry["numerator"], entry["denominator"]) for entry in row)
        for row in boundary["matrix"]
    )
    return RationalMatrixSystem(matrix, variable_name=boundary["variable"], name=name)


def _initial_vectors(
    data: dict[str, Any], digits: int
) -> tuple[Any, Any | None, dict[str, Any]]:
    """分别构造主链与参考链初值；有限精确边界只需保存一份。"""

    from flint import acb, acb_mat
    from flintnde import (
        NamedPoint, build_adaptive_path, frobenius_boundary,
        gaussian_rational, transport_frobenius_boundaries_refined,
    )

    boundary = data["boundary"]
    dimension = data["dimension"]
    if boundary["kind"] == "finite":
        vector = acb_mat([[acb(item["real"], item["imag"])] for item in boundary["values"]])
        return vector, None, {"kind": "finite", "path": [], "seconds": 0.0}
    system = _rational_system(boundary, dimension, data["masterDigest"] + "-boundary")
    start = gaussian_rational(boundary["start"]).to_acb()
    target = gaussian_rational(boundary["target"]).to_acb()
    path = build_adaptive_path(
        system, NamedPoint("madstree_boundary", start), NamedPoint("madstree_anchor", target),
        path_name="madstree_boundary_to_anchor", max_step_over_radius=0.45,
        singularity_mode="avoid", message_language=data["messageLanguage"],
    )
    branches = [frobenius_boundary([record]) for record in boundary["branches"]]
    clock = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = transport_frobenius_boundaries_refined(
            system, branches, path,
            primary_order=data["primaryOrder"], reference_order=data["referenceOrder"],
            radius_fraction=0.60, target_relative_error=data["targetRelativeError"],
            certification_mode=data["certificationMode"],
        )
    weights = [acb(item["real"], item["imag"]) for item in boundary["weights"]]

    def weighted_vector(matrix: Any) -> Any:
        """按同一物理边界权重合并一条输运链的基解矩阵。"""

        vector = acb_mat(dimension, 1)
        for row in range(dimension):
            vector[row, 0] = sum(
                (weights[column] * matrix[row, column] for column in range(len(weights))),
                acb(0),
            )
        return vector

    primary_vector = weighted_vector(result["primary_snapshots"][-1])
    reference_vector = weighted_vector(result["reference_snapshots"][-1])
    return primary_vector, reference_vector, {
        "kind": "regular_singular",
        "path": [_point_record(point, digits) for point in path],
        "pathPointCount": len(path),
        "seconds": time.perf_counter() - clock,
        "warnings": [str(item.message) for item in caught],
    }


def _point_assignments(plan: Any, points: list[Any], user_indices: list[int]) -> list[dict[str, Any]]:
    """把每个用户点映射到节点快照或 dense sample。"""

    sample_by_index = {
        int(item["user_point_index"]): item for item in plan.sample_assignments
        if item["user_point_index"] is not None
    }
    assignments = []
    for local_index, (point, user_index) in enumerate(zip(points, user_indices)):
        if local_index in sample_by_index:
            sample = sample_by_index[local_index]
            assignments.append({
                "userIndex": user_index, "localPointIndex": local_index,
                "source": sample["source"], "nodeIndex": None,
                "segmentIndex": int(sample["segment_index"]),
            })
            continue
        nodes = [index for index, node in enumerate(plan.nodes) if abs(node - point).contains(0)]
        if not nodes:
            raise ValueError(f"user point {user_index} has no FlintNDE value assignment")
        assignments.append({
            "userIndex": user_index, "localPointIndex": local_index,
            "source": "node_snapshot", "nodeIndex": nodes[-1], "segmentIndex": None,
        })
    return assignments


def _point_values(result: dict[str, Any], assignments: list[dict[str, Any]], digits: int, dimension: int) -> list[dict[str, Any]]:
    """按映射合并节点快照与批量 dense 输出。"""

    dense = {
        int(item["user_point_index"]): item for item in result["sample_results"]
        if item["user_point_index"] is not None
    }
    singular = {
        int(item["user_point_index"]): item
        for item in result.get("singular_target_results", [])
        if item["user_point_index"] is not None
    }
    records = []
    for assignment in assignments:
        local_index = assignment["localPointIndex"]
        node_index = assignment["nodeIndex"]
        singular_item = singular.get(local_index)
        if singular_item is not None:
            records.append({
                "userIndex": assignment["userIndex"],
                "localPointIndex": local_index,
                "source": assignment["source"],
                "evaluationAlgorithm": None,
                "values": [
                    _value_record(value, digits)
                    for value in singular_item["values"]
                ],
            })
            continue
        item = dense.get(local_index) if node_index is None else None
        vector = item["value"] if item is not None else result["primary_snapshots"][node_index]
        records.append({
            "userIndex": assignment["userIndex"], "localPointIndex": local_index,
            "source": assignment["source"],
            "evaluationAlgorithm": None if item is None else item.get("evaluation_algorithm"),
            "values": [_acb_record(vector[row, 0], digits) for row in range(dimension)],
        })
    return records


def _value_record(value: Any, digits: int) -> dict[str, Any]:
    """把有限 Acb 或奇点状态转换为统一的用户结果记录。"""

    if isinstance(value, str):
        return {"text": value}
    return _acb_record(value, digits)


def _singular_result_difference(primary: Any, reference: Any) -> tuple[Any, str | None]:
    """比较末端奇点主/参考结果；发散状态只允许逐分量完全一致。"""

    from flint import arb

    if primary.component_classifications != reference.component_classifications:
        return arb(0), "primary/reference singular-target classifications do not agree"
    difference = arb(0)
    for primary_value, reference_value in zip(primary.values, reference.values):
        if isinstance(primary_value, str) or isinstance(reference_value, str):
            if primary_value != reference_value:
                return difference, "primary/reference singular-target divergence states do not agree"
            continue
        denominator = max(
            arb(1), abs(reference_value), key=lambda item: float(item.mid())
        )
        difference = max(difference, abs(primary_value - reference_value) / denominator)
    return difference, None


def _letters_to_rational_system(records: list[dict[str, Any]], dimension: int, name: str) -> Any:
    """从 exact dlog letters 重建局部幂对数分析所需的有理矩阵系统。"""

    from flintnde import RationalMatrixSystem, gaussian_rational, rational_function

    rows = [[rational_function(0) for _ in range(dimension)] for _ in range(dimension)]
    for record in records:
        alpha = gaussian_rational(record["alpha"])
        beta = gaussian_rational(record["beta"])
        if beta.is_zero:
            continue
        denominator = (alpha, beta)
        for row in range(dimension):
            for column in range(dimension):
                item = record["residue"][row][column]
                residue = gaussian_rational(item)
                # record 表示 R dlog(alpha+beta*s)，局部系统必须保留导数因子 beta。
                rows[row][column] += rational_function(beta * residue, denominator)
    return RationalMatrixSystem(tuple(tuple(row) for row in rows), variable_name="s", name=name)


def _run(data: dict[str, Any]) -> dict[str, Any]:
    """在单一 Python 进程内完成全部边界和复仿射段。"""

    backend_path = _resolve_backend_path(data)
    from flint import acb, acb_mat, arb
    from flintnde import (
        SingularPathError, configure_working_precision, direct_user_point_path,
        evaluate_singular_target, gaussian_rational, plan_singular_target_match,
        plan_transport_path, transport_planned_path_refined,
    )

    digits = data["workingPrecisionDigits"]
    configure_working_precision(digits)
    primary_vector, reference_vector, boundary_report = _initial_vectors(data, digits)
    if reference_vector is None:
        reference_vector = acb_mat(primary_vector)
    segments_out = []
    all_met = True
    primary_seconds = 0.0
    reference_seconds = 0.0
    planning_seconds = 0.0
    singularity_classifications = []
    for segment_index, segment in enumerate(data["segments"]):
        system = _letters_to_system(
            segment["letters"], data["dimension"],
            f"MadStree-{data['masterDigest'][:12]}-{segment_index}",
        )
        start = gaussian_rational(segment["start"]).to_acb()
        points = [gaussian_rational(point).to_acb() for point in segment["points"]]
        user_indices = list(segment["userIndices"])
        singular_positions = [
            index for index, point in enumerate(points)
            if any((point - pole).contains(0) for pole in system.singularities)
        ]
        terminal_singular_position = (
            len(points) - 1
            if singular_positions and singular_positions[-1] == len(points) - 1
            else None
        )
        intermediate_singular_positions = [
            position for position in singular_positions
            if position != terminal_singular_position
        ]
        if intermediate_singular_positions and (
            not data["pathPlanning"] or data["singularityMode"] != "singularity_jump"
        ):
            return {
                "status": "intermediateSingularPointRequiresPlannedJump", "schema": data["schema"],
                "segmentIndex": segment_index,
                "message": (
                    "intermediate singular user points require FlintNDE path planning "
                    "with singularity-jump mode"
                ),
                "singularUserIndices": [
                    user_indices[index] for index in intermediate_singular_positions
                ],
            }
        singular_position = terminal_singular_position
        regular_points = points if singular_position is None else points[:singular_position]
        regular_indices = user_indices if singular_position is None else user_indices[:singular_position]
        clock = time.perf_counter()
        point_records = []
        plan = None
        result = None
        terminal_plan = None
        terminal_result = None
        assignments = []
        met = True
        if regular_points:
            try:
                plan = (
                    plan_transport_path(
                        system, start, regular_points, singularity_mode=data["singularityMode"],
                        message_language=data["messageLanguage"],
                    )
                    if data["pathPlanning"]
                    else direct_user_point_path(
                        system, start, regular_points, message_language=data["messageLanguage"]
                    )
                )
            except SingularPathError as error:
                return {
                    "status": "singularPathRefused", "schema": data["schema"],
                    "segmentIndex": segment_index, "message": str(error),
                    "singularPathPairs": [[_point_record(a, digits), _point_record(b, digits)] for a, b in error.singular_path_pairs],
                }
            planning_seconds += time.perf_counter() - clock
            assignments = _point_assignments(plan, regular_points, regular_indices)
            result = transport_planned_path_refined(
                system, primary_vector, plan,
                reference_initial_vector=reference_vector,
                primary_order=data["primaryOrder"], reference_order=data["referenceOrder"],
                radius_fraction=0.60, target_relative_error=data["targetRelativeError"],
                certification_mode=data["certificationMode"],
            )
            primary_vector = result["primary_snapshots"][-1]
            reference_vector = result["reference_snapshots"][-1]
            primary_seconds += float(result["primary_seconds"])
            reference_seconds += float(result["reference_seconds"])
            met = result["target_relative_error_met"]
            all_met = all_met and (met is None or bool(met))
            point_records = _point_values(result, assignments, digits, data["dimension"])
            for singular_item in result.get("singular_target_results", []):
                local_index = int(singular_item["user_point_index"])
                singularity_classifications.append({
                    "segmentIndex": segment_index,
                    "userIndex": regular_indices[local_index],
                    "coordinate": _point_record(regular_points[local_index], digits),
                    **singular_item["report"],
                })
        if singular_position is not None:
            if segment_index != len(data["segments"]) - 1:
                return {
                    "status": "singularTurnContinuationUnsupported",
                    "schema": data["schema"],
                    "segmentIndex": segment_index,
                    "message": (
                        "a singular point at an affine-segment turn does not fix a unique "
                        "continuation branch"
                    ),
                    "singularUserIndex": user_indices[singular_position],
                }
            target = points[singular_position]
            target_exact = gaussian_rational(segment["points"][singular_position])
            rational_system = _letters_to_rational_system(
                segment["letters"], data["dimension"],
                f"MadStree-{data['masterDigest'][:12]}-{segment_index}-target",
            )
            start_for_target = (
                gaussian_rational(segment["points"][singular_position - 1])
                if regular_points else gaussian_rational(segment["start"])
            )
            match_clock = time.perf_counter()
            match_plan = plan_singular_target_match(
                rational_system, start_for_target, target_exact,
            )
            if match_plan.inserted:
                if not data["pathPlanning"]:
                    return {
                        "status": "singularTargetMatchPlanningRequired",
                        "schema": data["schema"],
                        "segmentIndex": segment_index,
                        "message": (
                            "the terminal singular target needs a hidden ordinary match point "
                            "inside its convergence disk; enable FlintNDE path planning"
                        ),
                        "singularUserIndex": user_indices[singular_position],
                        "terminalSingularMatch": match_plan.report,
                    }
                terminal_plan = plan_transport_path(
                    system,
                    start_for_target.to_acb(),
                    [match_plan.match_point],
                    singularity_mode=data["singularityMode"],
                    message_language=data["messageLanguage"],
                )
                planning_seconds += time.perf_counter() - match_clock
                terminal_result = transport_planned_path_refined(
                    system,
                    primary_vector,
                    terminal_plan,
                    reference_initial_vector=reference_vector,
                    primary_order=data["primaryOrder"],
                    reference_order=data["referenceOrder"],
                    radius_fraction=0.60,
                    target_relative_error=data["targetRelativeError"],
                    certification_mode=data["certificationMode"],
                )
                primary_vector = terminal_result["primary_snapshots"][-1]
                reference_vector = terminal_result["reference_snapshots"][-1]
                primary_seconds += float(terminal_result["primary_seconds"])
                reference_seconds += float(terminal_result["reference_seconds"])
                terminal_transport_met = terminal_result["target_relative_error_met"]
                all_met = all_met and (
                    terminal_transport_met is None or bool(terminal_transport_met)
                )
            else:
                planning_seconds += time.perf_counter() - match_clock
            primary_singular = evaluate_singular_target(
                rational_system, match_plan.match_point, primary_vector, target_exact,
                order=data["primaryOrder"], target_relative_error=data["targetRelativeError"],
            )
            reference_singular = evaluate_singular_target(
                rational_system, match_plan.match_point, reference_vector, target_exact,
                order=data["referenceOrder"], target_relative_error=data["targetRelativeError"],
            )
            singular_difference, singular_error = _singular_result_difference(
                primary_singular, reference_singular
            )
            if singular_error is not None:
                return {
                    "status": "singularTargetReferenceMismatch",
                    "schema": data["schema"],
                    "segmentIndex": segment_index,
                    "message": singular_error,
                    "singularUserIndex": user_indices[singular_position],
                    "primaryClassification": primary_singular.classification,
                    "referenceClassification": reference_singular.classification,
                    "terminalSingularMatch": match_plan.report,
                }
            singular_met = bool(singular_difference < arb(data["targetRelativeError"]))
            all_met = all_met and singular_met
            met = singular_met if met is None else bool(met) and singular_met
            point_records.append({
                "userIndex": user_indices[singular_position],
                "localPointIndex": singular_position,
                "source": "singular_target",
                "evaluationAlgorithm": None,
                "values": [_value_record(value, digits) for value in primary_singular.values],
            })
            singularity_classifications.append({
                "segmentIndex": segment_index,
                "userIndex": user_indices[singular_position],
                "coordinate": _point_record(target, digits),
                **primary_singular.report,
                "referenceClassification": reference_singular.classification,
                "primaryReferenceRelativeDifference": singular_difference.str(digits),
                "targetRelativeErrorMet": singular_met,
                "terminalSingularMatch": match_plan.report,
            })
        actual_nodes = [] if plan is None else list(plan.nodes)
        if singular_position is not None and terminal_plan is not None:
            terminal_nodes = list(terminal_plan.nodes)
            if actual_nodes and terminal_nodes and abs(actual_nodes[-1] - terminal_nodes[0]).contains(0):
                actual_nodes.extend(terminal_nodes[1:])
            else:
                actual_nodes.extend(terminal_nodes)
        plan_report = {} if plan is None else dict(plan.report)
        if singular_position is not None:
            plan_report["terminal_singular_match"] = match_plan.report
        difference_candidates = []
        if result is not None:
            difference_candidates.append(result["relative_difference_inf"])
        if singular_position is not None:
            difference_candidates.append(singular_difference)
        if terminal_result is not None:
            difference_candidates.append(terminal_result["relative_difference_inf"])
        segment_difference = max(difference_candidates) if difference_candidates else None
        point_records.sort(key=lambda item: item["userIndex"])
        segments_out.append({
            "segmentIndex": segment_index,
            "pathPlanning": data["pathPlanning"],
            "actualNodes": [_point_record(node, digits) for node in actual_nodes],
            "nodeCount": len(actual_nodes),
            "coveredSampleCount": 0 if plan is None else len(plan.sample_assignments),
            "planReport": plan_report,
            "pointAssignments": assignments,
            "pointValues": point_records,
            "relativeDifferenceInf": None if segment_difference is None else segment_difference.str(digits),
            "relativeDifferenceMidpoint": None if segment_difference is None else float(segment_difference.mid()),
            "targetRelativeErrorMet": True if met is None else bool(met),
            "certificationMode": None if result is None else result["certification_mode"],
            "continuationReady": singular_position is None,
            "endpointValues": (
                [_value_record(value, digits) for value in primary_singular.values]
                if singular_position is not None
                else [
                    _acb_record(primary_vector[row, 0], digits)
                    for row in range(data["dimension"])
                ]
            ),
        })
    message = (
        "FlintNDE planned and evaluated every MadStree affine segment in one process."
        if data["pathPlanning"] else
        "FlintNDE evaluated every MadStree affine segment using the supplied user points as nodes."
    ) if data["messageLanguage"] == "EN" else (
        "FlintNDE 已在同一进程内规划并计算全部 MadStree 复仿射段。"
        if data["pathPlanning"] else
        "FlintNDE 已把用户给定点严格作为节点，在同一进程内计算全部 MadStree 复仿射段。"
    )
    return {
        "status": "success", "schema": data["schema"],
        "workerPid": os.getpid(),
        "executionAction": "plan_and_execute" if data["pathPlanning"] else "execute_user_nodes",
        "masterDigest": data["masterDigest"], "dimension": data["dimension"],
        "backendPackagePath": str(backend_path), "messageLanguage": data["messageLanguage"],
        "message": message, "pathPlanning": data["pathPlanning"],
        "segmentCount": len(segments_out), "segments": segments_out,
        "singularityClassifications": singularity_classifications,
        "boundary": boundary_report, "planningSeconds": planning_seconds,
        "primarySeconds": primary_seconds, "referenceSeconds": reference_seconds,
        "targetRelativeErrorMet": all_met,
        "finalValues": [
            _acb_record(primary_vector[row, 0], digits)
            for row in range(data["dimension"])
        ],
    }


def _run_indexed_task(item: tuple[int, dict[str, Any]]) -> tuple[int, dict[str, Any]]:
    """在独立进程执行一个固定 ep 请求，并保留输入序号。"""

    index, data = item
    return index, _run(data)


def _validate_batch_request(value: Any) -> dict[str, Any]:
    """验证 MadStree 固定 ep 批量请求和缺省 12 并行合同。"""

    data = _exact_keys(
        value, {"schema", "parallelTaskCount", "tasks", "messageLanguage"},
        "ep batch request",
    )
    if data["schema"] != EVALUATE_BATCH_SCHEMA:
        raise ValueError("unsupported MadStree-FlintNDE ep batch schema")
    parallel_count = _integer(data["parallelTaskCount"], "parallelTaskCount")
    if data["messageLanguage"] not in {"EN", "CN"}:
        raise ValueError('messageLanguage must be "EN" or "CN"')
    if not isinstance(data["tasks"], list) or not data["tasks"]:
        raise ValueError("ep batch tasks must be a nonempty list")
    tasks = []
    for index, task in enumerate(data["tasks"]):
        item = _exact_keys(task, {"ep", "request"}, f"tasks[{index}]")
        tasks.append({"ep": item["ep"], "request": _validate_request(item["request"])})
    return {
        "schema": data["schema"],
        "parallelTaskCount": parallel_count,
        "tasks": tasks,
        "messageLanguage": data["messageLanguage"],
    }


def _run_batch(data: dict[str, Any]) -> dict[str, Any]:
    """用有界进程池运行不同 ep；每项隔离 FLINT 上下文并按输入顺序返回。"""

    tasks = data["tasks"]
    requested = data["parallelTaskCount"]
    effective = min(requested, len(tasks))
    payloads = [(index, task["request"]) for index, task in enumerate(tasks)]
    with get_context("spawn").Pool(
        processes=effective, maxtasksperchild=1
    ) as pool:
        completed = pool.map(_run_indexed_task, payloads, chunksize=1)
    results = [
        {"ep": tasks[index]["ep"], "result": result}
        for index, result in completed
    ]
    success = all(item["result"].get("status") == "success" for item in results)
    language = data["messageLanguage"]
    message = (
        f"MadStree/FlintNDE ep task pool used {effective} of requested {requested} "
        "workers (default 12); queued tasks started automatically as workers finished."
        if language == "EN" else
        f"MadStree/FlintNDE 不同 ep 任务池：请求并行数 {requested}，实际并行数 {effective}"
        "（缺省 12）；任务完成后已自动续交队列。"
    )
    return {
        "status": "success" if success else "failed",
        "schema": data["schema"],
        "message": message,
        "messageLanguage": language,
        "parallelTaskCountRequested": requested,
        "parallelTaskCountEffective": effective,
        "taskCount": len(results),
        "results": results,
    }


def _validate_series_control_request(value: Any) -> dict[str, Any]:
    """验证自适应正规化控制请求；每个阶段只接受自己的严格字段集合。"""

    if not isinstance(value, dict) or value.get("schema") != SERIES_CONTROL_SCHEMA:
        raise ValueError("unsupported MadStree-FlintNDE ep series control schema")
    action = value.get("action")
    common = {"schema", "action", "backendPackagePath", "maximumPower", "goalDigits"}
    production_base = common | {
            "leadingPower", "sampleSpacing", "validationSampleCount", "validationScale",
            "maximumSamples", "extraWorkingPrecision", "productionRound",
            "fitExtraOrder", "fitOrderIncrement", "fitMaximumRounds",
        }
    expected = {
        "production_plan": production_base,
        "fit": common | {
            "leadingPower", "workingPrecisionDigits", "points", "values",
            "validationPoints", "validationValues", "validationTolerance",
        },
    }.get(action)
    if expected is None:
        raise ValueError(f"unsupported ep series control action: {action}")
    if action == "production_plan":
        optional = {
            "samplePoints", "sampleAngleRange", "validationPoints",
            "initialInternalMaximumPower"
        }
        actual = set(value)
        if not production_base <= actual or not actual <= production_base | optional:
            raise ValueError(
                "ep series production_plan request has unexpected or missing fields"
            )
        data = dict(value)
    else:
        data = _exact_keys(value, expected, f"ep series {action} request")
    _string(data["backendPackagePath"], "backendPackagePath")
    _integer(data["maximumPower"], "maximumPower", -10**9)
    _integer(data["goalDigits"], "goalDigits")
    return data


def _fraction_text(value: Fraction) -> str:
    """把自动网格点写成 Wolfram 和 FLINT 都能精确读取的有理数字符串。"""

    return f"{value.numerator}/{value.denominator}"


def _wolfram_sample_text(value: Any, digits: int) -> str:
    """把实点或 Acb 复点写成 MadStree 可精确拉回的 Gaussian rational。"""

    from flint import acb, fmpq

    if isinstance(value, fmpq):
        return str(value)
    point = acb(value)
    real = _fraction_text(Fraction(point.real.mid().str(digits, radius=False)))
    imag = _fraction_text(Fraction(point.imag.mid().str(digits, radius=False)))
    return f"({real})+({imag})*I"


def _sample_text_acb(value: Any) -> Any:
    """读取 adapter 自己生成的 Wolfram 复数文本或原有实标量文本。"""

    from flint import acb

    if isinstance(value, str):
        match = re.fullmatch(r"\(([^()]*)\)\+\(([^()]*)\)\*I", value)
        if match is not None:
            return acb(
                match.group(1).replace("*^", "e"),
                match.group(2).replace("*^", "e"),
            )
    return acb(value)


def _acb_from_ball_record(record: Any, label: str) -> Any:
    """从 MadStree 求值结果恢复完整 Acb 球，不使用仅供显示的中点。"""

    from flint import acb
    from flintnde.core import arb_ball_from_json

    item = _exact_keys(
        record, {"real", "imag", "real_ball", "imag_ball"}, label
    )
    return acb(
        arb_ball_from_json(item["real_ball"], f"{label}.real_ball"),
        arb_ball_from_json(item["imag_ball"], f"{label}.imag_ball"),
    )


def _series_vectors(records: Any, label: str) -> tuple[Any, ...]:
    """恢复一组同维 Acb 列矢量。"""

    from flintnde.core import column_vector

    if not isinstance(records, list) or not records:
        raise ValueError(f"{label} must be a nonempty list")
    vectors = []
    for sample_index, vector in enumerate(records):
        if not isinstance(vector, list) or not vector:
            raise ValueError(f"{label}[{sample_index}] must be a nonempty vector")
        vectors.append(column_vector([
            _acb_from_ball_record(item, f"{label}[{sample_index}][{row}]")
            for row, item in enumerate(vector)
        ]))
    dimension = vectors[0].nrows()
    if any(vector.nrows() != dimension for vector in vectors):
        raise ValueError(f"{label} vectors have inconsistent dimensions")
    return tuple(vectors)


def _run_series_control(data: dict[str, Any]) -> dict[str, Any]:
    """按 MadStree 符号证书规划 ep 网格，或重构并验证 Laurent 系数。"""

    _resolve_backend_path(data)
    from flint import acb
    from flintnde import configure_working_precision, fit_sampled_series
    from flintnde.regularization import (
        AUTOMATIC,
        _resolve_plan,
    )

    action = data["action"]
    maximum_power = data["maximumPower"]
    goal_digits = data["goalDigits"]
    if action == "production_plan":
        leading = _integer(data["leadingPower"], "leadingPower", -10**9)
        production_round = _integer(data["productionRound"], "productionRound")
        fit_extra_order = _integer(data["fitExtraOrder"], "fitExtraOrder", 0)
        fit_order_increment = _integer(
            data["fitOrderIncrement"], "fitOrderIncrement", 1
        )
        fit_maximum_rounds = _integer(
            data["fitMaximumRounds"], "fitMaximumRounds", 1
        )
        if production_round > fit_maximum_rounds:
            raise ValueError("productionRound exceeds fitMaximumRounds")
        pole_depth = max(0, -leading)
        reconstruction_depth = maximum_power - leading
        automatic_initial_count = max(
            math.ceil(2.5 * reconstruction_depth + pole_depth),
            reconstruction_depth + 1 + fit_extra_order,
        )
        initial_internal_maximum = data.get("initialInternalMaximumPower")
        if initial_internal_maximum is None:
            initial_count = automatic_initial_count
        else:
            initial_internal_maximum = _integer(
                initial_internal_maximum, "initialInternalMaximumPower", maximum_power
            )
            initial_count = initial_internal_maximum - leading + 1
        target_count = initial_count + fit_order_increment * (production_round - 1)
        explicit_points = data.get("samplePoints")
        explicit_validation_points = data.get("validationPoints")
        if explicit_points is not None:
            if not isinstance(explicit_points, list) or not explicit_points:
                raise ValueError("samplePoints must be a nonempty list")
            capacity_count = len(explicit_points)
            if initial_internal_maximum is not None and initial_count > capacity_count:
                raise ValueError(
                    f"samplePoints has {capacity_count} candidates but first-round internal "
                    f"power {initial_internal_maximum} requires {initial_count}"
                )
            initial_count = min(initial_count, capacity_count)
            if initial_count < reconstruction_depth + 1:
                raise ValueError("samplePoints do not cover the requested regulator powers")
            target_count = min(
                capacity_count,
                initial_count + fit_order_increment * (production_round - 1),
            )
        else:
            capacity_count = initial_count + fit_order_increment * (fit_maximum_rounds - 1)
        maximum_samples = _integer(data["maximumSamples"], "maximumSamples")
        if target_count > maximum_samples:
            raise ValueError(
                f"incremental production needs {target_count} samples, "
                f"above maximumSamples={maximum_samples}"
            )
        capacity_count = min(capacity_count, maximum_samples)
        alpha = Decimal(pole_depth) / Decimal(4) + Decimal(goal_digits) / Decimal(
            reconstruction_depth + 1
        )
        plan = _resolve_plan(
            leading_power=leading, maximum_power=maximum_power,
            goal_digits=goal_digits,
            sample_points=(explicit_points if explicit_points is not None else AUTOMATIC),
            sample_angle_range=data.get("sampleAngleRange", AUTOMATIC),
            sample_count=capacity_count,
            base_sample=AUTOMATIC, sample_spacing=data["sampleSpacing"],
            working_precision_digits=AUTOMATIC,
            extra_working_precision=float(data["extraWorkingPrecision"]),
            transport_order=AUTOMATIC, transport_extra_order=AUTOMATIC,
            transport_sample_count=AUTOMATIC, transport_extra_sample_count=AUTOMATIC,
            validation_sample_count=_integer(
                data["validationSampleCount"], "validationSampleCount"
            ),
            validation_points=(
                explicit_validation_points
                if explicit_validation_points is not None
                else AUTOMATIC
            ),
            validation_scale=data["validationScale"],
            maximum_samples=maximum_samples,
            rationalize_sample_points=True,
        )
        return {
            "status": "success", "schema": data["schema"], "action": action,
            "points": [
                _wolfram_sample_text(item, plan.working_precision_digits)
                for item in plan.sample_arguments[:target_count]
            ],
            "validationPoints": [
                _wolfram_sample_text(item, plan.working_precision_digits)
                for item in plan.validation_arguments
            ],
            "sampleCount": target_count,
            "initialSampleCount": initial_count,
            "capacitySampleCount": capacity_count,
            "unusedCandidateCount": capacity_count - target_count,
            "internalMaximumPower": leading + target_count - 1,
            "workingPrecisionDigits": plan.working_precision_digits,
            "primaryOrder": plan.transport_order,
            "referenceOrder": plan.transport_order + plan.transport_extra_order,
            "targetRelativeError": f"1e-{goal_digits}",
            "baseSample": plan.base_sample,
            "alphaEpsilon": plan.alpha_epsilon,
            "productionRound": production_round,
            "sampleSource": plan.source,
        }
    leading = _integer(data["leadingPower"], "leadingPower", -10**9)
    precision = _integer(data["workingPrecisionDigits"], "workingPrecisionDigits")
    configure_working_precision(precision)
    values = _series_vectors(data["values"], "values")
    validation_values = _series_vectors(data["validationValues"], "validationValues")
    internal_maximum = leading + len(values) - 1
    result = fit_sampled_series(
        sample_points=[_sample_text_acb(value) for value in data["points"]],
        sample_values=values,
        maximum_power=internal_maximum, leading_power=leading,
        validation_points=[
            _sample_text_acb(value) for value in data["validationPoints"]
        ],
        validation_values=validation_values,
        validation_tolerance=None, series_parameter="ep",
    )
    tolerance = acb(data["validationTolerance"]).real
    maximum_residual = acb(
        result.diagnostics["maximum_validation_relative_residual"]
    ).real
    coefficients = {
        str(power): [_acb_record(vector[row, 0], precision) for row in range(vector.nrows())]
        for power, vector in zip(result.powers, result.coefficients)
        if power <= maximum_power
    }
    if not maximum_residual < tolerance:
        return {
            "status": "success", "schema": data["schema"], "action": action,
            "fitStatus": "retry",
            "reason": (
                f"series validation residual {maximum_residual.str(20)} is not "
                f"below {tolerance.str(20)}"
            ),
            "internalMaximumPower": internal_maximum,
            "leadingPower": result.leading_power,
            "maximumPower": result.maximum_power,
            "coefficients": coefficients,
            "diagnostics": result.diagnostics,
            "effectiveParameters": result.effective_parameters,
        }
    return {
        "status": "success", "schema": data["schema"], "action": action,
        "fitStatus": "accepted",
        "leadingPower": result.leading_power,
        "maximumPower": result.maximum_power,
        "internalMaximumPower": internal_maximum,
        "coefficients": coefficients,
        "diagnostics": result.diagnostics,
        "effectiveParameters": result.effective_parameters,
    }


def main() -> int:
    """命令行入口；失败也以 UTF-8 JSON 返回完整原因。"""

    if len(sys.argv) == 1:
        print(_citation_text())
        print("usage: flintnde_transport.py INPUT.json OUTPUT.json")
        return 0
    if len(sys.argv) != 3:
        raise SystemExit("usage: flintnde_transport.py INPUT.json OUTPUT.json")
    # 文件协议把 MadStree 引用写入响应 metadata，且不让依赖包的 import 提示进入日志。
    os.environ["FLINTNDE_SUPPRESS_CITATION_NOTICE"] = "1"
    input_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])
    try:
        raw = json.loads(input_path.read_text(encoding="utf-8"))
        if raw.get("schema") == EVALUATE_BATCH_SCHEMA:
            result = _run_batch(_validate_batch_request(raw))
        elif raw.get("schema") == SERIES_CONTROL_SCHEMA:
            result = _run_series_control(_validate_series_control_request(raw))
        else:
            result = _run(_validate_request(raw))
        exit_code = 0
    except Exception as error:  # noqa: BLE001 - CLI 边界必须序列化所有错误
        result = {
            "status": "failed", "errorType": type(error).__name__,
            "error": str(error), "traceback": traceback.format_exc(),
        }
        exit_code = 1
    result = _with_citation_metadata(result)
    output_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
