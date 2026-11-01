"""正规化参数的幂级数解重构。

本模块在多个固定、非零 regulator 样本上复用 ``transport_path_refined``，随后用
FLINT Acb 方阵插值重构最终解矢量的连续 Laurent 系数。自动样本数、样本尺度和
工作精度采用 AMFlow 2.0 的经验公式；最低幂必须由调用方根据符号边界与 DE 先行认证，
独立验证和失败门禁由本程序包补充。模块不对奇点局部解做多项式回归。

实现思路：先核对调用方提供的最低整数幂，再生成带高阶截断缓冲的生产网格；每个
样本完整调用同一 NDE 主线，最后在未参与插值的更小尺度上验证内部级数。
验证失败时只追加新生产点并提高内部拟合阶数，已有生产值和验证值全程复用。
"""

from __future__ import annotations

import math
import pickle
import warnings
from dataclasses import dataclass, replace
from decimal import Decimal, localcontext
from fractions import Fraction
from multiprocessing import get_context
from pathlib import Path
from typing import Any, Callable

from flint import acb, acb_mat, arb, ctx, fmpq

from .boundary import FrobeniusBoundary, frobenius_boundary
from .core import (
    DEFAULT_WORKING_PRECISION_DIGITS,
    arb_ball_from_json,
    arb_ball_to_json,
    column_vector,
    configure_working_precision,
    exact_rational,
    vector_norm_inf,
)
from .routing import AdaptivePath, adaptive_path_from_json, adaptive_path_to_json
from .singularities import RationalMatrixSystem, rational_function
from .systems import AnalyticMatrixSystem
from .transport import transport_path_refined
from .parallel import DEFAULT_PARALLEL_TASK_COUNT


AUTOMATIC = "automatic"
MatrixSystem = AnalyticMatrixSystem | RationalMatrixSystem


class SeriesValidationError(ValueError):
    """表示 NDE refinement 或独立 regulator 样本未达到目标精度。"""


@dataclass(frozen=True)
class SeriesReconstructionResult:
    """保存返回系数、实际样本、有效参数和全部验证诊断。

    ``coefficients[index]`` 对应 ``powers[index]``，每个系数都是与 NDE 解同维的
    Acb 列向量。``sample_values`` 是生产样本处的完整 NDE 终点值。
    """

    series_parameter: str
    powers: tuple[int, ...]
    coefficients: tuple[acb_mat, ...]
    leading_power: int
    maximum_power: int
    internal_maximum_power: int
    sample_points: tuple[acb, ...]
    sample_values: tuple[acb_mat, ...]
    validation_points: tuple[acb, ...]
    validation_values: tuple[acb_mat, ...]
    effective_parameters: dict[str, Any]
    diagnostics: dict[str, Any]
    output_file: str | None = None
    regularization_points: tuple[acb, ...] | None = None
    validation_regularization_points: tuple[acb, ...] | None = None

    def coefficient(self, power: int) -> acb_mat:
        """返回指定幂的系数列向量；超出用户请求范围时拒绝。"""

        if power < self.leading_power or power > self.maximum_power:
            raise KeyError(
                f"power {power} is outside [{self.leading_power}, {self.maximum_power}]"
            )
        return acb_mat(self.coefficients[power - self.leading_power])

    def evaluate(self, point: Any) -> acb_mat:
        """在非零 regulator 点计算已经截断到 ``maximum_power`` 的返回级数。"""

        regulator = point if isinstance(point, acb) else acb(point)
        if regulator.is_zero() or abs(regulator).contains(0):
            raise ValueError("series evaluation requires a definitely nonzero regulator")
        return _evaluate_coefficients(self.coefficients, self.leading_power, regulator)

    def to_summary(self) -> dict[str, Any]:
        """返回可直接写入 JSON 的稳定摘要。"""

        summary = {
            "schema": "flintnde_series_reconstruction_v1",
            "series_parameter": self.series_parameter,
            "leading_power": self.leading_power,
            "maximum_power": self.maximum_power,
            "internal_maximum_power": self.internal_maximum_power,
            "sample_points": [point.str(50) for point in self.sample_points],
            "validation_points": [point.str(50) for point in self.validation_points],
            "coefficients": {
                str(power): _vector_to_strings(coefficient)
                for power, coefficient in zip(self.powers, self.coefficients)
            },
            "sample_values": [
                _vector_to_strings(value) for value in self.sample_values
            ],
            "validation_values": [
                _vector_to_strings(value) for value in self.validation_values
            ],
            "effective_parameters": self.effective_parameters,
            "diagnostics": self.diagnostics,
            "output_file": self.output_file,
        }
        if self.regularization_points is not None:
            summary["regularization_points"] = [
                point.str(50) for point in self.regularization_points
            ]
        if self.validation_regularization_points is not None:
            summary["validation_regularization_points"] = [
                point.str(50) for point in self.validation_regularization_points
            ]
        return summary

@dataclass(frozen=True)
class _ReconstructionPlan:
    """保存内部生产/验证点以及已经解析的 NDE 数值设置。"""

    sample_arguments: tuple[Any, ...]
    sample_points: tuple[acb, ...]
    validation_arguments: tuple[Any, ...]
    validation_points: tuple[acb, ...]
    sample_count: int
    base_sample: str
    alpha_epsilon: float
    base_precision_digits: int
    working_precision_digits: int
    transport_order: int
    transport_extra_order: int
    transport_sample_count: int | None
    transport_extra_sample_count: int | None
    source: str
    validation_source: str


def _vector_to_strings(vector: acb_mat) -> list[str]:
    """把 Acb 列向量转换为不降到 binary64 的字符串列表。"""

    return [vector[row, 0].str(50) for row in range(vector.nrows())]


def _require_integer(value: Any, name: str, *, minimum: int) -> int:
    """验证非布尔整数及其下界。"""

    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _require_automatic_or_integer(value: Any, name: str, *, minimum: int) -> int | None:
    """把字面量 ``automatic`` 解析为 ``None``，否则验证整数。"""

    if value == AUTOMATIC:
        return None
    return _require_integer(value, name, minimum=minimum)


def _fraction_from_decimal_control(value: Any, name: str) -> Fraction:
    """把采样比例控制量转为精确十进制分数。"""

    try:
        result = Fraction(str(value))
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError(f"{name} must be a finite real number") from error
    return result


def _sample_pair(value: Any, *, rationalize: bool, name: str) -> tuple[Any, acb]:
    """同时保留工厂所需 exact 参数和插值所需 Acb 点。

    自动或显式有理实数以 ``fmpq`` 传给 ``DEmatrix(ep)``；Acb 点保持数值输入。
    Python float/complex 会丢失用户声称的高精度，因此显式拒绝并要求字符串或 Acb。
    """

    if isinstance(value, bool) or isinstance(value, (float, complex)):
        raise TypeError(f"{name} must use an exact string/fmpq or an acb value, not float/complex")
    if isinstance(value, acb):
        point = acb(value)
        if point.is_zero() or abs(point).contains(0):
            raise ValueError(f"{name} must be definitely nonzero")
        return point, point
    if rationalize:
        try:
            exact = exact_rational(value)
        except ValueError as error:
            raise TypeError(f"{name} must be an exact real scalar or acb") from error
        if exact == 0:
            raise ValueError(f"{name} must be nonzero")
        return exact, acb(exact)
    try:
        point = acb(str(value))
    except (TypeError, ValueError) as error:
        raise TypeError(f"{name} cannot be converted to acb") from error
    if point.is_zero() or abs(point).contains(0):
        raise ValueError(f"{name} must be definitely nonzero")
    return point, point


def _power_of_ten_fraction(exponent: Decimal, precision: int) -> Fraction:
    """在受控十进制精度下计算 ``10**exponent`` 并转为 exact rational。"""

    with localcontext() as decimal_context:
        decimal_context.prec = max(50, precision)
        value = Decimal(10) ** exponent
    return Fraction(str(value))


def _automatic_grid(
    base: Fraction,
    count: int,
    spacing: Fraction,
    *,
    rationalize: bool,
) -> tuple[tuple[Any, ...], tuple[acb, ...]]:
    """生成 AMFlow 式网格；按选项以 fmpq 或 Acb 传入用户工厂。"""

    arguments: list[fmpq] = []
    points: list[acb] = []
    for index in range(1, count + 1):
        value = base * (Fraction(1) + index * spacing)
        if rationalize:
            argument = fmpq(value.numerator, value.denominator)
        else:
            argument = acb(value.numerator) / acb(value.denominator)
        arguments.append(argument)
        points.append(acb(argument))
    return tuple(arguments), tuple(points)


def _angle_range_texts(value: Any) -> tuple[str, str] | None:
    """解析自动采样的复角域；普通浮点只作为区域边界控制量。"""

    if value == AUTOMATIC:
        return None
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise TypeError('sample_angle_range must be "automatic" or a two-item sequence')
    if any(isinstance(item, bool) or isinstance(item, complex) for item in value):
        raise TypeError("sample_angle_range endpoints must be finite real angles in radians")
    texts = tuple(str(item) for item in value)
    try:
        bounds = tuple(arb(text) for text in texts)
    except (TypeError, ValueError) as error:
        raise TypeError(
            "sample_angle_range endpoints must be finite real angles in radians"
        ) from error
    if not all(bound.is_finite() for bound in bounds):
        raise ValueError("sample_angle_range endpoints must be finite")
    if not bounds[0] < bounds[1]:
        raise ValueError("sample_angle_range must define a nonempty open angle interval")
    return texts


def _rotate_automatic_grid(
    radii: tuple[acb, ...],
    angle_range: tuple[str, str],
) -> tuple[tuple[acb, ...], tuple[acb, ...]]:
    """在开角域内部均匀选至多三条射线，并循环分配自动模长。"""

    lower = arb(angle_range[0])
    upper = arb(angle_range[1])
    width = upper - lower
    angle_count = min(3, len(radii))
    angles = tuple(
        lower + width * arb(index) / arb(angle_count + 1)
        for index in range(1, angle_count + 1)
    )
    points: list[acb] = []
    for index, radius in enumerate(radii, start=1):
        theta = angles[(index - 1) % angle_count]
        point = acb(radius) * acb(0, theta).exp()
        points.append(point)
    resolved = tuple(points)
    _require_distinct_points(resolved, "automatic angle-range sample_points")
    return resolved, resolved


def _point_alpha(points: tuple[acb, ...]) -> float:
    """由最大生产样本模估计 ``-log10(r)``，供自定义网格精度规划。"""

    log_magnitudes = [float(abs(point).log().mid()) / math.log(10) for point in points]
    return max(0.0, -max(log_magnitudes))


def _scaled_validation_points(
    sample_arguments: tuple[Any, ...],
    count: int,
    scale: Fraction,
) -> tuple[tuple[Any, ...], tuple[acb, ...]]:
    """从生产网格前部生成未参与插值的更小尺度验证点。"""

    arguments: list[Any] = []
    points: list[acb] = []
    exact_scale = fmpq(scale.numerator, scale.denominator)
    for argument in sample_arguments[:count]:
        if isinstance(argument, fmpq):
            value = argument * exact_scale
        else:
            value = acb(argument) * acb(exact_scale)
        if any(abs(acb(value) - existing).contains(0) for existing in points):
            raise ValueError("validation grid contains duplicate points")
        arguments.append(value)
        points.append(acb(value))
    return tuple(arguments), tuple(points)


def _resolve_plan(
    *,
    leading_power: int,
    maximum_power: int,
    goal_digits: int,
    sample_points: Any,
    sample_angle_range: Any = AUTOMATIC,
    sample_count: Any,
    base_sample: Any,
    sample_spacing: Any,
    working_precision_digits: Any,
    extra_working_precision: float,
    transport_order: Any,
    transport_extra_order: Any,
    transport_sample_count: Any,
    transport_extra_sample_count: Any,
    validation_sample_count: int,
    validation_points: Any,
    validation_scale: Any,
    maximum_samples: int,
    rationalize_sample_points: bool,
    fit_extra_order: int = 0,
) -> _ReconstructionPlan:
    """解析生产网格、AMFlow 式精度以及底层 NDE 阶数。"""

    pole_depth = max(0, -leading_power)
    reconstruction_depth = maximum_power - leading_power
    if reconstruction_depth < 0:
        raise ValueError("maximum_power must not be below the detected leading_power")
    automatic_count = max(
        math.ceil(2.5 * reconstruction_depth + pole_depth),
        reconstruction_depth + 1 + fit_extra_order,
    )
    resolved_count = _require_automatic_or_integer(
        sample_count, "sample_count", minimum=reconstruction_depth + 1
    )

    angle_range = _angle_range_texts(sample_angle_range)
    if sample_points != AUTOMATIC and angle_range is not None:
        raise ValueError("sample_angle_range cannot be combined with explicit sample_points")

    if sample_points == AUTOMATIC:
        count = automatic_count if resolved_count is None else resolved_count
        if count > maximum_samples:
            raise ValueError(
                f"automatic production needs {count} samples, above maximum_samples={maximum_samples}"
            )
        alpha_decimal = Decimal(pole_depth) / Decimal(4) + Decimal(goal_digits) / Decimal(
            reconstruction_depth + 1
        )
        if base_sample == AUTOMATIC:
            base_fraction = _power_of_ten_fraction(
                -alpha_decimal, goal_digits + 30
            )
        else:
            base_fraction = _fraction_from_decimal_control(base_sample, "base_sample")
        if base_fraction <= 0:
            raise ValueError("base_sample must be positive")
        spacing = _fraction_from_decimal_control(sample_spacing, "sample_spacing")
        if spacing <= 0:
            raise ValueError("sample_spacing must be positive")
        arguments, points = _automatic_grid(
            base_fraction,
            count,
            spacing,
            rationalize=rationalize_sample_points,
        )
        alpha = float(alpha_decimal) if base_sample == AUTOMATIC else _point_alpha(points)
        source = "automatic"
        base_label = str(base_fraction)
    else:
        if base_sample != AUTOMATIC:
            raise ValueError("base_sample cannot be combined with explicit sample_points")
        if not isinstance(sample_points, (list, tuple)) or not sample_points:
            raise TypeError('sample_points must be "automatic" or a nonempty sequence')
        pairs = tuple(
            _sample_pair(
                value,
                rationalize=rationalize_sample_points,
                name=f"sample_points[{index}]",
            )
            for index, value in enumerate(sample_points)
        )
        arguments = tuple(pair[0] for pair in pairs)
        points = tuple(pair[1] for pair in pairs)
        count = len(points)
        if resolved_count is not None and resolved_count != count:
            raise ValueError("sample_count must equal the explicit sample_points length")
        if count < reconstruction_depth + 1:
            raise ValueError("explicit sample_points do not cover the requested powers")
        if count > maximum_samples:
            raise ValueError("explicit sample_points exceed maximum_samples")
        _require_distinct_points(points, "sample_points")
        alpha = _point_alpha(points)
        source = "user"
        base_label = "user-supplied"

    base_precision = max(math.ceil((count + pole_depth) * alpha), 30)
    resolved_precision = _require_automatic_or_integer(
        working_precision_digits, "working_precision_digits", minimum=1
    )
    precision = resolved_precision or max(
        DEFAULT_WORKING_PRECISION_DIGITS,
        math.ceil(2 * (1 + extra_working_precision) * base_precision),
    )
    if sample_points == AUTOMATIC and angle_range is not None:
        old_precision = ctx.prec
        try:
            ctx.prec = max(old_precision, math.ceil(precision * math.log2(10)) + 32)
            arguments, points = _rotate_automatic_grid(points, angle_range)
        finally:
            ctx.prec = old_precision
        source = "automatic-angle-range"
    resolved_order = _require_automatic_or_integer(
        transport_order, "transport_order", minimum=1
    )
    primary_order = resolved_order or 4 * base_precision
    resolved_extra_order = _require_automatic_or_integer(
        transport_extra_order, "transport_extra_order", minimum=1
    )
    extra_order = resolved_extra_order or math.ceil(max(50, base_precision / 5))
    primary_samples = _require_automatic_or_integer(
        transport_sample_count, "transport_sample_count", minimum=primary_order + 1
    )
    reference_order = primary_order + extra_order
    reference_samples = _require_automatic_or_integer(
        transport_extra_sample_count,
        "transport_extra_sample_count",
        minimum=reference_order + 1,
    )
    if validation_points == AUTOMATIC:
        scale = _fraction_from_decimal_control(validation_scale, "validation_scale")
        if not 0 < scale < 1:
            raise ValueError("validation_scale must lie strictly between zero and one")
        validation_arguments, resolved_validation_points = _scaled_validation_points(
            arguments, validation_sample_count, scale
        )
        validation_source = "scaled-production-grid"
    else:
        if not isinstance(validation_points, (list, tuple)) or not validation_points:
            raise TypeError(
                'validation_points must be "automatic" or a nonempty sequence'
            )
        validation_pairs = tuple(
            _sample_pair(
                value,
                rationalize=rationalize_sample_points,
                name=f"validation_points[{index}]",
            )
            for index, value in enumerate(validation_points)
        )
        validation_arguments = tuple(pair[0] for pair in validation_pairs)
        resolved_validation_points = tuple(pair[1] for pair in validation_pairs)
        _require_distinct_points(resolved_validation_points, "validation_points")
        validation_source = "user"
    production_set = points
    if any(
        abs(validation_point - production_point).contains(0)
        for validation_point in resolved_validation_points
        for production_point in production_set
    ):
        raise ValueError("validation points must be distinct from production points")

    return _ReconstructionPlan(
        arguments,
        points,
        validation_arguments,
        resolved_validation_points,
        count,
        base_label,
        alpha,
        base_precision,
        precision,
        primary_order,
        extra_order,
        primary_samples,
        reference_samples,
        source,
        validation_source,
    )


def _resolve_matrix_system(DEmatrix: Any, sample: Any) -> MatrixSystem:
    """在一个 regulator 样本上解析固定系统或 ``DEmatrix(ep)`` 工厂。"""

    system = DEmatrix(sample) if callable(DEmatrix) else DEmatrix
    if not isinstance(system, (AnalyticMatrixSystem, RationalMatrixSystem)):
        raise TypeError(
            "DEmatrix must be an AnalyticMatrixSystem/RationalMatrixSystem or ep -> system"
        )
    return system


def _resolve_boundary(
    boundary: Any, sample: Any, path: list[acb]
) -> acb_mat | FrobeniusBoundary:
    """按路径起点解析普通列向量或 exact ``{a,b,C}`` 边界。"""

    value = boundary(sample) if callable(boundary) else boundary
    start_classification = getattr(path, "start_classification", None)
    if start_classification is not None and start_classification.kind != "ordinary":
        return frobenius_boundary(value)
    if isinstance(value, acb_mat):
        vector = acb_mat(value)
    elif isinstance(value, (list, tuple)):
        if value and isinstance(value[0], (dict, FrobeniusBoundary)):
            raise TypeError("{a,b,C} boundary requires a regular-singular path start")
        vector = column_vector(list(value))
    else:
        raise TypeError(
            "ordinary boundary must be an acb_mat/list; singular boundary must use {a,b,C}"
        )
    if vector.ncols() != 1:
        raise ValueError("boundary must be a column vector")
    return vector


def _resolve_path(path: Any, sample: Any, system: MatrixSystem) -> list[acb]:
    """解析固定路径或 ``path(ep, DEmatrix_at_ep)`` 工厂。"""

    value = path(sample, system) if callable(path) else path
    if not isinstance(value, list):
        raise TypeError("path must be a list[acb] or (ep, system) -> list[acb]")
    if type(value) is list:
        return [point if isinstance(point, acb) else acb(point) for point in value]
    return value


def _solve_sample(
    *,
    DEmatrix: Any,
    boundary: Any,
    path: Any,
    sample_argument: Any,
    sample_point: acb,
    transport_order: int,
    transport_extra_order: int,
    transport_sample_count: int | None,
    transport_extra_sample_count: int | None,
    radius_fraction: float,
    nde_tolerance: arb,
) -> tuple[acb_mat, dict[str, Any]]:
    """在固定 regulator 上运行基础 NDE refinement，并执行样本级精度门禁。"""

    system = _resolve_matrix_system(DEmatrix, sample_argument)
    resolved_path = _resolve_path(path, sample_argument, system)
    initial_vector = _resolve_boundary(boundary, sample_argument, resolved_path)
    if isinstance(initial_vector, acb_mat) and initial_vector.nrows() != system.dimension:
        raise ValueError("boundary dimension does not match DEmatrix")
    result = transport_path_refined(
        system,
        initial_vector,
        resolved_path,
        primary_order=transport_order,
        reference_order=transport_order + transport_extra_order,
        primary_sample_count=transport_sample_count,
        reference_sample_count=transport_extra_sample_count,
        radius_fraction=radius_fraction,
        target_relative_error=nde_tolerance,
    )
    nde_error = result["relative_difference_inf"]
    if not nde_error < nde_tolerance:
        raise SeriesValidationError(
            f"NDE refinement failed at {sample_point.str(30)}: "
            f"relative difference {nde_error.str(20)} is not below {nde_tolerance.str(20)}"
        )
    endpoint = acb_mat(result["primary_snapshots"][-1])
    boundary_initialization = next(
        (
            report
            for report in result["primary_segments"]
            if report.get("method", "").endswith("boundary_initialization")
        ),
        None,
    )
    return endpoint, {
        "sample_point": sample_point.str(50),
        "endpoint": _vector_to_strings(endpoint),
        "nde_relative_difference": nde_error.str(30),
        "primary_seconds": result["primary_seconds"],
        "reference_seconds": result["reference_seconds"],
        "boundary_initialization": boundary_initialization,
    }


def _solve_sample_job(payload: dict[str, Any]) -> tuple[list[Any], dict[str, Any]]:
    """在独立进程重建 FLINT 标量并执行一个固定 regulator 样本。"""

    ctx.prec = payload["working_precision_bits"]
    DEmatrix, boundary, path = _restore_parallel_inputs(payload["inputs"])
    sample_argument = payload["sample_argument"]
    if isinstance(sample_argument, dict) and sample_argument.get("kind") == "acb-ball":
        sample_argument = _acb_ball_from_record(
            sample_argument["value"], "sample_argument"
        )
    sample_point = payload["sample_point"]
    if isinstance(sample_point, dict) and sample_point.get("kind") == "acb-ball":
        sample_point = _acb_ball_from_record(sample_point["value"], "sample_point")
    else:
        sample_point = acb(sample_point)
    value, report = _solve_sample(
        DEmatrix=DEmatrix,
        boundary=boundary,
        path=path,
        sample_argument=sample_argument,
        sample_point=sample_point,
        transport_order=payload["transport_order"],
        transport_extra_order=payload["transport_extra_order"],
        transport_sample_count=payload["transport_sample_count"],
        transport_extra_sample_count=payload["transport_extra_sample_count"],
        radius_fraction=payload["radius_fraction"],
        nde_tolerance=arb(payload["nde_tolerance"]),
    )
    record_digits = math.ceil(int(ctx.prec) / math.log2(10)) + 10
    values = []
    for row in range(value.nrows()):
        entry = value[row, 0]
        values.append(
            entry.str(50)
            if entry.imag.is_zero()
            else {"kind": "acb-ball", "value": _acb_ball_record(entry, record_digits)}
        )
    return values, report


def _rational_system_record(system: RationalMatrixSystem) -> dict[str, Any]:
    """把固定 exact 有理矩阵转换为不含 FLINT 扩展对象的进程记录。"""

    entries = []
    for row in system.entries:
        output_row = []
        for entry in row:
            numerator, denominator = entry.exact_polynomials()
            output_row.append(
                {
                    "numerator": numerator.records(),
                    "denominator": denominator.records(),
                }
            )
        entries.append(output_row)
    return {
        "kind": "fixed-rational-system",
        "entries": entries,
        "variable_name": system.variable_name,
        "name": system.name,
    }


def _rational_system_from_record(record: dict[str, Any]) -> RationalMatrixSystem:
    """在 worker 内恢复固定 exact 有理矩阵。"""

    return RationalMatrixSystem(
        tuple(
            tuple(
                rational_function(entry["numerator"], entry["denominator"])
                for entry in row
            )
            for row in record["entries"]
        ),
        variable_name=record["variable_name"],
        name=record["name"],
    )


def _acb_ball_record(value: acb, decimal_digits: int) -> dict[str, Any]:
    """把 Acb 球拆成两个可严格恢复的 Arb 球记录。"""

    return {
        "real": arb_ball_to_json(value.real, decimal_digits),
        "imag": arb_ball_to_json(value.imag, decimal_digits),
    }


def _acb_ball_from_record(record: dict[str, Any], field_name: str) -> acb:
    """从进程记录恢复 Acb 球，不把固定边界截断到有限十进制字符串。"""

    return acb(
        arb_ball_from_json(record["real"], f"{field_name}.real"),
        arb_ball_from_json(record["imag"], f"{field_name}.imag"),
    )


def _parallel_input_records(DEmatrix: Any, boundary: Any, path: Any) -> dict[str, Any]:
    """序列化固定 FLINT 输入；顶层工厂保持原对象交给 pickle。"""

    record_digits = math.ceil(int(ctx.prec) / math.log2(10)) + 10
    system_record = (
        _rational_system_record(DEmatrix)
        if isinstance(DEmatrix, RationalMatrixSystem)
        else {"kind": "factory", "value": DEmatrix}
    )
    if isinstance(boundary, acb_mat):
        boundary_record = {
            "kind": "fixed-vector",
            "values": [
                _acb_ball_record(boundary[row, 0], record_digits)
                for row in range(boundary.nrows())
            ],
        }
    elif callable(boundary):
        boundary_record = {"kind": "factory", "value": boundary}
    elif isinstance(boundary, (list, tuple)) and not any(
        isinstance(item, dict) for item in boundary
    ):
        boundary_record = {
            "kind": "fixed-vector",
            "values": [
                _acb_ball_record(acb(item), record_digits) for item in boundary
            ],
        }
    else:
        boundary_record = {"kind": "fixed-value", "value": boundary}
    if isinstance(path, AdaptivePath):
        if not isinstance(DEmatrix, RationalMatrixSystem):
            raise TypeError(
                "a fixed AdaptivePath requires its fixed RationalMatrixSystem in parallel mode"
            )
        path_record = {
            "kind": "fixed-adaptive-path",
            "value": adaptive_path_to_json(path, digits=80),
        }
    elif callable(path):
        path_record = {"kind": "factory", "value": path}
    elif isinstance(path, list):
        path_record = {
            "kind": "fixed-path",
            "values": [acb(item).str(80) for item in path],
        }
    else:
        path_record = {"kind": "fixed-value", "value": path}
    return {"system": system_record, "boundary": boundary_record, "path": path_record}


def _restore_parallel_inputs(records: dict[str, Any]) -> tuple[Any, Any, Any]:
    """恢复 worker 所需系统、边界和路径，不在子进程重新规划路径。"""

    system_record = records["system"]
    DEmatrix = (
        _rational_system_from_record(system_record)
        if system_record["kind"] == "fixed-rational-system"
        else system_record["value"]
    )
    boundary_record = records["boundary"]
    boundary = (
        column_vector(
            [
                _acb_ball_from_record(value, f"boundary[{index}]")
                for index, value in enumerate(boundary_record["values"])
            ]
        )
        if boundary_record["kind"] == "fixed-vector"
        else boundary_record["value"]
    )
    path_record = records["path"]
    if path_record["kind"] == "fixed-adaptive-path":
        path = adaptive_path_from_json(DEmatrix, path_record["value"])
    elif path_record["kind"] == "fixed-path":
        path = [acb(value) for value in path_record["values"]]
    else:
        path = path_record["value"]
    return DEmatrix, boundary, path


def _solve_samples(
    *,
    DEmatrix: Any,
    boundary: Any,
    path: Any,
    sample_arguments: tuple[Any, ...],
    sample_points: tuple[acb, ...],
    transport_order: int,
    transport_extra_order: int,
    transport_sample_count: int | None,
    transport_extra_sample_count: int | None,
    radius_fraction: float,
    nde_tolerance: arb,
    parallel_task_count: int,
) -> tuple[list[acb_mat], list[dict[str, Any]], int]:
    """按输入顺序求解一组 ep，并返回程序实际使用的 worker 数。"""

    effective_count = min(parallel_task_count, len(sample_arguments))
    input_records = _parallel_input_records(DEmatrix, boundary, path)
    record_digits = math.ceil(int(ctx.prec) / math.log2(10)) + 10
    payloads = [
        {
            "inputs": input_records,
            "sample_argument": (
                {"kind": "acb-ball", "value": _acb_ball_record(argument, record_digits)}
                if isinstance(argument, acb)
                else argument
            ),
            "sample_point": (
                point.str(record_digits)
                if point.imag.is_zero()
                else {"kind": "acb-ball", "value": _acb_ball_record(point, record_digits)}
            ),
            "transport_order": transport_order,
            "transport_extra_order": transport_extra_order,
            "transport_sample_count": transport_sample_count,
            "transport_extra_sample_count": transport_extra_sample_count,
            "radius_fraction": radius_fraction,
            "nde_tolerance": nde_tolerance.str(50),
            "working_precision_bits": int(ctx.prec),
        }
        for argument, point in zip(sample_arguments, sample_points)
    ]
    if effective_count == 1:
        completed = [_solve_sample_job(payload) for payload in payloads]
        values = [
            column_vector([
                (
                    _acb_ball_from_record(entry["value"], f"worker_value[{row}]")
                    if isinstance(entry, dict) and entry.get("kind") == "acb-ball"
                    else acb(entry)
                )
                for row, entry in enumerate(result[0])
            ])
            for result in completed
        ]
        reports = [result[1] for result in completed]
        return values, reports, effective_count
    try:
        pickle.dumps(payloads[0])
    except Exception as error:
        raise TypeError(
            "parallel regulator solves require non-fixed DEmatrix/boundary/path "
            "factories defined at module top level with pickleable inputs; fixed "
            "RationalMatrixSystem, ordinary boundary vectors, and paths are serialized "
            "automatically; use parallel_task_count=1 for local closures or an "
            "AnalyticMatrixSystem instance"
        ) from error
    with get_context("spawn").Pool(
        processes=effective_count, maxtasksperchild=1
    ) as pool:
        completed = pool.map(_solve_sample_job, payloads, chunksize=1)
    values = [
        column_vector([
            (
                _acb_ball_from_record(entry["value"], f"worker_value[{row}]")
                if isinstance(entry, dict) and entry.get("kind") == "acb-ball"
                else acb(entry)
            )
            for row, entry in enumerate(result[0])
        ])
        for result in completed
    ]
    reports = [result[1] for result in completed]
    return values, reports, effective_count


def _fit_coefficient_vectors(
    points: tuple[acb, ...],
    values: tuple[acb_mat, ...],
    leading_power: int,
) -> tuple[acb_mat, ...]:
    """一次求解多右端 Vandermonde 方阵，返回全部内部系数矢量。"""

    count = len(points)
    if count != len(values):
        raise ValueError("production point/value counts do not match")
    if not values:
        raise ValueError("production samples must not be empty")
    dimension = values[0].nrows()
    if any(value.nrows() != dimension or value.ncols() != 1 for value in values):
        raise ValueError("production NDE values have inconsistent dimensions")
    powers = range(leading_power, leading_power + count)
    matrix = acb_mat([[point**power for power in powers] for point in points])
    right_hand_sides = acb_mat(
        [[value[row, 0] for row in range(dimension)] for value in values]
    )
    try:
        solved = matrix.solve(right_hand_sides)
    except RuntimeError as error:
        raise SeriesValidationError(
            "Laurent interpolation matrix is singular or unresolved at the working precision"
        ) from error
    return tuple(
        acb_mat([[solved[degree, component]] for component in range(dimension)])
        for degree in range(count)
    )


def _evaluate_coefficients(
    coefficients: tuple[acb_mat, ...], leading_power: int, point: acb
) -> acb_mat:
    """用 Horner 法计算从 ``leading_power`` 开始的连续系数矢量。"""

    value = acb_mat(coefficients[-1])
    for coefficient in reversed(coefficients[:-1]):
        value = value * point + coefficient
    return value * (point**leading_power)


def _relative_vector_residual(predicted: acb_mat, actual: acb_mat) -> arb:
    """计算验证矢量残差；全零参考值退化为绝对无穷范数。"""

    difference = vector_norm_inf(actual - predicted)
    scale = vector_norm_inf(actual)
    if scale.contains(0):
        return difference
    return difference / scale


def _sample_value_vector(value: Any, name: str) -> acb_mat:
    """把外部标量或一维序列规范为 Acb 列矢量，并保留现有球精度。"""

    if isinstance(value, acb_mat):
        vector = acb_mat(value)
    elif isinstance(value, (list, tuple)):
        vector = column_vector(list(value))
    else:
        try:
            vector = acb_mat([[value if isinstance(value, acb) else acb(value)]])
        except (TypeError, ValueError) as error:
            raise TypeError(
                f"{name} must be an acb scalar, a one-dimensional sequence, or an acb_mat"
            ) from error
    if vector.ncols() != 1 or vector.nrows() < 1:
        raise ValueError(f"{name} must be a nonempty column vector")
    return vector


def _require_distinct_points(points: tuple[acb, ...], name: str) -> None:
    """拒绝互异的 Acb 点；零距离按 contains 判定（球重叠也算重复）。"""

    for index, point in enumerate(points):
        if any(abs(point - previous).contains(0) for previous in points[:index]):
            raise ValueError(f"{name} must contain distinct points")


def _explicit_series_points(values: Any, name: str) -> tuple[acb, ...]:
    """解析外部拟合点；拒绝 binary64、零点、重复点和隐式精度损失。"""

    if not isinstance(values, (list, tuple)) or not values:
        raise TypeError(f"{name} must be a nonempty sequence")
    points = tuple(
        _sample_pair(value, rationalize=False, name=f"{name}[{index}]")[1]
        for index, value in enumerate(values)
    )
    _require_distinct_points(points, name)
    return points


def fit_sampled_series(
    *,
    sample_points: Any,
    sample_values: Any,
    regularization_points: Any | None = None,
    maximum_power: int,
    leading_power: int = 0,
    validation_points: Any = (),
    validation_values: Any = (),
    validation_regularization_points: Any | None = None,
    validation_tolerance: Any | None = None,
    series_parameter: str = "ep",
) -> SeriesReconstructionResult:
    """从调用方已计算的 regulator 样本重构连续 Laurent/Taylor 系数。

    sample_points 与 sample_values 一一对应。缺省时 sample_points 同时作为插值变量；
    显式给出 regularization_points 时，前者只保存源采样标签，后者才代入
    Vandermonde。例如同一组正 delta 标签可分别配对 +delta 和 -delta，以重构实轴
    奇点两侧的 lateral 级数。validation 使用相同的可选配对规则。

    标量样本会提升为一维列矢量，多分量输入必须保持同一维数。全部生产点一次进入
    方阵求解，不使用最小二乘或伪逆；验证点不参与拟合，残差超限时 fail closed。
    本接口不改变 FLINT 全局精度，调用方必须先配置工作精度并在该精度下生成样本。
    """

    resolved_leading = _require_integer(
        leading_power, "leading_power", minimum=-10**9
    )
    resolved_maximum = _require_integer(
        maximum_power, "maximum_power", minimum=-10**9
    )
    if resolved_maximum < resolved_leading:
        raise ValueError("maximum_power must not be below leading_power")
    if not isinstance(series_parameter, str) or not series_parameter.strip():
        raise ValueError("series_parameter must be a nonempty string")

    points = _explicit_series_points(sample_points, "sample_points")
    explicit_regularization_points = regularization_points is not None
    fit_points = (
        _explicit_series_points(regularization_points, "regularization_points")
        if explicit_regularization_points
        else points
    )
    if len(points) != len(fit_points):
        raise ValueError("sample and regularization point counts do not match")
    if not isinstance(sample_values, (list, tuple)):
        raise TypeError("sample_values must be a sequence")
    values = tuple(
        _sample_value_vector(value, f"sample_values[{index}]")
        for index, value in enumerate(sample_values)
    )
    if len(points) != len(values):
        raise ValueError("sample point/value counts do not match")
    internal_maximum = resolved_leading + len(points) - 1
    if resolved_maximum > internal_maximum:
        raise ValueError(
            "sample count does not cover the requested leading_power..maximum_power range"
        )

    internal_coefficients = _fit_coefficient_vectors(
        fit_points, values, resolved_leading
    )
    returned_count = resolved_maximum - resolved_leading + 1
    returned_coefficients = internal_coefficients[:returned_count]

    if bool(validation_points) != bool(validation_values):
        raise ValueError(
            "validation_points and validation_values must be supplied together"
        )
    resolved_validation_points: tuple[acb, ...] = ()
    resolved_validation_fit_points: tuple[acb, ...] = ()
    resolved_validation_values: tuple[acb_mat, ...] = ()
    validation_reports: list[dict[str, Any]] = []
    maximum_residual: arb | None = None
    tolerance = None
    if validation_regularization_points is not None and not validation_points:
        raise ValueError(
            "validation_regularization_points require validation_points and validation_values"
        )
    if validation_points:
        resolved_validation_points = _explicit_series_points(
            validation_points, "validation_points"
        )
        if any(
            abs(validation_point - production_point).contains(0)
            for validation_point in resolved_validation_points
            for production_point in points
        ):
            raise ValueError("validation points must be distinct from production points")
        if not isinstance(validation_values, (list, tuple)):
            raise TypeError("validation_values must be a sequence")
        resolved_validation_values = tuple(
            _sample_value_vector(value, f"validation_values[{index}]")
            for index, value in enumerate(validation_values)
        )
        if len(resolved_validation_points) != len(resolved_validation_values):
            raise ValueError("validation point/value counts do not match")
        resolved_validation_fit_points = (
            _explicit_series_points(
                validation_regularization_points,
                "validation_regularization_points",
            )
            if validation_regularization_points is not None
            else resolved_validation_points
        )
        if len(resolved_validation_points) != len(resolved_validation_fit_points):
            raise ValueError(
                "validation and validation regularization point counts do not match"
            )
        if any(
            abs(validation_point - production_point).contains(0)
            for validation_point in resolved_validation_fit_points
            for production_point in fit_points
        ):
            raise ValueError(
                "validation regularization points must be distinct from production regularization points"
            )
        if validation_tolerance is not None:
            tolerance = arb(str(validation_tolerance))
            if tolerance <= 0:
                raise ValueError("validation_tolerance must be positive")

        for source_point, fit_point, actual in zip(
            resolved_validation_points,
            resolved_validation_fit_points,
            resolved_validation_values,
        ):
            predicted = _evaluate_coefficients(
                returned_coefficients, resolved_leading, fit_point
            )
            residual = _relative_vector_residual(predicted, actual)
            maximum_residual = (
                residual
                if maximum_residual is None or residual > maximum_residual
                else maximum_residual
            )
            passed = tolerance is None or bool(residual < tolerance)
            validation_reports.append(
                {
                    "point": source_point.str(50),
                    **(
                        {"regularization_point": fit_point.str(50)}
                        if validation_regularization_points is not None
                        else {}
                    ),
                    "predicted": _vector_to_strings(predicted),
                    "actual": _vector_to_strings(actual),
                    "series_relative_residual": residual.str(30),
                    "passed": passed,
                }
            )
            if not passed:
                raise SeriesValidationError(
                    f"series validation failed at {fit_point.str(30)}: residual "
                    f"{residual.str(20)} is not below {tolerance.str(20)}"
                )

    powers = tuple(range(resolved_leading, resolved_maximum + 1))
    return SeriesReconstructionResult(
        series_parameter,
        powers,
        tuple(acb_mat(value) for value in returned_coefficients),
        resolved_leading,
        resolved_maximum,
        internal_maximum,
        points,
        tuple(acb_mat(value) for value in values),
        resolved_validation_points,
        tuple(acb_mat(value) for value in resolved_validation_values),
        {
            "sample_source": "user-supplied-values",
            "sample_count": len(points),
            "validation_source": (
                "user-supplied-values" if resolved_validation_points else "not-run"
            ),
            "validation_sample_count": len(resolved_validation_points),
            "validation_tolerance": (
                tolerance.str(30) if tolerance is not None else None
            ),
            **(
                {"regularization_point_source": "user-supplied-values"}
                if explicit_regularization_points
                else {}
            ),
            **(
                {
                    "validation_regularization_point_source":
                    "user-supplied-values"
                }
                if validation_regularization_points is not None
                else {}
            ),
        },
        {
            "validation_samples": validation_reports,
            "maximum_validation_relative_residual": (
                maximum_residual.str(30) if maximum_residual is not None else None
            ),
            "internal_buffer_powers": internal_maximum - resolved_maximum,
            "interpolation": "square Acb Vandermonde with all production samples",
            "least_squares_used": False,
        },
        regularization_points=(
            tuple(acb(point) for point in fit_points)
            if explicit_regularization_points
            else None
        ),
        validation_regularization_points=(
            tuple(acb(point) for point in resolved_validation_fit_points)
            if validation_regularization_points is not None
            else None
        ),
    )

def reconstruct_series_solution(
    *,
    DEmatrix: MatrixSystem | Callable[[Any], MatrixSystem],
    boundary: Any,
    path: list[acb] | Callable[[Any, MatrixSystem], list[acb]],
    maximum_power: int,
    leading_power: int,
    leading_power_certificate: dict[str, Any],
    series_parameter: str = "ep",
    goal_digits: int = 30,
    sample_points: Any = AUTOMATIC,
    sample_angle_range: Any = AUTOMATIC,
    sample_count: Any = AUTOMATIC,
    base_sample: Any = AUTOMATIC,
    sample_spacing: Any = 0.01,
    working_precision_digits: Any = AUTOMATIC,
    extra_working_precision: float = 0.0,
    transport_order: Any = AUTOMATIC,
    transport_extra_order: Any = AUTOMATIC,
    transport_sample_count: Any = AUTOMATIC,
    transport_extra_sample_count: Any = AUTOMATIC,
    radius_fraction: float = 0.60,
    guard_bits: int = 32,
    validation_sample_count: int = 2,
    validation_points: Any = AUTOMATIC,
    validation_scale: Any = 0.5,
    validation_tolerance: Any = AUTOMATIC,
    fit_extra_order: int = 2,
    initial_internal_maximum_power: Any = AUTOMATIC,
    fit_order_increment: int = 2,
    fit_max_rounds: int = 3,
    maximum_samples: int = 100,
    rationalize_sample_points: bool = True,
    parallel_task_count: int = DEFAULT_PARALLEL_TASK_COUNT,
    output_layout: Any | None = None,
    result_name: str = "series_reconstruction",
) -> SeriesReconstructionResult:
    """用固定 regulator 的完整 NDE 解重构最终解矢量的 Laurent 系数。

    参数 ``DEmatrix``、``boundary`` 和 ``path`` 可直接给固定对象；依赖 regulator 时，
    分别给出 ``DEmatrix(ep)``、``boundary(ep)`` 和 ``path(ep, system)``。自动生成的
    exact-rational 样本以 ``fmpq`` 传给这些工厂。普通点使用列向量；正则奇点起点使用
    ``frobenius_boundary`` 生成的 exact ``{a,b,C}`` 记录，并由基础输运验证 indicial root、
    最高 log 次数和领头向量相容性。

    ``leading_power`` 及其证书必须由调用方在数值求解前从符号边界与 DE 得到；仅有
    regulator callable 时 FlintNDE 不尝试用数值点猜测结构最低阶。返回值包含用户请求到
    ``maximum_power`` 的系数；内部拟合缺省至少多两阶，并用
    这些高阶项在独立点估计截断误差。验证失败时每轮缺省再增加两阶，只求解新增生产
    点；已有生产点、验证点、容差和数值设置均复用。达到轮数或样本上限仍不满足
    精度时返回当前最佳系数并显式标记未达到目标精度，不使用最小二乘、伪逆或静默
    降精度。

    显式 ``sample_points`` 是有序候选池：首轮只消费当前内部最高幂所需的前缀，验证
    失败后再从剩余候选点增量取用，绝不生成池外点。``initial_internal_maximum_power``
    可显式指定首轮内部最高幂；缺省保持原自动规则。候选池比缺省首轮点数更短时，仍按
    其完整长度拟合，只要足以覆盖用户请求的系数。

    ``sample_angle_range=(theta_min, theta_max)`` 以弧度限制自动样本的复角度。模长仍由
    原 AMFlow 式精度策略自动决定且不设置模长上限；角域视为开区间，程序在内部均匀
    选择至多三条射线并循环分配样本，因此不会刻意贴近边界。缺省 ``automatic`` 严格
    保留原正实轴网格。
    """

    maximum_power = _require_integer(maximum_power, "maximum_power", minimum=-10**9)
    goal_digits = _require_integer(goal_digits, "goal_digits", minimum=1)
    validation_sample_count = _require_integer(
        validation_sample_count, "validation_sample_count", minimum=1
    )
    fit_extra_order = _require_integer(
        fit_extra_order, "fit_extra_order", minimum=0
    )
    resolved_initial_internal_maximum = _require_automatic_or_integer(
        initial_internal_maximum_power,
        "initial_internal_maximum_power",
        minimum=maximum_power,
    )
    fit_order_increment = _require_integer(
        fit_order_increment, "fit_order_increment", minimum=1
    )
    fit_max_rounds = _require_integer(
        fit_max_rounds, "fit_max_rounds", minimum=1
    )
    maximum_samples = _require_integer(maximum_samples, "maximum_samples", minimum=1)
    parallel_task_count = _require_integer(
        parallel_task_count, "parallel_task_count", minimum=1
    )
    guard_bits = _require_integer(guard_bits, "guard_bits", minimum=0)
    if not isinstance(series_parameter, str) or not series_parameter.strip():
        raise ValueError("series_parameter must be a nonempty string")
    if extra_working_precision < 0:
        raise ValueError("extra_working_precision must be nonnegative")
    if not 0 < radius_fraction < 1:
        raise ValueError("radius_fraction must lie strictly between zero and one")
    if not isinstance(rationalize_sample_points, bool):
        raise TypeError("rationalize_sample_points must be bool")
    if sample_points != AUTOMATIC and sample_angle_range != AUTOMATIC:
        raise ValueError("sample_angle_range cannot be combined with explicit sample_points")

    resolved_leading = _require_integer(
        leading_power, "leading_power", minimum=-10**9
    )
    if not isinstance(leading_power_certificate, dict) or set(leading_power_certificate) != {
        "status", "leading_power", "method"
    }:
        raise ValueError(
            "leading_power_certificate must contain exactly status, leading_power and method"
        )
    if leading_power_certificate["status"] != "certified":
        raise ValueError("leading_power_certificate.status must be 'certified'")
    certified_power = _require_integer(
        leading_power_certificate["leading_power"],
        "leading_power_certificate.leading_power",
        minimum=-10**9,
    )
    if certified_power != resolved_leading:
        raise ValueError("leading_power does not match leading_power_certificate")
    if not isinstance(leading_power_certificate["method"], str) or not leading_power_certificate[
        "method"
    ].strip():
        raise ValueError("leading_power_certificate.method must be a nonempty string")

    requested_coefficient_count = maximum_power - resolved_leading + 1
    empirical_count = max(
        math.ceil(
            2.5 * (maximum_power - resolved_leading)
            + max(0, -resolved_leading)
        ),
        requested_coefficient_count + fit_extra_order,
    )
    requested_initial_count = (
        empirical_count
        if resolved_initial_internal_maximum is None
        else resolved_initial_internal_maximum - resolved_leading + 1
    )
    if sample_points == AUTOMATIC:
        if (
            resolved_initial_internal_maximum is not None
            and sample_count != AUTOMATIC
            and sample_count != requested_initial_count
        ):
            raise ValueError(
                "sample_count must match the count required by "
                "initial_internal_maximum_power"
            )
        initial_count = (
            requested_initial_count
            if sample_count == AUTOMATIC
            else _require_integer(
                sample_count,
                "sample_count",
                minimum=requested_coefficient_count + fit_extra_order,
            )
        )
        capacity_count = min(
            maximum_samples,
            initial_count + fit_order_increment * (fit_max_rounds - 1),
        )
        if initial_count > maximum_samples:
            raise ValueError(
                f"automatic production needs {initial_count} samples, "
                f"above maximum_samples={maximum_samples}"
            )
    else:
        if not isinstance(sample_points, (list, tuple)) or not sample_points:
            raise TypeError('sample_points must be "automatic" or a nonempty sequence')
        candidate_count = len(sample_points)
        if sample_count != AUTOMATIC:
            explicit_count = _require_integer(
                sample_count, "sample_count", minimum=requested_coefficient_count
            )
            if explicit_count != candidate_count:
                raise ValueError("sample_count must equal the explicit sample_points length")
        if (
            resolved_initial_internal_maximum is not None
            and requested_initial_count > candidate_count
        ):
            raise ValueError(
                "explicit sample_points candidate pool has "
                f"{candidate_count} points but initial_internal_maximum_power="
                f"{resolved_initial_internal_maximum} requires {requested_initial_count}"
            )
        initial_count = min(requested_initial_count, candidate_count)
        if initial_count < requested_coefficient_count:
            raise ValueError(
                "explicit sample_points candidate pool does not cover the requested powers"
            )
        capacity_count = candidate_count

    capacity_plan = _resolve_plan(
        leading_power=resolved_leading,
        maximum_power=maximum_power,
        goal_digits=goal_digits,
        sample_points=sample_points,
        sample_angle_range=sample_angle_range,
        sample_count=capacity_count,
        base_sample=base_sample,
        sample_spacing=sample_spacing,
        working_precision_digits=working_precision_digits,
        extra_working_precision=extra_working_precision,
        transport_order=transport_order,
        transport_extra_order=transport_extra_order,
        transport_sample_count=transport_sample_count,
        transport_extra_sample_count=transport_extra_sample_count,
        validation_sample_count=validation_sample_count,
        validation_points=validation_points,
        validation_scale=validation_scale,
        maximum_samples=maximum_samples,
        rationalize_sample_points=rationalize_sample_points,
        fit_extra_order=fit_extra_order,
    )
    plan = replace(
        capacity_plan,
        sample_arguments=capacity_plan.sample_arguments[:initial_count],
        sample_points=capacity_plan.sample_points[:initial_count],
        sample_count=initial_count,
    )
    configure_working_precision(plan.working_precision_digits, guard_bits)
    tolerance = (
        arb(f"1e-{goal_digits}")
        if validation_tolerance == AUTOMATIC
        else arb(str(validation_tolerance))
    )
    if tolerance <= 0:
        raise ValueError("validation_tolerance must be positive")

    # 验证点只求解一次并保存在独立缓存；它们从不进入 Vandermonde 拟合矩阵。
    validation_actual_values, validation_solve_reports, validation_parallel_count = _solve_samples(
        DEmatrix=DEmatrix,
        boundary=boundary,
        path=path,
        sample_arguments=plan.validation_arguments,
        sample_points=plan.validation_points,
        transport_order=plan.transport_order,
        transport_extra_order=plan.transport_extra_order,
        transport_sample_count=plan.transport_sample_count,
        transport_extra_sample_count=plan.transport_extra_sample_count,
        radius_fraction=radius_fraction,
        nde_tolerance=tolerance,
        parallel_task_count=parallel_task_count,
    )
    production_values: list[acb_mat] = []
    production_reports: list[dict[str, Any]] = []
    expansion_history: list[dict[str, Any]] = []
    production_parallel_counts: list[int] = []
    validation_values = [acb_mat(value) for value in validation_actual_values]
    accepted = False
    validation_reports: list[dict[str, Any]] = []
    internal_coefficients: tuple[acb_mat, ...] = ()
    current_count = 0
    final_plan = plan
    for round_index in range(1, fit_max_rounds + 1):
        target_count = min(
            capacity_plan.sample_count,
            plan.sample_count + fit_order_increment * (round_index - 1),
        )
        new_arguments = capacity_plan.sample_arguments[current_count:target_count]
        new_points = capacity_plan.sample_points[current_count:target_count]
        if new_arguments:
            new_values, new_reports, effective_parallel = _solve_samples(
                DEmatrix=DEmatrix,
                boundary=boundary,
                path=path,
                sample_arguments=new_arguments,
                sample_points=new_points,
                transport_order=capacity_plan.transport_order,
                transport_extra_order=capacity_plan.transport_extra_order,
                transport_sample_count=capacity_plan.transport_sample_count,
                transport_extra_sample_count=capacity_plan.transport_extra_sample_count,
                radius_fraction=radius_fraction,
                nde_tolerance=tolerance,
                parallel_task_count=parallel_task_count,
            )
            production_values.extend(new_values)
            production_reports.extend(new_reports)
            production_parallel_counts.append(effective_parallel)
        current_count = target_count
        current_points = capacity_plan.sample_points[:current_count]
        internal_coefficients = _fit_coefficient_vectors(
            current_points, tuple(production_values), resolved_leading
        )
        validation_reports = []
        residuals: list[arb] = []
        for point, actual, solve_report in zip(
            capacity_plan.validation_points,
            validation_values,
            validation_solve_reports,
        ):
            predicted = _evaluate_coefficients(
                internal_coefficients, resolved_leading, point
            )
            residual = _relative_vector_residual(predicted, actual)
            residuals.append(residual)
            validation_reports.append(
                {
                    **solve_report,
                    "predicted": _vector_to_strings(predicted),
                    "series_relative_residual": residual.str(30),
                    "passed": bool(residual < tolerance),
                }
            )
        accepted = all(residual < tolerance for residual in residuals)
        expansion_history.append(
            {
                "round": round_index,
                "internal_maximum_power": resolved_leading + current_count - 1,
                "new_production_sample_count": len(new_arguments),
                "reused_production_sample_count": current_count - len(new_arguments),
                "reused_validation_sample_count": (
                    0 if round_index == 1 else len(validation_values)
                ),
                "maximum_validation_relative_residual": max(residuals).str(30),
                "passed": accepted,
            }
        )
        final_plan = replace(
            capacity_plan,
            sample_arguments=capacity_plan.sample_arguments[:current_count],
            sample_points=current_points,
            sample_count=current_count,
        )
        if accepted:
            break
        if current_count >= capacity_plan.sample_count:
            break

    precision_failure_reason = None
    precision_warning = None
    if not accepted:
        last_residual = expansion_history[-1]["maximum_validation_relative_residual"]
        if sample_points != AUTOMATIC and current_count >= capacity_plan.sample_count:
            precision_failure_reason = "candidate_pool_exhausted"
            detail = (
                "the explicit sample_points candidate pool was exhausted; no points "
                "outside the user-supplied pool were generated"
            )
        elif current_count >= maximum_samples:
            precision_failure_reason = "maximum_samples_reached"
            detail = f"maximum_samples={maximum_samples} was reached"
        else:
            precision_failure_reason = "fit_round_limit_reached"
            detail = f"fit_max_rounds={fit_max_rounds} was reached"
        precision_warning = (
            "series validation remained above tolerance after incremental fitting "
            f"through internal power {resolved_leading + current_count - 1}: "
            f"maximum residual {last_residual}, tolerance {tolerance.str(20)}; {detail}. "
            "Returning the current best coefficients without precision certification."
        )
        warnings.warn(precision_warning, RuntimeWarning, stacklevel=2)

    returned_count = requested_coefficient_count
    returned_coefficients = internal_coefficients[:returned_count]
    internal_maximum = resolved_leading + final_plan.sample_count - 1
    effective_parameters = {
        "goal_digits": goal_digits,
        "sample_source": final_plan.source,
        "sample_count": final_plan.sample_count,
        "initial_sample_count": plan.sample_count,
        "sample_candidate_count": capacity_plan.sample_count,
        "unused_sample_candidate_count": capacity_plan.sample_count - final_plan.sample_count,
        "base_sample": final_plan.base_sample,
        "sample_spacing": str(sample_spacing),
        "sample_angle_range": (
            AUTOMATIC if sample_angle_range == AUTOMATIC else list(_angle_range_texts(sample_angle_range))
        ),
        "alpha_epsilon": final_plan.alpha_epsilon,
        "base_precision_digits": final_plan.base_precision_digits,
        "working_precision_digits": final_plan.working_precision_digits,
        "extra_working_precision": extra_working_precision,
        "transport_order": final_plan.transport_order,
        "transport_extra_order": final_plan.transport_extra_order,
        "transport_sample_count": final_plan.transport_sample_count or AUTOMATIC,
        "transport_extra_sample_count": final_plan.transport_extra_sample_count or AUTOMATIC,
        "radius_fraction": radius_fraction,
        "guard_bits": guard_bits,
        "validation_source": final_plan.validation_source,
        "validation_sample_count": len(final_plan.validation_points),
        "validation_scale": (
            str(validation_scale) if validation_points == AUTOMATIC else None
        ),
        "validation_tolerance": tolerance.str(30),
        "maximum_samples": maximum_samples,
        "fit_extra_order": fit_extra_order,
        "initial_internal_maximum_power_requested": (
            initial_internal_maximum_power
            if initial_internal_maximum_power != AUTOMATIC
            else AUTOMATIC
        ),
        "fit_order_increment": fit_order_increment,
        "fit_max_rounds": fit_max_rounds,
        "fit_rounds_used": len(expansion_history),
        "rationalize_sample_points": rationalize_sample_points,
        "parallel_task_count_requested": parallel_task_count,
        "production_parallel_task_count_effective": max(production_parallel_counts),
        "validation_parallel_task_count_effective": validation_parallel_count,
        "precision_target_met": accepted,
        "precision_failure_reason": precision_failure_reason,
    }
    diagnostics = {
        "leading_power_certificate": dict(leading_power_certificate),
        "production_solves": production_reports,
        "validation_solves": validation_reports,
        "fit_expansion_history": expansion_history,
        "production_values_reused": sum(
            item["reused_production_sample_count"] for item in expansion_history[1:]
        ),
        "validation_values_reused": (
            len(validation_values) * max(0, len(expansion_history) - 1)
        ),
        "internal_buffer_powers": internal_maximum - maximum_power,
        "interpolation": "square Acb Vandermonde with all production samples",
        "least_squares_used": False,
        "precision_target_met": accepted,
        "precision_failure_reason": precision_failure_reason,
        "precision_warning": precision_warning,
    }
    powers = tuple(range(resolved_leading, maximum_power + 1))
    result = SeriesReconstructionResult(
        series_parameter,
        powers,
        tuple(acb_mat(value) for value in returned_coefficients),
        resolved_leading,
        maximum_power,
        internal_maximum,
        final_plan.sample_points,
        tuple(acb_mat(value) for value in production_values),
        final_plan.validation_points,
        tuple(acb_mat(value) for value in validation_values),
        effective_parameters,
        diagnostics,
    )
    if output_layout is not None:
        filename = f"{result_name}_series_reconstruction.json"
        relative_output = str(Path("regularization") / filename)
        result = replace(result, output_file=relative_output)
        output_layout.write_json("regularization", filename, result.to_summary())
    return result
