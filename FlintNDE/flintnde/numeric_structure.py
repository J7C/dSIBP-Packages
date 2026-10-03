"""浮点 residue 的精度感知 indicial、Jordan 与 resonance 诊断。

本模块只在调用者没有 exact rational 矩阵时使用。它按可靠输入精度和矩阵尺度
生成可覆盖的零斩杀线，记录最大被判零量；无法稳定判定的谱结构会 fail closed，
不会把不确定的奇点当成无 log 的普通幂级数。
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from typing import Any

from flint import acb, acb_mat, arb


_BINARY64_RELIABLE_DECIMAL_DIGITS = 15


def _scalar_precision_digits(value: Any) -> tuple[int, str] | None:
    """返回浮点标量可证明的可靠十进制位数；exact/无元数据记录返回 ``None``。"""

    if isinstance(value, (float, complex)):
        return _BINARY64_RELIABLE_DECIMAL_DIGITS, "python_binary64_input"
    if isinstance(value, (acb, arb)):
        accuracy_bits = int(value.rel_accuracy_bits())
        if 0 < accuracy_bits < 1_000_000_000:
            return max(1, math.floor(accuracy_bits / math.log2(10))), "flint_ball_input"
        return None
    if isinstance(value, (tuple, list)) and len(value) == 2:
        candidates = [
            candidate
            for item in value
            if (candidate := _scalar_precision_digits(item)) is not None
        ]
        return min(candidates, default=None, key=lambda item: item[0])
    return None


def _numeric_scalar(value: Any) -> acb:
    """把用户数值转换为 Acb；Python float 先按其十进制打印值转换。"""

    if isinstance(value, acb):
        return acb(value)
    if isinstance(value, arb):
        return acb(value)
    if isinstance(value, complex):
        return acb(str(value.real), str(value.imag))
    if isinstance(value, (tuple, list)) and len(value) == 2:
        return acb(str(value[0]), str(value[1]))
    return acb(str(value))


def _numeric_matrix(records: list[list[Any]]) -> acb_mat:
    """把二维数值记录转换为 Acb 矩阵。"""

    return acb_mat([[_numeric_scalar(value) for value in row] for row in records])


def _midpoint(value: acb) -> acb:
    """提取 Acb 高精度中点，避免把结构诊断降到 binary64。"""

    return acb(value.real.mid(), value.imag.mid())


def _matrix_midpoint(matrix: acb_mat) -> acb_mat:
    """逐元提取 Acb 矩阵中点。"""

    return acb_mat(
        [
            [_midpoint(matrix[row, column]) for column in range(matrix.ncols())]
            for row in range(matrix.nrows())
        ]
    )


def _identity(dimension: int) -> acb_mat:
    """构造 Acb 单位矩阵。"""

    return acb_mat(
        [[1 if row == column else 0 for column in range(dimension)] for row in range(dimension)]
    )


def _entry_scale(matrices: list[acb_mat]) -> arb:
    """返回一组矩阵的最大元素绝对值中点。"""

    magnitudes = [
        abs(matrix[row, column]).mid()
        for matrix in matrices
        for row in range(matrix.nrows())
        for column in range(matrix.ncols())
    ]
    return max(magnitudes) if magnitudes else arb(0)


def _matrix_max_abs(matrix: acb_mat) -> arb:
    """返回矩阵最大元素绝对值中点。"""

    return _entry_scale([matrix])


def _acb_record(value: acb) -> dict[str, str]:
    """把 Acb 中点保存为 JSON 友好的实部、虚部十进制字符串。"""

    midpoint = _midpoint(value)
    return {
        "real": midpoint.real.mid().str(radius=False),
        "imag": midpoint.imag.mid().str(radius=False),
    }


def _matrix_records(matrix: acb_mat) -> list[list[dict[str, str]]]:
    """把 Acb 矩阵保存为数值 manifest 记录。"""

    return [
        [_acb_record(matrix[row, column]) for column in range(matrix.ncols())]
        for row in range(matrix.nrows())
    ]


def _vector_records(vector: acb_mat) -> list[dict[str, str]]:
    """把 Acb 列向量保存为数值 manifest 记录。"""

    return [_acb_record(vector[row, 0]) for row in range(vector.nrows())]


@dataclass(frozen=True)
class NumericalFrobeniusOptions:
    """配置浮点结构诊断使用的输入精度和相对零斩杀线。"""

    precision_digits: int | None = None
    relative_zero_tolerance: str | None = None

    def resolve(
        self,
        dimension: int,
        input_precision: tuple[int, str] | None,
    ) -> tuple[int, arb, str]:
        """按输入可靠精度返回十进制精度、相对阈值及其来源。"""

        declared_digits = self.precision_digits
        if declared_digits is not None and declared_digits <= 0:
            raise ValueError("precision_digits must be positive")
        if declared_digits is None and input_precision is None:
            raise ValueError(
                "numerical input precision is unknown; set input_precision_digits on "
                "NumericalRegularSingularSystem or precision_digits in NumericalFrobeniusOptions"
            )
        if declared_digits is None:
            digits, precision_source = input_precision
        elif input_precision is None:
            digits, precision_source = declared_digits, "user"
        else:
            inferred_digits, inferred_source = input_precision
            digits = min(declared_digits, inferred_digits)
            precision_source = (
                "user" if declared_digits <= inferred_digits else inferred_source
            )
        if self.relative_zero_tolerance is None:
            relative = arb(dimension) * arb(10) ** (-digits)
            tolerance_source = f"dimension_times_10^-precision ({precision_source})"
        else:
            relative = arb(self.relative_zero_tolerance)
            if relative <= 0:
                raise ValueError("relative_zero_tolerance must be positive")
            tolerance_source = "user"
        return digits, relative, tolerance_source


@dataclass(frozen=True)
class NumericalRegularSingularSystem:
    """保存浮点 ``Y'=(R/z+sum A_m z^m)Y`` 的局部矩阵数据。"""

    residue: list[list[Any]]
    regular_coefficients: tuple[list[list[Any]], ...]
    name: str = "numerical-regular-singular-system"
    input_precision_digits: int | None = None

    @property
    def dimension(self) -> int:
        """返回矩阵系统维数。"""

        return len(self.residue)

    def acb_data(self) -> tuple[acb_mat, list[acb_mat]]:
        """返回 Acb residue 和解析部分系数并验证维数。"""

        residue = _numeric_matrix(self.residue)
        if residue.nrows() != residue.ncols():
            raise ValueError(f"{self.name}: residue matrix must be square")
        regular = [_numeric_matrix(records) for records in self.regular_coefficients]
        if any(
            matrix.nrows() != residue.nrows() or matrix.ncols() != residue.ncols()
            for matrix in regular
        ):
            raise ValueError(f"{self.name}: regular coefficient dimension mismatch")
        return residue, regular

    def input_precision(self) -> tuple[int, str] | None:
        """返回输入可靠精度；显式元数据优先，机器浮点或 ball 精度作为上限。"""

        if self.input_precision_digits is not None and self.input_precision_digits <= 0:
            raise ValueError("input_precision_digits must be positive")
        candidates = [
            candidate
            for matrix in (self.residue, *self.regular_coefficients)
            for row in matrix
            for value in row
            if (candidate := _scalar_precision_digits(value)) is not None
        ]
        inferred = min(candidates, default=None, key=lambda item: item[0])
        if self.input_precision_digits is None:
            return inferred
        if inferred is None or self.input_precision_digits <= inferred[0]:
            return self.input_precision_digits, "system_input_precision"
        return inferred


class _ZeroAudit:
    """集中记录所有按阈值判零或保留的结构量。"""

    def __init__(self, cutoff: arb) -> None:
        self.cutoff = cutoff
        self.zeroed_count = 0
        self.largest_zeroed = arb(0)
        self.smallest_retained: arb | None = None

    def is_zero(self, value: acb | arb) -> bool:
        """按统一绝对阈值判零并更新诊断统计。"""

        magnitude = (abs(value) if isinstance(value, acb) else abs(value)).mid()
        if magnitude <= self.cutoff:
            self.zeroed_count += 1
            self.largest_zeroed = max(self.largest_zeroed, magnitude)
            return True
        if self.smallest_retained is None or magnitude < self.smallest_retained:
            self.smallest_retained = magnitude
        return False

    def records(self) -> dict[str, Any]:
        """返回可写入 manifest 的零斩杀诊断。"""

        return {
            "zeroed_quantity_count": self.zeroed_count,
            "largest_zeroed_absolute_value": self.largest_zeroed.str(radius=False),
            "smallest_retained_absolute_value": (
                None
                if self.smallest_retained is None
                else self.smallest_retained.str(radius=False)
            ),
        }


def _warn_threshold_classification(system_name: str, diagnostics: dict[str, Any]) -> None:
    """提醒用户本次 log 结构依赖数值斩杀，并完整报告绝对/相对尺度。"""

    largest_relative = diagnostics["largest_zeroed_relative_to_matrix_scale"]
    smallest_relative = diagnostics["smallest_retained_relative_to_matrix_scale"]
    if largest_relative is None:
        largest_relative = "undefined (matrix scale is zero)"
    if smallest_relative is None:
        smallest_relative = "undefined (matrix scale is zero or no quantity was retained)"

    warnings.warn(
        (
            f"{system_name}: numerical Frobenius classification used relative cutoff "
            f"{diagnostics['relative_zero_tolerance']} and absolute cutoff "
            f"{diagnostics['absolute_zero_cutoff']} at matrix scale "
            f"{diagnostics['matrix_scale']}; largest zeroed value: absolute "
            f"{diagnostics['largest_zeroed_absolute_value']}, relative to matrix scale "
            f"{largest_relative}; smallest retained value: absolute "
            f"{diagnostics['smallest_retained_absolute_value']}, relative to matrix scale "
            f"{smallest_relative}. Use exact rational input for a certified "
            "log-structure decision."
        ),
        UserWarning,
        stacklevel=2,
    )


def _threshold_matrix(matrix: acb_mat, audit: _ZeroAudit) -> acb_mat:
    """把小于阈值的输入元素显式设为零。"""

    result = _matrix_midpoint(matrix)
    for row in range(result.nrows()):
        for column in range(result.ncols()):
            if audit.is_zero(result[row, column]):
                result[row, column] = 0
    return result


def _rref(matrix: acb_mat, audit: _ZeroAudit) -> tuple[acb_mat, list[int]]:
    """用精度感知 pivot 斩杀构造数值 RREF 与 pivot 列。"""

    reduced = _matrix_midpoint(matrix)
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(reduced.ncols()):
        candidates = list(range(pivot_row, reduced.nrows()))
        if not candidates:
            break
        best_row = max(candidates, key=lambda row: abs(reduced[row, column]).mid())
        pivot = reduced[best_row, column]
        if audit.is_zero(pivot):
            continue
        if best_row != pivot_row:
            for current_column in range(reduced.ncols()):
                reduced[pivot_row, current_column], reduced[best_row, current_column] = (
                    reduced[best_row, current_column],
                    reduced[pivot_row, current_column],
                )
        pivot = reduced[pivot_row, column]
        for current_column in range(reduced.ncols()):
            reduced[pivot_row, current_column] /= pivot
        for row in range(reduced.nrows()):
            if row == pivot_row:
                continue
            factor = reduced[row, column]
            if audit.is_zero(factor):
                reduced[row, column] = 0
                continue
            for current_column in range(reduced.ncols()):
                reduced[row, current_column] -= factor * reduced[pivot_row, current_column]
        pivot_columns.append(column)
        pivot_row += 1
    return reduced, pivot_columns


def _nullspace(matrix: acb_mat, audit: _ZeroAudit) -> list[acb_mat]:
    """由数值 RREF 构造右零空间基。"""

    reduced, pivot_columns = _rref(matrix, audit)
    free_columns = [
        column for column in range(matrix.ncols()) if column not in pivot_columns
    ]
    basis: list[acb_mat] = []
    for free_column in free_columns:
        vector = acb_mat(matrix.ncols(), 1)
        vector[free_column, 0] = 1
        for row, column in enumerate(pivot_columns):
            vector[column, 0] = -reduced[row, free_column]
        basis.append(vector)
    return basis


def _matrix_is_zero(matrix: acb_mat, audit: _ZeroAudit) -> bool:
    """判断矩阵所有元素是否都落在斩杀线内。"""

    maximum = _matrix_max_abs(matrix)
    return audit.is_zero(maximum)


def _cluster_eigenvalues(eigenvalues: list[acb], audit: _ZeroAudit) -> list[list[int]]:
    """按同一绝对阈值聚类近重复特征值。"""

    groups: list[list[int]] = []
    centers: list[acb] = []
    for index, eigenvalue in enumerate(eigenvalues):
        value = _midpoint(eigenvalue)
        target = next(
            (
                position
                for position, center in enumerate(centers)
                if audit.is_zero(value - center)
            ),
            None,
        )
        if target is None:
            groups.append([index])
            centers.append(value)
        else:
            groups[target].append(index)
            centers[target] = sum(
                (_midpoint(eigenvalues[item]) for item in groups[target]), acb(0)
            ) / len(groups[target])
    return groups


def _diagnostic_block(
    *,
    digits: int,
    relative: arb,
    cutoff: arb,
    scale: arb,
    tolerance_source: str,
    audit: _ZeroAudit,
) -> dict[str, Any]:
    """统一构造浮点判别 provenance。"""

    audit_records = audit.records()
    largest_relative = None
    smallest_relative = None
    if scale > 0:
        largest_relative = (audit.largest_zeroed / scale).str(radius=False)
        if audit.smallest_retained is not None:
            smallest_relative = (audit.smallest_retained / scale).str(radius=False)

    return {
        "classification": "numerical_threshold",
        "precision_digits": digits,
        "relative_zero_tolerance": relative.str(radius=False),
        "absolute_zero_cutoff": cutoff.str(radius=False),
        "matrix_scale": scale.str(radius=False),
        "matrix_scale_is_zero": bool(scale == 0),
        "tolerance_source": tolerance_source,
        **audit_records,
        "largest_zeroed_relative_to_matrix_scale": largest_relative,
        "smallest_retained_relative_to_matrix_scale": smallest_relative,
    }


def build_numerical_frobenius_manifest(
    system: NumericalRegularSingularSystem,
    options: NumericalFrobeniusOptions | None = None,
) -> dict[str, Any]:
    """按输入精度斩杀线构造浮点 Jordan/projector/resonance manifest。

    该结果是阈值依赖的数值判别，不等同于 exact 证明。无法稳定得到完整局部基时抛出
    ``ValueError``，调用者必须提高精度或改用 exact rational 输入。
    """

    options = options or NumericalFrobeniusOptions()
    residue_raw, regular_raw = system.acb_data()
    dimension = residue_raw.nrows()
    digits, relative, tolerance_source = options.resolve(
        dimension,
        system.input_precision(),
    )
    scale = _entry_scale([residue_raw, *regular_raw])
    reference_scale = scale if scale > 0 else arb(1)
    cutoff = relative * reference_scale
    audit = _ZeroAudit(cutoff)
    residue = _threshold_matrix(residue_raw, audit)
    regular = [_threshold_matrix(matrix, audit) for matrix in regular_raw]
    identity = _identity(dimension)

    eigenvalues, _ = residue.eig(right=True, algorithm="approx")
    if any("nan" in str(value).lower() for value in eigenvalues):
        raise ValueError(
            f"{system.name}: numerical indicial spectrum is uncertain; increase precision or use exact input"
        )
    groups = _cluster_eigenvalues(eigenvalues, audit)
    roots = [
        sum((_midpoint(eigenvalues[index]) for index in group), acb(0)) / len(group)
        for group in groups
    ]
    ordering = sorted(
        range(len(roots)),
        key=lambda position: (
            float(roots[position].real.mid()),
            float(roots[position].imag.mid()),
        ),
    )
    roots = [roots[position] for position in ordering]
    groups = [groups[position] for position in ordering]
    nullspaces = [_nullspace(residue - identity * root, audit) for root in roots]
    multiplicities = [len(group) for group in groups]
    geometric_multiplicities = [len(vectors) for vectors in nullspaces]
    common = {
        "schema": "flintnde_numerical_frobenius_manifest_v1",
        "status": "passed",
        "system_name": system.name,
        "dimension": dimension,
        "roots_numeric": [_acb_record(root) for root in roots],
        "root_multiplicities": multiplicities,
        "geometric_multiplicities": geometric_multiplicities,
    }

    if len(roots) == 1:
        root = _midpoint(residue.trace() / dimension)
        nilpotent = residue - identity * root
        nilpotency_index = None
        power = acb_mat(identity)
        for degree in range(1, dimension + 1):
            power = power * nilpotent
            if _matrix_is_zero(power, audit):
                nilpotency_index = degree
                break
        if nilpotency_index is None:
            raise ValueError(
                f"{system.name}: repeated numerical root did not pass the nilpotency threshold; use exact input"
            )
        diagnostics = _diagnostic_block(
            digits=digits,
            relative=relative,
            cutoff=cutoff,
            scale=scale,
            tolerance_source=tolerance_source,
            audit=audit,
        )
        _warn_threshold_classification(system.name, diagnostics)
        return {
            **common,
            "route": "single_root_jordan_numeric_gate",
            "root_numeric": _acb_record(root),
            "nilpotent_numeric": _matrix_records(nilpotent),
            "maximum_log_degree": nilpotency_index - 1,
            "projectors_numeric": None,
            "solution_roots_numeric": None,
            "initial_vectors_numeric": None,
            "resonance_gates": {},
            "resonance_targets": {},
            "gate_series_order": 0,
            "numeric_diagnostics": diagnostics,
        }

    if sum(geometric_multiplicities) != dimension:
        raise ValueError(
            f"{system.name}: mixed-root numerical defective structure is uncertain; use exact input"
        )
    projectors: list[acb_mat] = []
    for root in roots:
        projector = acb_mat(identity)
        for other in roots:
            if not audit.is_zero(root - other):
                projector = projector * (residue - identity * other) / (root - other)
        projectors.append(projector)
    completeness = acb_mat(dimension, dimension)
    for projector in projectors:
        completeness += projector
    if not _matrix_is_zero(completeness - identity, audit):
        raise ValueError(
            f"{system.name}: numerical spectral projectors failed the configured threshold"
        )

    integer_differences: list[int] = []
    resonance_by_power: dict[tuple[int, int], int] = {}
    for left_position, left in enumerate(roots):
        for right_position, right in enumerate(roots):
            difference = left - right
            nearest = round(float(difference.real.mid()))
            if nearest > 0 and audit.is_zero(difference - nearest):
                integer_differences.append(nearest)
                resonance_by_power[(right_position, nearest)] = left_position
    gate_order = max(integer_differences, default=0)
    regular_gate = list(regular[:gate_order])
    regular_gate.extend(
        acb_mat(dimension, dimension) for _ in range(gate_order - len(regular_gate))
    )
    solution_roots: list[acb] = []
    initial_vectors: list[acb_mat] = []
    solution_root_positions: list[int] = []
    for root_position, (root, vectors) in enumerate(zip(roots, nullspaces)):
        for vector in vectors:
            solution_roots.append(root)
            initial_vectors.append(vector)
            solution_root_positions.append(root_position)
    gates: dict[str, bool] = {}
    resonance_targets: dict[str, int] = {}
    actual_maximum_log_degree = 0

    for solution_index, (root, root_position, initial) in enumerate(
        zip(solution_roots, solution_root_positions, initial_vectors), start=1
    ):
        series: list[list[acb_mat]] = [[initial]]
        for degree in range(1, gate_order + 1):
            previous_log_degree = max(len(item) for item in series) - 1
            right_by_log = [acb_mat(dimension, 1) for _ in range(previous_log_degree + 1)]
            for regular_degree in range(degree):
                earlier = series[degree - 1 - regular_degree]
                for log_degree, coefficient in enumerate(earlier):
                    right_by_log[log_degree] += regular_gate[regular_degree] * coefficient
            absolute_power = root + degree
            target_position = resonance_by_power.get((root_position, degree))
            if target_position is None:
                inverse = acb_mat(dimension, dimension)
                for projector, spectral_root in zip(projectors, roots):
                    inverse += projector / (absolute_power - spectral_root)
                coefficients = [
                    acb_mat(dimension, 1) for _ in range(previous_log_degree + 1)
                ]
                for log_degree in range(previous_log_degree, -1, -1):
                    corrected = acb_mat(right_by_log[log_degree])
                    if log_degree < previous_log_degree:
                        corrected -= coefficients[log_degree + 1] * (log_degree + 1)
                    coefficients[log_degree] = inverse * corrected
            else:
                reduced_inverse = acb_mat(dimension, dimension)
                for position, (projector, spectral_root) in enumerate(zip(projectors, roots)):
                    if position != target_position:
                        reduced_inverse += projector / (absolute_power - spectral_root)
                coefficients = [
                    acb_mat(dimension, 1) for _ in range(previous_log_degree + 2)
                ]
                retained_log_degree = previous_log_degree
                for log_degree in range(previous_log_degree, -1, -1):
                    corrected = (
                        right_by_log[log_degree]
                        - coefficients[log_degree + 1] * (log_degree + 1)
                    )
                    defect = projectors[target_position] * corrected
                    active = not _matrix_is_zero(defect, audit)
                    gate_key = f"{solution_index}:{degree}:{log_degree}"
                    gates[gate_key] = active
                    resonance_targets[gate_key] = target_position
                    if active:
                        coefficients[log_degree + 1] += defect / (log_degree + 1)
                        corrected -= defect
                        retained_log_degree = max(retained_log_degree, log_degree + 1)
                    coefficients[log_degree] = reduced_inverse * corrected
                coefficients = coefficients[: retained_log_degree + 1]
            actual_maximum_log_degree = max(
                actual_maximum_log_degree,
                len(coefficients) - 1,
            )
            series.append(coefficients)

    diagnostics = _diagnostic_block(
        digits=digits,
        relative=relative,
        cutoff=cutoff,
        scale=scale,
        tolerance_source=tolerance_source,
        audit=audit,
    )
    _warn_threshold_classification(system.name, diagnostics)
    return {
        **common,
        "route": "diagonalizable_roots_numeric_gate",
        "root_numeric": None,
        "nilpotent_numeric": None,
        "maximum_log_degree": actual_maximum_log_degree,
        "projectors_numeric": [_matrix_records(projector) for projector in projectors],
        "solution_roots_numeric": [_acb_record(root) for root in solution_roots],
        "solution_root_positions": solution_root_positions,
        "initial_vectors_numeric": [_vector_records(vector) for vector in initial_vectors],
        "resonance_gates": gates,
        "resonance_targets": resonance_targets,
        "gate_series_order": gate_order,
        "numeric_diagnostics": diagnostics,
    }


def acb_from_numeric_record(record: dict[str, str]) -> acb:
    """从数值 manifest 的实部、虚部字符串恢复 Acb。"""

    return acb(arb(record["real"]), arb(record["imag"]))


def acb_matrix_from_numeric_records(records: list[list[dict[str, str]]]) -> acb_mat:
    """从数值 manifest 恢复 Acb 矩阵。"""

    return acb_mat([[acb_from_numeric_record(value) for value in row] for row in records])
