"""正则奇点的 matrix-indicial 分析和 generalized power-log 局部基。

本模块用 FLINT exact Q(i) 矩阵决定 indicial roots、Jordan nilpotency 和整数差共振 gate，
再用 FLINT/Acb 递推数值局部基。重根本身不会自动制造 log；只有非平凡 Jordan 链
或 exact 共振相容性缺陷会提高 log 次数。
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any, Callable

from flint import acb, acb_mat

from .core import identity_matrix
from .exact_gaussian import (
    GaussianMatrix,
    GaussianRational,
    canonical_exact_solve,
    exact_column,
    exact_matrix_from_columns,
    gaussian_rational,
)
from .numeric_structure import (
    NumericalFrobeniusOptions,
    NumericalRegularSingularSystem,
    acb_from_numeric_record,
    acb_matrix_from_numeric_records,
    build_numerical_frobenius_manifest,
)


class IndicialSpectrumError(ValueError):
    """表示 exact 谱门禁失败：indicium 谱不在 Q(i) 上分裂或声称 roots 无法 exact 验证。

    小数化输入（根式参数、非分裂特征多项式）经此异常把调用方显式引导到
    ``NumericalRegularSingularSystem`` 的数值 Frobenius 路线；本模块自身
    绝不暗中降级。
    """


def _exact_matrix(records: list[list[Any]]) -> GaussianMatrix:
    """把二维 exact 记录恢复为原维 FLINT Q(i) 矩阵。"""

    try:
        return GaussianMatrix.from_records(records)
    except TypeError as error:
        raise TypeError(
            "exact Frobenius gate received floating input; use NumericalRegularSingularSystem"
        ) from error


def _matrix_records(matrix: GaussianMatrix) -> list[list[str]]:
    """把 FLINT 有理矩阵保存为可移植 exact 字符串。"""

    return matrix.to_records()


def _vector_records(vector: GaussianMatrix) -> list[str]:
    """把 FLINT 有理列向量保存为 exact 字符串。"""

    return [str(vector.scalar(row, 0)) for row in range(vector.nrows)]


def _exact_zero_vector(vector: GaussianMatrix) -> bool:
    """严格判断一个 exact 向量是否为零。"""

    return vector.is_zero


def _exact_identity(dimension: int) -> GaussianMatrix:
    """构造 FLINT exact 有理单位矩阵。"""

    return GaussianMatrix.identity(dimension)


def _exact_nullspace(matrix: GaussianMatrix) -> list[GaussianMatrix]:
    """由 FLINT RREF 构造 exact 右零空间基。"""

    return matrix.nullspace()


@dataclass(frozen=True)
class RegularSingularSystem:
    """保存 ``Y'=(R/z+sum A_m z^m)Y`` 的 exact 局部数据。"""

    residue_exact: list[list[Any]]
    regular_coefficients_exact: tuple[list[list[Any]], ...]
    name: str = "regular-singular-system"

    @property
    def dimension(self) -> int:
        """返回矩阵系统维数。"""

        return len(self.residue_exact)

    def exact_data(self) -> tuple[GaussianMatrix, list[GaussianMatrix]]:
        """返回 exact residue 和解析部分系数。"""

        residue = _exact_matrix(self.residue_exact)
        if residue.nrows != residue.ncols:
            raise ValueError(f"{self.name}: residue matrix must be square")
        regular = [_exact_matrix(records) for records in self.regular_coefficients_exact]
        if any(
            matrix.nrows != residue.nrows or matrix.ncols != residue.ncols
            for matrix in regular
        ):
            raise ValueError(f"{self.name}: regular coefficient dimension mismatch")
        return residue, regular

    def acb_data(self, order: int) -> tuple[acb_mat, list[acb_mat]]:
        """返回 Acb residue 和补零到指定阶数的解析部分系数。"""

        if order <= 0:
            raise ValueError("power-log series order must be positive")
        residue = _exact_matrix(self.residue_exact).to_acb()
        regular = [_exact_matrix(records).to_acb() for records in self.regular_coefficients_exact]
        regular.extend(acb_mat(self.dimension, self.dimension) for _ in range(order - len(regular)))
        return residue, regular[:order]


def _manifest_common(
    system: RegularSingularSystem,
    dimension: int,
    roots: list[GaussianRational],
    root_multiplicities: list[int],
    geometric_multiplicities: list[int],
    spectrum_source: str | None,
) -> dict[str, Any]:
    """拼出两种 exact 谱来源共用的 manifest 公共字段。"""

    common = {
        "schema": "flintnde_exact_frobenius_manifest_v1",
        "status": "passed",
        "system_name": system.name,
        "dimension": dimension,
        "exact_field": "Q(i)",
        "roots_exact": [str(root) for root in roots],
        "root_multiplicities": root_multiplicities,
        "geometric_multiplicities": geometric_multiplicities,
    }
    if spectrum_source is not None:
        common["spectrum_source"] = spectrum_source
    return common


def _single_root_manifest(
    system: RegularSingularSystem,
    dimension: int,
    root: GaussianRational,
    nilpotent: GaussianMatrix,
    nilpotency_index: int,
    spectrum_source: str | None,
) -> dict[str, Any]:
    """单根任意 Jordan 结构的 manifest；n=0 对数阶 = nilpotency index - 1。"""

    return {
        **_manifest_common(system, dimension, [root], [dimension], [dimension], spectrum_source),
        "route": "single_root_jordan_exact_gate",
        "root_exact": str(root),
        "nilpotent_exact": _matrix_records(nilpotent),
        "maximum_log_degree": nilpotency_index - 1,
        "projectors_exact": None,
        "solution_roots_exact": None,
        "initial_vectors_exact": None,
        "resonance_gates": {},
        "gate_series_order": 0,
    }


def _diagonalizable_manifest(
    system: RegularSingularSystem,
    dimension: int,
    residue: GaussianMatrix,
    regular: list[GaussianMatrix],
    roots: list[GaussianRational],
    root_multiplicities: list[int],
    geometric_multiplicities: list[int],
    nullspaces: list[list[GaussianMatrix]],
    spectrum_source: str | None,
) -> dict[str, Any]:
    """可对角化 residue 的 Lagrange projector 与整数差共振 gate manifest。"""

    identity = _exact_identity(dimension)
    if sum(geometric_multiplicities) != dimension:
        raise ValueError(
            f"{system.name}: mixed-root defective residue needs explicit analytic Jordan-chain data"
        )
    projectors: list[GaussianMatrix] = []
    for root in roots:
        projector = identity
        for other in roots:
            if other != root:
                projector = projector * (residue - identity * other) / (root - other)
        projectors.append(projector)
    projector_sum = GaussianMatrix.zero(dimension, dimension)
    for projector in projectors:
        projector_sum += projector
    if not (projector_sum - identity).is_zero:
        raise ValueError(f"{system.name}: exact spectral projectors failed completeness")

    positive_integer_differences = [
        int((left - right).real)
        for left in roots
        for right in roots
        if (left - right).imag == 0
        and (left - right).real.q == 1
        and (left - right).real > 0
    ]
    gate_order = max(positive_integer_differences, default=0)
    regular_gate = list(regular[:gate_order])
    regular_gate.extend(
        GaussianMatrix.zero(dimension, dimension)
        for _ in range(gate_order - len(regular_gate))
    )
    solution_roots: list[GaussianRational] = []
    initial_vectors: list[GaussianMatrix] = []
    for root, vectors in zip(roots, nullspaces):
        for vector in vectors:
            solution_roots.append(root)
            initial_vectors.append(vector)
    gates: dict[str, bool] = {}
    actual_maximum_log_degree = 0

    for solution_index, (root, initial) in enumerate(zip(solution_roots, initial_vectors), start=1):
        series: list[list[GaussianMatrix]] = [[initial]]
        for degree in range(1, gate_order + 1):
            previous_log_degree = max(len(item) for item in series) - 1
            right_by_log = [
                GaussianMatrix.zero(dimension, 1)
                for _ in range(previous_log_degree + 1)
            ]
            for regular_degree in range(degree):
                earlier = series[degree - 1 - regular_degree]
                for log_degree, coefficient in enumerate(earlier):
                    right_by_log[log_degree] += regular_gate[regular_degree] * coefficient
            absolute_power = root + degree
            if absolute_power not in roots:
                inverse = GaussianMatrix.zero(dimension, dimension)
                for projector, spectral_root in zip(projectors, roots):
                    inverse += projector / (absolute_power - spectral_root)
                coefficients = [
                    GaussianMatrix.zero(dimension, 1)
                    for _ in range(previous_log_degree + 1)
                ]
                for log_degree in range(previous_log_degree, -1, -1):
                    corrected = right_by_log[log_degree]
                    if log_degree < previous_log_degree:
                        corrected -= (log_degree + 1) * coefficients[log_degree + 1]
                    coefficients[log_degree] = inverse * corrected
            else:
                resonant_position = roots.index(absolute_power)
                reduced_inverse = GaussianMatrix.zero(dimension, dimension)
                for projector, spectral_root in zip(projectors, roots):
                    if spectral_root != absolute_power:
                        reduced_inverse += projector / (absolute_power - spectral_root)
                coefficients = [
                    GaussianMatrix.zero(dimension, 1)
                    for _ in range(previous_log_degree + 2)
                ]
                retained_log_degree = previous_log_degree
                for log_degree in range(previous_log_degree, -1, -1):
                    corrected = (
                        right_by_log[log_degree]
                        - (log_degree + 1) * coefficients[log_degree + 1]
                    )
                    defect = projectors[resonant_position] * corrected
                    active = not _exact_zero_vector(defect)
                    gates[f"{solution_index}:{degree}:{log_degree}"] = active
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

    return {
        **_manifest_common(
            system, dimension, roots, root_multiplicities,
            geometric_multiplicities, spectrum_source,
        ),
        "route": "diagonalizable_roots_exact_gate",
        "root_exact": None,
        "nilpotent_exact": None,
        "maximum_log_degree": actual_maximum_log_degree,
        "projectors_exact": [_matrix_records(projector) for projector in projectors],
        "solution_roots_exact": [str(root) for root in solution_roots],
        "initial_vectors_exact": [_vector_records(vector) for vector in initial_vectors],
        "resonance_gates": gates,
        "gate_series_order": gate_order,
    }


def build_exact_frobenius_manifest(
    system: RegularSingularSystem,
) -> dict[str, Any]:
    """由 exact 局部矩阵构造 Jordan/projector/resonance manifest。"""

    residue, regular = system.exact_data()
    dimension = residue.nrows
    identity = _exact_identity(dimension)
    try:
        root_records = residue.charpoly().gaussian_rational_roots()
    except ValueError as error:
        raise IndicialSpectrumError(
            f"{system.name}: exact Q(i) gate currently requires an indicial spectrum "
            "that splits completely over Q(i)"
        ) from error
    roots = [root for root, _ in root_records]
    multiplicity_by_root = {root: int(multiplicity) for root, multiplicity in root_records}
    root_multiplicities = [multiplicity_by_root[root] for root in roots]
    nullspaces = [_exact_nullspace(residue - identity * root) for root in roots]
    geometric_multiplicities = [len(vectors) for vectors in nullspaces]

    if len(roots) == 1:
        root = roots[0]
        nilpotent = residue - identity * root
        nilpotency_index = next(
            (power for power in range(1, dimension + 1) if (nilpotent**power).is_zero),
            None,
        )
        if nilpotency_index is None:
            raise ValueError(f"{system.name}: single-root residue failed the nilpotency gate")
        return _single_root_manifest(
            system, dimension, root, nilpotent, nilpotency_index, None,
        )
    return _diagonalizable_manifest(
        system, dimension, residue, regular, roots, root_multiplicities,
        geometric_multiplicities, nullspaces, None,
    )


def _ascending_kernel_chain(nilpotent: GaussianMatrix) -> tuple[list[list[GaussianMatrix]], int]:
    """升幂零空间链 K_k = NullSpace(N^k) 及其稳定幂；不解特征多项式。"""

    kernels = [nilpotent.nullspace()]
    index = 1
    while len(kernels[-1]) < nilpotent.nrows:
        extended = (nilpotent ** (index + 1)).nullspace()
        if len(extended) == len(kernels[-1]):
            break
        kernels.append(extended)
        index += 1
    return kernels, index


def build_known_roots_frobenius_manifest(
    system: RegularSingularSystem,
    known_roots: Sequence[Any],
) -> dict[str, Any]:
    """用已知 indicial roots（含重复，重复次数即代数重数）构造 manifest。

    只做验证性线性代数：每个声称 root 的广义本征子空间维数必须等于其声称
    重数，且重数总和等于系统维数。二者合起来唯一确定
    ``det(lam*I-R) = prod (lam-alpha)^m_alpha``，因此无需 charpoly 或求根。
    """

    residue, regular = system.exact_data()
    dimension = residue.nrows
    identity = _exact_identity(dimension)
    roots: list[GaussianRational] = []
    multiplicity: dict[GaussianRational, int] = {}
    for value in known_roots:
        root = gaussian_rational(value)
        if root not in multiplicity:
            roots.append(root)
            multiplicity[root] = 0
        multiplicity[root] += 1
    root_multiplicities = [multiplicity[root] for root in roots]
    if sum(root_multiplicities) != dimension:
        raise IndicialSpectrumError(
            f"{system.name}: known indicial roots claim total multiplicity "
            f"{sum(root_multiplicities)} but the system dimension is {dimension}"
        )
    kernels_by_root: dict[GaussianRational, list[list[GaussianMatrix]]] = {}
    nilpotency_by_root: dict[GaussianRational, int] = {}
    nullspaces: list[list[GaussianMatrix]] = []
    geometric_multiplicities: list[int] = []
    for root in roots:
        nilpotent = residue - identity * root
        kernels, index = _ascending_kernel_chain(nilpotent)
        generalized = len(kernels[-1])
        if generalized != multiplicity[root]:
            raise IndicialSpectrumError(
                f"{system.name}: known root {root} claims multiplicity "
                f"{multiplicity[root]} but its exact generalized eigenspace has "
                f"dimension {generalized}"
            )
        kernels_by_root[root] = kernels
        nilpotency_by_root[root] = index
        nullspaces.append(kernels[0])
        geometric_multiplicities.append(len(kernels[0]))
    if len(roots) == 1:
        root = roots[0]
        return _single_root_manifest(
            system, dimension, root, residue - identity * root,
            nilpotency_by_root[root], "known_roots",
        )
    if geometric_multiplicities == root_multiplicities:
        return _diagonalizable_manifest(
            system, dimension, residue, regular, roots, root_multiplicities,
            geometric_multiplicities, nullspaces, "known_roots",
        )
    return _known_roots_defective_manifest(
        system, dimension, residue, roots, root_multiplicities,
        geometric_multiplicities, kernels_by_root, nilpotency_by_root,
    )


def _known_roots_defective_manifest(
    system: RegularSingularSystem,
    dimension: int,
    residue: GaussianMatrix,
    roots: list[GaussianRational],
    root_multiplicities: list[int],
    geometric_multiplicities: list[int],
    kernels_by_root: dict[GaussianRational, list[list[GaussianMatrix]]],
    nilpotency_by_root: dict[GaussianRational, int],
) -> dict[str, Any]:
    """已知谱、部分 root 有 Jordan 缺陷时的 exact manifest。

    spectral projector 由广义本征基拼成的相似变换直接构造，不做多项式插值。
    canonical 解按每个广义本征基向量给出，其 n=0 领头 log 多项式为
    ``exp(N_alpha log z) w``；log 的实际首现阶不在此预估，由 exact
    power-log 递推按需给出。
    """

    columns: list[GaussianMatrix] = []
    for root in roots:
        columns.extend(kernels_by_root[root][-1])
    similarity = exact_matrix_from_columns(columns, dimension)
    inverse = similarity.inverse()
    identity = _exact_identity(dimension)
    projectors: list[GaussianMatrix] = []
    nilpotents: list[GaussianMatrix] = []
    offset = 0
    for root in roots:
        count = len(kernels_by_root[root][-1])
        block = GaussianMatrix.zero(dimension, dimension)
        for row in range(offset, offset + count):
            for column in range(dimension):
                block.real[row, column] = inverse.real[row, column]
                block.imag[row, column] = inverse.imag[row, column]
        projector = similarity * block
        if not (projector * projector - projector).is_zero:
            raise ArithmeticError(
                f"{system.name}: known-roots spectral projector for root {root} "
                "is not idempotent"
            )
        projectors.append(projector)
        nilpotents.append((residue - identity * root) * projector)
        offset += count
    projector_sum = GaussianMatrix.zero(dimension, dimension)
    for projector in projectors:
        projector_sum += projector
    if not (projector_sum - identity).is_zero:
        raise ArithmeticError(
            f"{system.name}: known-roots spectral projectors failed completeness"
        )
    solution_roots: list[str] = []
    initial_vectors: list[list[str]] = []
    for root in roots:
        for vector in kernels_by_root[root][-1]:
            solution_roots.append(str(root))
            initial_vectors.append(_vector_records(vector))
    return {
        **_manifest_common(
            system, dimension, roots, root_multiplicities,
            geometric_multiplicities, "known_roots",
        ),
        "route": "known_roots_defective_exact_gate",
        "root_exact": None,
        "nilpotent_exact": None,
        "maximum_log_degree": max(nilpotency_by_root[root] for root in roots) - 1,
        "projectors_exact": [_matrix_records(projector) for projector in projectors],
        "nilpotents_exact": [_matrix_records(nilpotent) for nilpotent in nilpotents],
        "nilpotency_indices": [nilpotency_by_root[root] for root in roots],
        "solution_roots_exact": solution_roots,
        "initial_vectors_exact": initial_vectors,
        "resonance_gates": {},
        "gate_series_order": 0,
    }


def build_frobenius_manifest(
    system: RegularSingularSystem | NumericalRegularSingularSystem,
    options: NumericalFrobeniusOptions | None = None,
    *,
    known_roots: Sequence[Any] | None = None,
) -> dict[str, Any]:
    """按输入系统类型选择结构判别；exact 失败时绝不自动退到阈值路线。

    ``known_roots`` 提供已知 indicial roots（exact 输入专用）时走验证式
    manifest，不解特征多项式；不提供则维持原有 exact/numeric 判别。
    """

    if isinstance(system, RegularSingularSystem):
        if options is not None:
            raise TypeError("exact Frobenius input does not accept numerical options")
        if known_roots is not None:
            return build_known_roots_frobenius_manifest(system, known_roots)
        return build_exact_frobenius_manifest(system)
    if isinstance(system, NumericalRegularSingularSystem):
        if known_roots is not None:
            raise TypeError("known indicial roots require exact RegularSingularSystem input")
        return build_numerical_frobenius_manifest(system, options)
    raise TypeError(
        "system must be RegularSingularSystem for exact input or "
        "NumericalRegularSingularSystem for floating input"
    )


@dataclass(frozen=True)
class PowerLogBasis:
    """保存可在指定分支上求值的局部基本矩阵。"""

    dimension: int
    maximum_log_degree: int
    manifest: dict[str, Any]
    _evaluator: Callable[[acb], acb_mat]

    def evaluate(self, point: acb) -> acb_mat:
        """在非零局部变量处计算 fundamental matrix。"""

        if abs(point).contains(0):
            raise ZeroDivisionError("power-log basis cannot be evaluated at the singular point")
        return self._evaluator(point)


@dataclass(frozen=True)
class ExactPowerLogSeries:
    """保存每个 canonical solution 的 exact ``z^(rho+n) log(z)^b`` 系数。"""

    roots: tuple[GaussianRational, ...]
    # solution -> degree -> log degree -> exact column vector
    coefficients: tuple[tuple[tuple[GaussianMatrix, ...], ...], ...]


def build_exact_power_log_series(
    system: RegularSingularSystem,
    manifest: dict[str, Any],
    series_order: int,
) -> ExactPowerLogSeries:
    """按已认证 manifest 重算 canonical power-log 基的 exact 系数。

    数值求值仍走 Acb；本接口专供一般 gauge 变换后的原基边界认证，避免浮点阈值
    参与 ``{a,b,C}`` 合法性判断。
    """

    if series_order <= 0:
        raise ValueError("exact power-log series order must be positive")
    if manifest.get("status") != "passed" or manifest.get("system_name") != system.name:
        raise ValueError(
            "the supplied exact power-log manifest does not belong to this system or "
            "has not passed validation; regenerate the manifest for the current "
            "system before solving"
        )
    _residue, regular = system.exact_data()
    dimension = system.dimension
    regular = list(regular[:series_order])
    regular.extend(
        GaussianMatrix.zero(dimension, dimension)
        for _ in range(series_order - len(regular))
    )
    identity = GaussianMatrix.identity(dimension)
    route = manifest.get("route")

    if route == "single_root_jordan_exact_gate":
        root = gaussian_rational(manifest["root_exact"])
        nilpotent = GaussianMatrix.from_records(manifest["nilpotent_exact"])
        maximum_log_degree = int(manifest["maximum_log_degree"])
        by_log: list[list[GaussianMatrix]] = [
            [nilpotent**log_degree / math.factorial(log_degree)]
            for log_degree in range(maximum_log_degree + 1)
        ]
        for degree in range(1, series_order + 1):
            inverse = (identity * degree - nilpotent).inverse()
            degree_coefficients = [
                GaussianMatrix.zero(dimension, dimension)
                for _ in range(maximum_log_degree + 1)
            ]
            for log_degree in range(maximum_log_degree, -1, -1):
                right = GaussianMatrix.zero(dimension, dimension)
                for regular_degree in range(degree):
                    right += regular[regular_degree] * by_log[log_degree][degree - 1 - regular_degree]
                if log_degree < maximum_log_degree:
                    right -= degree_coefficients[log_degree + 1] * (log_degree + 1)
                degree_coefficients[log_degree] = inverse * right
            for log_degree, coefficient in enumerate(degree_coefficients):
                by_log[log_degree].append(coefficient)
        solutions: list[tuple[tuple[GaussianMatrix, ...], ...]] = []
        for column in range(dimension):
            solutions.append(
                tuple(
                    tuple(
                        by_log[log_degree][degree].column(column)
                        for log_degree in range(maximum_log_degree + 1)
                    )
                    for degree in range(series_order + 1)
                )
            )
        return ExactPowerLogSeries(tuple(root for _ in range(dimension)), tuple(solutions))

    if route == "known_roots_defective_exact_gate":
        return _defective_exact_series(system, manifest, series_order)[0]

    if route != "diagonalizable_roots_exact_gate":
        raise ValueError(
            f"the exact power-log manifest declares route '{route}', which this "
            f"solver does not recognize; regenerate the manifest with a supported "
            f"route (a single Jordan root or diagonalizable roots)"
        )
    roots = [gaussian_rational(value) for value in manifest["roots_exact"]]
    projectors = {
        root: GaussianMatrix.from_records(records)
        for root, records in zip(roots, manifest["projectors_exact"])
    }
    solution_roots = [
        gaussian_rational(value) for value in manifest["solution_roots_exact"]
    ]
    initial_vectors = [
        GaussianMatrix.from_records([[gaussian_rational(value)] for value in records])
        for records in manifest["initial_vectors_exact"]
    ]
    gates = {str(key): bool(value) for key, value in manifest["resonance_gates"].items()}
    all_series: list[tuple[tuple[GaussianMatrix, ...], ...]] = []
    for solution_index, (root, initial) in enumerate(
        zip(solution_roots, initial_vectors), start=1
    ):
        series: list[list[GaussianMatrix]] = [[initial]]
        for degree in range(1, series_order + 1):
            previous_log_degree = max(len(item) for item in series) - 1
            right_by_log = [
                GaussianMatrix.zero(dimension, 1)
                for _ in range(previous_log_degree + 1)
            ]
            for regular_degree in range(degree):
                earlier = series[degree - 1 - regular_degree]
                for log_degree, coefficient in enumerate(earlier):
                    right_by_log[log_degree] += regular[regular_degree] * coefficient
            absolute_power = root + degree
            if absolute_power not in projectors:
                inverse = GaussianMatrix.zero(dimension, dimension)
                for spectral_root, projector in projectors.items():
                    inverse += projector / (absolute_power - spectral_root)
                coefficients = [
                    GaussianMatrix.zero(dimension, 1)
                    for _ in range(previous_log_degree + 1)
                ]
                for log_degree in range(previous_log_degree, -1, -1):
                    corrected = right_by_log[log_degree]
                    if log_degree < previous_log_degree:
                        corrected -= coefficients[log_degree + 1] * (log_degree + 1)
                    coefficients[log_degree] = inverse * corrected
            else:
                resonant_projector = projectors[absolute_power]
                reduced_inverse = GaussianMatrix.zero(dimension, dimension)
                for spectral_root, projector in projectors.items():
                    if spectral_root != absolute_power:
                        reduced_inverse += projector / (absolute_power - spectral_root)
                coefficients = [
                    GaussianMatrix.zero(dimension, 1)
                    for _ in range(previous_log_degree + 2)
                ]
                for log_degree in range(previous_log_degree, -1, -1):
                    corrected = (
                        right_by_log[log_degree]
                        - coefficients[log_degree + 1] * (log_degree + 1)
                    )
                    gate_key = f"{solution_index}:{degree}:{log_degree}"
                    if gates.get(gate_key, False):
                        defect = resonant_projector * corrected
                        coefficients[log_degree + 1] += defect / (log_degree + 1)
                        corrected -= defect
                    coefficients[log_degree] = reduced_inverse * corrected
                while len(coefficients) > 1 and coefficients[-1].is_zero:
                    coefficients.pop()
            series.append(coefficients)
        all_series.append(tuple(tuple(item) for item in series))
    return ExactPowerLogSeries(tuple(solution_roots), tuple(all_series))


class _ResonanceLogRaise(ValueError):
    """共振阶在当前 log 次数下不相容；提高一档重解的内部信号。"""


def _stack_columns(vectors: list[GaussianMatrix]) -> GaussianMatrix:
    """把若干同维列向量按顺序竖向叠成一个长列向量。"""

    return GaussianMatrix.from_records(
        [[vector.scalar(row, 0)] for vector in vectors for row in range(vector.nrows)]
    )


def _resonance_block_matrix(
    matrix: GaussianMatrix, identity: GaussianMatrix, size: int
) -> GaussianMatrix:
    """联立 ``A c_b + (b+1) c_{b+1} = rhs_b``（b=0..size-1, c_size=0）的块矩阵。"""

    dimension = matrix.nrows

    def entry(block_row: int, block_column: int, row: int, column: int) -> GaussianRational:
        if block_column == block_row:
            return matrix.scalar(row, column)
        if block_column == block_row + 1:
            return identity.scalar(row, column) * (block_row + 1)
        return GaussianRational()

    return GaussianMatrix.from_records(
        [
            [
                entry(block_row, block_column, row, column)
                for block_column in range(size)
                for column in range(dimension)
            ]
            for block_row in range(size)
            for row in range(dimension)
        ]
    )


def _defective_exact_series(
    system: RegularSingularSystem,
    manifest: dict[str, Any],
    series_order: int,
) -> tuple[ExactPowerLogSeries, dict[str, bool], int]:
    """已知谱 mixed-defective 路线的 exact power-log 递推。

    非共振阶直接用 ``(gamma*I-R)^{-1}`` 从高 log 次数向下解；gamma 命中已知
    root 时把该阶全部 log 次数联立成块方程，不相容则提高一档 log 次数重解，
    新 log 的首现阶由此自然给出，全程不求特征多项式。
    """

    residue, regular_full = system.exact_data()
    dimension = system.dimension
    regular = list(regular_full[:series_order])
    regular.extend(
        GaussianMatrix.zero(dimension, dimension)
        for _ in range(series_order - len(regular))
    )
    identity = GaussianMatrix.identity(dimension)
    roots = [gaussian_rational(value) for value in manifest["roots_exact"]]
    root_positions = {root: position for position, root in enumerate(roots)}
    nilpotents = [
        GaussianMatrix.from_records(records) for records in manifest["nilpotents_exact"]
    ]
    solution_roots = [
        gaussian_rational(value) for value in manifest["solution_roots_exact"]
    ]
    initial_vectors = [
        GaussianMatrix.from_records([[gaussian_rational(value)] for value in records])
        for records in manifest["initial_vectors_exact"]
    ]
    zero_column = GaussianMatrix.zero(dimension, 1)
    inverse_cache: dict[GaussianRational, GaussianMatrix] = {}
    gates: dict[str, bool] = {}
    maximum_log_degree = 0
    all_series: list[tuple[tuple[GaussianMatrix, ...], ...]] = []
    for solution_index, (root, initial) in enumerate(
        zip(solution_roots, initial_vectors), start=1
    ):
        nilpotent = nilpotents[root_positions[root]]
        degree_zero = [initial]
        image = initial
        power = 0
        while True:
            power += 1
            image = nilpotent * image
            if image.is_zero:
                break
            degree_zero.append(image / math.factorial(power))
        maximum_log_degree = max(maximum_log_degree, len(degree_zero) - 1)
        series: list[list[GaussianMatrix]] = [degree_zero]
        for degree in range(1, series_order + 1):
            gamma = root + degree
            right_top = max(len(series[degree - k]) - 1 for k in range(1, degree + 1))
            right = []
            for log_degree in range(right_top + 1):
                total = GaussianMatrix.zero(dimension, 1)
                for k in range(1, degree + 1):
                    earlier = series[degree - k]
                    if log_degree < len(earlier):
                        total += regular[k - 1] * earlier[log_degree]
                right.append(total)
            if gamma not in root_positions:
                inverse = inverse_cache.get(gamma)
                if inverse is None:
                    inverse = (identity * gamma - residue).inverse()
                    inverse_cache[gamma] = inverse
                coefficients: list[GaussianMatrix] = [zero_column] * (right_top + 1)
                for log_degree in range(right_top, -1, -1):
                    corrected = right[log_degree]
                    if log_degree < right_top:
                        corrected = corrected - coefficients[log_degree + 1] * (log_degree + 1)
                    coefficients[log_degree] = inverse * corrected
            else:
                matrix = identity * gamma - residue
                budget = right_top
                while True:
                    stacked = _resonance_block_matrix(matrix, identity, budget + 1)
                    rhs = [
                        right[block] if block <= right_top else zero_column
                        for block in range(budget + 1)
                    ]
                    try:
                        solution = canonical_exact_solve(
                            stacked,
                            _stack_columns(rhs),
                            incompatible_error=_ResonanceLogRaise,
                            incompatible_message="resonant log-degree budget exhausted",
                            solve_message="known-roots resonant solve failed substitution",
                        )
                    except _ResonanceLogRaise:
                        budget += 1
                        if budget > right_top + dimension:
                            raise ValueError(
                                f"{system.name}: resonant degree {degree} of solution "
                                f"{solution_index} did not admit a power-log solution "
                                f"within log degree {budget}"
                            )
                        continue
                    break
                gates[f"{solution_index}:{degree}"] = budget > right_top
                coefficients = [
                    GaussianMatrix.from_records(
                        [
                            [solution.scalar(block * dimension + row, 0)]
                            for row in range(dimension)
                        ]
                    )
                    for block in range(budget + 1)
                ]
            while len(coefficients) > 1 and coefficients[-1].is_zero:
                coefficients.pop()
            maximum_log_degree = max(maximum_log_degree, len(coefficients) - 1)
            series.append(coefficients)
        all_series.append(tuple(tuple(item) for item in series))
    return (
        ExactPowerLogSeries(tuple(solution_roots), tuple(all_series)),
        gates,
        maximum_log_degree,
    )


def _build_known_roots_defective_basis(
    system: RegularSingularSystem,
    manifest: dict[str, Any],
    series_order: int,
) -> PowerLogBasis:
    """已知谱 mixed-defective 残差的 exact power-log 数值基。"""

    exact_series, gates, maximum_log_degree = _defective_exact_series(
        system, manifest, series_order
    )
    dimension = system.dimension
    solution_roots = [value.to_acb() for value in exact_series.roots]
    numeric_series = [
        [
            [vector.to_acb() for vector in log_coefficients]
            for log_coefficients in degrees
        ]
        for degrees in exact_series.coefficients
    ]

    def evaluate(point: acb) -> acb_mat:
        logarithm = point.log()
        columns: list[acb_mat] = []
        for root, series in zip(solution_roots, numeric_series):
            value = acb_mat(dimension, 1)
            for degree_coefficients in reversed(series):
                log_polynomial = acb_mat(dimension, 1)
                for coefficient in reversed(degree_coefficients):
                    log_polynomial = log_polynomial * logarithm + coefficient
                value = value * point + log_polynomial
            columns.append(value * point**root)
        return acb_mat(
            [[columns[column][row, 0] for column in range(dimension)] for row in range(dimension)]
        )

    basis_manifest = {
        **manifest,
        "maximum_log_degree": maximum_log_degree,
        "resonance_gates": gates,
    }
    return PowerLogBasis(dimension, maximum_log_degree, basis_manifest, evaluate)


def _build_single_root_basis_from_data(
    *,
    dimension: int,
    regular: list[acb_mat],
    root: acb,
    nilpotent: acb_mat,
    manifest: dict[str, Any],
    series_order: int,
) -> PowerLogBasis:
    """由已数值化的单根 Jordan 数据构造 power-log 基。"""

    identity = identity_matrix(dimension)
    maximum_log_degree = int(manifest["maximum_log_degree"])
    coefficient_by_log: list[list[acb_mat]] = []
    nilpotent_power = acb_mat(identity)
    for log_degree in range(maximum_log_degree + 1):
        if log_degree > 0:
            nilpotent_power = nilpotent_power * nilpotent
        coefficient_by_log.append([nilpotent_power / acb(math.factorial(log_degree))])
    for degree in range(1, series_order + 1):
        operator_inverse = acb_mat(dimension, dimension)
        nilpotent_power = acb_mat(identity)
        for nilpotent_degree in range(maximum_log_degree + 1):
            operator_inverse += nilpotent_power / acb(degree ** (nilpotent_degree + 1))
            nilpotent_power = nilpotent_power * nilpotent
        degree_coefficients = [acb_mat(dimension, dimension) for _ in range(maximum_log_degree + 1)]
        for log_degree in range(maximum_log_degree, -1, -1):
            right = acb_mat(dimension, dimension)
            for regular_degree in range(degree):
                right += regular[regular_degree] * coefficient_by_log[log_degree][degree - 1 - regular_degree]
            if log_degree < maximum_log_degree:
                right -= degree_coefficients[log_degree + 1] * acb(log_degree + 1)
            degree_coefficients[log_degree] = operator_inverse * right
        for log_degree, coefficient in enumerate(degree_coefficients):
            coefficient_by_log[log_degree].append(coefficient)

    def evaluate(point: acb) -> acb_mat:
        logarithm = point.log()
        value = acb_mat(dimension, dimension)
        for log_degree, series in enumerate(coefficient_by_log):
            polynomial = acb_mat(dimension, dimension)
            for coefficient in reversed(series):
                polynomial = polynomial * point + coefficient
            value += polynomial * logarithm**log_degree
        return value * point**root

    return PowerLogBasis(dimension, maximum_log_degree, manifest, evaluate)


def _build_single_root_basis(
    system: RegularSingularSystem,
    manifest: dict[str, Any],
    series_order: int,
) -> PowerLogBasis:
    """构造 exact 单根任意 Jordan 链的 power-log 基。"""

    _, regular = system.acb_data(series_order)
    return _build_single_root_basis_from_data(
        dimension=system.dimension,
        regular=regular,
        root=gaussian_rational(manifest["root_exact"]).to_acb(),
        nilpotent=_exact_matrix(manifest["nilpotent_exact"]).to_acb(),
        manifest=manifest,
        series_order=series_order,
    )


def _build_diagonalizable_basis(
    system: RegularSingularSystem,
    manifest: dict[str, Any],
    series_order: int,
) -> PowerLogBasis:
    """构造可对角化 residue，并按 exact gate 动态加入共振 log。"""

    _, regular = system.acb_data(series_order)
    dimension = system.dimension
    roots_exact = [gaussian_rational(value) for value in manifest["roots_exact"]]
    projectors = {
        root: _exact_matrix(records).to_acb()
        for root, records in zip(roots_exact, manifest["projectors_exact"])
    }
    solution_roots_exact = [
        gaussian_rational(value) for value in manifest["solution_roots_exact"]
    ]
    solution_roots = [value.to_acb() for value in solution_roots_exact]
    initial_vectors = [
        acb_mat([[gaussian_rational(value).to_acb()] for value in records])
        for records in manifest["initial_vectors_exact"]
    ]
    gates = {str(key): bool(value) for key, value in manifest["resonance_gates"].items()}
    series_by_solution: list[list[list[acb_mat]]] = []
    actual_maximum_log_degree = 0
    for solution_index, (root_exact, initial) in enumerate(
        zip(solution_roots_exact, initial_vectors), start=1
    ):
        series: list[list[acb_mat]] = [[initial]]
        for degree in range(1, series_order + 1):
            previous_log_degree = max(len(item) for item in series) - 1
            right_by_log = [acb_mat(dimension, 1) for _ in range(previous_log_degree + 1)]
            for regular_degree in range(degree):
                earlier = series[degree - 1 - regular_degree]
                for log_degree, coefficient in enumerate(earlier):
                    right_by_log[log_degree] += regular[regular_degree] * coefficient
            absolute_power = root_exact + degree
            if absolute_power not in projectors:
                inverse = acb_mat(dimension, dimension)
                for root, projector in projectors.items():
                    inverse += projector / (absolute_power - root).to_acb()
                coefficients = [acb_mat(dimension, 1) for _ in range(previous_log_degree + 1)]
                for log_degree in range(previous_log_degree, -1, -1):
                    corrected = acb_mat(right_by_log[log_degree])
                    if log_degree < previous_log_degree:
                        corrected -= coefficients[log_degree + 1] * acb(log_degree + 1)
                    coefficients[log_degree] = inverse * corrected
            else:
                resonant_projector = projectors[absolute_power]
                reduced_inverse = acb_mat(dimension, dimension)
                for root, projector in projectors.items():
                    if root != absolute_power:
                        reduced_inverse += projector / (absolute_power - root).to_acb()
                coefficients = [acb_mat(dimension, 1) for _ in range(previous_log_degree + 2)]
                for log_degree in range(previous_log_degree, -1, -1):
                    corrected = acb_mat(right_by_log[log_degree])
                    corrected -= coefficients[log_degree + 1] * acb(log_degree + 1)
                    gate_key = f"{solution_index}:{degree}:{log_degree}"
                    if gate_key not in gates:
                        raise ValueError(
                            f"the exact Frobenius manifest is missing the resonance "
                            f"decision encoded as '{gate_key}' (solution index, "
                            f"series degree, log degree); regenerate the manifest so "
                            f"every resonance gate is recorded"
                        )
                    if gates[gate_key]:
                        defect = resonant_projector * corrected
                        coefficients[log_degree + 1] += defect / acb(log_degree + 1)
                        corrected -= defect
                    coefficients[log_degree] = reduced_inverse * corrected
                while len(coefficients) > 1 and all(
                    abs(coefficients[-1][row, 0]).contains(0) for row in range(dimension)
                ):
                    coefficients.pop()
            actual_maximum_log_degree = max(actual_maximum_log_degree, len(coefficients) - 1)
            series.append(coefficients)
        series_by_solution.append(series)

    def evaluate(point: acb) -> acb_mat:
        logarithm = point.log()
        columns: list[acb_mat] = []
        for root, series in zip(solution_roots, series_by_solution):
            value = acb_mat(dimension, 1)
            for degree_coefficients in reversed(series):
                log_polynomial = acb_mat(dimension, 1)
                for coefficient in reversed(degree_coefficients):
                    log_polynomial = log_polynomial * logarithm + coefficient
                value = value * point + log_polynomial
            columns.append(value * point**root)
        return acb_mat(
            [[columns[column][row, 0] for column in range(dimension)] for row in range(dimension)]
        )

    return PowerLogBasis(dimension, actual_maximum_log_degree, manifest, evaluate)


def _build_numeric_diagonalizable_basis(
    system: NumericalRegularSingularSystem,
    manifest: dict[str, Any],
    series_order: int,
) -> PowerLogBasis:
    """按浮点 manifest 构造可对角化 residue 的 threshold-gated power-log 基。"""

    _, regular = system.acb_data()
    dimension = system.dimension
    regular = [_matrix for _matrix in regular[:series_order]]
    regular.extend(acb_mat(dimension, dimension) for _ in range(series_order - len(regular)))
    roots = [acb_from_numeric_record(record) for record in manifest["roots_numeric"]]
    projectors = [
        acb_matrix_from_numeric_records(records)
        for records in manifest["projectors_numeric"]
    ]
    solution_roots = [
        acb_from_numeric_record(record) for record in manifest["solution_roots_numeric"]
    ]
    solution_root_positions = [int(value) for value in manifest["solution_root_positions"]]
    initial_vectors = [
        acb_mat([[acb_from_numeric_record(value)] for value in records])
        for records in manifest["initial_vectors_numeric"]
    ]
    gates = {str(key): bool(value) for key, value in manifest["resonance_gates"].items()}
    resonance_targets = {
        str(key): int(value) for key, value in manifest["resonance_targets"].items()
    }
    series_by_solution: list[list[list[acb_mat]]] = []
    actual_maximum_log_degree = 0
    for solution_index, (root, root_position, initial) in enumerate(
        zip(solution_roots, solution_root_positions, initial_vectors), start=1
    ):
        series: list[list[acb_mat]] = [[initial]]
        for degree in range(1, series_order + 1):
            previous_log_degree = max(len(item) for item in series) - 1
            right_by_log = [acb_mat(dimension, 1) for _ in range(previous_log_degree + 1)]
            for regular_degree in range(degree):
                earlier = series[degree - 1 - regular_degree]
                for log_degree, coefficient in enumerate(earlier):
                    right_by_log[log_degree] += regular[regular_degree] * coefficient
            absolute_power = root + degree
            gate_prefix = f"{solution_index}:{degree}:"
            target_positions = {
                target
                for key, target in resonance_targets.items()
                if key.startswith(gate_prefix)
            }
            if not target_positions:
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
                if len(target_positions) != 1:
                    raise ValueError(
                        f"numeric Frobenius manifest has ambiguous resonance targets at solution {solution_index}, degree {degree}"
                    )
                target_position = next(iter(target_positions))
                reduced_inverse = acb_mat(dimension, dimension)
                for position, (projector, spectral_root) in enumerate(zip(projectors, roots)):
                    if position != target_position:
                        reduced_inverse += projector / (absolute_power - spectral_root)
                coefficients = [
                    acb_mat(dimension, 1) for _ in range(previous_log_degree + 2)
                ]
                for log_degree in range(previous_log_degree, -1, -1):
                    corrected = acb_mat(right_by_log[log_degree])
                    corrected -= coefficients[log_degree + 1] * (log_degree + 1)
                    gate_key = f"{solution_index}:{degree}:{log_degree}"
                    if gate_key not in gates:
                        raise ValueError(
                            f"the numerical Frobenius manifest is missing the "
                            f"resonance decision encoded as '{gate_key}' (solution "
                            f"index, series degree, log degree); regenerate the "
                            f"manifest so every resonance gate is recorded"
                        )
                    if gates[gate_key]:
                        defect = projectors[target_position] * corrected
                        coefficients[log_degree + 1] += defect / (log_degree + 1)
                        corrected -= defect
                    coefficients[log_degree] = reduced_inverse * corrected
                while len(coefficients) > 1 and all(
                    abs(coefficients[-1][row, 0]).contains(0) for row in range(dimension)
                ):
                    coefficients.pop()
            actual_maximum_log_degree = max(actual_maximum_log_degree, len(coefficients) - 1)
            series.append(coefficients)
        series_by_solution.append(series)

    def evaluate(point: acb) -> acb_mat:
        logarithm = point.log()
        columns: list[acb_mat] = []
        for root, series in zip(solution_roots, series_by_solution):
            value = acb_mat(dimension, 1)
            for degree_coefficients in reversed(series):
                log_polynomial = acb_mat(dimension, 1)
                for coefficient in reversed(degree_coefficients):
                    log_polynomial = log_polynomial * logarithm + coefficient
                value = value * point + log_polynomial
            columns.append(value * point**root)
        return acb_mat(
            [[columns[column][row, 0] for column in range(dimension)] for row in range(dimension)]
        )

    return PowerLogBasis(dimension, actual_maximum_log_degree, manifest, evaluate)


def build_power_log_basis(
    system: RegularSingularSystem | NumericalRegularSingularSystem,
    manifest: dict[str, Any],
    *,
    series_order: int,
) -> PowerLogBasis:
    """按 exact manifest 构造正则奇点的 FLINT power-log fundamental matrix。"""

    if manifest.get("status") != "passed" or manifest.get("system_name") != system.name:
        raise ValueError(
            "the supplied Frobenius manifest does not belong to this system or has "
            "not passed validation; regenerate it for the current system before "
            "solving"
        )
    route = manifest.get("route")
    if route == "single_root_jordan_exact_gate":
        return _build_single_root_basis(system, manifest, series_order)
    if route == "diagonalizable_roots_exact_gate":
        return _build_diagonalizable_basis(system, manifest, series_order)
    if route == "known_roots_defective_exact_gate":
        if not isinstance(system, RegularSingularSystem):
            raise TypeError("known-roots defective manifest requires RegularSingularSystem")
        return _build_known_roots_defective_basis(system, manifest, series_order)
    if route == "single_root_jordan_numeric_gate":
        if not isinstance(system, NumericalRegularSingularSystem):
            raise TypeError("numeric Frobenius manifest requires NumericalRegularSingularSystem")
        _, regular = system.acb_data()
        regular = regular[:series_order]
        regular.extend(
            acb_mat(system.dimension, system.dimension)
            for _ in range(series_order - len(regular))
        )
        return _build_single_root_basis_from_data(
            dimension=system.dimension,
            regular=regular,
            root=acb_from_numeric_record(manifest["root_numeric"]),
            nilpotent=acb_matrix_from_numeric_records(manifest["nilpotent_numeric"]),
            manifest=manifest,
            series_order=series_order,
        )
    if route == "diagonalizable_roots_numeric_gate":
        if not isinstance(system, NumericalRegularSingularSystem):
            raise TypeError("numeric Frobenius manifest requires NumericalRegularSingularSystem")
        return _build_numeric_diagonalizable_basis(system, manifest, series_order)
    raise ValueError(
        f"the Frobenius manifest declares route '{route}', which this solver does "
        f"not recognize; regenerate the manifest with a supported route"
    )
