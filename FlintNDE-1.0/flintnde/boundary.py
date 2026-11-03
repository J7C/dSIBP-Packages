"""奇点起点的 power-log/指数边界输入、精确验证与基常数恢复。

用户以 ``{a,b,C}`` 记录 ``z**a log(z)**b C`` 的最高 log 领头项；``a`` 和
``C`` 必须属于 Q(i)，``b`` 为非负整数。正则奇点直接根据 exact Frobenius manifest
验证 indicial root、log 次数和向量相容性；高阶 pole 则由局部调度器先完成 shearing
或指数 sector 分解，再复用本协议验证各 sector 的 power-log 部分。可复用的本性奇点
边界另用 ``{phi,a,b,C}``，表示 ``exp(phi(z)) z**a log(z)**b C``，其中 ``phi``
只含局部变量的负整数幂及 exact Q(i) 系数。
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from flint import acb, acb_mat, arb

from .exact_gaussian import (
    GaussianMatrix,
    GaussianRational,
    canonical_exact_solve,
    exact_column as _column,
    exact_matrix_from_columns as _matrix_from_columns,
    gaussian_rational,
)
from .frobenius import RegularSingularSystem
from .numeric_structure import (
    acb_from_numeric_record,
    numerical_rref_solve,
)


@dataclass(frozen=True)
class FrobeniusBoundaryTerm:
    """保存一个经过 exact 解析的 ``z**a log(z)**b C`` 领头项。"""

    a: GaussianRational
    b: int
    C: tuple[GaussianRational, ...]

    def __post_init__(self) -> None:
        """即使用户直接构造 dataclass，也执行与工厂相同的 exact 门禁。"""

        if isinstance(self.b, bool) or not isinstance(self.b, int) or self.b < 0:
            raise TypeError("Frobenius boundary b must be a nonnegative integer")
        if not isinstance(self.C, (list, tuple)) or not self.C:
            raise TypeError("Frobenius boundary C must be a nonempty vector")
        power = gaussian_rational(self.a)
        vector = tuple(gaussian_rational(value) for value in self.C)
        if all(value.is_zero for value in vector):
            raise ValueError("Frobenius boundary C must not be the zero vector")
        object.__setattr__(self, "a", power)
        object.__setattr__(self, "C", vector)

    def to_json(self) -> dict[str, Any]:
        """返回不损失 Q(i) 精度的用户边界记录。"""

        return {"a": str(self.a), "b": self.b, "C": [str(value) for value in self.C]}


@dataclass(frozen=True)
class FrobeniusBoundary:
    """保存一个或多个可叠加的奇点 power-log 领头分支。"""

    terms: tuple[FrobeniusBoundaryTerm, ...]

    def __post_init__(self) -> None:
        """验证直接构造时的 term 类型、非空性和共同维数。"""

        parsed = tuple(_parse_term(term, position) for position, term in enumerate(self.terms, 1))
        if not parsed:
            raise ValueError("Frobenius boundary must contain at least one term")
        if len({len(term.C) for term in parsed}) != 1:
            raise ValueError("all Frobenius boundary C vectors must have the same dimension")
        object.__setattr__(self, "terms", parsed)

    def to_json(self) -> list[dict[str, Any]]:
        """按用户输入顺序序列化全部分支。"""

        return [term.to_json() for term in self.terms]


def _parse_phi(records: Any, position: int) -> tuple[tuple[int, GaussianRational], ...]:
    """解析并 canonicalize 一个无常数项的 exact Laurent 指数 ``phi``。"""

    if not isinstance(records, (list, tuple)):
        raise TypeError(
            f"exponential boundary term {position}: phi must be a sequence"
        )
    coefficients: dict[int, GaussianRational] = {}
    for phi_position, record in enumerate(records, 1):
        if not isinstance(record, Mapping) or set(record) != {"power", "coefficient"}:
            raise ValueError(
                f"exponential boundary term {position}: phi term {phi_position} must contain "
                "exactly power and coefficient"
            )
        power = record["power"]
        if isinstance(power, bool) or not isinstance(power, int) or power >= 0:
            raise TypeError(
                f"exponential boundary term {position}: phi power must be a negative integer"
            )
        try:
            coefficient = gaussian_rational(record["coefficient"])
        except (TypeError, ValueError) as error:
            raise type(error)(
                f"exponential boundary term {position}: phi coefficient must be exact Q(i)"
            ) from error
        coefficients[power] = coefficients.get(power, GaussianRational()) + coefficient
    canonical = tuple(
        (power, coefficient)
        for power, coefficient in sorted(coefficients.items())
        if not coefficient.is_zero
    )
    return canonical


@dataclass(frozen=True)
class ExponentialBoundaryTerm:
    """保存 ``exp(phi(z)) z**a log(z)**b C`` 的一个 exact 领头项。"""

    phi: tuple[tuple[int, GaussianRational], ...]
    a: GaussianRational
    b: int
    C: tuple[GaussianRational, ...]

    def __post_init__(self) -> None:
        """对直接构造的对象执行与公开工厂一致的 exact 门禁。"""

        phi_records = [
            record
            if isinstance(record, Mapping)
            else {"power": record[0], "coefficient": record[1]}
            for record in self.phi
        ]
        parsed_phi = _parse_phi(phi_records, 1)
        power_log = FrobeniusBoundaryTerm(self.a, self.b, self.C)
        object.__setattr__(self, "phi", parsed_phi)
        object.__setattr__(self, "a", power_log.a)
        object.__setattr__(self, "b", power_log.b)
        object.__setattr__(self, "C", power_log.C)

    def to_frobenius_term(self) -> FrobeniusBoundaryTerm:
        """丢开已单独验证的指数签名，返回同一 power-log 领头项。"""

        return FrobeniusBoundaryTerm(self.a, self.b, self.C)

    def to_json(self) -> dict[str, Any]:
        """返回可直接重新输入且不损失 Q(i) 精度的记录。"""

        return {
            "phi": [
                {"power": power, "coefficient": str(coefficient)}
                for power, coefficient in self.phi
            ],
            "a": str(self.a),
            "b": self.b,
            "C": [str(value) for value in self.C],
        }


@dataclass(frozen=True)
class ExponentialBoundary:
    """保存一个或多个带明确指数签名的本性奇点边界分支。"""

    terms: tuple[ExponentialBoundaryTerm, ...]

    def __post_init__(self) -> None:
        """验证 term 类型、非空性和共同原积分基维数。"""

        parsed = tuple(
            _parse_exponential_term(term, position)
            for position, term in enumerate(self.terms, 1)
        )
        if not parsed:
            raise ValueError("exponential boundary must contain at least one term")
        if len({len(term.C) for term in parsed}) != 1:
            raise ValueError("all exponential boundary C vectors must have the same dimension")
        object.__setattr__(self, "terms", parsed)

    def to_frobenius_boundary(self) -> FrobeniusBoundary:
        """返回供各指数 sector 的 power-log 引擎复用的边界。"""

        return FrobeniusBoundary(tuple(term.to_frobenius_term() for term in self.terms))

    def to_json(self) -> list[dict[str, Any]]:
        """按输入次序序列化所有指数分支。"""

        return [term.to_json() for term in self.terms]


def _parse_exponential_term(record: Any, position: int) -> ExponentialBoundaryTerm:
    """解析一个 ``{phi,a,b,C}`` 映射并拒绝任何含混或浮点字段。"""

    if isinstance(record, ExponentialBoundaryTerm):
        return record
    if not isinstance(record, Mapping) or set(record) not in (
        {"phi", "a", "b", "C"},
        {"phi", "a", "b", "C", "CBalls"},
    ):
        raise ValueError(
            f"exponential boundary term {position} must contain phi, a, b, C and optional CBalls"
        )
    phi = _parse_phi(record["phi"], position)
    power_log = _parse_term(
        {"a": record["a"], "b": record["b"], "C": record["C"]},
        position,
    )
    return ExponentialBoundaryTerm(phi, power_log.a, power_log.b, power_log.C)


def exponential_boundary(
    records: Iterable[Any] | ExponentialBoundary,
) -> ExponentialBoundary:
    """验证并冻结一列 ``{phi,a,b,C}`` 本性奇点边界记录。"""

    if isinstance(records, ExponentialBoundary):
        return records
    if isinstance(records, (str, bytes, Mapping)):
        raise TypeError(
            "exponential boundary must be a nonempty sequence of {phi,a,b,C} terms"
        )
    try:
        terms = tuple(
            _parse_exponential_term(record, position)
            for position, record in enumerate(records, 1)
        )
    except TypeError as error:
        if "not iterable" in str(error):
            raise TypeError(
                "exponential boundary must be a nonempty sequence of {phi,a,b,C} terms"
            ) from error
        raise
    return ExponentialBoundary(terms)


def _parse_term(record: Any, position: int) -> FrobeniusBoundaryTerm:
    """解析单个字典或三元组，并严格拒绝多余字段和浮点 exact 输入。"""

    if isinstance(record, FrobeniusBoundaryTerm):
        return record
    if isinstance(record, Mapping):
        if set(record) not in ({"a", "b", "C"}, {"a", "b", "C", "CBalls"}):
            raise ValueError(
                f"Frobenius boundary term {position} must contain a, b, C and optional CBalls"
            )
        a_value, b_value, vector_value = record["a"], record["b"], record["C"]
    elif isinstance(record, (list, tuple)) and len(record) == 3:
        a_value, b_value, vector_value = record
    else:
        raise TypeError(
            f"Frobenius boundary term {position} must be a {{a,b,C}} mapping or (a,b,C)"
        )
    if isinstance(b_value, bool) or not isinstance(b_value, int) or b_value < 0:
        raise TypeError(f"Frobenius boundary term {position}: b must be a nonnegative integer")
    if not isinstance(vector_value, (list, tuple)) or not vector_value:
        raise TypeError(f"Frobenius boundary term {position}: C must be a nonempty vector")
    try:
        power = gaussian_rational(a_value)
        vector = tuple(gaussian_rational(value) for value in vector_value)
    except (TypeError, ValueError) as error:
        raise type(error)(
            f"Frobenius boundary term {position}: a and every C entry must be exact Q(i)"
        ) from error
    if all(value.is_zero for value in vector):
        raise ValueError(f"Frobenius boundary term {position}: C must not be the zero vector")
    return FrobeniusBoundaryTerm(power, b_value, vector)


def frobenius_boundary(records: Iterable[Any] | FrobeniusBoundary) -> FrobeniusBoundary:
    """验证并冻结一列 ``{a,b,C}`` 记录。

    推荐输入为 ``[{"a": ..., "b": ..., "C": [...]}, ...]``；每项也可写成
    ``(a,b,C)``。这里只验证记录结构和 exact 数域，indicial/log 相容性在获得具体
    微分方程与路径起点后验证。指数型局部解不要求用户另填指数：程序从 exact 高阶
    Laurent 矩阵推断指数 sector；若一个 ``C`` 跨越多个 sector，用户必须拆成多项。
    """

    if isinstance(records, FrobeniusBoundary):
        return records
    if isinstance(records, (str, bytes, Mapping)):
        raise TypeError("Frobenius boundary must be a nonempty sequence of {a,b,C} terms")
    try:
        terms = tuple(_parse_term(record, position) for position, record in enumerate(records, 1))
    except TypeError as error:
        if "not iterable" in str(error):
            raise TypeError(
                "Frobenius boundary must be a nonempty sequence of {a,b,C} terms"
            ) from error
        raise
    if not terms:
        raise ValueError("Frobenius boundary must contain at least one term")
    dimensions = {len(term.C) for term in terms}
    if len(dimensions) != 1:
        raise ValueError("all Frobenius boundary C vectors must have the same dimension")
    return FrobeniusBoundary(terms)


def _canonical_exact_solve(left: GaussianMatrix, right: GaussianMatrix) -> GaussianMatrix:
    """本模块错误消息合同的薄适配；RREF 与代入验证共用 exact 公共实现。"""

    return canonical_exact_solve(
        left,
        right,
        incompatible_message="Frobenius leading vector is incompatible with the exact local basis",
        solve_message="exact Frobenius boundary solve failed its substitution check",
        dimension_message="exact boundary linear solve dimension mismatch",
    )


def _single_root_constants(
    boundary: FrobeniusBoundary,
    manifest: dict[str, Any],
    dimension: int,
) -> tuple[GaussianMatrix, list[dict[str, Any]]]:
    """由最高 log 领头项恢复单根 Jordan 基常数。"""

    root = gaussian_rational(manifest["root_exact"])
    nilpotent = GaussianMatrix.from_records(manifest["nilpotent_exact"])
    maximum_log_degree = int(manifest["maximum_log_degree"])
    constants = GaussianMatrix.zero(dimension, 1)
    details: list[dict[str, Any]] = []
    for position, term in enumerate(boundary.terms, 1):
        if term.a != root:
            raise ValueError(
                f"Frobenius boundary term {position}: a={term.a} is not the indicial root {root}"
            )
        if term.b > maximum_log_degree:
            raise ValueError(
                f"Frobenius boundary term {position}: b={term.b} exceeds exact maximum log degree "
                f"{maximum_log_degree}"
            )
        leading = _column(term.C)
        if not (nilpotent * leading).is_zero:
            raise ValueError(
                f"Frobenius boundary term {position}: C is not a highest-log coefficient "
                "because the nilpotent residue still acts nontrivially"
            )
        preimage = _canonical_exact_solve(
            nilpotent ** term.b,
            leading * math.factorial(term.b),
        )
        if not (nilpotent ** (term.b + 1) * preimage).is_zero:
            raise ArithmeticError("canonical Jordan preimage has a higher unexpected log power")
        constants += preimage
        details.append(
            {
                "term_index": position,
                "input": term.to_json(),
                "canonical_basis_constants": [
                    str(preimage.scalar(row, 0)) for row in range(dimension)
                ],
                "completion": "lower log powers are generated by exp(N log(z))",
            }
        )
    return constants, details


def _diagonalizable_constants(
    boundary: FrobeniusBoundary,
    manifest: dict[str, Any],
    dimension: int,
) -> tuple[GaussianMatrix, list[dict[str, Any]]]:
    """把半单 residue 的任意 exact 本征向量展开到 canonical solution columns。"""

    roots = [gaussian_rational(value) for value in manifest["solution_roots_exact"]]
    columns = [
        _column(tuple(gaussian_rational(value) for value in records))
        for records in manifest["initial_vectors_exact"]
    ]
    constants = GaussianMatrix.zero(dimension, 1)
    details: list[dict[str, Any]] = []
    for position, term in enumerate(boundary.terms, 1):
        if term.b != 0:
            raise ValueError(
                f"Frobenius boundary term {position}: a semisimple leading branch requires b=0; "
                "resonance logs generated at higher series order are automatic"
            )
        positions = [index for index, root in enumerate(roots) if root == term.a]
        if not positions:
            spectrum = sorted({str(root) for root in roots})
            raise ValueError(
                f"Frobenius boundary term {position}: a={term.a} is not in the exact indicial "
                f"spectrum {spectrum}"
            )
        eigenbasis = _matrix_from_columns([columns[index] for index in positions], dimension)
        local_constants = _canonical_exact_solve(eigenbasis, _column(term.C))
        full_constants = GaussianMatrix.zero(dimension, 1)
        for local_index, solution_index in enumerate(positions):
            full_constants.real[solution_index, 0] = local_constants.real[local_index, 0]
            full_constants.imag[solution_index, 0] = local_constants.imag[local_index, 0]
        constants += full_constants
        details.append(
            {
                "term_index": position,
                "input": term.to_json(),
                "canonical_basis_constants": [
                    str(full_constants.scalar(row, 0)) for row in range(dimension)
                ],
                "completion": "higher-order resonance log terms are generated automatically",
            }
        )
    return constants, details


def resolve_frobenius_boundary(
    system: RegularSingularSystem,
    manifest: dict[str, Any],
    value: Any,
) -> tuple[FrobeniusBoundary, GaussianMatrix, dict[str, Any]]:
    """验证具体局部系统的边界，并返回 power-log 基常数及审计记录。"""

    boundary = frobenius_boundary(value)
    dimension = system.dimension
    if any(len(term.C) != dimension for term in boundary.terms):
        raise ValueError(
            f"Frobenius boundary C dimension must equal differential-equation dimension {dimension}"
        )
    route = manifest.get("route")
    if route == "single_root_jordan_exact_gate":
        constants, details = _single_root_constants(boundary, manifest, dimension)
    elif route == "diagonalizable_roots_exact_gate":
        constants, details = _diagonalizable_constants(boundary, manifest, dimension)
    else:
        raise ValueError(
            "Frobenius boundary requires a supported exact manifest; route "
            f"'{route}' cannot certify {{a,b,C}} exactly; numerical manifests use "
            "resolve_numerical_frobenius_boundary"
        )
    if constants.is_zero:
        raise ValueError("combined Frobenius boundary selects the zero solution")
    report = {
        "schema": "flintnde_frobenius_boundary_v1",
        "system_name": system.name,
        "exact_field": "Q(i)",
        "terms": boundary.to_json(),
        "canonical_basis_constants": [
            str(constants.scalar(row, 0)) for row in range(dimension)
        ],
        "term_resolutions": details,
    }
    return boundary, constants, report


def _acb_column(vector: tuple[GaussianRational, ...]) -> acb_mat:
    """把 exact Q(i) 列向量转换为 Acb 列。"""

    column = acb_mat(len(vector), 1)
    for row, value in enumerate(vector):
        column[row, 0] = value.to_acb()
    return column


def _acb_column_from_records(records: list[dict[str, str]]) -> acb_mat:
    """由数值 manifest 的向量记录恢复 Acb 列。"""

    column = acb_mat(len(records), 1)
    for row, record in enumerate(records):
        column[row, 0] = acb_from_numeric_record(record)
    return column


def _acb_matrix_from_columns(columns: list[acb_mat], dimension: int) -> acb_mat:
    """按列拼接 Acb 矩阵。"""

    if not columns:
        raise ValueError("matrix requires at least one column")
    matrix = acb_mat(dimension, len(columns))
    for column_index, column in enumerate(columns):
        for row in range(dimension):
            matrix[row, column_index] = column[row, 0]
    return matrix


def _acb_matrix_max_abs(matrix: acb_mat) -> arb:
    """返回 Acb 矩阵元素绝对值中点的最大值。"""

    magnitudes = [
        abs(matrix[row, column]).mid()
        for row in range(matrix.nrows())
        for column in range(matrix.ncols())
    ]
    return max(magnitudes) if magnitudes else arb(0)


def _acb_matrix_power(matrix: acb_mat, exponent: int) -> acb_mat:
    """非负整数次幂；零次幂为单位矩阵。"""

    if exponent < 0:
        raise ValueError("matrix exponent must be nonnegative")
    result = acb_mat(matrix.nrows(), matrix.ncols())
    for index in range(matrix.nrows()):
        result[index, index] = 1
    for _ in range(exponent):
        result = result * matrix
    return result


def _acb_string(value: acb) -> str:
    """记录 Acb 中点的无损十进制串。"""

    return acb(value.real.mid(), value.imag.mid()).str(radius=False)


def _numerical_cutoff(manifest: dict[str, Any]) -> arb:
    """从数值 manifest 诊断恢复统一的绝对零斩杀线。"""

    return arb(manifest["numeric_diagnostics"]["absolute_zero_cutoff"])


def _numeric_single_root_constants(
    boundary: FrobeniusBoundary,
    manifest: dict[str, Any],
    dimension: int,
    cutoff: arb,
) -> tuple[acb_mat, list[dict[str, Any]]]:
    """数值单根 Jordan 路线的 ``{a,b,C}`` 常数恢复，镜像 exact 版本。

    谱匹配、``N C=0`` 最高 log 判定和 ``N^b w = C b!`` 原像都在同一条
    绝对斩杀线内执行；原像由阈值 RREF 求解并附带残差记录。
    """

    root = acb_from_numeric_record(manifest["root_numeric"])
    nilpotent = _acb_matrix_from_records(manifest["nilpotent_numeric"])
    maximum_log_degree = int(manifest["maximum_log_degree"])
    constants = acb_mat(dimension, 1)
    details: list[dict[str, Any]] = []
    for position, term in enumerate(boundary.terms, 1):
        claimed = term.a.to_acb()
        if abs(root - claimed).mid() > cutoff:
            raise ValueError(
                f"Frobenius boundary term {position}: a={term.a} does not match the "
                f"numerical indicial root {_acb_string(root)} within the "
                "classification cutoff"
            )
        if term.b > maximum_log_degree:
            raise ValueError(
                f"Frobenius boundary term {position}: b={term.b} exceeds numerical "
                f"maximum log degree {maximum_log_degree}"
            )
        leading = _acb_column(term.C)
        nilpotent_leading = _acb_matrix_max_abs(nilpotent * leading)
        if nilpotent_leading > cutoff:
            raise ValueError(
                f"Frobenius boundary term {position}: C is not a highest-log "
                "coefficient because the nilpotent residue still acts beyond the "
                "classification cutoff"
            )
        factorial = acb(math.factorial(term.b))
        preimage = numerical_rref_solve(
            _acb_matrix_power(nilpotent, term.b),
            leading * factorial,
            cutoff,
        )
        higher_residual = _acb_matrix_max_abs(
            _acb_matrix_power(nilpotent, term.b + 1) * preimage
        )
        if higher_residual > cutoff:
            raise ArithmeticError(
                "numerical Jordan preimage retains a higher log power beyond the "
                "classification cutoff"
            )
        constants += preimage
        details.append(
            {
                "term_index": position,
                "input": term.to_json(),
                "canonical_basis_constants": [
                    _acb_string(preimage[row, 0]) for row in range(dimension)
                ],
                "nilpotent_leading_residual": nilpotent_leading.str(radius=False),
                "higher_log_residual": higher_residual.str(radius=False),
                "completion": "lower log powers are generated by exp(N log(z))",
            }
        )
    return constants, details


def _numeric_diagonalizable_constants(
    boundary: FrobeniusBoundary,
    manifest: dict[str, Any],
    dimension: int,
    cutoff: arb,
) -> tuple[acb_mat, list[dict[str, Any]]]:
    """数值对角化路线的 ``{a,b,C}`` 常数恢复，镜像 exact 半单版本。

    每个 term 要求 ``b=0``；``a`` 与 manifest 谱根的差必须落在斩杀线内，
    ``C`` 再由对应本征子空间列组成的矩形系统做阈值 RREF 求解。
    """

    roots = [acb_from_numeric_record(record) for record in manifest["solution_roots_numeric"]]
    columns = [
        _acb_column_from_records(records) for records in manifest["initial_vectors_numeric"]
    ]
    constants = acb_mat(dimension, 1)
    details: list[dict[str, Any]] = []
    for position, term in enumerate(boundary.terms, 1):
        if term.b != 0:
            raise ValueError(
                f"Frobenius boundary term {position}: a semisimple leading branch requires b=0; "
                "resonance logs generated at higher series order are automatic"
            )
        claimed = term.a.to_acb()
        positions = [
            index
            for index, root in enumerate(roots)
            if abs(root - claimed).mid() <= cutoff
        ]
        if not positions:
            spectrum = sorted(
                {
                    _acb_string(acb_from_numeric_record(record))
                    for record in manifest["roots_numeric"]
                }
            )
            raise ValueError(
                f"Frobenius boundary term {position}: a={term.a} does not match the "
                f"numerical indicial spectrum {spectrum} within the classification "
                "cutoff"
            )
        eigenbasis = _acb_matrix_from_columns(
            [columns[index] for index in positions], dimension
        )
        leading = _acb_column(term.C)
        local_constants = numerical_rref_solve(eigenbasis, leading, cutoff)
        residual = _acb_matrix_max_abs(eigenbasis * local_constants - leading)
        for local_index, solution_index in enumerate(positions):
            constants[solution_index, 0] += local_constants[local_index, 0]
        details.append(
            {
                "term_index": position,
                "input": term.to_json(),
                "canonical_basis_constants": [
                    _acb_string(constants[row, 0]) for row in range(dimension)
                ],
                "eigen_expansion_residual": residual.str(radius=False),
                "completion": "higher-order resonance log terms are generated automatically",
            }
        )
    return constants, details


def _acb_matrix_from_records(records: list[list[dict[str, str]]]) -> acb_mat:
    """由数值 manifest 的矩阵记录恢复 Acb 矩阵。"""

    matrix = acb_mat(len(records), len(records[0]))
    for row, matrix_row in enumerate(records):
        for column, record in enumerate(matrix_row):
            matrix[row, column] = acb_from_numeric_record(record)
    return matrix


def resolve_numerical_frobenius_boundary(
    manifest: dict[str, Any],
    value: Any,
) -> tuple[FrobeniusBoundary, acb_mat, dict[str, Any]]:
    """验证数值 Frobenius manifest 的边界，并返回 power-log 基常数及审计记录。

    谱匹配、最高 log 系数判定和 Jordan 原像都按 manifest 自带的绝对零斩杀线
    执行；报告记录阈值、残差和每项常数，供调用方与证书谱交叉核对。
    只支持数值单根 Jordan 与对角化两条路线。
    """

    boundary = frobenius_boundary(value)
    dimension = int(manifest["dimension"])
    if any(len(term.C) != dimension for term in boundary.terms):
        raise ValueError(
            f"Frobenius boundary C dimension must equal differential-equation dimension {dimension}"
        )
    diagnostics = manifest["numeric_diagnostics"]
    cutoff = _numerical_cutoff(manifest)
    route = manifest.get("route")
    if route == "single_root_jordan_numeric_gate":
        constants, details = _numeric_single_root_constants(
            boundary, manifest, dimension, cutoff
        )
    elif route == "diagonalizable_roots_numeric_gate":
        constants, details = _numeric_diagonalizable_constants(
            boundary, manifest, dimension, cutoff
        )
    else:
        raise ValueError(
            "numerical Frobenius boundary matching is only implemented for the "
            f"single-root Jordan and diagonalizable routes; manifest declares '{route}'"
        )
    if _acb_matrix_max_abs(constants) <= cutoff:
        raise ValueError("combined Frobenius boundary selects the zero solution")
    report = {
        "schema": "flintnde_numerical_frobenius_boundary_v1",
        "system_name": manifest["system_name"],
        "classification": "numerical_threshold",
        "absolute_zero_cutoff": diagnostics["absolute_zero_cutoff"],
        "terms": boundary.to_json(),
        "canonical_basis_constants": [
            _acb_string(constants[row, 0]) for row in range(dimension)
        ],
        "term_resolutions": details,
    }
    return boundary, constants, report
