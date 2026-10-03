"""奇点局部基的统一构造、边界验证和 generalized-series 调度。

高阶 pole 先调用 exact Fuchsian reduction。若单项式 shearing 成功，模块在降阶基中
复用现有 power-log 引擎，并在普通匹配点乘回原始积分基。若仍有高阶 pole，只接受
能够 exact 常数对角化为互不耦合 sector、且每个 sector 的高阶部分为标量的系统；此时
严格抽出 ``exp(Phi(z))`` 后再调用同一 power-log 引擎。对二阶 pole、主导特征值单重且
互异的耦合系统，还可构造逐 sector 的形式渐近级数；求值固定保留用户要求的 ``N`` 阶，
并额外计算五阶用于收敛诊断。该路线只供奇点起点初始化，不提供跨 Stokes sector 的连接矩阵。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

from flint import acb, acb_mat, arb

from .boundary import (
    ExponentialBoundary,
    FrobeniusBoundary,
    FrobeniusBoundaryTerm,
    frobenius_boundary,
    resolve_frobenius_boundary,
)
from .exact_gaussian import (
    GaussianMatrix,
    GaussianRational,
    canonical_exact_solve,
    exact_column as _column,
    exact_matrix_from_columns,
    gaussian_rational,
)
from .frobenius import (
    build_exact_power_log_series,
    build_frobenius_manifest,
    build_power_log_basis,
)
from .fuchsian import (
    FuchsianReductionResult,
    MoserBalanceTransformation,
    attempt_fuchsian_reduction,
)
from .singularities import RationalMatrixSystem, rational_function


class LocalReductionError(NotImplementedError):
    """表示当前 exact 降阶与指数 sector 门禁都不能认证局部基。"""


@dataclass(frozen=True)
class LocalSolutionBasis:
    """统一封装普通 Frobenius、shearing 和指数广义幂级数基。"""

    dimension: int
    method: str
    maximum_log_degree: int
    manifest: dict[str, Any]
    _evaluate: Callable[[acb], acb_mat] = field(repr=False, compare=False)
    _resolve_boundary: Callable[[Any], tuple[GaussianMatrix, dict[str, Any]]] = field(
        repr=False,
        compare=False,
    )
    continuation_ready: bool = True
    _evaluation_report: Callable[[acb], dict[str, Any]] | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    _match_distance_report: Callable[[int, float], tuple[arb, dict[str, Any]]] | None = field(
        default=None,
        repr=False,
        compare=False,
    )

    def evaluate(self, local_point: acb) -> acb_mat:
        """在非零局部点计算原始积分基中的 fundamental matrix。"""

        if local_point.contains(0):
            raise ZeroDivisionError("a singular local basis cannot be evaluated at its center")
        return self._evaluate(acb(local_point))

    def resolve_boundary(self, value: Any) -> tuple[GaussianMatrix, dict[str, Any]]:
        """验证用户 ``{a,b,C}`` 记录并返回本局部基的常数列。"""

        return self._resolve_boundary(value)

    def evaluation_report(self, local_point: acb) -> dict[str, Any] | None:
        """返回局部求值诊断；收敛 power-log 基无需额外截断报告。"""

        if self._evaluation_report is None:
            return None
        return self._evaluation_report(acb(local_point))

    def suggest_match_distance(
        self,
        target_order: int,
        minimum_order_factor: float = 3.0,
    ) -> tuple[arb, dict[str, Any]]:
        """由最近指数根差估计形式渐近起点的保守匹配距离。"""

        if self._match_distance_report is None:
            raise ValueError(f"{self.method} does not use a formal-asymptotic match estimate")
        if target_order <= 0:
            raise ValueError("target_order must be positive")
        if minimum_order_factor <= 1:
            raise ValueError("minimum_order_factor must exceed one")
        return self._match_distance_report(target_order, minimum_order_factor)


def _exponential_power_log_boundary(
    value: Any,
) -> tuple[FrobeniusBoundary, tuple[tuple[tuple[int, GaussianRational], ...], ...]]:
    """要求显式指数签名，并把 power-log 字段交给 Frobenius 引擎。"""

    if not isinstance(value, ExponentialBoundary):
        raise TypeError(
            "an exponential singularity requires ExponentialBoundary {phi,a,b,C} data"
        )
    return value.to_frobenius_boundary(), tuple(term.phi for term in value.terms)


def _matrix_from_columns(columns: list[GaussianMatrix], dimension: int) -> GaussianMatrix:
    """本模块防御合同的薄适配：要求完整列集，再交给 exact 公共实现。"""

    if len(columns) != dimension:
        raise ValueError("exact eigenbasis does not contain a full set of columns")
    return exact_matrix_from_columns(columns, dimension)


def _transpose(matrix: GaussianMatrix) -> GaussianMatrix:
    """返回 exact Q(i) 转置；左本征向量按双线性而非 Hermitian 约定。"""

    return GaussianMatrix.from_records(
        [
            [matrix.scalar(column, row) for column in range(matrix.nrows)]
            for row in range(matrix.ncols)
        ]
    )


def _bilinear(left: GaussianMatrix, right: GaussianMatrix) -> GaussianRational:
    """计算两个 exact 列向量的无共轭双线性配对。"""

    if left.ncols != 1 or right.ncols != 1 or left.nrows != right.nrows:
        raise ValueError("exact bilinear pairing requires equal-length column vectors")
    return sum(
        (
            left.scalar(index, 0) * right.scalar(index, 0)
            for index in range(left.nrows)
        ),
        GaussianRational(),
    )


def _canonical_exact_solve(left: GaussianMatrix, right: GaussianMatrix) -> GaussianMatrix:
    """本模块错误消息合同的薄适配；RREF 与代入验证共用 exact 公共实现。"""

    return canonical_exact_solve(
        left,
        right,
        incompatible_error=LocalReductionError,
        incompatible_message="formal exponential recurrence is exactly incompatible",
        solve_message="exact formal recurrence solve failed substitution",
        dimension_message="exact formal recurrence solve dimension mismatch",
    )


def _power_log_local_basis(
    system: RationalMatrixSystem,
    center: GaussianRational,
    series_order: int,
    *,
    allow_ordinary: bool,
    method: str,
    extra_manifest: dict[str, Any] | None = None,
) -> LocalSolutionBasis:
    """把一个至多 simple-pole 的系统包装为统一局部基。"""

    local_system = system.regular_singular_system_at(
        center,
        series_order,
        allow_ordinary=allow_ordinary,
    )
    frobenius_manifest = build_frobenius_manifest(local_system)
    basis = build_power_log_basis(
        local_system,
        frobenius_manifest,
        series_order=series_order,
    )

    def resolve(value: Any) -> tuple[GaussianMatrix, dict[str, Any]]:
        _parsed, constants, report = resolve_frobenius_boundary(
            local_system,
            frobenius_manifest,
            value,
        )
        return constants, report

    manifest = {
        "schema": "flintnde_local_solution_basis_v1",
        "method": method,
        "center_exact": str(center),
        "frobenius": basis.manifest,
    }
    if extra_manifest:
        manifest.update(extra_manifest)
    return LocalSolutionBasis(
        system.dimension,
        method,
        basis.maximum_log_degree,
        manifest,
        basis.evaluate,
        resolve,
    )


def _integer_shift(left: GaussianRational, right: GaussianRational) -> int | None:
    """若 ``left-right`` 是实整数则返回该整数，否则返回 ``None``。"""

    difference = left - right
    if difference.imag != 0 or difference.real.q != 1:
        return None
    return int(difference.real)


def _resolve_moser_boundary(
    boundary: FrobeniusBoundary,
    transformation: MoserBalanceTransformation,
    local_system: Any,
    frobenius_manifest: dict[str, Any],
    series_order: int,
) -> tuple[GaussianMatrix, dict[str, Any]]:
    """在原积分基上 exact 匹配一般 Moser 变换后的 ``{a,b,C}``。

    对每个 canonical 降阶解生成有限 exact power-log jet，乘回累计 Laurent 变换，
    再约束目标幂次之前的系数及同幂更高 log 系数为零。这样相邻阶分量对原基
    leading vector 的贡献不会被遗漏。
    """

    exact_series = build_exact_power_log_series(
        local_system, frobenius_manifest, series_order
    )
    laurent = transformation.laurent_coefficients_exact()
    dimension = transformation.dimension
    solution_maps: list[dict[tuple[GaussianRational, int], GaussianMatrix]] = []
    for root, degrees in zip(exact_series.roots, exact_series.coefficients):
        coefficient_map: dict[tuple[GaussianRational, int], GaussianMatrix] = {}
        for degree, log_coefficients in enumerate(degrees):
            for log_degree, vector in enumerate(log_coefficients):
                for gauge_degree, gauge_matrix in laurent.items():
                    key = (root + degree + gauge_degree, log_degree)
                    coefficient_map[key] = coefficient_map.get(
                        key, GaussianMatrix.zero(dimension, 1)
                    ) + gauge_matrix * vector
        solution_maps.append(
            {key: vector for key, vector in coefficient_map.items() if not vector.is_zero}
        )

    constants = GaussianMatrix.zero(dimension, 1)
    details: list[dict[str, Any]] = []
    for position, term in enumerate(boundary.terms, 1):
        if len(term.C) != dimension:
            raise ValueError(
                f"Frobenius boundary C dimension must equal differential-equation dimension {dimension}"
            )
        candidates = [
            index
            for index, root in enumerate(exact_series.roots)
            if _integer_shift(term.a, root) is not None
        ]
        if not candidates:
            raise ValueError(
                f"Frobenius boundary term {position}: a={term.a} is not integer-related "
                "to any reduced indicial root"
            )
        earlier_keys: set[tuple[GaussianRational, int]] = set()
        higher_log_keys: set[tuple[GaussianRational, int]] = set()
        for candidate in candidates:
            for power, log_degree in solution_maps[candidate]:
                shift = _integer_shift(term.a, power)
                if shift is not None and shift > 0:
                    earlier_keys.add((power, log_degree))
                elif power == term.a and log_degree > term.b:
                    higher_log_keys.add((power, log_degree))
        ordered_constraints = sorted(
            earlier_keys,
            key=lambda item: (-int((term.a - item[0]).real), item[1]),
        ) + sorted(higher_log_keys, key=lambda item: item[1], reverse=True)
        ordered_constraints.append((term.a, term.b))

        rows: list[list[GaussianRational]] = []
        right_values: list[GaussianRational] = []
        target_key = (term.a, term.b)
        for key in ordered_constraints:
            for component in range(dimension):
                row = [
                    solution_maps[candidate].get(
                        key, GaussianMatrix.zero(dimension, 1)
                    ).scalar(component, 0)
                    for candidate in candidates
                ]
                right = term.C[component] if key == target_key else GaussianRational()
                if all(value.is_zero for value in row):
                    if not right.is_zero:
                        raise ValueError(
                            f"Frobenius boundary term {position}: C is incompatible with the "
                            "exact original-basis Moser jet"
                        )
                    continue
                rows.append(row)
                right_values.append(right)
        if not rows:
            raise ValueError(
                f"Frobenius boundary term {position}: exact Moser jet produced no constraints"
            )
        try:
            local_constants = _canonical_exact_solve(
                GaussianMatrix.from_records(rows),
                GaussianMatrix.from_records([[value] for value in right_values]),
            )
        except LocalReductionError as error:
            raise ValueError(
                f"Frobenius boundary term {position}: leading data are incompatible with "
                "the exact original-basis local solutions"
            ) from error
        full_constants = GaussianMatrix.zero(dimension, 1)
        for local_index, solution_index in enumerate(candidates):
            full_constants.real[solution_index, 0] = local_constants.real[local_index, 0]
            full_constants.imag[solution_index, 0] = local_constants.imag[local_index, 0]
        constants += full_constants
        details.append(
            {
                "term_index": position,
                "input": term.to_json(),
                "candidate_solution_indices": [index + 1 for index in candidates],
                "earlier_power_constraints": len(earlier_keys),
                "higher_log_constraints": len(higher_log_keys),
                "canonical_basis_constants": [
                    str(full_constants.scalar(row, 0)) for row in range(dimension)
                ],
                "verification": "exact_original_basis_power_log_jet",
            }
        )
    if constants.is_zero:
        raise ValueError("combined Frobenius boundary selects the zero solution")
    return constants, {
        "schema": "flintnde_moser_frobenius_boundary_v1",
        "exact_field": "Q(i)",
        "terms": boundary.to_json(),
        "canonical_basis_constants": [
            str(constants.scalar(row, 0)) for row in range(dimension)
        ],
        "term_resolutions": details,
        "jet_series_order": series_order,
        "fuchsian_transformation": transformation.to_json(),
    }


def _sheared_power_log_basis(
    reduction: FuchsianReductionResult,
    series_order: int,
) -> LocalSolutionBasis:
    """在 Fuchsian 降阶基中求 power-log，并映射回原始积分基。"""

    transformation = reduction.transformation
    local_system = reduction.transformed_system.regular_singular_system_at(
        transformation.center,
        series_order,
        allow_ordinary=True,
    )
    frobenius_manifest = build_frobenius_manifest(local_system)
    reduced_basis = build_power_log_basis(
        local_system,
        frobenius_manifest,
        series_order=series_order,
    )

    def evaluate(local_point: acb) -> acb_mat:
        return transformation.to_original(local_point, reduced_basis.evaluate(local_point))

    def resolve(value: Any) -> tuple[GaussianMatrix, dict[str, Any]]:
        boundary = frobenius_boundary(value)
        if isinstance(transformation, MoserBalanceTransformation):
            constants, report = _resolve_moser_boundary(
                boundary,
                transformation,
                local_system,
                frobenius_manifest,
                series_order,
            )
            report["fuchsian_reduction"] = reduction.to_json()
            return constants, report
        transformed_terms: list[FrobeniusBoundaryTerm] = []
        details: list[dict[str, Any]] = []
        for position, term in enumerate(boundary.terms, 1):
            if len(term.C) != transformation.dimension:
                raise ValueError(
                    "Frobenius boundary C dimension must equal the reduction dimension"
                )
            transformed, detail = transformation.transform_boundary_term(term, position)
            transformed_terms.append(transformed)
            details.append(detail)
        _parsed, constants, report = resolve_frobenius_boundary(
            local_system,
            frobenius_manifest,
            FrobeniusBoundary(tuple(transformed_terms)),
        )
        report = {
            **report,
            "original_basis_terms": boundary.to_json(),
            "shearing_boundary_map": details,
            "fuchsian_reduction": reduction.to_json(),
        }
        return constants, report

    manifest = {
        "schema": "flintnde_local_solution_basis_v1",
        "method": "fuchsian_reduced_power_log",
        "center_exact": str(transformation.center),
        "frobenius": frobenius_manifest,
        "fuchsian_reduction": reduction.to_json(),
    }
    return LocalSolutionBasis(
        reduction.transformed_system.dimension,
        "fuchsian_reduced_power_log",
        reduced_basis.maximum_log_degree,
        manifest,
        evaluate,
        resolve,
    )


@dataclass(frozen=True)
class _ExponentialSector:
    """保存一个互不耦合指数 sector 及抽掉指数后的 power-log 基。"""

    indices: tuple[int, ...]
    derivative_coefficients: tuple[tuple[int, GaussianRational], ...]
    local_basis: LocalSolutionBasis

    def exponent_records(self) -> list[dict[str, Any]]:
        """返回 ``Phi'(z)`` 与 ``Phi(z)`` 的 exact 系数。"""

        return [
            {
                "matrix_power": -pole_order,
                "phi_power": -(pole_order - 1),
                "derivative_coefficient": str(coefficient),
                "phi_coefficient": str(-coefficient / (pole_order - 1)),
            }
            for pole_order, coefficient in self.derivative_coefficients
            if not coefficient.is_zero
        ]

    def exponential(self, local_point: acb) -> acb:
        """计算 ``exp(Phi(z))``；其中 ``Phi'=`` 高阶标量部分。"""

        phi = acb(0)
        for pole_order, coefficient in self.derivative_coefficients:
            if coefficient.is_zero:
                continue
            phi += (
                (-coefficient / (pole_order - 1)).to_acb()
                * local_point ** (-(pole_order - 1))
            )
        return phi.exp()

    def phi_signature(self) -> tuple[tuple[int, GaussianRational], ...]:
        """返回 ``phi(z)`` 的 canonical exact 负幂系数。"""

        return tuple(
            (-(pole_order - 1), -coefficient / (pole_order - 1))
            for pole_order, coefficient in self.derivative_coefficients
            if not coefficient.is_zero
        )


def _diagonalizing_basis(leading: GaussianMatrix) -> tuple[GaussianMatrix, list[GaussianRational]]:
    """在 Q(i) 中构造最高阶矩阵的完整 exact 本征基。"""

    try:
        roots = leading.charpoly().gaussian_rational_roots()
    except ValueError as error:
        raise LocalReductionError(
            "leading high-pole matrix does not split over Q(i); algebraic field extension is required"
        ) from error
    identity = GaussianMatrix.identity(leading.nrows)
    columns: list[GaussianMatrix] = []
    column_roots: list[GaussianRational] = []
    for root, multiplicity in roots:
        vectors = (leading - identity * root).nullspace()
        if len(vectors) != int(multiplicity):
            raise LocalReductionError(
                "leading high-pole matrix is defective; a Jordan/shearing or ramified formal gauge is required"
            )
        columns.extend(vectors)
        column_roots.extend(root for _ in vectors)
    transformation = _matrix_from_columns(columns, leading.nrows)
    transformation.inverse()
    return transformation, column_roots


def _exponential_power_log_basis(
    system: RationalMatrixSystem,
    center: GaussianRational,
    series_order: int,
) -> LocalSolutionBasis:
    """构造已认证互不耦合 sector 的 ``exp(Phi)`` 乘 power-log 基。"""

    pole_order = system.pole_order_at(center)
    if pole_order < 2:
        raise ValueError("exponential local basis requires a higher-order pole")
    leading = system.local_laurent_matrices_at(center, -pole_order, -pole_order)[
        -pole_order
    ]
    constant_basis, _leading_roots = _diagonalizing_basis(leading)
    inverse_basis = constant_basis.inverse()
    transformed = system.constant_basis_transform(constant_basis)
    high_coefficients = transformed.local_laurent_matrices_at(center, -pole_order, -2)

    for power, matrix in high_coefficients.items():
        for row in range(system.dimension):
            for column in range(system.dimension):
                if row != column and not matrix.scalar(row, column).is_zero:
                    raise LocalReductionError(
                        f"high-order coefficient A_{power} is not diagonal in the exact leading eigenbasis; "
                        "formal block decoupling is required"
                    )

    signatures = [
        tuple(
            high_coefficients[-order].scalar(index, index)
            for order in range(pole_order, 1, -1)
        )
        for index in range(system.dimension)
    ]
    groups: list[tuple[int, ...]] = []
    for index, signature in enumerate(signatures):
        match = next(
            (group for group in groups if signatures[group[0]] == signature),
            None,
        )
        if match is None:
            groups.append((index,))
        else:
            groups[groups.index(match)] = (*match, index)
    group_index = {
        component: group_position
        for group_position, group in enumerate(groups)
        for component in group
    }
    for row in range(system.dimension):
        for column in range(system.dimension):
            if group_index[row] != group_index[column] and not transformed.entries[row][column].is_zero:
                raise LocalReductionError(
                    "different exponential sectors remain coupled by lower-order terms; "
                    "Levelt-Turrittin block decoupling and Stokes data are required"
                )

    coordinate = rational_function((-center, 1))
    sectors: list[_ExponentialSector] = []
    for group in groups:
        derivative_coefficients = tuple(
            (
                order,
                high_coefficients[-order].scalar(group[0], group[0]),
            )
            for order in range(pole_order, 1, -1)
        )
        block_rows = []
        for local_row, row in enumerate(group):
            block_row = []
            for local_column, column in enumerate(group):
                value = transformed.entries[row][column]
                if local_row == local_column:
                    for order, coefficient in derivative_coefficients:
                        if not coefficient.is_zero:
                            value -= rational_function(coefficient) / (coordinate**order)
                block_row.append(value)
            block_rows.append(tuple(block_row))
        residual = RationalMatrixSystem(
            tuple(block_rows),
            variable_name=system.variable_name,
            name=f"{system.name}-exponential-sector-{len(sectors) + 1:03d}",
        )
        residual_order = residual.pole_order_at(center)
        if residual_order > 1:
            raise LocalReductionError(
                f"exponential sector {len(sectors) + 1} retains pole order {residual_order}; "
                "a further meromorphic gauge is required"
            )
        local_basis = _power_log_local_basis(
            residual,
            center,
            series_order,
            allow_ordinary=True,
            method="exponential_sector_power_log",
        )
        sectors.append(_ExponentialSector(group, derivative_coefficients, local_basis))

    constant_acb = constant_basis.to_acb()

    def evaluate(local_point: acb) -> acb_mat:
        transformed_fundamental = acb_mat(system.dimension, system.dimension)
        for sector in sectors:
            block = sector.local_basis.evaluate(local_point)
            exponential = sector.exponential(local_point)
            for local_row, row in enumerate(sector.indices):
                for local_column, column in enumerate(sector.indices):
                    transformed_fundamental[row, column] = exponential * block[
                        local_row, local_column
                    ]
        return constant_acb * transformed_fundamental

    def resolve(value: Any) -> tuple[GaussianMatrix, dict[str, Any]]:
        boundary, phi_signatures = _exponential_power_log_boundary(value)
        if any(len(term.C) != system.dimension for term in boundary.terms):
            raise ValueError(
                f"generalized boundary C dimension must equal differential-equation dimension {system.dimension}"
            )
        constants = GaussianMatrix.zero(system.dimension, 1)
        term_details: list[dict[str, Any]] = []
        for position, term in enumerate(boundary.terms, 1):
            transformed_vector = inverse_basis * _column(term.C)
            active_groups = {
                group_index[index]
                for index in range(system.dimension)
                if not transformed_vector.scalar(index, 0).is_zero
            }
            if len(active_groups) != 1:
                raise ValueError(
                    f"generalized boundary term {position}: C contains multiple exponential sectors; "
                    "split it into one {phi,a,b,C} term per exact sector"
                )
            sector_index = next(iter(active_groups))
            sector = sectors[sector_index]
            if (
                phi_signatures[position - 1] != sector.phi_signature()
            ):
                raise ValueError(
                    f"exponential boundary term {position}: phi does not match exact sector "
                    f"{sector_index + 1}"
                )
            block_vector = tuple(
                transformed_vector.scalar(index, 0) for index in sector.indices
            )
            block_boundary = FrobeniusBoundary(
                (FrobeniusBoundaryTerm(term.a, term.b, block_vector),)
            )
            block_constants, block_report = sector.local_basis.resolve_boundary(block_boundary)
            for local_index, global_index in enumerate(sector.indices):
                constants.real[global_index, 0] += block_constants.real[local_index, 0]
                constants.imag[global_index, 0] += block_constants.imag[local_index, 0]
            term_details.append(
                {
                    "term_index": position,
                    "input": term.to_json(),
                    "verified_sector": sector_index + 1,
                    "component_indices_zero_based": list(sector.indices),
                    "exponential": sector.exponent_records(),
                    "power_log_resolution": block_report,
                }
            )
        if constants.is_zero:
            raise ValueError("combined generalized boundary selects the zero solution")
        report = {
            "schema": "flintnde_generalized_boundary_v1",
            "terms": value.to_json(),
            "constant_basis_exact": constant_basis.to_records(),
            "canonical_basis_constants": [
                str(constants.scalar(row, 0)) for row in range(system.dimension)
            ],
            "term_resolutions": term_details,
        }
        return constants, report

    manifest = {
        "schema": "flintnde_exponential_power_log_v1",
        "method": "exponential_power_log",
        "center_exact": str(center),
        "original_pole_order": pole_order,
        "constant_basis_convention": "I=S J",
        "constant_basis_exact": constant_basis.to_records(),
        "sectors": [
            {
                "sector_index": index,
                "component_indices_zero_based": list(sector.indices),
                "exponential": sector.exponent_records(),
                "maximum_log_degree": sector.local_basis.maximum_log_degree,
                "frobenius": sector.local_basis.manifest["frobenius"],
            }
            for index, sector in enumerate(sectors, 1)
        ],
        "scope": (
            "exactly decoupled Q(i) exponential sectors; general formal block decoupling "
            "and Stokes matching are not inferred"
        ),
    }
    return LocalSolutionBasis(
        system.dimension,
        "exponential_power_log",
        max(sector.local_basis.maximum_log_degree for sector in sectors),
        manifest,
        evaluate,
        resolve,
    )


@dataclass(frozen=True)
class _FormalExponentialBranch:
    """保存一个单重指数根的 exact 形式系数及五阶前瞻诊断。"""

    k: GaussianRational
    rho: GaussianRational
    leading: GaussianMatrix
    coefficients: tuple[GaussianMatrix, ...]
    evaluation_degree: int

    def _term_magnitudes(self, local_point: acb) -> list[arb]:
        """返回各向量项的无穷范数，用于前瞻比值与最小项诊断。"""

        magnitudes: list[arb] = []
        power = acb(1)
        for coefficient in self.coefficients:
            values = [abs(coefficient.scalar(row, 0).to_acb() * power) for row in range(coefficient.nrows)]
            magnitudes.append(max(values, key=lambda value: float(value.mid())))
            power *= local_point
        return magnitudes

    def _partial_series(self, local_point: acb, degree: int) -> acb_mat:
        """计算不含公共指数/幂因子的指定阶局部向量和。"""

        series = acb_mat(self.leading.nrows, 1)
        power = acb(1)
        for coefficient in self.coefficients[: degree + 1]:
            coefficient_acb = coefficient.to_acb()
            for row in range(self.leading.nrows):
                series[row, 0] += coefficient_acb[row, 0] * power
            power *= local_point
        return series

    @staticmethod
    def _vector_norm(vector: acb_mat) -> arb:
        """返回列向量无穷范数。"""

        return max(
            (abs(vector[row, 0]) for row in range(vector.nrows())),
            key=lambda value: float(value.mid()),
        )

    def truncation(self, local_point: acb) -> tuple[int, dict[str, Any]]:
        """固定使用正式阶数，并报告随后五阶的下降趋势与相对 refinement。"""

        magnitudes = self._term_magnitudes(local_point)
        nonzero = [index for index, value in enumerate(magnitudes) if value > 0]
        if not nonzero:
            raise ArithmeticError("formal exponential branch has no nonzero numerical term")
        least_index = min(nonzero, key=lambda index: float(magnitudes[index].mid()))
        maximum_degree = len(self.coefficients) - 1
        terminated = least_index < maximum_degree and all(
            coefficient.is_zero for coefficient in self.coefficients[least_index + 1 :]
        )
        least_term_reached = terminated or least_index < maximum_degree
        primary = self._partial_series(local_point, self.evaluation_degree)
        reference = self._partial_series(local_point, self.evaluation_degree + 5)
        previous_base_degree = self.evaluation_degree - 5
        previous_base = (
            acb_mat(self.leading.nrows, 1)
            if previous_base_degree < 0
            else self._partial_series(local_point, previous_base_degree)
        )
        previous_five_sum = primary - previous_base
        next_five_sum = reference - primary
        previous_five_norm = self._vector_norm(previous_five_sum)
        next_five_norm = self._vector_norm(next_five_sum)
        if previous_five_norm.contains(0):
            five_order_block_ratio = arb(0) if next_five_norm.contains(0) else None
        else:
            five_order_block_ratio = next_five_norm / previous_five_norm
        reference_norm = self._vector_norm(reference)
        relative_refinement = (
            None
            if reference_norm.contains(0)
            else self._vector_norm(reference - primary) / reference_norm
        )
        report = {
            "k_exact": str(self.k),
            "rho_exact": str(self.rho),
            "selected_truncation_degree": self.evaluation_degree,
            "available_maximum_degree": maximum_degree,
            "least_term_degree": least_index,
            "least_term_magnitude_midpoint": f"{float(magnitudes[least_index].mid()):.17e}",
            "last_available_term_magnitude_midpoint": f"{float(magnitudes[-1].mid()):.17e}",
            "least_term_reached": least_term_reached,
            "series_terminated_exactly": terminated,
            "lookahead_order_count": 5,
            "previous_five_sum_norm": previous_five_norm.str(30),
            "next_five_sum_norm": next_five_norm.str(30),
            "five_order_block_ratio_definition": (
                "norm_inf(sum(T[N+1:N+6]))/norm_inf(sum(T[N-4:N+1]))"
            ),
            "next_five_over_previous_five": (
                None if five_order_block_ratio is None else five_order_block_ratio.str(30)
            ),
            "next_five_over_previous_five_below_one": (
                None if five_order_block_ratio is None else bool(five_order_block_ratio < 1)
            ),
            "five_order_relative_refinement": (
                None if relative_refinement is None else relative_refinement.str(30)
            ),
            "trend": (
                "undefined_zero_previous_five_sum"
                if five_order_block_ratio is None
                else (
                    "decreasing_five_order_blocks"
                    if five_order_block_ratio < 1
                    else "nondecreasing_five_order_blocks"
                )
            ),
        }
        return self.evaluation_degree, report

    def evaluate(self, local_point: acb) -> acb_mat:
        """固定按用户要求的 ``N`` 阶计算 ``exp(-k/z) z^rho sum(v_n z^n)``。"""

        selected_degree, _report = self.truncation(local_point)
        series = self._partial_series(local_point, selected_degree)
        common = (-self.k.to_acb() / local_point).exp() * local_point ** self.rho.to_acb()
        return series * common


def _formal_branch_coefficients(
    leading_matrix: GaussianMatrix,
    laurent: dict[int, GaussianMatrix],
    k: GaussianRational,
    leading_vector: GaussianMatrix,
    series_order: int,
    lookahead_order: int = 5,
) -> _FormalExponentialBranch:
    """用左本征相容条件递推单重 rank-one irregular 形式分支。"""

    dimension = leading_matrix.nrows
    identity = GaussianMatrix.identity(dimension)
    recurrence_matrix = identity * k - leading_matrix
    left_vectors = _transpose(recurrence_matrix).nullspace()
    if len(left_vectors) != 1:
        raise LocalReductionError("formal exponential branch does not have a unique left eigenvector")
    left = left_vectors[0]
    pairing = _bilinear(left, leading_vector)
    if pairing.is_zero:
        raise LocalReductionError("left/right leading eigenvectors have zero exact pairing")
    residue = laurent[-1]
    rho = _bilinear(left, residue * leading_vector) / pairing
    coefficients = [leading_vector]

    for degree in range(1, series_order + lookahead_order + 1):
        right = (residue - identity * (rho + degree - 1)) * coefficients[degree - 1]
        for regular_power in range(degree - 1):
            right += laurent[regular_power] * coefficients[degree - 2 - regular_power]
        if not _bilinear(left, right).is_zero:
            raise LocalReductionError(
                f"formal exponential recurrence compatibility fails at degree {degree}"
            )
        particular = _canonical_exact_solve(recurrence_matrix, right)

        # 下一阶相容条件唯一固定当前阶沿 v0 的自由分量。
        lookahead = (residue - identity * (rho + degree)) * particular
        for regular_power in range(degree):
            lookahead += laurent[regular_power] * coefficients[degree - 1 - regular_power]
        alpha = _bilinear(left, lookahead) / (pairing * degree)
        coefficient = particular + leading_vector * alpha
        coefficients.append(coefficient)

        next_right = (residue - identity * (rho + degree)) * coefficient
        for regular_power in range(degree):
            next_right += laurent[regular_power] * coefficients[degree - 1 - regular_power]
        if not _bilinear(left, next_right).is_zero:
            raise ArithmeticError(
                f"formal exponential look-ahead compatibility failed at degree {degree + 1}"
            )

    return _FormalExponentialBranch(
        k,
        rho,
        leading_vector,
        tuple(coefficients),
        series_order,
    )


def _formal_exponential_asymptotic_basis(
    system: RationalMatrixSystem,
    center: GaussianRational,
    series_order: int,
) -> LocalSolutionBasis:
    """构造二阶 pole、单重不同主导根的 start-only 形式渐近基。"""

    if series_order < 4:
        raise ValueError(
            "formal asymptotic evaluation requires series_order at least four for its five-order diagnostic"
        )
    if system.pole_order_at(center) != 2:
        raise LocalReductionError("formal rank-one route requires pole order exactly two")
    leading_matrix = system.local_laurent_matrices_at(center, -2, -2)[-2]
    try:
        roots = leading_matrix.charpoly().gaussian_rational_roots()
    except ValueError as error:
        raise LocalReductionError(
            "formal rank-one leading spectrum does not split over Q(i)"
        ) from error
    if len(roots) != system.dimension or any(multiplicity != 1 for _root, multiplicity in roots):
        raise LocalReductionError(
            "formal rank-one route requires simple pairwise-distinct leading eigenvalues"
        )

    laurent = system.local_laurent_matrices_at(center, -2, series_order + 4)
    branches: list[_FormalExponentialBranch] = []
    for root, _multiplicity in roots:
        vectors = (leading_matrix - GaussianMatrix.identity(system.dimension) * root).nullspace()
        if len(vectors) != 1:
            raise LocalReductionError("simple leading root does not have one exact eigenvector")
        branches.append(
            _formal_branch_coefficients(
                leading_matrix,
                laurent,
                root,
                vectors[0],
                series_order,
            )
        )
    leading_basis = _matrix_from_columns(
        [branch.leading for branch in branches], system.dimension
    )
    leading_basis.inverse()

    def evaluate(local_point: acb) -> acb_mat:
        columns = [branch.evaluate(local_point) for branch in branches]
        return acb_mat(
            [
                [columns[column][row, 0] for column in range(system.dimension)]
                for row in range(system.dimension)
            ]
        )

    def evaluation_report(local_point: acb) -> dict[str, Any]:
        return {
            "schema": "flintnde_formal_asymptotic_evaluation_v1",
            "local_point": local_point.str(40),
            "branch_diagnostics": [branch.truncation(local_point)[1] for branch in branches],
            "interpretation": (
                "formal rank-one series evaluated at the requested order with five-order "
                "look-ahead diagnostics; no convergence radius or Stokes connection is asserted"
            ),
        }

    def match_distance_report(
        target_order: int,
        minimum_order_factor: float,
    ) -> tuple[arb, dict[str, Any]]:
        target_minimum_degree = minimum_order_factor * target_order
        branch_results: list[tuple[arb, dict[str, Any]]] = []
        for branch_index, branch in enumerate(branches, 1):
            gaps = []
            for other_index, other in enumerate(branches, 1):
                if other_index == branch_index:
                    continue
                difference = branch.k - other.k
                gaps.append((abs(difference.to_acb()), other_index, difference))
            if not gaps:
                raise LocalReductionError(
                    "formal root-gap matching requires at least two exponential roots"
                )
            nearest_gap, nearest_index, nearest_difference = min(
                gaps, key=lambda item: float(item[0].mid())
            )
            distance = nearest_gap / target_minimum_degree
            branch_results.append(
                (
                    distance,
                    {
                        "sector_index": branch_index,
                        "k_exact": str(branch.k),
                        "nearest_other_sector_index": nearest_index,
                        "nearest_root_difference_exact": str(nearest_difference),
                        "nearest_root_gap": nearest_gap.str(30),
                        "target_order": target_order,
                        "minimum_order_factor": minimum_order_factor,
                        "estimated_least_term_degree": target_minimum_degree,
                        "suggested_match_distance": distance.str(30),
                    },
                )
            )
        distance = min(branch_results, key=lambda item: float(item[0].mid()))[0]
        return distance, {
            "schema": "flintnde_formal_asymptotic_path_estimate_v1",
            "selection_rule": (
                "minimum nearest-exponential-root gap divided by "
                "minimum_order_factor*target_order"
            ),
            "late_term_model": "n_min approximately abs(k_i-k_j)/abs(z)",
            "target_order": target_order,
            "minimum_order_factor": minimum_order_factor,
            "estimated_least_term_degree": target_minimum_degree,
            "suggested_match_distance": distance.str(30),
            "branches": [report for _branch_distance, report in branch_results],
        }

    def resolve(value: Any) -> tuple[GaussianMatrix, dict[str, Any]]:
        boundary, phi_signatures = _exponential_power_log_boundary(value)
        if any(len(term.C) != system.dimension for term in boundary.terms):
            raise ValueError(
                f"generalized boundary C dimension must equal differential-equation dimension {system.dimension}"
            )
        constants = GaussianMatrix.zero(system.dimension, 1)
        details: list[dict[str, Any]] = []
        for position, term in enumerate(boundary.terms, 1):
            if term.b != 0:
                raise ValueError(
                    f"formal exponential boundary term {position}: current simple sector requires b=0"
                )
            input_vector = _column(term.C)
            matches: list[tuple[int, GaussianRational]] = []
            for branch_index, branch in enumerate(branches):
                pivot = next(
                    (
                        row
                        for row in range(system.dimension)
                        if not branch.leading.scalar(row, 0).is_zero
                    ),
                    None,
                )
                if pivot is None:
                    raise ArithmeticError("formal branch has a zero leading vector")
                scale = input_vector.scalar(pivot, 0) / branch.leading.scalar(pivot, 0)
                if scale.is_zero:
                    continue
                if (branch.leading * scale - input_vector).is_zero and term.a == branch.rho:
                    matches.append((branch_index, scale))
            if len(matches) != 1:
                raise ValueError(
                    f"formal exponential boundary term {position}: {{a,C}} must select exactly one "
                    "exactly one exponential sector"
                )
            branch_index, scale = matches[0]
            expected_phi = () if branches[branch_index].k.is_zero else ((-1, -branches[branch_index].k),)
            if (
                phi_signatures[position - 1] != expected_phi
            ):
                raise ValueError(
                    f"exponential boundary term {position}: phi does not match exact formal "
                    f"sector {branch_index + 1}"
                )
            constants.real[branch_index, 0] += scale.real
            constants.imag[branch_index, 0] += scale.imag
            branch = branches[branch_index]
            details.append(
                {
                    "term_index": position,
                    "input": term.to_json(),
                    "verified_sector": branch_index + 1,
                    "k_exact": str(branch.k),
                    "rho_exact": str(branch.rho),
                    "canonical_scale": str(scale),
                }
            )
        if constants.is_zero:
            raise ValueError("combined formal exponential boundary selects the zero solution")
        return constants, {
            "schema": "flintnde_formal_exponential_boundary_v1",
            "terms": value.to_json(),
            "canonical_basis_constants": [
                str(constants.scalar(row, 0)) for row in range(system.dimension)
            ],
            "term_resolutions": details,
        }

    manifest = {
        "schema": "flintnde_formal_exponential_asymptotic_v1",
        "method": "formal_exponential_asymptotic",
        "center_exact": str(center),
        "original_pole_order": 2,
        "series_order": series_order,
        "leading_matrix_exact": leading_matrix.to_records(),
        "branches": [
            {
                "sector_index": index,
                "k_exact": str(branch.k),
                "rho_exact": str(branch.rho),
                "leading_vector_exact": branch.leading.to_records(),
            }
            for index, branch in enumerate(branches, 1)
        ],
        "scope": (
            "simple distinct Q(i) leading eigenvalues at a rank-one irregular point; "
            "start-only fixed-order evaluation with five-order diagnostics and no Stokes connection"
        ),
    }
    return LocalSolutionBasis(
        system.dimension,
        "formal_exponential_asymptotic",
        0,
        manifest,
        evaluate,
        resolve,
        False,
        evaluation_report,
        match_distance_report,
    )


def build_local_solution_basis(
    system: RationalMatrixSystem,
    point: Any,
    series_order: int,
) -> LocalSolutionBasis:
    """按 regular、Lee--Moser、exponential 的顺序构造可认证局部基。"""

    if series_order <= 0:
        raise ValueError("local solution series order must be positive")
    center = gaussian_rational(point)
    reduction = attempt_fuchsian_reduction(system, center)
    if reduction.status == "already_fuchsian":
        if reduction.original_pole_order != 1:
            raise ValueError("local singular basis requires a singular point")
        return _power_log_local_basis(
            system,
            center,
            series_order,
            allow_ordinary=False,
            method="regular_singular_power_log",
            extra_manifest={"fuchsian_reduction": reduction.to_json()},
        )
    if reduction.status == "reduced_to_fuchsian":
        return _sheared_power_log_basis(reduction, series_order)
    try:
        basis = _exponential_power_log_basis(system, center, series_order)
    except LocalReductionError as exact_sector_error:
        try:
            basis = _formal_exponential_asymptotic_basis(system, center, series_order)
        except LocalReductionError as formal_error:
            raise LocalReductionError(
                f"high-order pole remains after exact Lee-Moser reduction ({reduction.reason}); "
                f"decoupled exponential route is not certified ({exact_sector_error}); "
                f"formal rank-one route is not certified ({formal_error})"
            ) from formal_error
    basis.manifest["fuchsian_reduction_attempt"] = reduction.to_json()
    return basis
