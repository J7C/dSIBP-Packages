"""高阶 pole 的 exact Lee--Moser projector-balance 降阶。

本模块对 ``I'=A(z)I`` 逐轮构造 ``I=T_Q J``，其中
``T_Q=1-Q+Q/z``、``Q^2=Q``。每一步都在 ``Q(i)`` 上精确计算最高两阶
Laurent 矩阵、nilpotent Jordan 链和 Moser projector，并对变换后的完整有理
矩阵重新检查 pole 阶。所有 projector、顺序和正逆映射都会进入 manifest；算法
适用域之外保持 fail closed，不把搜索失败冒充真实不规则奇点证明。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from flint import acb, acb_mat

from .boundary import FrobeniusBoundaryTerm
from .exact_gaussian import (
    GaussianMatrix,
    GaussianRational,
    exact_matrix_from_columns as _matrix_from_columns,
    gaussian_rational,
)
from .singularities import RationalFunction, RationalMatrixSystem, rational_function


def _local_coordinate(center: GaussianRational) -> RationalFunction:
    """返回 exact 有理函数 ``z=x-center``。"""

    return rational_function((-center, 1))


def _column_rank(columns: list[GaussianMatrix], dimension: int) -> int:
    """返回一组 ``Q(i)`` 列向量的 exact rank。"""

    if not columns:
        return 0
    return _matrix_from_columns(columns, dimension).rref()[1]


def _submatrix(
    matrix: GaussianMatrix,
    rows: list[int],
    columns: list[int],
) -> GaussianMatrix:
    """按给定行列抽取 exact 子矩阵；空列只用于 rank=0 的内部判定。"""

    if not columns:
        return GaussianMatrix.zero(len(rows), 1)
    return GaussianMatrix.from_records(
        [[matrix.scalar(row, column) for column in columns] for row in rows]
    )


def _outer_product(column: GaussianMatrix, row: GaussianMatrix) -> GaussianMatrix:
    """构造 exact 列向量与行向量的外积。"""

    if column.ncols != 1 or row.nrows != 1:
        raise ValueError("outer product requires a column and a row")
    return GaussianMatrix.from_records(
        [
            [column.scalar(i, 0) * row.scalar(0, j) for j in range(row.ncols)]
            for i in range(column.nrows)
        ]
    )


def _row(matrix: GaussianMatrix, index: int) -> GaussianMatrix:
    """返回 exact 单行矩阵。"""

    return GaussianMatrix.from_records(
        [[matrix.scalar(index, column) for column in range(matrix.ncols)]]
    )


def _rational_multiply(
    left: list[list[RationalFunction]],
    right: list[list[RationalFunction]],
) -> list[list[RationalFunction]]:
    """乘两个同维 exact 有理函数矩阵。"""

    dimension = len(left)
    return [
        [
            sum(
                (left[row][inner] * right[inner][column] for inner in range(dimension)),
                rational_function(0),
            )
            for column in range(dimension)
        ]
        for row in range(dimension)
    ]


def _balance_entries(
    projector: GaussianMatrix,
    coordinate: RationalFunction,
    *,
    inverse: bool,
) -> list[list[RationalFunction]]:
    """构造 ``1-Q+z Q`` 或 ``1-Q+Q/z`` 的 exact 矩阵元。"""

    dimension = projector.nrows
    identity = GaussianMatrix.identity(dimension)
    regular = identity - projector
    scale = coordinate if inverse else 1 / coordinate
    return [
        [rational_function(regular.scalar(row, column))
         + scale * projector.scalar(row, column)
         for column in range(dimension)]
        for row in range(dimension)
    ]


def _apply_projector_balance(
    system: RationalMatrixSystem,
    center: GaussianRational,
    projector: GaussianMatrix,
) -> RationalMatrixSystem:
    """施加 ``I=(1-Q+Q/z)J`` 并返回 ``J'`` 的 exact 系统。"""

    coordinate = _local_coordinate(center)
    transformation = _balance_entries(projector, coordinate, inverse=False)
    inverse = _balance_entries(projector, coordinate, inverse=True)
    product = _rational_multiply(
        [list(row) for row in system.entries], transformation
    )
    # ``-T' = Q/z^2``；先加入括号再左乘 T^-1，与 AMFlow Balance 公式一致。
    for row in range(system.dimension):
        for column in range(system.dimension):
            product[row][column] += (
                rational_function(projector.scalar(row, column)) / (coordinate**2)
            )
    transformed = _rational_multiply(inverse, product)
    return RationalMatrixSystem(
        tuple(tuple(row) for row in transformed),
        variable_name=system.variable_name,
        name=f"{system.name}-moser-balance",
    )


def _restore_projector_balance(
    system: RationalMatrixSystem,
    center: GaussianRational,
    projector: GaussianMatrix,
) -> RationalMatrixSystem:
    """把降阶系统反变换回应用当前 balance 之前的原系统。"""

    coordinate = _local_coordinate(center)
    transformation = _balance_entries(projector, coordinate, inverse=False)
    inverse = _balance_entries(projector, coordinate, inverse=True)
    restored = _rational_multiply(
        _rational_multiply(transformation, [list(row) for row in system.entries]),
        inverse,
    )
    derivative_times_inverse = _rational_multiply(
        [
            [
                -rational_function(projector.scalar(row, column)) / (coordinate**2)
                for column in range(system.dimension)
            ]
            for row in range(system.dimension)
        ],
        inverse,
    )
    for row in range(system.dimension):
        for column in range(system.dimension):
            restored[row][column] += derivative_times_inverse[row][column]
    return RationalMatrixSystem(
        tuple(tuple(row) for row in restored),
        variable_name=system.variable_name,
        name=f"{system.name}-inverse-moser-balance",
    )


@dataclass(frozen=True)
class MoserBalanceTransformation:
    """保存有序 Lee--Moser projector balances 及全部正逆映射。"""

    center: GaussianRational
    dimension: int
    projectors: tuple[GaussianMatrix, ...] = ()

    def __post_init__(self) -> None:
        """验证每个 projector 都同维且严格满足 ``Q^2=Q``。"""

        center = gaussian_rational(self.center)
        if self.dimension <= 0:
            raise ValueError("Moser balance dimension must be positive")
        projectors = tuple(self.projectors)
        for projector in projectors:
            if projector.nrows != self.dimension or projector.ncols != self.dimension:
                raise ValueError("Moser projector dimension mismatch")
            if not (projector * projector - projector).is_zero:
                raise ValueError("Moser projector must satisfy Q^2=Q exactly")
        object.__setattr__(self, "center", center)
        object.__setattr__(self, "projectors", projectors)

    @staticmethod
    def _evaluate_factor(local_point: acb, projector: GaussianMatrix, inverse: bool) -> acb_mat:
        """在非零局部点计算一个 balance 或其逆矩阵。"""

        identity = acb_mat(projector.nrows, projector.ncols)
        for index in range(projector.nrows):
            identity[index, index] = 1
        q = projector.to_acb()
        scale = local_point if inverse else 1 / local_point
        return identity - q + q * scale

    def apply_to_system(self, system: RationalMatrixSystem) -> RationalMatrixSystem:
        """按记录顺序施加全部 projector balances。"""

        if system.dimension != self.dimension:
            raise ValueError("Moser balance system dimension mismatch")
        transformed = system
        for projector in self.projectors:
            transformed = _apply_projector_balance(transformed, self.center, projector)
        return transformed

    def restore_system(self, reduced: RationalMatrixSystem) -> RationalMatrixSystem:
        """从最终降阶系统精确恢复原系统，用于独立 round-trip 认证。"""

        restored = reduced
        for projector in reversed(self.projectors):
            restored = _restore_projector_balance(restored, self.center, projector)
        return restored

    def to_reduced(self, local_point: acb, value: acb_mat) -> acb_mat:
        """计算 ``J=T_m^-1 ... T_1^-1 I``。"""

        if local_point.contains(0):
            raise ZeroDivisionError("Moser balance cannot be evaluated at its center")
        result = acb_mat(value)
        for projector in self.projectors:
            result = self._evaluate_factor(local_point, projector, True) * result
        return result

    def to_original(self, local_point: acb, value: acb_mat) -> acb_mat:
        """计算 ``I=T_1 ... T_m J``。"""

        if local_point.contains(0):
            raise ZeroDivisionError("Moser balance cannot be evaluated at its center")
        result = acb_mat(value)
        for projector in reversed(self.projectors):
            result = self._evaluate_factor(local_point, projector, False) * result
        return result

    def transform_boundary_term(
        self, term: FrobeniusBoundaryTerm, position: int
    ) -> tuple[FrobeniusBoundaryTerm, dict[str, Any]]:
        """由 ``T^-1`` 的最低非零幂确定降阶基 leading term。

        每个逆 balance 是 ``1-Q+zQ``，故其乘积是有限矩阵多项式。这里只保留作用
        于用户 ``C`` 后的首个非零系数；随后由降阶系统的 indicial 方程独立验证。
        """

        if len(term.C) != self.dimension:
            raise ValueError("Frobenius boundary C dimension must equal Moser dimension")
        coefficients = [GaussianMatrix.from_records([[value] for value in term.C])]
        identity = GaussianMatrix.identity(self.dimension)
        for projector in self.projectors:
            previous = coefficients
            coefficients = [GaussianMatrix.zero(self.dimension, 1) for _ in range(len(previous) + 1)]
            for degree, vector in enumerate(previous):
                coefficients[degree] = coefficients[degree] + (identity - projector) * vector
                coefficients[degree + 1] = coefficients[degree + 1] + projector * vector
        leading_degree = next(
            (degree for degree, vector in enumerate(coefficients) if not vector.is_zero),
            None,
        )
        if leading_degree is None:
            raise ValueError(f"Frobenius boundary term {position}: Moser inverse annihilated C")
        leading = coefficients[leading_degree]
        transformed = FrobeniusBoundaryTerm(
            term.a + leading_degree,
            term.b,
            tuple(leading.scalar(index, 0) for index in range(self.dimension)),
        )
        return transformed, {
            "term_index": position,
            "original": term.to_json(),
            "inverse_balance_leading_degree": leading_degree,
            "inverse_balance_vector_coefficients": [
                [str(vector.scalar(index, 0)) for index in range(self.dimension)]
                for vector in coefficients
            ],
            "reduced_basis": transformed.to_json(),
        }

    def laurent_coefficients_exact(self) -> dict[int, GaussianMatrix]:
        """返回累计 ``T_1...T_m`` 的有限 exact Laurent 系数。

        该表示用于把降阶 Frobenius 级数完整乘回原基后验证用户 ``{a,b,C}``，
        避免只变换 leading vector 时漏掉相邻阶贡献。
        """

        coefficients = {0: GaussianMatrix.identity(self.dimension)}
        identity = GaussianMatrix.identity(self.dimension)
        for projector in self.projectors:
            updated: dict[int, GaussianMatrix] = {}
            for degree, coefficient in coefficients.items():
                updated[degree] = updated.get(
                    degree, GaussianMatrix.zero(self.dimension, self.dimension)
                ) + coefficient * (identity - projector)
                updated[degree - 1] = updated.get(
                    degree - 1, GaussianMatrix.zero(self.dimension, self.dimension)
                ) + coefficient * projector
            coefficients = {degree: value for degree, value in updated.items() if not value.is_zero}
        return coefficients

    def to_json(self) -> dict[str, Any]:
        """返回有序 projectors 和 exact idempotency 认证。"""

        return {
            "route": "exact_lee_moser_projector_balance",
            "convention": "I=Product_k(1-Q_k+Q_k/z) J",
            "center_exact": str(self.center),
            "balance_count": len(self.projectors),
            "projectors_exact": [projector.to_records() for projector in self.projectors],
            "projector_idempotency_exact": [True for _ in self.projectors],
        }


@dataclass(frozen=True)
class FuchsianReductionResult:
    """保存降阶状态、变换后的系统和完整验证信息。"""

    status: str
    original_pole_order: int
    reduced_pole_order: int
    transformed_system: RationalMatrixSystem
    transformation: MoserBalanceTransformation
    reason: str
    pole_order_history: tuple[int, ...] = ()

    @property
    def succeeded(self) -> bool:
        """返回变换后是否已经至多为 simple pole。"""

        return self.status in {"already_fuchsian", "reduced_to_fuchsian"}

    def to_json(self) -> dict[str, Any]:
        """序列化路线、逐步 pole 阶和失败原因。"""

        return {
            "schema": "flintnde_fuchsian_reduction_v2",
            "status": self.status,
            "original_pole_order": self.original_pole_order,
            "reduced_pole_order": self.reduced_pole_order,
            "pole_order_history": list(self.pole_order_history),
            "reason": self.reason,
            "transformation": self.transformation.to_json(),
        }


def _nilpotent_jordan_chains(
    leading: GaussianMatrix,
) -> tuple[GaussianMatrix, tuple[int, ...]] | None:
    """构造最高阶 nilpotent 矩阵的 exact Jordan-chain 基及块长。

    链内列按 ``N c_j=c_(j-1)`` 排列，块长从大到小。返回前严格检查
    ``P^-1 N P`` 是否等于对应的零特征值 Jordan 矩阵。
    """

    if leading.is_zero or leading.nrows != leading.ncols:
        return None
    dimension = leading.nrows
    nilpotency_index = next(
        (power for power in range(1, dimension + 1) if (leading**power).is_zero),
        None,
    )
    if nilpotency_index is None:
        return None

    columns: list[GaussianMatrix] = []
    lengths: list[int] = []
    rank = 0
    for length in range(nilpotency_index, 0, -1):
        for head in (leading**length).nullspace():
            chain = [leading**power * head for power in range(length - 1, -1, -1)]
            candidate_rank = _column_rank([*columns, *chain], dimension)
            if candidate_rank != rank + length:
                continue
            columns.extend(chain)
            lengths.append(length)
            rank = candidate_rank
            if rank == dimension:
                basis = _matrix_from_columns(columns, dimension)
                inverse = basis.inverse()
                jordan = GaussianMatrix.zero(dimension, dimension)
                offset = 0
                records = jordan.to_records()
                for block_length in lengths:
                    for index in range(1, block_length):
                        records[offset + index - 1][offset + index] = "1"
                    offset += block_length
                jordan = GaussianMatrix.from_records(records)
                if not (inverse * leading * basis - jordan).is_zero:
                    raise ArithmeticError("nilpotent Jordan-chain substitution check failed")
                return basis, tuple(lengths)
    return None


def _reduce_l0(
    l0: GaussianMatrix,
    long_block_count: int,
    block_lengths: tuple[int, ...],
) -> tuple[int, tuple[int, ...], GaussianMatrix]:
    """Exact 移植 AMFlow ``ReduceL0`` 的 block-chain 消元。

    该步骤只在首个可选链不属于长 Jordan 链时使用 ``A1`` 调整链基。每次 nullspace
    必须严格一维，否则当前 Moser 实现返回 inconclusive。
    """

    block_count = len(block_lengths)
    selected: list[int] = []
    delta = GaussianMatrix.zero(block_count, block_count)
    current = l0
    candidate = 0
    guard = 0
    while candidate >= long_block_count:
        remaining = [index for index in range(block_count) if index not in selected]
        dependent = None
        for index in remaining:
            before_columns = list(range(index))
            through_columns = list(range(index + 1))
            before_rank = 0 if not before_columns else _submatrix(current, remaining, before_columns).rref()[1]
            through_rank = _submatrix(current, remaining, through_columns).rref()[1]
            if before_rank == through_rank:
                dependent = index
                break
        if dependent is None:
            raise ArithmeticError("Moser L0 reduction found no dependent block column")
        candidate = dependent
        kernel = _submatrix(current, remaining, list(range(candidate + 1))).nullspace()
        if len(kernel) != 1:
            raise ArithmeticError("Moser L0 reduction requires a one-dimensional nullspace")
        relation = kernel[0]
        last = relation.scalar(candidate, 0)
        if last.is_zero:
            raise ArithmeticError("Moser L0 null vector has a zero normalization entry")
        relation = relation / (-last)
        delta0_records = [[0 for _ in range(block_count)] for _ in range(block_count)]
        delta0t_records = [[0 for _ in range(block_count)] for _ in range(block_count)]
        for index in range(candidate):
            value = -relation.scalar(index, 0)
            delta0_records[index][candidate] = value
            if block_lengths[index] == block_lengths[candidate]:
                delta0t_records[index][candidate] = value
        delta0 = GaussianMatrix.from_records(delta0_records)
        delta0t = GaussianMatrix.from_records(delta0t_records)
        identity = GaussianMatrix.identity(block_count)
        current = (identity - delta0t) * current * (identity + delta0)
        delta = delta + delta0 + delta * delta0
        selected.append(candidate)
        guard += 1
        if guard > block_count:
            raise ArithmeticError("Moser L0 reduction exceeded its finite block bound")
    return candidate, tuple(index for index in selected if index != candidate), delta


def _find_moser_projector(
    leading: GaussianMatrix,
    subleading: GaussianMatrix,
) -> GaussianMatrix:
    """由 ``A0,A1`` exact 构造一次 Lee--Moser projector。

    实现对应 AMFlow ``FindProjector``：先取零特征值 Jordan chains，再在 block 空间
    计算 ``L0`` 和必要的链间调整，最后回到原基构造 idempotent projector。
    """

    decomposition = _nilpotent_jordan_chains(leading)
    if decomposition is None:
        raise ArithmeticError("highest Laurent coefficient is not nonzero nilpotent")
    basis, block_lengths = decomposition
    inverse = basis.inverse()
    starts: list[int] = []
    offset = 0
    for length in block_lengths:
        starts.append(offset)
        offset += length
    block_count = len(block_lengths)
    l0_records: list[list[GaussianRational]] = []
    l1_records: list[list[GaussianRational]] = []
    for left_block, left_start in enumerate(starts):
        left_index = left_start + block_lengths[left_block] - 1
        l0_row: list[GaussianRational] = []
        l1_row: list[GaussianRational] = []
        left_vector = _row(inverse, left_index)
        for right_start in starts:
            right_vector = basis.column(right_start)
            l0_row.append((left_vector * subleading * right_vector).scalar(0, 0))
            l1_row.append((left_vector * right_vector).scalar(0, 0))
        l0_records.append(l0_row)
        l1_records.append(l1_row)
    l0 = GaussianMatrix.from_records(l0_records)
    l1 = GaussianMatrix.from_records(l1_records)
    long_block_count = sum(
        1 for index in range(block_count) if l1.scalar(index, index).is_zero
    )
    k0, selected, delta = _reduce_l0(l0, long_block_count, block_lengths)

    et_records = [[0 for _ in range(leading.nrows)] for _ in range(leading.nrows)]
    for left_block in range(block_count):
        for right_block in range(left_block + 1, block_count):
            coefficient = delta.scalar(left_block, right_block)
            if coefficient.is_zero:
                continue
            common = min(block_lengths[left_block], block_lengths[right_block])
            for position in range(common):
                et_records[starts[left_block] + position][starts[right_block] + position] = coefficient
    chain_adjustment = GaussianMatrix.identity(leading.nrows) + GaussianMatrix.from_records(et_records)
    adjusted_basis = basis * chain_adjustment
    adjusted_inverse = chain_adjustment.inverse() * inverse
    projector = GaussianMatrix.zero(leading.nrows, leading.nrows)
    for block in (*selected, k0):
        projector = projector + _outer_product(
            adjusted_basis.column(starts[block]),
            _row(adjusted_inverse, starts[block]),
        )
    if projector.is_zero or not (projector * projector - projector).is_zero:
        raise ArithmeticError("constructed Moser projector failed exact idempotency")
    return projector


def _systems_equal(left: RationalMatrixSystem, right: RationalMatrixSystem) -> bool:
    """逐元检查两个 exact 有理矩阵系统完全相等。"""

    if left.dimension != right.dimension:
        return False
    return all(
        (left.entries[row][column] - right.entries[row][column]).is_zero
        for row in range(left.dimension)
        for column in range(left.dimension)
    )


def attempt_fuchsian_reduction(
    system: RationalMatrixSystem,
    point: Any,
) -> FuchsianReductionResult:
    """执行 exact Lee--Moser projector-balance 降阶并完整认证结果。

    只有每一步 pole 阶严格下降且最终至多单极点时返回 ``reduced_to_fuchsian``。
    非 nilpotent 最高阶项返回 exact ``moser_irreducible``；构造假设未满足、projector
    失败或 round-trip 不成立则返回 ``reduction_inconclusive``。
    """

    center = gaussian_rational(point)
    original_order = system.pole_order_at(center)
    identity = MoserBalanceTransformation(center, system.dimension)
    if original_order <= 1:
        return FuchsianReductionResult(
            "already_fuchsian", original_order, original_order, system, identity,
            "input system already has at most a simple pole", (original_order,),
        )

    current = system
    projectors: list[GaussianMatrix] = []
    history = [original_order]
    maximum_steps = max(1, system.dimension * original_order)
    reason = ""
    for _step in range(maximum_steps):
        order = current.pole_order_at(center)
        if order <= 1:
            break
        coefficients = current.local_laurent_matrices_at(
            center, -order, -order + 1
        )
        leading = coefficients[-order]
        subleading = coefficients[-order + 1]
        if _nilpotent_jordan_chains(leading) is None:
            transformation = MoserBalanceTransformation(
                center, system.dimension, tuple(projectors)
            )
            return FuchsianReductionResult(
                "moser_irreducible", original_order, order, current, transformation,
                "highest Laurent coefficient is not nilpotent; a rational Fuchsian gauge "
                "does not exist in the unramified meromorphic Lee-Moser class",
                tuple(history),
            )
        try:
            projector = _find_moser_projector(leading, subleading)
            candidate = _apply_projector_balance(current, center, projector)
        except (ArithmeticError, ValueError, ZeroDivisionError) as error:
            transformation = MoserBalanceTransformation(
                center, system.dimension, tuple(projectors)
            )
            return FuchsianReductionResult(
                "reduction_inconclusive", original_order, order, current, transformation,
                f"exact Lee-Moser projector construction failed: {error}", tuple(history),
            )
        candidate_order = candidate.pole_order_at(center)
        if candidate_order >= order:
            transformation = MoserBalanceTransformation(
                center, system.dimension, tuple(projectors)
            )
            return FuchsianReductionResult(
                "reduction_inconclusive", original_order, order, current, transformation,
                "exact projector balance did not strictly lower the complete matrix pole order",
                tuple(history),
            )
        projectors.append(projector)
        current = candidate
        history.append(candidate_order)
    else:
        reason = "exact Lee-Moser iteration exceeded its finite safety bound"

    transformation = MoserBalanceTransformation(center, system.dimension, tuple(projectors))
    reduced_order = current.pole_order_at(center)
    if reduced_order <= 1:
        restored = transformation.restore_system(current)
        if not _systems_equal(restored, system):
            return FuchsianReductionResult(
                "reduction_inconclusive", original_order, reduced_order, current,
                transformation, "exact Moser transformation failed system round-trip check",
                tuple(history),
            )
        return FuchsianReductionResult(
            "reduced_to_fuchsian", original_order, reduced_order, current, transformation,
            "exact Lee-Moser projector balances reduced and round-trip certified the system",
            tuple(history),
        )
    return FuchsianReductionResult(
        "reduction_inconclusive", original_order, reduced_order, current, transformation,
        reason or "exact Lee-Moser reduction did not reach Fuchsian form", tuple(history),
    )
