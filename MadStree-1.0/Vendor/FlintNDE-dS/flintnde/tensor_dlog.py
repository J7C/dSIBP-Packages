"""MadStree 张量积 dlog 留数的精确重建。

MadStree 只提交每个 sector 的局部基标签、quotient 稀疏映射和按字母的结构化
贡献；本模块在 FlintNDE 内部按 slot 轴作用重建 U 与 U^{-1}，再计算

    R = postFactor * U^{-1} diag(indicator) U * (preVector 或 middle)

并写回 master 基。局部基的公共 1/sqrt(2) 因子在 U^{-1} D U 中整体抵消，所以
全程只用整数高斯矩阵，留数保持精确 Q(i)，不会引入任何代数数近似。

约定与 MadStree 的 msTensorApplyU 完全一致：state 编号是
``IntegerDigits[m-1, 2, nSlots]``，slot 1 为最高位；U = P (kron locals) S。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable

from .exact_gaussian import GaussianRational


_ZERO = GaussianRational(0)
_ONE = GaussianRational(1)
_MINUS_I = GaussianRational(0, -1)
_HALF = _ONE / 2
_I_HALF = GaussianRational(0, 1) / 2

# msPaperT / msHadamard / msIdentity2 去掉公共 1/sqrt(2) 因子后的整数高斯矩阵。
LOCAL_BASES: dict[str, tuple[tuple[GaussianRational, ...], ...]] = {
    "massiveT": (
        (_ONE, _MINUS_I),
        (_MINUS_I, _ONE),
    ),
    "masslessH": (
        (_ONE, _ONE),
        (_ONE, -_ONE),
    ),
    "identity": (
        (_ONE, _ZERO),
        (_ZERO, _ONE),
    ),
}

# 同一组基的逐元素逆；整数矩阵的逆带 1/2 因子，仍落在 Q(i) 内。
LOCAL_BASES_INVERSE: dict[str, tuple[tuple[GaussianRational, ...], ...]] = {
    "massiveT": (
        (_HALF, _I_HALF),
        (_I_HALF, _HALF),
    ),
    "masslessH": (
        (_HALF, _HALF),
        (_HALF, -_HALF),
    ),
    "identity": (
        (_ONE, _ZERO),
        (_ZERO, _ONE),
    ),
}


@dataclass(frozen=True)
class TensorSector:
    """一个 MadStree sector 的张量因子数据。"""

    sector_key: str
    dimension: int
    offset: int
    raw_count: int
    slot_bases: tuple[str, ...]
    raw_to_master_column: tuple[int, ...]
    master_representative_raw: tuple[int, ...]

    def __post_init__(self) -> None:
        """检查标签、维数和 quotient 映射的自洽性。"""

        if not isinstance(self.sector_key, str):
            raise ValueError("tensor sector requires a string sectorKey")
        # 无根传播子的 sector 键是空串，属于合法取值。
        if self.dimension < 1 or self.raw_count < 1 or self.offset < 1:
            raise ValueError(f"{self.sector_key}: dimension/rawCount/offset must be positive")
        if self.raw_count != 1 << len(self.slot_bases):
            raise ValueError(
                f"{self.sector_key}: rawCount must equal 2**slotCount"
            )
        unknown = set(self.slot_bases) - set(LOCAL_BASES)
        if unknown:
            raise ValueError(f"{self.sector_key}: unknown local bases {sorted(unknown)}")
        if len(self.raw_to_master_column) != self.raw_count:
            raise ValueError(f"{self.sector_key}: rawToMasterColumn length mismatch")
        for column in self.raw_to_master_column:
            if not 0 <= column <= self.dimension:
                raise ValueError(f"{self.sector_key}: rawToMasterColumn out of range")
        if len(self.master_representative_raw) != self.dimension:
            raise ValueError(f"{self.sector_key}: masterRepresentativeRaw length mismatch")
        for raw in self.master_representative_raw:
            if not 1 <= raw <= self.raw_count:
                raise ValueError(f"{self.sector_key}: masterRepresentativeRaw out of range")
        if any(
            self.raw_to_master_column[raw - 1] != column
            for column, raw in enumerate(self.master_representative_raw, start=1)
        ):
            raise ValueError(
                f"{self.sector_key}: representative raw state does not map to its column"
            )

    @property
    def first_global_index(self) -> int:
        """返回该 sector 在全局 master 顺序中的首个 0-based 下标。"""

        return self.offset - 1


@dataclass(frozen=True)
class TensorContribution:
    """一条结构化 dlog 贡献。"""

    kind: str
    scope: str
    sector_key: str | None = None
    vector: tuple[GaussianRational, ...] | None = None
    pre_vector: tuple[GaussianRational, ...] | None = None
    letter_indicator: tuple[int, ...] | None = None
    middle_rows: tuple[int, ...] = ()
    middle_columns: tuple[int, ...] = ()
    middle_values: tuple[GaussianRational, ...] = ()
    post_factor: GaussianRational = _ZERO


def apply_local_chain(
    locals_: Iterable[tuple[tuple[GaussianRational, ...], ...]],
    vector: list[GaussianRational],
) -> list[GaussianRational]:
    """按 slot 轴依次作用局部 2x2 矩阵，复杂度 O(n 2^n)。"""

    matrices = list(locals_)
    rank = len(matrices)
    if rank == 0:
        return vector
    if len(vector) != 1 << rank:
        raise ValueError("tensor chain dimension mismatch")
    for axis in range(rank):
        step = 1 << (rank - 1 - axis)
        block = step << 1
        matrix = matrices[axis]
        m00, m01 = matrix[0]
        m10, m11 = matrix[1]
        for base in range(0, len(vector), block):
            for offset in range(step):
                first = base + offset
                second = first + step
                lower = vector[first]
                upper = vector[second]
                vector[first] = m00 * lower + m01 * upper
                vector[second] = m10 * lower + m11 * upper
    return vector


@dataclass
class TensorDLogSystem:
    """把结构化贡献组装成 master 基上的精确留数。

    每个 sector 的 U 只在第一次使用时 materialize，并且只到该 sector 自身的
    master 维数（MadStree 中不超过 2^slots）；全局 2^n 维矩阵从不出现。
    """

    sectors: tuple[TensorSector, ...]
    dimension: int
    _by_key: dict[str, TensorSector] = field(init=False, repr=False)
    _u: dict[str, list[list[GaussianRational]]] = field(init=False, repr=False)
    _u_inverse: dict[str, list[list[GaussianRational]]] = field(init=False, repr=False)
    _blocks: dict[tuple[str, tuple[int, ...]], list[list[GaussianRational]]] = field(
        init=False, repr=False
    )

    def __post_init__(self) -> None:
        """检查 sector 覆盖全局 master 顺序且互不重叠。"""

        self._by_key = {}
        for sector in self.sectors:
            if sector.sector_key in self._by_key:
                raise ValueError(f"duplicate tensor sector {sector.sector_key}")
            self._by_key[sector.sector_key] = sector
        ordered = sorted(self.sectors, key=lambda item: item.offset)
        expected = 1
        for sector in ordered:
            if sector.offset != expected:
                raise ValueError("tensor sectors must tile the global master order")
            expected += sector.dimension
        if expected - 1 != self.dimension:
            raise ValueError("tensor sector dimensions do not match the master count")
        self._u = {}
        self._u_inverse = {}
        self._blocks = {}

    def _sector(self, key: str) -> TensorSector:
        try:
            return self._by_key[key]
        except KeyError as error:
            raise ValueError(f"unknown tensor sector {key}") from error

    def _apply(self, sector: TensorSector, vector: list[GaussianRational], inverse: bool) -> list[GaussianRational]:
        """U 或 U^{-1} 作用于 sector 内 master 基向量。"""

        raw = [
            0 if column == 0 else vector[column - 1]
            for column in sector.raw_to_master_column
        ]
        raw = [value if isinstance(value, GaussianRational) else _ZERO for value in raw]
        bases = LOCAL_BASES_INVERSE if inverse else LOCAL_BASES
        raw = apply_local_chain((bases[tag] for tag in sector.slot_bases), raw)
        return [raw[raw_index - 1] for raw_index in sector.master_representative_raw]

    def sector_u(self, key: str) -> list[list[GaussianRational]]:
        """返回 U 的 dense 表示（仅该 sector 自身维数）。"""

        if key not in self._u:
            sector = self._sector(key)
            columns = [
                self._apply(sector, [_ONE if row == column else _ZERO for row in range(sector.dimension)], False)
                for column in range(sector.dimension)
            ]
            self._u[key] = [
                [columns[column][row] for column in range(sector.dimension)]
                for row in range(sector.dimension)
            ]
        return self._u[key]

    def sector_u_inverse(self, key: str) -> list[list[GaussianRational]]:
        """返回 U^{-1} 的 dense 表示（仅该 sector 自身维数）。"""

        if key not in self._u_inverse:
            sector = self._sector(key)
            columns = [
                self._apply(sector, [_ONE if row == column else _ZERO for row in range(sector.dimension)], True)
                for column in range(sector.dimension)
            ]
            self._u_inverse[key] = [
                [columns[column][row] for column in range(sector.dimension)]
                for row in range(sector.dimension)
            ]
        return self._u_inverse[key]

    def conjugated_block(self, key: str, indicator: tuple[int, ...]) -> list[list[GaussianRational]]:
        """返回 K = U^{-1} diag(indicator) U。"""

        cache_key = (key, indicator)
        if cache_key not in self._blocks:
            u = self.sector_u(key)
            inverse = self.sector_u_inverse(key)
            result = [[_ZERO] * len(u) for _ in range(len(u))]
            # 只有 indicator 为 1 的行参与，先做 U 的行筛选再乘 U^{-1}。
            active = [index for index, bit in enumerate(indicator) if bit]
            for row in range(len(u)):
                inverse_row = inverse[row]
                accumulator = [_ZERO] * len(u)
                for source in active:
                    coefficient = inverse_row[source]
                    if coefficient.is_zero:
                        continue
                    basis_row = u[source]
                    for column in range(len(u)):
                        value = basis_row[column]
                        if not value.is_zero:
                            accumulator[column] = accumulator[column] + coefficient * value
                result[row] = accumulator
            self._blocks[cache_key] = result
        return self._blocks[cache_key]

    def residue_entries(
        self,
        contributions: Iterable[TensorContribution],
    ) -> dict[tuple[int, int], GaussianRational]:
        """把一组结构化贡献累加成 master 基上的稀疏留数。"""

        entries: dict[tuple[int, int], GaussianRational] = {}
        for contribution in contributions:
            if contribution.kind == "diagonal":
                if contribution.scope == "sector":
                    sector = self._sector(contribution.sector_key or "")
                    offset = sector.first_global_index
                else:
                    sector = None
                    offset = 0
                for index, value in enumerate(contribution.vector or ()):
                    if not value.is_zero:
                        position = (offset + index, offset + index)
                        entries[position] = entries.get(position, _ZERO) + value
                continue
            if contribution.kind != "conjugated":
                raise ValueError(f"unsupported tensor contribution kind {contribution.kind}")
            sector = self._sector(contribution.sector_key or "")
            block = self.conjugated_block(
                sector.sector_key, tuple(contribution.letter_indicator or ())
            )
            offset = sector.first_global_index
            post = contribution.post_factor
            if contribution.scope == "sector":
                pre_vector = contribution.pre_vector or ()
                if len(pre_vector) != sector.dimension:
                    raise ValueError(f"{sector.sector_key}: preVector dimension mismatch")
                for column, pre in enumerate(pre_vector):
                    factor = post * pre
                    if factor.is_zero:
                        continue
                    for row in range(sector.dimension):
                        value = block[row][column]
                        if value.is_zero:
                            continue
                        position = (offset + row, offset + column)
                        entries[position] = entries.get(position, _ZERO) + factor * value
            elif contribution.scope == "global":
                for local_row, column, value in zip(
                    contribution.middle_rows,
                    contribution.middle_columns,
                    contribution.middle_values,
                ):
                    factor = post * value
                    if factor.is_zero:
                        continue
                    for row in range(sector.dimension):
                        block_value = block[row][local_row]
                        if block_value.is_zero:
                            continue
                        position = (offset + row, column)
                        entries[position] = entries.get(position, _ZERO) + factor * block_value
            else:
                raise ValueError(f"unsupported tensor contribution scope {contribution.scope}")
        return {position: value for position, value in entries.items() if not value.is_zero}
