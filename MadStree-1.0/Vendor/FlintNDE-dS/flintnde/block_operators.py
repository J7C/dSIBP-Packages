"""将精确 sector/source 留数编译为整块直积算子，并输运完整的多列解。"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from flint import acb, acb_mat

from .exact_gaussian import GaussianRational


_ZERO = GaussianRational(0)
Layout = tuple[tuple[int, ...], ...]
Entries = dict[tuple[int, int], GaussianRational]


def _remove_bit(index: int, bit: int) -> int:
    return (index & ((1 << bit) - 1)) | ((index >> (bit + 1)) << bit)


def _identity_layout(
    entries: Entries, row_offset: int, row_count: int,
    column_offset: int, column_count: int,
) -> tuple[Layout, Layout, Entries]:
    """逐次证明配对 bit 为单位因子，返回 active 行列与 spectator 布局。"""

    row_bases = list(range(row_count))
    column_bases = list(range(column_count))
    row_bits = list(range(row_count.bit_length() - 1))
    column_bits = list(range(column_count.bit_length() - 1))
    spectators: list[tuple[int, int]] = []
    binary = row_count & (row_count - 1) == 0 and column_count & (column_count - 1) == 0
    while binary and row_bits and column_bits:
        pair = None
        for row_bit in range(len(row_bits)):
            for column_bit in range(len(column_bits)):
                if all(
                    ((row >> row_bit) & 1) == ((column >> column_bit) & 1)
                    and entries.get((row ^ (1 << row_bit), column ^ (1 << column_bit))) == value
                    for (row, column), value in entries.items()
                ):
                    pair = row_bit, column_bit
                    break
            if pair is not None:
                break
        if pair is None:
            break
        row_bit, column_bit = pair
        spectators.append((row_bits.pop(row_bit), column_bits.pop(column_bit)))
        entries = {
            (_remove_bit(row, row_bit), _remove_bit(column, column_bit)): value
            for (row, column), value in entries.items()
            if not (row & (1 << row_bit))
        }
        row_bases = [base for index, base in enumerate(row_bases) if not (index & (1 << row_bit))]
        column_bases = [base for index, base in enumerate(column_bases) if not (index & (1 << column_bit))]

    active_rows = sorted({row for row, _column in entries})
    active_columns = sorted({column for _row, column in entries})
    row_index = {row: index for index, row in enumerate(active_rows)}
    column_index = {column: index for index, column in enumerate(active_columns)}
    row_shifts = tuple(
        sum(((spectator >> bit) & 1) << pair[0] for bit, pair in enumerate(spectators))
        for spectator in range(1 << len(spectators))
    )
    column_shifts = tuple(
        sum(((spectator >> bit) & 1) << pair[1] for bit, pair in enumerate(spectators))
        for spectator in range(1 << len(spectators))
    )
    rows = tuple(tuple(row_offset + row_bases[row] + shift for shift in row_shifts) for row in active_rows)
    columns = tuple(tuple(column_offset + column_bases[column] + shift for shift in column_shifts) for column in active_columns)
    reduced = {(row_index[row], column_index[column]): value for (row, column), value in entries.items()}
    return rows, columns, reduced


@dataclass(frozen=True)
class ActiveBlock:
    """矩阵仅保存 active 下标；布局的第二轴为不混合的 spectator。"""

    rows: Layout
    columns: Layout
    matrix: acb_mat


@dataclass(frozen=True)
class BlockResidue:
    dimension: int
    blocks: tuple[ActiveBlock, ...]

    def nrows(self) -> int:
        return self.dimension

    def ncols(self) -> int:
        return self.dimension

    def to_dense(self) -> acb_mat:
        """仅供显式矩阵接口和奇点处理重建，不进入普通段递推。"""

        dense = acb_mat(self.dimension, self.dimension)
        for block in self.blocks:
            for row, destinations in enumerate(block.rows):
                for column, sources in enumerate(block.columns):
                    value = block.matrix[row, column]
                    for destination, source in zip(destinations, sources):
                        dense[destination, source] += value
        return dense


def compile_block_residue(
    dimension: int, partition: tuple[tuple[int, int], ...], terms: Sequence[Entries],
) -> BlockResidue:
    """在实际 master 分区内编译精确贡献，保留矩形 source 及重叠支撑。"""

    owners = [0] * dimension
    for sector, (offset, count) in enumerate(partition):
        owners[offset:offset + count] = [sector] * count
    merged: dict[Layout, dict[tuple[int, tuple[int, ...]], GaussianRational]] = {}
    for term in terms:
        sector_entries: dict[tuple[int, int], Entries] = {}
        for (row, column), value in term.items():
            if value.is_zero:
                continue
            row_sector, column_sector = owners[row], owners[column]
            row_offset, _ = partition[row_sector]
            column_offset, _ = partition[column_sector]
            sector_entries.setdefault((row_sector, column_sector), {})[
                row - row_offset, column - column_offset
            ] = value
        for (row_sector, column_sector), entries in sector_entries.items():
            rows, columns, reduced = _identity_layout(
                entries, *partition[row_sector], *partition[column_sector]
            )
            bucket = merged.setdefault(rows, {})
            for (row, column), value in reduced.items():
                position = row, columns[column]
                bucket[position] = bucket.get(position, _ZERO) + value
    blocks = []
    for rows, entries in merged.items():
        entries = {position: value for position, value in entries.items() if not value.is_zero}
        if not entries:
            continue
        columns = tuple(sorted({column for _row, column in entries}))
        column_index = {column: index for index, column in enumerate(columns)}
        matrix = acb_mat(len(rows), len(columns))
        for (row, column), value in entries.items():
            matrix[row, column_index[column]] = value.to_acb()
        blocks.append(ActiveBlock(rows, columns, matrix))
    return BlockResidue(dimension, tuple(blocks))


class BlockResidues(Sequence):
    """普通段读取编译块；显式 residues 索引按需返回完整 Acb 矩阵。"""

    def __init__(self, operators: Sequence[BlockResidue]):
        self.operators = tuple(operators)
        grouped: dict[Layout, dict[Layout, list[tuple[int, acb_mat]]]] = {}
        for pole_index, operator in enumerate(self.operators):
            for block in operator.blocks:
                grouped.setdefault(block.rows, {}).setdefault(block.columns, []).append(
                    (pole_index, block.matrix)
                )
        self.groups = tuple(
            (rows, tuple((columns, tuple(actions)) for columns, actions in sources.items()))
            for rows, sources in grouped.items()
        )
        self.source_layouts = tuple(dict.fromkeys(
            columns for _rows, sources in self.groups for columns, _actions in sources
        ))
        outputs = [
            destination
            for rows, _sources in self.groups
            for destinations in rows
            for destination in destinations
        ]
        dimension = self.operators[0].dimension if self.operators else 0
        self.tiled_outputs = sorted(outputs) == list(range(dimension))
        row_plans: list[list[tuple[int, int, int]]] = [[] for _ in range(dimension)]
        for group_index, (rows, _sources) in enumerate(self.groups):
            for row, destinations in enumerate(rows):
                for spectator, destination in enumerate(destinations):
                    row_plans[destination].append((group_index, row, spectator))
        self.row_plans = tuple(tuple(plan) for plan in row_plans)
        self._sum_vectors: dict[int, acb_mat] = {}

    def __len__(self) -> int:
        return len(self.operators)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return tuple(operator.to_dense() for operator in self.operators[index])
        return self.operators[index].to_dense()

    def solution_coefficients(
        self, initial: acb_mat, order: int, offsets: Sequence[acb],
        polynomial: Sequence[acb_mat],
    ) -> list[acb_mat]:
        """按输出状态递推，整块共享打包并在累计完成后写回完整解。"""

        dimension, rhs_count = initial.nrows(), initial.ncols()
        solution = [acb_mat(initial)]
        polynomial_terms = tuple(
            (degree, matrix) for degree, matrix in enumerate(polynomial)
            if any(matrix[row, column] != 0 for row in range(dimension) for column in range(dimension))
        )
        states = [
            [
                [acb_mat(len(rows), len(rows[0]) * rhs_count) for _action in actions]
                for _columns, actions in sources
            ]
            for rows, sources in self.groups
        ]
        full_layout = tuple((row,) for row in range(dimension))
        if not polynomial_terms:
            return self._assembled_solution(
                initial, order, offsets, states, full_layout,
            )
        for degree in range(order):
            current = solution[degree]
            packed = {
                layout: current if layout == full_layout else acb_mat([
                    [current[source, rhs] for source in sources for rhs in range(rhs_count)]
                    for sources in layout
                ])
                for layout in self.source_layouts
            }
            total = acb_mat(dimension, rhs_count)
            for power, matrix in polynomial_terms:
                if power <= degree:
                    total += matrix * solution[degree - power]
            for group_index, (rows, sources) in enumerate(self.groups):
                accumulated = acb_mat(len(rows), len(rows[0]) * rhs_count)
                for source_index, (columns, actions) in enumerate(sources):
                    source = packed[columns]
                    source_states = states[group_index][source_index]
                    for action_index, (pole_index, matrix) in enumerate(actions):
                        # w[n]=R*u[n]，故 w[n]=(w[n-1]-R*c[n])/(p-center)。
                        state = (source_states[action_index] - matrix * source) / offsets[pole_index]
                        source_states[action_index] = state
                        accumulated += state
                if rows == full_layout:
                    total += accumulated
                else:
                    for row, destinations in enumerate(rows):
                        for spectator, destination in enumerate(destinations):
                            for rhs in range(rhs_count):
                                total[destination, rhs] += accumulated[row, spectator * rhs_count + rhs]
            solution.append(total / acb(degree + 1))
        return solution

    def _assembled_solution(
        self, initial: acb_mat, order: int, offsets: Sequence[acb],
        states: list, full_layout: Layout,
    ) -> list[acb_mat]:
        """按全局行汇总各输出块，避免 Python 逐元素 gather/scatter。

        source 打包用整行拼接；每个输出块在 C 层完成乘法和缩放后，单一贡献行
        直接引用，多贡献行用一次 C 层小矩阵求和。acb 算术不经过逐元素接口。
        """

        dimension, rhs_count = initial.nrows(), initial.ncols()
        solution = [acb_mat(initial)]
        for degree in range(order):
            current = solution[degree]
            current_table = None
            packed_cache: dict[Layout, acb_mat] = {}
            divisor = acb(degree + 1)
            scaled_tables = []
            for group_index, (rows, sources) in enumerate(self.groups):
                accumulated = None
                for source_index, (columns, actions) in enumerate(sources):
                    if columns == full_layout:
                        source = current
                    else:
                        source = packed_cache.get(columns)
                        if source is None:
                            if current_table is None:
                                current_table = current.table()
                            source = acb_mat([
                                [entry for index in indices for entry in current_table[index]]
                                for indices in columns
                            ])
                            packed_cache[columns] = source
                    source_states = states[group_index][source_index]
                    for action_index, (pole_index, matrix) in enumerate(actions):
                        state = (source_states[action_index] - matrix * source) / offsets[pole_index]
                        source_states[action_index] = state
                        accumulated = state if accumulated is None else accumulated + state
                scaled_tables.append((accumulated / divisor).table())
            next_rows = []
            for plan in self.row_plans:
                if not plan:
                    # 该输出行在当前留数算子里没有贡献；无多项式时系数严格为零。
                    next_rows.append([acb(0)] * rhs_count)
                    continue
                if len(plan) == 1:
                    group_index, row, spectator = plan[0]
                    next_rows.append(scaled_tables[group_index][row][
                        spectator * rhs_count:(spectator + 1) * rhs_count
                    ])
                    continue
                stack = acb_mat([
                    scaled_tables[group_index][row][
                        spectator * rhs_count:(spectator + 1) * rhs_count
                    ]
                    for group_index, row, spectator in plan
                ])
                ones = self._sum_vectors.get(len(plan))
                if ones is None:
                    ones = acb_mat([[acb(1)] * len(plan)])
                    self._sum_vectors[len(plan)] = ones
                next_rows.append((ones * stack).table()[0])
            solution.append(acb_mat(next_rows))
        return solution
