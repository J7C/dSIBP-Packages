#!/usr/bin/env python3
"""用 python-flint/Acb 求解 pole/residue 形式的 Step3 NDEs。

输入系统必须已经写成
``Y'(k)=(C+Sum_p R_p/(k-p))Y(k)``，并完成物理参数与 convention 数值代入；
模块完成 Watson 边界、普通点 Taylor 输运、exact-manifest-gated 的 ``k=0``
generalized Frobenius 匹配和 finite part。数值矩阵只使用 FLINT AcbMat；SymPy
仅解析 exact InputForm 和控制严格根/共振标签，不参与数值矩阵运算。

公开接口：``FlintPoleResidueSystem.from_payload`` 读取系统，
``build_positive_real_fractional_path`` 生成固定 ``step/R`` 的正实轴路径，
``transport_path`` 完成指定阶数的逐段输运，``solve_k0_finite_part`` 读取 exact
log manifest 后提取 ``k=0`` finite part。
"""

from __future__ import annotations

import math
import re
import time
from dataclasses import dataclass
from typing import Any, Iterable

import sympy as sp
from flint import acb, acb_mat, acb_series, arb, fmpq

from .core import (
    acb_midpoint_matrix,
    column_vector,
    configure_working_precision,
    identity_matrix,
    matrix_norm_inf,
    relative_difference_inf,
    vector_norm_inf,
)


def _rational_to_fmpq(value: sp.Expr) -> fmpq:
    """把有理实数或 InputForm 十进制转换为 FLINT fmpq。"""

    rational = sp.cancel(value)
    if isinstance(rational, sp.Float):
        rational = sp.Rational(str(rational))
    if not isinstance(rational, sp.Rational):
        raise ValueError(f"expected exact rational scalar, got {rational}")
    return fmpq(int(rational.p), int(rational.q))


def _normalize_mathematica_numeric_inputform(value: Any) -> str:
    """删除 Mathematica 数值精度标记，并转换科学计数法供 SymPy 读取。"""

    text = str(value)
    text = re.sub(r"`(?:\d+(?:\.\d*)?|\.\d+)?", "", text)
    return text.replace("*^", "e")


def exact_inputform_to_acb(value: Any) -> acb:
    """把 exporter 的有理实数或复数 InputForm 记录精确转换为 Acb。"""

    expression = exact_inputform_expression(value)
    return acb(
        _rational_to_fmpq(sp.re(expression)),
        _rational_to_fmpq(sp.im(expression)),
    )


def exact_inputform_expression(value: Any) -> sp.Expr:
    """恢复 InputForm 数值表达式，只用于严格结构和标签控制。"""

    normalized = _normalize_mathematica_numeric_inputform(value)
    expression = sp.expand(sp.sympify(normalized, locals={"I": sp.I}))
    if not bool(expression.is_number):
        raise ValueError(f"non-numeric exact scalar: {value}")
    return expression


def arb_midpoint_float(value: arb) -> float:
    """仅为排序和 JSON 诊断读取 Arb 中点；正式计算保留完整 ball。"""

    return float(value.mid())


def acb_record(value: acb, digits: int = 35) -> dict[str, str]:
    """把 Acb 数值保存为包含误差半径的可读记录。"""

    return {
        "re_ball": str(value.real.str(digits)),
        "im_ball": str(value.imag.str(digits)),
        "abs_ball": str(abs(value).str(digits)),
    }


def acb_value_record(value: acb, digits: int = 70) -> dict[str, str]:
    """保存 Acb 中点与 ball，供正式结果和精度审计同时使用。"""

    return {
        "re": value.real.str(digits, radius=False),
        "im": value.imag.str(digits, radius=False),
        "abs": abs(value).str(digits, radius=False),
        "re_ball": value.real.str(digits),
        "im_ball": value.imag.str(digits),
    }


def arb_record(value: arb, digits: int = 35) -> str:
    """把 Arb 数值保存为包含误差半径的字符串。"""

    return str(value.str(digits))


def vector_value_records(vector: acb_mat, digits: int = 70) -> list[dict[str, str]]:
    """序列化 Acb 列向量的全部分量。"""

    if vector.ncols() != 1:
        raise ValueError("vector_value_records requires a column vector")
    return [acb_value_record(vector[row, 0], digits) for row in range(vector.nrows())]


def minimum_vector_accuracy_digits(vector: acb_mat) -> float:
    """返回非零分量中最小 Acb 相对准确位数。"""

    bits = [
        vector[row, 0].rel_accuracy_bits()
        for row in range(vector.nrows())
        if not abs(vector[row, 0]).contains(0)
    ]
    return math.inf if not bits else min(bits) / math.log2(10)


@dataclass(frozen=True)
class FractionalPathSegment:
    """记录一段按当前最近 pole 半径取步长的正实轴路径。"""

    index: int
    start: acb
    target: acb
    nearest_pole: acb
    radius: arb
    step: arb
    step_over_radius: arb


@dataclass(frozen=True)
class K0LeadingLogCoefficient:
    """保存会影响 ``k -> 0`` 极限的实际 ``k^0 log(k)^q`` 系数向量。"""

    absolute_power_exact: str
    log_degree: int
    coefficient_vector: acb_mat


class FlintPoleResidueSystem:
    """数值矩阵 ``Omega(k)=C+Sum_p R_p/(k-p)`` 的 Acb 表示。"""

    def __init__(
        self,
        *,
        request_id: str,
        dimension: int,
        poles: list[acb],
        constant: acb_mat,
        residues: list[acb_mat],
        zero_pole_position: int,
    ) -> None:
        if constant.nrows() != dimension or constant.ncols() != dimension:
            raise ValueError(f"{request_id}: constant matrix dimension mismatch")
        if len(poles) != len(residues):
            raise ValueError(f"{request_id}: pole/residue count mismatch")
        if any(matrix.nrows() != dimension or matrix.ncols() != dimension for matrix in residues):
            raise ValueError(f"{request_id}: residue matrix dimension mismatch")
        self.request_id = request_id
        self.dimension = dimension
        self.poles = poles
        self.constant = constant
        self.residues = residues
        self.zero_pole_position = zero_pole_position

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "FlintPoleResidueSystem":
        """读取含精确 pole、常数矩阵和 residue 矩阵的系统记录。"""

        required = {"request_id", "dimension", "poles_exact", "constant_matrix_exact", "residue_matrices_exact"}
        missing = sorted(required.difference(payload))
        if missing:
            raise ValueError(f"pole/residue payload is missing fields: {missing}")
        cache: dict[str, acb] = {}

        def convert(value: Any) -> acb:
            """缓存重复 exact 标量，减少大矩阵解析开销。"""

            key = str(value)
            if key not in cache:
                cache[key] = exact_inputform_to_acb(value)
            return cache[key]

        def convert_matrix(records: list[list[Any]]) -> acb_mat:
            """把二维 exact 记录转换为 Acb 稠密矩阵。"""

            return acb_mat([[convert(value) for value in row] for row in records])

        pole_expressions = [
            exact_inputform_expression(value) for value in payload["poles_exact"]
        ]
        zero_positions = [
            index for index, value in enumerate(pole_expressions) if sp.simplify(value) == 0
        ]
        if len(zero_positions) != 1:
            raise ValueError(f"{payload['request_id']}: expected exactly one k=0 pole")
        return cls(
            request_id=str(payload["request_id"]),
            dimension=int(payload["dimension"]),
            poles=[convert(value) for value in payload["poles_exact"]],
            constant=convert_matrix(payload["constant_matrix_exact"]),
            residues=[convert_matrix(matrix) for matrix in payload["residue_matrices_exact"]],
            zero_pole_position=zero_positions[0],
        )

    def taylor_matrix_coefficients(self, center: acb, solution_order: int) -> list[acb_mat]:
        """生成递推到 ``Y_solution_order`` 所需的 ``A_0..A_(order-1)``。"""

        if solution_order <= 0:
            raise ValueError("solution order must be positive")
        inverse_offsets = [acb(1) / (center - pole) for pole in self.poles]
        factors = list(inverse_offsets)
        coefficients: list[acb_mat] = []
        for degree in range(solution_order):
            coefficient = acb_mat(self.constant) if degree == 0 else acb_mat(
                self.dimension, self.dimension
            )
            for residue, factor in zip(self.residues, factors):
                coefficient += residue * factor
            coefficients.append(coefficient)
            factors = [-factor * inverse for factor, inverse in zip(factors, inverse_offsets)]
        return coefficients

    def laurent_at_zero(self, order: int) -> tuple[acb_mat, list[acb_mat]]:
        """生成 ``R0/k + Sum_m A_m k^m`` 的 Acb Laurent/Taylor 系数。"""

        if order <= 0:
            raise ValueError("Laurent regular order must be positive")
        residue_zero = acb_mat(self.residues[self.zero_pole_position])
        regular: list[acb_mat] = []
        for degree in range(order):
            coefficient = acb_mat(self.constant) if degree == 0 else acb_mat(
                self.dimension, self.dimension
            )
            for index, (pole, residue) in enumerate(zip(self.poles, self.residues)):
                if index == self.zero_pole_position:
                    continue
                coefficient -= residue / (pole ** (degree + 1))
            regular.append(coefficient)
        return residue_zero, regular


@dataclass(frozen=True)
class RegulatorRationalMatrixTerm:
    """保存一个 ``matrix/(k-pole)^pole_order`` 精确数值项。"""

    pole: acb
    pole_order: int
    matrix: acb_mat


class FlintRegulatorRationalSystem:
    """保存 ``A(k,rho)=Sum_p rho^p A_p(k)`` 的截断有理矩阵系统。

    每个 ``A_p`` 允许常数矩阵和任意阶固定 pole；因此能够表示移动 pole
    ``1/(k+c*rho)`` 在 ``rho=0`` 展开后产生的高阶 ``k`` pole。
    """

    def __init__(
        self,
        *,
        request_id: str,
        dimension: int,
        constants_by_regulator: list[acb_mat],
        terms_by_regulator: list[list[RegulatorRationalMatrixTerm]],
    ) -> None:
        if not constants_by_regulator or len(constants_by_regulator) != len(terms_by_regulator):
            raise ValueError("regulator system requires equally sized nonempty coefficient lists")
        if any(
            matrix.nrows() != dimension or matrix.ncols() != dimension
            for matrix in constants_by_regulator
        ):
            raise ValueError(f"{request_id}: regulator constant matrix dimension mismatch")
        for terms in terms_by_regulator:
            for term in terms:
                if term.pole_order <= 0:
                    raise ValueError(f"{request_id}: pole order must be positive")
                if term.matrix.nrows() != dimension or term.matrix.ncols() != dimension:
                    raise ValueError(f"{request_id}: regulator pole matrix dimension mismatch")
        self.request_id = request_id
        self.dimension = dimension
        self.constants_by_regulator = constants_by_regulator
        self.terms_by_regulator = terms_by_regulator

    @property
    def regulator_order(self) -> int:
        """返回当前截断的最高 regulator 幂次。"""

        return len(self.constants_by_regulator) - 1

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "FlintRegulatorRationalSystem":
        """读取 exact 字符串 payload，并转换为 FLINT Acb 系数矩阵。"""

        required = {
            "request_id", "dimension", "regulator_order", "regulator_coefficients"
        }
        missing = sorted(required.difference(payload))
        if missing:
            raise ValueError(f"regulator payload is missing fields: {missing}")
        dimension = int(payload["dimension"])
        regulator_order = int(payload["regulator_order"])
        records = payload["regulator_coefficients"]
        if len(records) != regulator_order + 1:
            raise ValueError(
                "regulator coefficient coverage does not match regulator_order"
            )
        cache: dict[str, acb] = {}

        def convert(value: Any) -> acb:
            """缓存 exact 标量转换，避免大矩阵重复解析。"""

            key = str(value)
            if key not in cache:
                cache[key] = exact_inputform_to_acb(value)
            return cache[key]

        def convert_matrix(values: list[list[Any]]) -> acb_mat:
            """转换一个二维 exact 矩阵记录。"""

            return acb_mat([[convert(value) for value in row] for row in values])

        constants: list[acb_mat] = []
        terms_by_regulator: list[list[RegulatorRationalMatrixTerm]] = []
        for expected_power, record in enumerate(records):
            if int(record["regulator_power"]) != expected_power:
                raise ValueError(
                    "regulator coefficient records are not contiguous and ordered"
                )
            constants.append(convert_matrix(record["constant_matrix_exact"]))
            terms_by_regulator.append(
                [
                    RegulatorRationalMatrixTerm(
                        pole=convert(term["pole_exact"]),
                        pole_order=int(term["pole_order"]),
                        matrix=convert_matrix(term["matrix_exact"]),
                    )
                    for term in record["rational_terms"]
                ]
            )
        return cls(
            request_id=str(payload["request_id"]),
            dimension=dimension,
            constants_by_regulator=constants,
            terms_by_regulator=terms_by_regulator,
        )

    def taylor_matrix_coefficients(
        self,
        center: acb,
        solution_order: int,
    ) -> list[list[acb_mat]]:
        """生成 ``A_{m,p}``，返回布局为 ``[regulator_power][taylor_degree]``。"""

        if solution_order <= 0:
            raise ValueError("solution order must be positive")
        result: list[list[acb_mat]] = []
        for constant, terms in zip(
            self.constants_by_regulator, self.terms_by_regulator
        ):
            coefficients = [
                acb_mat(constant) if degree == 0 else acb_mat(
                    self.dimension, self.dimension
                )
                for degree in range(solution_order)
            ]
            for term in terms:
                offset = center - term.pole
                if abs(offset).contains(0):
                    raise ZeroDivisionError(
                        f"{self.request_id}: Taylor center overlaps a pole"
                    )
                for degree in range(solution_order):
                    factor = acb(
                        (-1) ** degree
                        * math.comb(term.pole_order + degree - 1, degree)
                    ) / (offset ** (term.pole_order + degree))
                    coefficients[degree] += term.matrix * factor
            result.append(coefficients)
        return result

    def laurent_matrix_coefficients_at_zero(
        self,
        regular_order: int,
    ) -> dict[tuple[int, int], acb_mat]:
        """生成 ``A_{j,p}``，键为 ``(regulator_power, k_power)``。

        零点 pole 保留其实际负幂；其它固定 pole 在 ``k=0`` 展成 Taylor 级数。
        ``regular_order`` 表示每个 regulator 阶需保留到的最高非负 ``k`` 次数。
        """

        if regular_order < 0:
            raise ValueError("k0 regular order must be nonnegative")
        result: dict[tuple[int, int], acb_mat] = {}
        for regulator_power, (constant, terms) in enumerate(
            zip(self.constants_by_regulator, self.terms_by_regulator)
        ):
            result[(regulator_power, 0)] = acb_mat(constant)
            for term in terms:
                if abs(term.pole).contains(0):
                    key = (regulator_power, -term.pole_order)
                    result[key] = result.get(
                        key, acb_mat(self.dimension, self.dimension)
                    ) + term.matrix
                    continue
                offset = -term.pole
                for degree in range(regular_order + 1):
                    factor = acb(
                        (-1) ** degree
                        * math.comb(term.pole_order + degree - 1, degree)
                    ) / (offset ** (term.pole_order + degree))
                    key = (regulator_power, degree)
                    result[key] = result.get(
                        key, acb_mat(self.dimension, self.dimension)
                    ) + term.matrix * factor
        return result


def regulator_vector_taylor_coefficients(
    matrix_coefficients_by_regulator: list[list[acb_mat]],
    initial_vectors_by_regulator: list[acb_mat],
) -> list[list[acb_mat]]:
    """递推普通点双级数，返回 ``[taylor_degree][regulator_power]``。

    第 ``q`` 个 regulator 系数只读取 ``0..q`` 阶矩阵和解，因此递推严格下三角；
    输入矩阵布局必须是 ``[regulator_power][taylor_degree]``。
    """

    if not matrix_coefficients_by_regulator or not initial_vectors_by_regulator:
        raise ValueError("regulator Taylor recurrence requires nonempty coefficient lists")
    regulator_count = len(initial_vectors_by_regulator)
    if len(matrix_coefficients_by_regulator) < regulator_count:
        raise ValueError(
            "matrix regulator order is lower than the requested solution order"
        )
    taylor_order = len(matrix_coefficients_by_regulator[0])
    if taylor_order <= 0 or any(
        len(coefficients) != taylor_order
        for coefficients in matrix_coefficients_by_regulator[:regulator_count]
    ):
        raise ValueError("matrix Taylor coefficient coverage is inconsistent")
    dimension = initial_vectors_by_regulator[0].nrows()
    if any(
        vector.ncols() != 1 or vector.nrows() != dimension
        for vector in initial_vectors_by_regulator
    ):
        raise ValueError(
            "initial regulator coefficients must be equal-sized column vectors"
        )
    solution: list[list[acb_mat]] = [
        [acb_mat(vector) for vector in initial_vectors_by_regulator]
    ]
    for degree in range(taylor_order):
        next_coefficients: list[acb_mat] = []
        for regulator_power in range(regulator_count):
            total = acb_mat(dimension, 1)
            for matrix_degree in range(degree + 1):
                for matrix_regulator_power in range(regulator_power + 1):
                    total += (
                        matrix_coefficients_by_regulator[matrix_regulator_power][matrix_degree]
                        * solution[degree - matrix_degree][
                            regulator_power - matrix_regulator_power
                        ]
                    )
            next_coefficients.append(total / acb(degree + 1))
        solution.append(next_coefficients)
    return solution


def evaluate_regulator_vector_taylor_series(
    coefficients_by_taylor: list[list[acb_mat]],
    delta: acb,
) -> list[acb_mat]:
    """在一个普通点步长上计算全部 regulator 系数向量。"""

    if not coefficients_by_taylor or not coefficients_by_taylor[0]:
        raise ValueError("regulator Taylor coefficient table must not be empty")
    regulator_count = len(coefficients_by_taylor[0])
    values = [
        acb_mat(coefficients_by_taylor[0][regulator_power].nrows(), 1)
        for regulator_power in range(regulator_count)
    ]
    for coefficients in reversed(coefficients_by_taylor):
        if len(coefficients) != regulator_count:
            raise ValueError("regulator Taylor coefficient table is ragged")
        values = [
            value * delta + coefficient
            for value, coefficient in zip(values, coefficients)
        ]
    return values


def transport_regulator_one_step(
    system: FlintRegulatorRationalSystem,
    initial_vectors_by_regulator: list[acb_mat],
    start: acb,
    target: acb,
    order: int,
) -> tuple[list[acb_mat], float]:
    """用普通点双级数完成一段解析 regulator 输运。"""

    clock = time.perf_counter()
    matrix_coefficients = system.taylor_matrix_coefficients(start, order)
    solution_coefficients = regulator_vector_taylor_coefficients(
        matrix_coefficients,
        initial_vectors_by_regulator,
    )
    values = evaluate_regulator_vector_taylor_series(
        solution_coefficients,
        target - start,
    )
    return values, time.perf_counter() - clock


def transport_regulator_path(
    system: FlintRegulatorRationalSystem,
    initial_vectors_by_regulator: list[acb_mat],
    path: list[acb],
    order: int,
) -> tuple[list[list[acb_mat]], list[dict[str, Any]], float]:
    """沿普通点路径输运全部 regulator 系数并记录逐段耗时。"""

    if len(path) < 2:
        raise ValueError("regulator transport path needs at least two points")
    vectors = [acb_mat(vector) for vector in initial_vectors_by_regulator]
    snapshots = [[acb_mat(vector) for vector in vectors]]
    reports: list[dict[str, Any]] = []
    total_clock = time.perf_counter()
    for index, (start, target) in enumerate(zip(path[:-1], path[1:]), start=1):
        vectors, elapsed = transport_regulator_one_step(
            system,
            vectors,
            start,
            target,
            order,
        )
        vectors = [acb_midpoint_matrix(vector) for vector in vectors]
        snapshots.append([acb_mat(vector) for vector in vectors])
        reports.append(
            {
                "segment_index": index,
                "start": acb_record(start),
                "target": acb_record(target),
                "order": order,
                "regulator_order": system.regulator_order,
                "elapsed_seconds": elapsed,
            }
        )
    return snapshots, reports, time.perf_counter() - total_clock


def _block_identity(dimension: int, block_count: int, block: int) -> acb_mat:
    """构造只在指定常数块含单位阵的宽矩阵。"""

    return acb_mat(
        [
            [
                acb(1 if column == block * dimension + row else 0)
                for column in range(block_count * dimension)
            ]
            for row in range(dimension)
        ]
    )


@dataclass
class RegulatorK0NilpotentLocalBasis:
    """保存单根、二阶 nilpotent residue 下的双变量局部基本矩阵。"""

    dimension: int
    regulator_order: int
    series_order: int
    maximum_log_degree: int
    coefficients: dict[tuple[int, int, int], acb_mat]
    compatibility_residuals: list[arb]

    @property
    def constant_count(self) -> int:
        """返回所有 regulator 阶齐次常数的总数。"""

        return self.dimension * (self.regulator_order + 1)

    def evaluate(self, point: acb) -> acb_mat:
        """在非零 ``k`` 点计算堆叠后的局部基本矩阵。"""

        if abs(point).contains(0):
            raise ZeroDivisionError("k0 local basis cannot be evaluated at k=0")
        logarithm = point.log()
        rows: list[list[acb]] = []
        for regulator_power in range(self.regulator_order + 1):
            value = acb_mat(self.dimension, self.constant_count)
            n_min = -regulator_power
            n_max = self.series_order + self.regulator_order - regulator_power
            for k_power in range(n_min, n_max + 1):
                for log_degree in range(self.maximum_log_degree + 1):
                    coefficient = self.coefficients.get(
                        (regulator_power, k_power, log_degree)
                    )
                    if coefficient is not None:
                        value += coefficient * (point**k_power) * (logarithm**log_degree)
            rows.extend(
                [value[row, column] for column in range(self.constant_count)]
                for row in range(self.dimension)
            )
        return acb_mat(rows)

    def finite_part_matrix(self, regulator_power: int) -> acb_mat:
        """返回指定 regulator 阶的 ``k^0 log(k)^0`` 系数矩阵。"""

        if not 0 <= regulator_power <= self.regulator_order:
            raise ValueError("finite-part regulator power is out of range")
        return acb_mat(
            self.coefficients[(regulator_power, 0, 0)]
        )


def build_regulator_k0_nilpotent_local_basis(
    system: FlintRegulatorRationalSystem,
    *,
    solution_regulator_order: int,
    series_order: int,
) -> tuple[RegulatorK0NilpotentLocalBasis, float]:
    """递推 ``rho^q k^n log(k)^r`` 局部基本矩阵。

    该接口针对真实 16 维奇异层的单 indicial 根与 ``R^2=0`` 情形。每个高阶
    regulator 分量包含由低阶驱动的特解，以及属于本阶的独立齐次常数；在 ``n=0``
    采用“特解常数项为零”的规范，缺陷自动送入更高 ``log(k)`` 次数。
    """

    if not 0 <= solution_regulator_order <= system.regulator_order:
        raise ValueError("requested regulator order exceeds the DE payload")
    if series_order <= 0:
        raise ValueError("k0 local series order must be positive")
    clock = time.perf_counter()
    maximum_log_degree = solution_regulator_order + 1
    dimension = system.dimension
    block_count = solution_regulator_order + 1
    width = dimension * block_count
    matrix_coefficients = system.laurent_matrix_coefficients_at_zero(
        series_order + solution_regulator_order
    )
    residue = matrix_coefficients.get((0, -1))
    if residue is None:
        raise ValueError("regulator k0 recurrence requires the p=0 simple-pole residue")
    nilpotency_error = matrix_norm_inf(residue * residue)
    if arb_midpoint_float(nilpotency_error) > 1.0e-60:
        raise ValueError("current regulator k0 recurrence requires an R^2=0 residue")
    identity = identity_matrix(dimension)
    coefficients: dict[tuple[int, int, int], acb_mat] = {}
    compatibility_residuals: list[arb] = []

    def zero_wide() -> acb_mat:
        """生成一个局部递推使用的零宽矩阵。"""

        return acb_mat(dimension, width)

    def right_hand_side(
        regulator_power: int,
        k_power: int,
        log_degree: int,
    ) -> acb_mat:
        """收集除 ``A_{-1,0}`` 外的已知卷积项。"""

        total = zero_wide()
        for (matrix_regulator_power, matrix_k_power), matrix in matrix_coefficients.items():
            if matrix_regulator_power > regulator_power or (
                matrix_regulator_power == 0 and matrix_k_power == -1
            ):
                continue
            source_key = (
                regulator_power - matrix_regulator_power,
                k_power - 1 - matrix_k_power,
                log_degree,
            )
            source = coefficients.get(source_key)
            if source is not None:
                total += matrix * source
        return total

    for regulator_power in range(block_count):
        n_min = -regulator_power
        n_max = series_order + solution_regulator_order - regulator_power
        for k_power in range(n_min, n_max + 1):
            if k_power == 0:
                coefficients[(regulator_power, 0, 0)] = _block_identity(
                    dimension, block_count, regulator_power
                )
                for log_degree in range(maximum_log_degree):
                    current = coefficients[(regulator_power, 0, log_degree)]
                    right = right_hand_side(regulator_power, 0, log_degree)
                    coefficients[(regulator_power, 0, log_degree + 1)] = (
                        right + residue * current
                    ) / acb(log_degree + 1)
                final = (
                    -residue * coefficients[(regulator_power, 0, maximum_log_degree)]
                    - right_hand_side(regulator_power, 0, maximum_log_degree)
                )
                compatibility_residuals.append(matrix_norm_inf(final))
                continue
            operator = identity * acb(k_power) - residue
            next_log = zero_wide()
            for log_degree in range(maximum_log_degree, -1, -1):
                right = right_hand_side(regulator_power, k_power, log_degree)
                current = operator.solve(
                    right - next_log * acb(log_degree + 1)
                )
                coefficients[(regulator_power, k_power, log_degree)] = current
                next_log = current
    basis = RegulatorK0NilpotentLocalBasis(
        dimension=dimension,
        regulator_order=solution_regulator_order,
        series_order=series_order,
        maximum_log_degree=maximum_log_degree,
        coefficients=coefficients,
        compatibility_residuals=compatibility_residuals,
    )
    return basis, time.perf_counter() - clock


def zero_series(order: int) -> list[acb]:
    """生成指定阶数的 Acb 零级数。"""

    return [acb(0) for _ in range(order + 1)]


def cu_qnm_bc_coefficients(
    iw: acb,
    lam: acb,
    order: int,
    *,
    it0: int,
    it1: int,
) -> list[acb]:
    """递推 QNM horizon-in 分支对应的统一 ``u`` 正则系数。"""

    ah = 1 + 2 * it1 * iw
    ax = -1 + 2 * it0
    bh = 1 - it1 * iw - 2 * iw * iw - 2 * lam + it0 + 2 * it1 * iw * it0
    bx = -1 + it1 * iw + 2 * lam - it0 - 2 * it1 * iw * it0
    a_regular = [acb(ax * ((-1) ** degree)) for degree in range(order + 1)]
    b_regular = [acb(bx * ((-1) ** degree)) for degree in range(order + 1)]
    b_regular[0] -= iw * iw
    coefficients = [acb(1)]
    for degree in range(order):
        known = bh * coefficients[degree]
        for lower in range(degree):
            known += a_regular[lower] * (degree - lower) * coefficients[degree - lower]
            known += b_regular[lower] * coefficients[degree - 1 - lower]
        denominator = (degree + 1) * (degree + ah)
        if denominator.contains(0):
            raise ArithmeticError(f"QNM boundary recurrence denominator contains zero at {degree}")
        coefficients.append(-known / denominator)
    return coefficients


def _derivative_series(coefficients: list[acb]) -> list[acb]:
    """把 ``u`` 正则级数转换为 ``u'`` 级数。"""

    return [degree * coefficients[degree] for degree in range(1, len(coefficients))]


def _integer_binomial_series(power: int, order: int) -> list[acb]:
    """用整数递推生成 ``(1+t)^power``。"""

    coefficients = [acb(1)]
    for degree in range(1, order + 1):
        coefficients.append(coefficients[-1] * (power - degree + 1) / degree)
    return coefficients


def _multiply_series(series_list: Iterable[list[acb]], order: int) -> list[acb]:
    """逐次卷积并截断多个 Acb 单变量级数。"""

    result = [acb(1)] + [acb(0) for _ in range(order)]
    for series in series_list:
        product = zero_series(order)
        for left in range(order + 1):
            for right in range(min(len(series), order + 1 - left)):
                product[left + right] += result[left] * series[right]
        result = product
    return result


def watson_boundary_terms(
    system_payload: dict[str, Any],
    numeric_point: dict[str, Any],
    *,
    k_start: acb,
    order: int,
    it0: int,
    it1: int,
) -> list[list[acb]]:
    """构造一个 request 全部分量的逐阶 Watson 边界项。"""

    iws = [exact_inputform_to_acb(numeric_point[f"iw{leg}"]) for leg in range(1, 4)]
    cls = [exact_inputform_to_acb(numeric_point[f"cl2{leg}"]) for leg in range(1, 4)]
    leg_u = [
        cu_qnm_bc_coefficients(iw, lam, order + 2, it0=it0, it1=it1)
        for iw, lam in zip(iws, cls)
    ]
    leg_du = [_derivative_series(coefficients) for coefficients in leg_u]
    z_value = acb(1) / k_start
    iw_sum = sum(iws, acb(0))
    all_terms: list[list[acb]] = []
    for key in system_payload["physical_basis_keys"]:
        n_values = key[:3]
        a_x, a_x_minus_one = key[3:5]
        z_powers = key[5:]
        active_legs = [index for index, power in enumerate(z_powers) if power != 0]
        if len(active_legs) > 1:
            raise ValueError(f"multiple active Z factors are unsupported: {key}")
        x_factor = _integer_binomial_series(-a_x, order)
        if not active_legs:
            z_factor = [acb(1)] + [acb(0) for _ in range(order)]
        else:
            active_leg = active_legs[0]
            z_power = z_powers[active_leg]
            z_horizon = 3 + 2 * cls[active_leg]
            ratio = 2 * cls[active_leg] / z_horizon
            z_factor = [
                coefficient * ratio**degree * z_horizon ** (-z_power)
                for degree, coefficient in enumerate(
                    _integer_binomial_series(-z_power, order)
                )
            ]
        modes = [
            leg_u[leg] if n_values[leg] == 0 else leg_du[leg]
            for leg in range(3)
        ]
        endpoint = _multiply_series([x_factor, z_factor, *modes], order)
        sigma = it1 * iw_sum - a_x_minus_one
        all_terms.append(
            [
                endpoint[degree]
                * (sigma + degree + 1).gamma()
                * z_value ** (sigma + degree + 1)
                for degree in range(order + 1)
            ]
        )
    return all_terms


def watson_boundary_pair(
    system_payload: dict[str, Any],
    numeric_point: dict[str, Any],
    *,
    k_start: acb,
    primary_order: int,
    reference_order: int,
    it0: int,
    it1: int,
) -> tuple[acb_mat, acb_mat, arb, float]:
    """一次构造主阶/参考阶边界向量并返回相对差和耗时。"""

    if not 0 <= primary_order < reference_order:
        raise ValueError("boundary primary/reference orders are invalid")
    clock = time.perf_counter()
    terms = watson_boundary_terms(
        system_payload,
        numeric_point,
        k_start=k_start,
        order=reference_order,
        it0=it0,
        it1=it1,
    )
    primary = column_vector([sum(row[: primary_order + 1], acb(0)) for row in terms])
    reference = column_vector([sum(row, acb(0)) for row in terms])
    return primary, reference, relative_difference_inf(primary, reference), time.perf_counter() - clock


def _constant_acb_series(value: acb, order: int) -> acb_series:
    """把 Acb 常数提升为指定截断阶数的 regulator 级数。"""

    return acb_series([value], order + 1)


def _cu_qnm_bc_regulator_series(
    iw: acb_series,
    lam: acb,
    taylor_order: int,
    regulator_order: int,
    *,
    it0: int,
    it1: int,
) -> list[acb_series]:
    """递推同时保留 horizon Taylor 阶与解析 regulator 阶的单腿 QNM 解。"""

    ah = 1 + 2 * it1 * iw
    ax = -1 + 2 * it0
    bh = 1 - it1 * iw - 2 * iw * iw - 2 * lam + it0 + 2 * it1 * iw * it0
    bx = -1 + it1 * iw + 2 * lam - it0 - 2 * it1 * iw * it0
    a_regular = [
        _constant_acb_series(acb(ax * ((-1) ** degree)), regulator_order)
        for degree in range(taylor_order + 1)
    ]
    b_regular = [bx * ((-1) ** degree) for degree in range(taylor_order + 1)]
    b_regular[0] -= iw * iw
    coefficients = [_constant_acb_series(acb(1), regulator_order)]
    for degree in range(taylor_order):
        known = bh * coefficients[degree]
        for lower in range(degree):
            known += a_regular[lower] * (degree - lower) * coefficients[degree - lower]
            known += b_regular[lower] * coefficients[degree - 1 - lower]
        denominator = (degree + 1) * (degree + ah)
        if abs(denominator[0]).contains(0):
            raise ArithmeticError(
                f"analytic-regulator QNM boundary denominator contains zero at {degree}"
            )
        coefficients.append(-known / denominator)
    return coefficients


def _multiply_taylor_series_over_regulator(
    factors: Iterable[list[acb_series]],
    taylor_order: int,
    regulator_order: int,
) -> list[acb_series]:
    """卷积 horizon Taylor 级数，同时由 ``acb_series`` 截断 regulator 阶。"""

    result = [_constant_acb_series(acb(1), regulator_order)] + [
        _constant_acb_series(acb(0), regulator_order) for _ in range(taylor_order)
    ]
    for factor in factors:
        product = [
            _constant_acb_series(acb(0), regulator_order)
            for _ in range(taylor_order + 1)
        ]
        for left in range(taylor_order + 1):
            for right in range(min(len(factor), taylor_order + 1 - left)):
                product[left + right] += result[left] * factor[right]
        result = product
    return result


def vector_taylor_coefficients(
    matrix_coefficients: list[acb_mat],
    initial_vector: acb_mat,
) -> list[acb_mat]:
    """按 ``Y'=A Y`` 递推一个当前解向量，不构造 fundamental matrix。"""

    if initial_vector.ncols() != 1:
        raise ValueError("initial value must be a column vector")
    coefficients = [acb_mat(initial_vector)]
    for degree in range(len(matrix_coefficients)):
        total = acb_mat(initial_vector.nrows(), 1)
        for matrix_degree in range(degree + 1):
            total += matrix_coefficients[matrix_degree] * coefficients[degree - matrix_degree]
        coefficients.append(total / acb(degree + 1))
    return coefficients


def evaluate_vector_series(coefficients: list[acb_mat], delta: acb) -> acb_mat:
    """用 Horner 法计算列向量 Taylor 级数。"""

    if not coefficients:
        raise ValueError("Taylor coefficient list must not be empty")
    value = acb_mat(coefficients[0].nrows(), 1)
    for coefficient in reversed(coefficients):
        value = value * delta + coefficient
    return value


def transport_one_step(
    system: FlintPoleResidueSystem,
    initial_vector: acb_mat,
    start: acb,
    target: acb,
    order: int,
) -> tuple[acb_mat, float]:
    """用指定阶数完成一段普通点 Taylor 输运并返回墙钟时间。"""

    start_clock = time.perf_counter()
    matrix_coefficients = system.taylor_matrix_coefficients(start, order)
    vector_coefficients = vector_taylor_coefficients(matrix_coefficients, initial_vector)
    result = evaluate_vector_series(vector_coefficients, target - start)
    return result, time.perf_counter() - start_clock


def transport_path(
    system: FlintPoleResidueSystem,
    initial_vector: acb_mat,
    path: list[acb],
    order: int,
) -> tuple[list[acb_mat], list[dict[str, Any]], float]:
    """沿给定普通点路径逐段输运，返回全部快照和逐段耗时。"""

    if len(path) < 2:
        raise ValueError("transport path needs at least two points")
    vector = acb_mat(initial_vector)
    snapshots = [acb_mat(vector)]
    reports: list[dict[str, Any]] = []
    total_clock = time.perf_counter()
    for index, (start, target) in enumerate(zip(path[:-1], path[1:]), start=1):
        vector, elapsed = transport_one_step(system, vector, start, target, order)
        # 这里执行任意精度点算术；累计截断误差由独立阶数链估计，不把 Acb 外包球跨段传播。
        vector = acb_midpoint_matrix(vector)
        snapshots.append(acb_mat(vector))
        reports.append(
            {
                "segment_index": index,
                "start": acb_record(start),
                "target": acb_record(target),
                "order": order,
                "elapsed_seconds": elapsed,
            }
        )
    return snapshots, reports, time.perf_counter() - total_clock


def transport_path_refined(
    system: FlintPoleResidueSystem,
    primary_initial: acb_mat,
    reference_initial: acb_mat,
    path: list[acb],
    *,
    primary_order: int,
    reference_order: int,
) -> tuple[list[acb_mat], list[acb_mat], list[dict[str, Any]], float]:
    """独立推进主阶和参考阶分支，记录完整累计 refinement 误差。"""

    if not 0 < primary_order < reference_order:
        raise ValueError("transport primary/reference orders are invalid")
    primary = acb_mat(primary_initial)
    reference = acb_mat(reference_initial)
    primary_snapshots = [acb_mat(primary)]
    reference_snapshots = [acb_mat(reference)]
    reports: list[dict[str, Any]] = []
    total_clock = time.perf_counter()
    for index, (start, target) in enumerate(zip(path[:-1], path[1:]), start=1):
        segment_clock = time.perf_counter()
        primary, primary_seconds = transport_one_step(
            system, primary, start, target, primary_order
        )
        reference, reference_seconds = transport_one_step(
            system, reference, start, target, reference_order
        )
        # 两条链分别重取中点，保持其独立性并消除区间包裹效应。
        primary = acb_midpoint_matrix(primary)
        reference = acb_midpoint_matrix(reference)
        try:
            relative = relative_difference_inf(primary, reference)
        except ZeroDivisionError as error:
            denominator = vector_norm_inf(reference)
            raise ArithmeticError(
                "transport refinement lost a nonzero reference norm at "
                f"segment {index}: denominator={denominator}, "
                f"primary_accuracy_digits={minimum_vector_accuracy_digits(primary)}, "
                f"reference_accuracy_digits={minimum_vector_accuracy_digits(reference)}"
            ) from error
        primary_snapshots.append(acb_mat(primary))
        reference_snapshots.append(acb_mat(reference))
        reports.append(
            {
                "segment_index": index,
                "start": acb_record(start),
                "target": acb_record(target),
                "primary_order": primary_order,
                "reference_order": reference_order,
                "refinement_relative_delta_inf": arb_record(relative),
                "refinement_relative_delta_midpoint": arb_midpoint_float(relative),
                "primary_seconds": primary_seconds,
                "reference_seconds": reference_seconds,
                "elapsed_seconds": time.perf_counter() - segment_clock,
            }
        )
    return (
        primary_snapshots,
        reference_snapshots,
        reports,
        time.perf_counter() - total_clock,
    )


def build_positive_real_fractional_path(
    poles: list[acb],
    first_point: acb,
    *,
    segment_count: int,
    step_fraction: fmpq,
) -> tuple[list[acb], list[FractionalPathSegment]]:
    """每段沿正实轴前进当前最近 pole 距离的固定比例。"""

    if segment_count <= 0:
        raise ValueError("segment count must be positive")
    fraction = arb(step_fraction)
    if not 0 < arb_midpoint_float(fraction) < 1:
        raise ValueError("step fraction must lie strictly between zero and one")
    if first_point.real.contains(0) or not first_point.imag.contains(0):
        raise ValueError("first point must be positive and real")
    points = [acb(first_point)]
    segments: list[FractionalPathSegment] = []
    for index in range(1, segment_count + 1):
        start = points[-1]
        nearest_pole = min(poles, key=lambda pole: arb_midpoint_float(abs(start - pole)))
        radius = abs(start - nearest_pole)
        step = radius * fraction
        target = start - acb(step)
        segments.append(
            FractionalPathSegment(
                index=index,
                start=start,
                target=target,
                nearest_pole=nearest_pole,
                radius=radius,
                step=step,
                step_over_radius=step / radius,
            )
        )
        points.append(target)
    return points, segments


def build_positive_real_path_through_targets(
    poles: list[acb],
    first_point: acb,
    targets: list[acb],
    *,
    step_fraction: fmpq,
) -> tuple[list[acb], list[FractionalPathSegment], list[int]]:
    """按最近 pole 半径限步，并依次精确命中所有正实轴目标点。"""

    fraction = arb(step_fraction)
    if not 0 < arb_midpoint_float(fraction) < 1:
        raise ValueError("step fraction must lie strictly between zero and one")
    points = [acb(first_point)]
    segments: list[FractionalPathSegment] = []
    target_positions: list[int] = []
    for target in targets:
        if arb_midpoint_float(target.real) >= arb_midpoint_float(points[-1].real):
            raise ValueError("required positive-real targets must decrease strictly")
        while arb_midpoint_float(points[-1].real - target.real) > 0:
            start = points[-1]
            nearest_pole = min(
                poles,
                key=lambda pole: arb_midpoint_float(abs(start - pole)),
            )
            radius = abs(start - nearest_pole)
            maximum_step = radius * fraction
            remaining = start.real - target.real
            if arb_midpoint_float(maximum_step) >= arb_midpoint_float(remaining):
                next_point = acb(target)
                step = abs(next_point - start)
            else:
                step = maximum_step
                next_point = start - acb(step)
            segment = FractionalPathSegment(
                index=len(segments) + 1,
                start=start,
                target=next_point,
                nearest_pole=nearest_pole,
                radius=radius,
                step=step,
                step_over_radius=step / radius,
            )
            segments.append(segment)
            points.append(next_point)
        target_positions.append(len(points) - 1)
    return points, segments, target_positions


def _append_fractional_line_segment(
    points: list[acb],
    segments: list[FractionalPathSegment],
    poles: list[acb],
    target: acb,
    *,
    step_fraction: arb,
    maximum_segment_count: int,
) -> None:
    """沿一条复直线按最近 pole 半径分步，并精确命中给定折跃点。"""

    while not abs(points[-1] - target).contains(0):
        if len(segments) >= maximum_segment_count:
            raise RuntimeError(
                "transport path exceeded MaximumSegmentCount; "
                "a pole may lie on an unguarded line segment"
            )
        start = points[-1]
        nearest_pole = min(
            poles,
            key=lambda pole: arb_midpoint_float(abs(start - pole)),
        )
        radius = abs(start - nearest_pole)
        if radius.contains(0):
            raise ValueError("transport path starts on a DE pole")
        remaining = abs(target - start)
        maximum_step = radius * arb(str(step_fraction))
        if arb_midpoint_float(maximum_step) >= arb_midpoint_float(remaining):
            next_point = acb(target)
        else:
            next_point = start + (target - start) * (maximum_step / remaining)
        step = abs(next_point - start)
        if step.is_zero():
            raise ArithmeticError("transport path subdivision made no progress")
        segments.append(
            FractionalPathSegment(
                index=len(segments) + 1,
                start=start,
                target=next_point,
                nearest_pole=nearest_pole,
                radius=radius,
                step=step,
                step_over_radius=step / radius,
            )
        )
        points.append(next_point)


def _midpoint_line_clearance(start: acb, target: acb, pole: acb) -> float:
    """用 Acb 中点计算 pole 到有限线段的距离，供折线路径预检查。"""

    left = complex(
        arb_midpoint_float(start.real), arb_midpoint_float(start.imag)
    )
    right = complex(
        arb_midpoint_float(target.real), arb_midpoint_float(target.imag)
    )
    singularity = complex(
        arb_midpoint_float(pole.real), arb_midpoint_float(pole.imag)
    )
    direction = right - left
    if direction == 0:
        return abs(singularity - left)
    position = ((singularity - left).real * direction.real + (
        singularity - left
    ).imag * direction.imag) / abs(direction) ** 2
    position = min(1.0, max(0.0, position))
    return abs(singularity - (left + position * direction))


def build_real_path_through_targets_with_detours(
    poles: list[acb],
    first_point: acb,
    targets: list[acb],
    *,
    step_fraction: fmpq,
    real_axis_pole_policy: str,
    detour_side: str,
    detour_clearance_fraction: fmpq,
    maximum_segment_count: int,
) -> tuple[
    list[acb],
    list[FractionalPathSegment],
    list[int],
    list[dict[str, Any]],
]:
    """命中递减正实目标点，并对路径上的非零实 pole 作指定半平面绕行。

    无交叉 pole 时复用正实轴路径构造器；存在交叉 pole 时加入指定半平面绕行。
    当前物理 prescription 使用 ``lower``，即从右向左按 ``k-i0`` 连续延拓。
    """

    fraction = arb(step_fraction)
    clearance_fraction = arb(detour_clearance_fraction)
    if not 0 < arb_midpoint_float(fraction) < 1:
        raise ValueError("step fraction must lie strictly between zero and one")
    if not 0 < arb_midpoint_float(clearance_fraction) < 0.5:
        raise ValueError("detour clearance fraction must lie between zero and one half")
    if maximum_segment_count <= 0:
        raise ValueError("MaximumSegmentCount must be positive")
    if real_axis_pole_policy != "detour":
        raise ValueError("RealAxisPolePolicy must be 'detour'")
    if detour_side not in {"lower", "upper"}:
        raise ValueError("DetourSide must be 'lower' or 'upper'")
    if first_point.real.contains(0) or not first_point.imag.is_zero():
        raise ValueError("first point must be positive and real")

    previous = arb_midpoint_float(first_point.real)
    for target in targets:
        current = arb_midpoint_float(target.real)
        if not target.imag.is_zero() or current <= 0 or current >= previous:
            raise ValueError("required targets must be positive real and decrease strictly")
        previous = current

    real_poles = sorted(
        (
            pole
            for pole in poles
            if pole.imag.is_zero()
            and not pole.real.contains(0)
            and arb_midpoint_float(pole.real) > 0
        ),
        key=lambda pole: arb_midpoint_float(pole.real),
        reverse=True,
    )
    lowest_target = arb_midpoint_float(targets[-1].real)
    first_real = arb_midpoint_float(first_point.real)
    crossed_poles = [
        pole
        for pole in real_poles
        if lowest_target < arb_midpoint_float(pole.real) < first_real
    ]
    if not crossed_poles:
        path, segments, target_positions = build_positive_real_path_through_targets(
            poles, first_point, targets, step_fraction=step_fraction
        )
        if len(segments) > maximum_segment_count:
            raise RuntimeError("transport path exceeded MaximumSegmentCount")
        return path, segments, target_positions, []

    anchors = [acb(first_point), *[acb(target) for target in targets]]
    for pole in crossed_poles:
        if any(abs(pole - anchor).contains(0) for anchor in anchors):
            raise ValueError("a transport endpoint coincides with a nonzero real pole")

    points = [acb(first_point)]
    segments: list[FractionalPathSegment] = []
    target_positions: list[int] = []
    detours: list[dict[str, Any]] = []
    remaining_poles = list(crossed_poles)
    current_real = first_real
    side_sign = -1 if detour_side == "lower" else 1

    for target in targets:
        target_real = arb_midpoint_float(target.real)
        pending = [
            pole
            for pole in remaining_poles
            if target_real < arb_midpoint_float(pole.real) < current_real
        ]
        for pole in pending:
            separation_candidates = [
                abs(pole - other)
                for other in [*poles, *anchors]
                if not abs(pole - other).contains(0)
            ]
            if not separation_candidates:
                raise ValueError("cannot determine a safe real-pole detour clearance")
            clearance_limit = min(
                separation_candidates, key=arb_midpoint_float
            )
            clearance = clearance_limit * clearance_fraction
            imaginary_offset = clearance if side_sign > 0 else -clearance
            waypoints = [
                acb(pole.real + clearance),
                acb(pole.real + clearance, imaginary_offset),
                acb(pole.real - clearance, imaginary_offset),
                acb(pole.real - clearance),
            ]
            for waypoint in waypoints:
                scale = max(
                    1.0,
                    abs(complex(
                        arb_midpoint_float(waypoint.real),
                        arb_midpoint_float(waypoint.imag),
                    )),
                )
                minimum_clearance = min(
                    _midpoint_line_clearance(points[-1], waypoint, candidate)
                    for candidate in poles
                )
                if minimum_clearance <= 1.0e-14 * scale:
                    raise ValueError("a proposed detour segment intersects another DE pole")
                _append_fractional_line_segment(
                    points,
                    segments,
                    poles,
                    waypoint,
                    step_fraction=fraction,
                    maximum_segment_count=maximum_segment_count,
                )
            detours.append(
                {
                    "pole": acb_record(pole),
                    "side": detour_side,
                    "clearance": arb_record(clearance),
                    "waypoints": [acb_record(waypoint) for waypoint in waypoints],
                }
            )
            remaining_poles.remove(pole)
            current_real = arb_midpoint_float(waypoints[-1].real)
        _append_fractional_line_segment(
            points,
            segments,
            poles,
            target,
            step_fraction=fraction,
            maximum_segment_count=maximum_segment_count,
        )
        target_positions.append(len(points) - 1)
        current_real = target_real

    if remaining_poles:
        raise RuntimeError("real-axis pole routing left an unhandled crossed pole")
    return points, segments, target_positions, detours


def _exact_matrix_to_acb(records: list[list[Any]]) -> acb_mat:
    """把 exact 二维记录转换为 Acb 稠密矩阵。"""

    cache: dict[str, acb] = {}

    def convert(value: Any) -> acb:
        key = str(value)
        if key not in cache:
            cache[key] = exact_inputform_to_acb(value)
        return cache[key]

    return acb_mat([[convert(value) for value in row] for row in records])


def _exact_vector_to_acb(records: list[Any]) -> acb_mat:
    """把 exact 一维记录转换为 Acb 列向量。"""

    return column_vector([exact_inputform_to_acb(value) for value in records])


def _validate_k0_sample_indices(
    sample_points: list[acb],
    fit_index_1based: int,
    stability_index_1based: int,
    validation_indices_1based: list[int],
) -> tuple[int, int, list[int]]:
    """统一检查 k0 匹配、稳定性和验证点索引。"""

    fit_index = fit_index_1based - 1
    stability_index = stability_index_1based - 1
    validation = sorted({index - 1 for index in validation_indices_1based})
    all_indices = [fit_index, stability_index, *validation]
    if any(index < 0 or index >= len(sample_points) for index in all_indices):
        raise ValueError("k0 matching or validation sample index is out of range")
    if fit_index in validation or stability_index in validation:
        raise ValueError("k0 matching/stability sample overlaps validation set")
    return fit_index, stability_index, validation


def _local_single_root_matrix(
    point: acb,
    root: acb,
    coefficient_by_log: list[list[acb_mat]],
) -> acb_mat:
    """计算单根 Jordan 路线的局部基本矩阵。"""

    dimension = coefficient_by_log[0][0].nrows()
    logarithm = point.log()
    value = acb_mat(dimension, dimension)
    for log_degree, series in enumerate(coefficient_by_log):
        polynomial = acb_mat(dimension, dimension)
        for coefficient in reversed(series):
            polynomial = polynomial * point + coefficient
        value += polynomial * (logarithm ** log_degree)
    return value * (point**root)


def _solve_single_root_jordan_k0(
    system: FlintPoleResidueSystem,
    manifest: dict[str, Any],
    sample_points: list[acb],
    sample_vectors: list[acb_mat],
    *,
    series_order: int,
    fit_index_1based: int,
    stability_index_1based: int,
    validation_indices_1based: list[int],
) -> tuple[acb_mat, dict[str, Any], list[K0LeadingLogCoefficient]]:
    """按 exact nilpotent gate 用 Acb 递推单根 Jordan 局部基本解。"""

    clock = time.perf_counter()
    residue, regular = system.laurent_at_zero(series_order)
    dimension = system.dimension
    identity = identity_matrix(dimension)
    root_expression = exact_inputform_expression(manifest["root_exact"])
    root = exact_inputform_to_acb(manifest["root_exact"])
    maximum_log_degree = int(manifest["maximum_log_degree"])
    nilpotent = _exact_matrix_to_acb(manifest["nilpotent_exact"])
    coefficient_by_log: list[list[acb_mat]] = []
    nilpotent_power = acb_mat(identity)
    for log_degree in range(maximum_log_degree + 1):
        if log_degree > 0:
            nilpotent_power = nilpotent_power * nilpotent
        coefficient_by_log.append(
            [nilpotent_power / acb(math.factorial(log_degree))]
        )
    for degree in range(1, series_order + 1):
        operator_inverse = acb_mat(dimension, dimension)
        nilpotent_power = acb_mat(identity)
        for nilpotent_degree in range(maximum_log_degree + 1):
            operator_inverse += nilpotent_power / acb(
                degree ** (nilpotent_degree + 1)
            )
            nilpotent_power = nilpotent_power * nilpotent
        degree_coefficients = [
            acb_mat(dimension, dimension) for _ in range(maximum_log_degree + 1)
        ]
        for log_degree in range(maximum_log_degree, -1, -1):
            right = acb_mat(dimension, dimension)
            for regular_degree in range(degree):
                right += regular[regular_degree] * coefficient_by_log[log_degree][
                    degree - 1 - regular_degree
                ]
            if log_degree < maximum_log_degree:
                right -= degree_coefficients[log_degree + 1] * acb(log_degree + 1)
            degree_coefficients[log_degree] = operator_inverse * right
        for log_degree, coefficient in enumerate(degree_coefficients):
            coefficient_by_log[log_degree].append(coefficient)
    fit_index, stability_index, validation_indices = _validate_k0_sample_indices(
        sample_points,
        fit_index_1based,
        stability_index_1based,
        validation_indices_1based,
    )
    fit_matrix = _local_single_root_matrix(
        sample_points[fit_index], root, coefficient_by_log
    )
    constants = fit_matrix.solve(sample_vectors[fit_index])
    stability_constants = _local_single_root_matrix(
        sample_points[stability_index], root, coefficient_by_log
    ).solve(sample_vectors[stability_index])
    if sp.simplify(root_expression) == 0:
        endpoint = coefficient_by_log[0][0] * constants
        stability_endpoint = coefficient_by_log[0][0] * stability_constants
    else:
        endpoint = acb_mat(dimension, 1)
        stability_endpoint = acb_mat(dimension, 1)
    validation = []
    for index in validation_indices:
        reconstructed = _local_single_root_matrix(
            sample_points[index], root, coefficient_by_log
        ) * constants
        error = relative_difference_inf(reconstructed, sample_vectors[index])
        validation.append(
            {
                "sample_index_1based": index + 1,
                "k": acb_record(sample_points[index]),
                "relative_error_inf": arb_record(error),
                "relative_error_midpoint": arb_midpoint_float(error),
            }
        )
    # 只有绝对幂为零的 log 项会阻碍 k -> 0 极限；正幂 log 项随 k 一同消失。
    leading_log_coefficients = []
    leading_logs = []
    for degree in range(series_order + 1):
        absolute_power = sp.simplify(root_expression + degree)
        if absolute_power != 0:
            continue
        for log_degree in range(1, maximum_log_degree + 1):
            coefficient = coefficient_by_log[log_degree][degree] * constants
            norm = vector_norm_inf(coefficient)
            leading_log_coefficients.append(
                K0LeadingLogCoefficient(
                    absolute_power_exact=str(absolute_power),
                    log_degree=log_degree,
                    coefficient_vector=coefficient,
                )
            )
            leading_logs.append(
                {
                    "absolute_power_exact": str(absolute_power),
                    "log_degree": log_degree,
                    "coefficient_norm_inf": arb_record(norm),
                    "coefficient_norm_midpoint": arb_midpoint_float(norm),
                }
            )
    stability = relative_difference_inf(stability_endpoint, endpoint)
    return endpoint, {
        "route": "single_root_jordan_exact_manifest_acb",
        "series_order": series_order,
        "root_exact": str(root_expression),
        "maximum_log_degree": maximum_log_degree,
        "matching_sample_index_1based": fit_index + 1,
        "stability_sample_index_1based": stability_index + 1,
        "finite_part_cross_match_relative_delta_inf": arb_record(stability),
        "finite_part_cross_match_relative_delta_midpoint": arb_midpoint_float(stability),
        "actual_log_power_norms": leading_logs,
        "validation": validation,
        "validation_max_relative_error_midpoint": max(
            record["relative_error_midpoint"] for record in validation
        ),
        "endpoint_minimum_accuracy_digits": minimum_vector_accuracy_digits(endpoint),
        "elapsed_seconds": time.perf_counter() - clock,
    }, leading_log_coefficients


def _solve_projected_degree(
    projectors: dict[sp.Expr, acb_mat],
    absolute_power: sp.Expr,
    right_by_log: list[acb_mat],
    *,
    solution_index: int,
    resonance_gates: dict[str, bool],
) -> list[acb_mat]:
    """按 exact resonance gate 解一个幂次的 Acb power-log 系数。"""

    dimension = right_by_log[0].nrows()
    existing_log_degree = len(right_by_log) - 1
    if absolute_power not in projectors:
        inverse = acb_mat(dimension, dimension)
        for root, projector in projectors.items():
            inverse += projector / exact_inputform_to_acb(absolute_power - root)
        coefficients = [
            acb_mat(dimension, 1) for _ in range(existing_log_degree + 1)
        ]
        for log_degree in range(existing_log_degree, -1, -1):
            corrected = acb_mat(right_by_log[log_degree])
            if log_degree < existing_log_degree:
                corrected -= coefficients[log_degree + 1] * acb(log_degree + 1)
            coefficients[log_degree] = inverse * corrected
        return coefficients
    resonant_projector = projectors[absolute_power]
    reduced_inverse = acb_mat(dimension, dimension)
    for root, projector in projectors.items():
        if root != absolute_power:
            reduced_inverse += projector / exact_inputform_to_acb(absolute_power - root)
    coefficients = [
        acb_mat(dimension, 1) for _ in range(existing_log_degree + 2)
    ]
    for log_degree in range(existing_log_degree, -1, -1):
        corrected = acb_mat(right_by_log[log_degree])
        corrected -= coefficients[log_degree + 1] * acb(log_degree + 1)
        gate_key = f"{solution_index}:{absolute_power}:{log_degree}"
        if gate_key not in resonance_gates:
            raise ValueError(f"exact resonance manifest misses gate {gate_key}")
        if bool(resonance_gates[gate_key]):
            defect = resonant_projector * corrected
            coefficients[log_degree + 1] += defect / acb(log_degree + 1)
            corrected -= defect
        coefficients[log_degree] = reduced_inverse * corrected
    while len(coefficients) > 1 and all(
        abs(coefficients[-1][row, 0]).contains(0)
        for row in range(dimension)
    ):
        coefficients.pop()
    return coefficients


def _local_diagonalizable_matrix(
    point: acb,
    roots: list[acb],
    series_by_solution: list[list[list[acb_mat]]],
) -> acb_mat:
    """计算可对角化根路线的局部基本矩阵。"""

    dimension = len(series_by_solution)
    logarithm = point.log()
    columns: list[acb_mat] = []
    for root, series in zip(roots, series_by_solution):
        value = acb_mat(dimension, 1)
        for degree_coefficients in reversed(series):
            log_polynomial = acb_mat(dimension, 1)
            for coefficient in reversed(degree_coefficients):
                log_polynomial = log_polynomial * logarithm + coefficient
            value = value * point + log_polynomial
        columns.append(value * (point**root))
    return acb_mat(
        [
            [columns[column][row, 0] for column in range(dimension)]
            for row in range(dimension)
        ]
    )


def _solve_diagonalizable_k0(
    system: FlintPoleResidueSystem,
    manifest: dict[str, Any],
    sample_points: list[acb],
    sample_vectors: list[acb_mat],
    *,
    series_order: int,
    fit_index_1based: int,
    stability_index_1based: int,
    validation_indices_1based: list[int],
) -> tuple[acb_mat, dict[str, Any], list[K0LeadingLogCoefficient]]:
    """按 exact projector/resonance gate 用 Acb 递推可对角化局部基。"""

    clock = time.perf_counter()
    _, regular = system.laurent_at_zero(series_order)
    dimension = system.dimension
    roots_exact = [exact_inputform_expression(value) for value in manifest["roots_exact"]]
    projectors = {
        root: _exact_matrix_to_acb(records)
        for root, records in zip(roots_exact, manifest["projectors_exact"])
    }
    solution_roots_exact = [
        exact_inputform_expression(value) for value in manifest["solution_roots_exact"]
    ]
    solution_roots = [exact_inputform_to_acb(value) for value in solution_roots_exact]
    initial_vectors = [
        _exact_vector_to_acb(records) for records in manifest["initial_vectors_exact"]
    ]
    resonance_gates = {
        str(key): bool(value) for key, value in manifest["resonance_gates"].items()
    }
    series_by_solution: list[list[list[acb_mat]]] = []
    log_records = []
    for solution_index, (root_exact, initial) in enumerate(
        zip(solution_roots_exact, initial_vectors), start=1
    ):
        series: list[list[acb_mat]] = [[initial]]
        for degree in range(1, series_order + 1):
            previous_log_degree = max(len(item) for item in series) - 1
            right_by_log = [
                acb_mat(dimension, 1) for _ in range(previous_log_degree + 1)
            ]
            for regular_degree in range(degree):
                earlier = series[degree - 1 - regular_degree]
                for log_degree, coefficient in enumerate(earlier):
                    right_by_log[log_degree] += regular[regular_degree] * coefficient
            coefficients = _solve_projected_degree(
                projectors,
                root_exact + degree,
                right_by_log,
                solution_index=solution_index,
                resonance_gates=resonance_gates,
            )
            if len(coefficients) - 1 > previous_log_degree:
                log_records.append(
                    {
                        "solution_index": solution_index,
                        "root_exact": str(root_exact),
                        "series_degree": degree,
                        "absolute_power_exact": str(root_exact + degree),
                        "new_maximum_log_degree": len(coefficients) - 1,
                    }
                )
            series.append(coefficients)
        series_by_solution.append(series)
    fit_index, stability_index, validation_indices = _validate_k0_sample_indices(
        sample_points,
        fit_index_1based,
        stability_index_1based,
        validation_indices_1based,
    )
    fit_matrix = _local_diagonalizable_matrix(
        sample_points[fit_index], solution_roots, series_by_solution
    )
    constants = fit_matrix.solve(sample_vectors[fit_index])
    stability_constants = _local_diagonalizable_matrix(
        sample_points[stability_index], solution_roots, series_by_solution
    ).solve(sample_vectors[stability_index])
    endpoint = acb_mat(dimension, 1)
    stability_endpoint = acb_mat(dimension, 1)
    for root_exact, initial, index in zip(
        solution_roots_exact, initial_vectors, range(dimension)
    ):
        if sp.simplify(root_exact) == 0:
            endpoint += initial * constants[index, 0]
            stability_endpoint += initial * stability_constants[index, 0]
    grouped_leading_logs: dict[tuple[str, int], acb_mat] = {}
    for root_exact, series, index in zip(
        solution_roots_exact, series_by_solution, range(dimension)
    ):
        for degree, coefficients in enumerate(series):
            absolute_power = sp.simplify(root_exact + degree)
            if absolute_power != 0:
                continue
            for log_degree, coefficient in enumerate(coefficients[1:], start=1):
                key = (str(absolute_power), log_degree)
                if key not in grouped_leading_logs:
                    grouped_leading_logs[key] = acb_mat(dimension, 1)
                grouped_leading_logs[key] += coefficient * constants[index, 0]
    leading_log_coefficients = [
        K0LeadingLogCoefficient(power, log_degree, coefficient)
        for (power, log_degree), coefficient in sorted(grouped_leading_logs.items())
    ]
    leading_log_norms = [
        {
            "absolute_power_exact": record.absolute_power_exact,
            "log_degree": record.log_degree,
            "coefficient_norm_inf": arb_record(vector_norm_inf(record.coefficient_vector)),
            "coefficient_norm_midpoint": arb_midpoint_float(
                vector_norm_inf(record.coefficient_vector)
            ),
        }
        for record in leading_log_coefficients
    ]
    validation = []
    for index in validation_indices:
        reconstructed = _local_diagonalizable_matrix(
            sample_points[index], solution_roots, series_by_solution
        ) * constants
        error = relative_difference_inf(reconstructed, sample_vectors[index])
        validation.append(
            {
                "sample_index_1based": index + 1,
                "k": acb_record(sample_points[index]),
                "relative_error_inf": arb_record(error),
                "relative_error_midpoint": arb_midpoint_float(error),
            }
        )
    stability = relative_difference_inf(stability_endpoint, endpoint)
    return endpoint, {
        "route": "diagonalizable_real_roots_exact_manifest_acb",
        "series_order": series_order,
        "roots_exact": [str(root) for root in roots_exact],
        "log_terms": log_records,
        "actual_log_power_norms": leading_log_norms,
        "matching_sample_index_1based": fit_index + 1,
        "stability_sample_index_1based": stability_index + 1,
        "finite_part_cross_match_relative_delta_inf": arb_record(stability),
        "finite_part_cross_match_relative_delta_midpoint": arb_midpoint_float(stability),
        "validation": validation,
        "validation_max_relative_error_midpoint": max(
            record["relative_error_midpoint"] for record in validation
        ),
        "endpoint_minimum_accuracy_digits": minimum_vector_accuracy_digits(endpoint),
        "elapsed_seconds": time.perf_counter() - clock,
    }, leading_log_coefficients


def solve_k0_finite_part(
    system: FlintPoleResidueSystem,
    manifest: dict[str, Any],
    sample_points: list[acb],
    sample_vectors: list[acb_mat],
    *,
    series_order: int,
    fit_index_analytic_leading: int,
    fit_index_leading_log: int,
    stability_index_analytic_leading: int,
    stability_index_leading_log: int,
    validation_indices_1based: list[int],
) -> tuple[acb_mat, dict[str, Any], list[K0LeadingLogCoefficient]]:
    """提取 ``k=0`` finite part，并返回实际 ``k^0 log(k)^q`` 审计向量。"""

    if manifest.get("status") != "passed" or manifest.get("request_id") != system.request_id:
        raise ValueError(f"{system.request_id}: exact k0 manifest identity/status mismatch")
    route = manifest.get("route")
    if route == "single_root_jordan_exact_gate":
        return _solve_single_root_jordan_k0(
            system,
            manifest,
            sample_points,
            sample_vectors,
            series_order=series_order,
            fit_index_1based=fit_index_leading_log,
            stability_index_1based=stability_index_leading_log,
            validation_indices_1based=validation_indices_1based,
        )
    if route == "diagonalizable_real_roots_exact_gate":
        return _solve_diagonalizable_k0(
            system,
            manifest,
            sample_points,
            sample_vectors,
            series_order=series_order,
            fit_index_1based=fit_index_analytic_leading,
            stability_index_1based=stability_index_analytic_leading,
            validation_indices_1based=validation_indices_1based,
        )
    raise ValueError(f"{system.request_id}: unsupported exact k0 manifest route {route}")
