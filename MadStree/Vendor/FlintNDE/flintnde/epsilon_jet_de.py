#!/usr/bin/env python3
"""求解 ``(x-1)^epsilon exp[-k(x-1)]`` 积分族的局部矩阵微分方程。

本模块在固定物理 QNM 数值点保持三条频率严格满足 EC，只把 horizon regulator
实现为作用点平移 ``a2 -> a2-epsilon``。普通点使用截断 epsilon jet 与 FLINT Acb
递推 ``dI/dk=A(k,epsilon)I``；奇点接口保留联合
``epsilon^q k^n log(k)^r`` 局部展开，供完整 scalar 装配后提取有限部。

通用 regulator jet、路径输运与局部基由同一 FlintNDE 包的 regular-point 模块提供；
本模块只把这些中性接口映射为 epsilon 语义，不改变任何物理频率 iw_i。"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any, Sequence

from flint import acb, acb_mat, acb_series


from .regular_point_de import (
    FlintRegulatorRationalSystem as FlintEpsilonRationalSystem,
    RegulatorK0NilpotentLocalBasis as EpsilonK0NilpotentLocalBasis,
    _cu_qnm_bc_regulator_series as _cu_qnm_bc_epsilon_series,
    _derivative_series,
    _integer_binomial_series,
    _multiply_taylor_series_over_regulator as _multiply_taylor_series_over_epsilon,
    arb_midpoint_float,
    build_regulator_k0_nilpotent_local_basis as build_epsilon_k0_nilpotent_local_basis,
    exact_inputform_to_acb,
    relative_difference_inf,
    transport_regulator_path as transport_epsilon_path,
)


@dataclass
class EpsilonRegularTransportResult:
    """保存一条解析 epsilon 普通点路径的边界、快照、报告和耗时。"""

    boundary_vectors: list[acb_mat]
    snapshots: list[list[acb_mat]]
    segment_reports: list[dict[str, Any]]
    boundary_seconds: float
    transport_seconds: float


@dataclass
class EpsilonK0MatchResult:
    """保存解析 epsilon 普通点输运与 ``k=0`` 广义级数匹配结果。

    ``constants`` 是堆叠后局部基的齐次常数；``finite_part_vectors`` 按
    epsilon 阶数保存 ``k^0 log(k)^0`` 的 16 维向量。
    """

    regular_transport: EpsilonRegularTransportResult
    local_basis: EpsilonK0NilpotentLocalBasis
    constants: acb_mat
    stability_constants: acb_mat
    finite_part_vectors: list[acb_mat]
    validation_errors: list[Any]
    stability_errors: list[Any]
    local_basis_seconds: float
    fit_seconds: float


@dataclass
class EpsilonK0FitResult:
    """保存同一局部 power-log 基在一组普通点上的拟合与稳定性结果。"""

    constants: acb_mat
    stability_constants: acb_mat
    finite_part_vectors: list[acb_mat]
    validation_errors: list[Any]
    stability_errors: list[Any]
    fit_seconds: float


class FlintLocalEpsilonPreinverseSystem:
    """按局部 ``(k-center,epsilon)`` 双级数数值构造 16/24 维 DE。

    三条频率固定为 metadata 中的 exact EC 值。epsilon 只实现
    ``a2 -> a2-epsilon``，因此 preinverse 的 `010` block 在线性 epsilon 阶
    增加单位矩阵；积分标签与 Taylor 次数始终保持整数。无 Z 输入使用
    ``[XMinusOne,X]`` 16 维 basis，active-Z 输入使用
    ``[Z,XMinusOne,X]`` 24 维 basis。
    """

    def __init__(
        self,
        metadata: dict[str, Any],
        *,
        epsilon_order: int,
        it0: int,
        it1: int,
        a_x_zero: int,
        regulator_name: str = "epsilon",
    ) -> None:
        if epsilon_order < 0:
            raise ValueError("epsilon_order must be nonnegative")
        if regulator_name != "epsilon":
            raise ValueError("FlintNDE epsilon-jet interface accepts only regulator_name='epsilon'")
        point = metadata["numeric_point"]
        self.request_id = f"local_numeric_{metadata['request_id']}"
        physical_basis_keys = metadata["physical_basis_keys"]
        self.dimension = len(physical_basis_keys)
        if self.dimension not in {16, 24}:
            raise ValueError(
                f"epsilon preinverse system requires dimension 16 or 24, got {self.dimension}"
            )
        self.basis_mode = "subsector16" if self.dimension == 16 else "top24"
        # 共享 jet 输运器使用 regulator-neutral 协议；本模块固定解释为 epsilon。
        self.epsilon_order = epsilon_order
        self.regulator_order = epsilon_order
        self.regulator_name = regulator_name
        self.iws = [
            exact_inputform_to_acb(point["iw1"]),
            exact_inputform_to_acb(point["iw2"]),
            exact_inputform_to_acb(point["iw3"]),
        ]
        if not abs(self.iws[2] - self.iws[1] - self.iws[0]).contains(0):
            raise ValueError("epsilon route requires exact EC frequencies")
        self.cls = [
            exact_inputform_to_acb(point[f"cl2{leg}"])
            for leg in range(1, 4)
        ]
        self.it0 = acb(it0)
        self.it1 = acb(it1)
        self.a_x_zero = acb(a_x_zero)
        labels = [
            tuple(int(value) for value in key[:3])
            for key in physical_basis_keys[:8]
        ]
        expected = {
            (n1, n2, n3)
            for n1 in range(2)
            for n2 in range(2)
            for n3 in range(2)
        }
        if len(labels) != 8 or set(labels) != expected:
            raise ValueError(
                "physical basis does not contain the expected binary derivative labels"
            )
        self.labels = labels
        self.label_positions = {label: index for index, label in enumerate(labels)}
        self.active_z_sector = int(metadata.get("active_z_sector", 0))
        if self.basis_mode == "subsector16":
            if self.active_z_sector != 0:
                raise ValueError("subsector16 metadata must use active_z_sector=0")
            self.active_z_coefficient = acb(1)
        else:
            if self.active_z_sector not in {1, 2, 3}:
                raise ValueError("top24 metadata requires active_z_sector in {1,2,3}")
            expected_z_powers = [0, 0, 0]
            expected_z_powers[self.active_z_sector - 1] = 1
            if any(
                [int(value) for value in key[5:]] != expected_z_powers
                for key in physical_basis_keys[:8]
            ):
                raise ValueError(
                    "top24 first block does not match the declared physical active-Z sector"
                )
            self.active_z_coefficient = self.cls[self.active_z_sector - 1]
        self.poles_at_epsilon0 = [
            acb(0),
            *[sign * 2 * iw for iw in self.iws for sign in (1, -1)],
        ]
        # `transport_epsilon_path` 读取上面的 regulator-neutral 协议属性。

    def _zero_biseries(self, k_taylor_order: int) -> list[list[acb_mat]]:
        """创建布局为 ``[epsilon power][k Taylor degree]`` 的 8x8 零矩阵表。"""

        return [
            [acb_mat(8, 8) for _ in range(k_taylor_order)]
            for _ in range(self.epsilon_order + 1)
        ]

    def _preinverse_block_series(
        self,
        center: acb,
        k_taylor_order: int,
    ) -> dict[str, list[list[acb_mat]]]:
        """直接数值生成 generic ``000/100/010/001`` preinverse 双级数。"""

        if k_taylor_order <= 0:
            raise ValueError("k_taylor_order must be positive")
        blocks = {
            axis: self._zero_biseries(k_taylor_order)
            for axis in ("000", "100", "010", "001")
        }
        identity8 = acb_mat(8, 8)
        for index in range(8):
            identity8[index, index] = acb(1)
        blocks["000"][0][0] -= identity8 * center
        if k_taylor_order > 1:
            blocks["000"][0][1] -= identity8

        if self.epsilon_order >= 1:
            blocks["010"][1][0] += identity8
        for row, label in enumerate(self.labels):
            blocks["100"][0][0][row, row] += -self.a_x_zero + sum(
                acb(value) * (1 - 2 * self.it0) for value in label
            )
            blocks["010"][0][0][row, row] += self.it1 * sum(
                self.iws,
                acb(0),
            )
            for leg, derivative_power in enumerate(label):
                iw0 = self.iws[leg]
                blocks["010"][0][0][row, row] += acb(derivative_power) * (
                    -1 - 2 * self.it1 * iw0
                )

                toggled = list(label)
                toggled[leg] = 1 - derivative_power
                column = self.label_positions[tuple(toggled)]
                if derivative_power == 0:
                    blocks["000"][0][0][row, column] += 1
                    continue

                blocks["000"][0][0][row, column] += iw0 * iw0

                x_constant = 1 - 2 * self.cls[leg] + self.it0
                x_linear = self.it1 * (2 * self.it0 - 1)
                blocks["100"][0][0][row, column] += x_constant + x_linear * iw0

                xm_constant = -1 + 2 * self.cls[leg] - self.it0
                xm_linear = self.it1 * (1 - 2 * self.it0)
                blocks["010"][0][0][row, column] += (
                    xm_constant + xm_linear * iw0 + 2 * iw0 * iw0
                )
        return blocks

    def _inverse_biseries_generic(
        self,
        block: list[list[acb_mat]],
        k_taylor_order: int,
    ) -> list[list[acb_mat]]:
        """在截断双级数环中递推矩阵逆，仅供独立回归检查。"""

        inverse = self._zero_biseries(k_taylor_order)
        identity8 = acb_mat(8, 8)
        for index in range(8):
            identity8[index, index] = acb(1)
        inverse[0][0] = block[0][0].solve(identity8)
        for k_degree in range(k_taylor_order):
            for epsilon_power in range(self.epsilon_order + 1):
                if k_degree == 0 and epsilon_power == 0:
                    continue
                total = acb_mat(8, 8)
                for left_k in range(k_degree + 1):
                    for left_epsilon in range(epsilon_power + 1):
                        if left_k == 0 and left_epsilon == 0:
                            continue
                        total += (
                            block[left_epsilon][left_k]
                            * inverse[epsilon_power - left_epsilon][k_degree - left_k]
                        )
                inverse[epsilon_power][k_degree] = -(inverse[0][0] * total)
        return inverse

    def _inverse_preinverse000_biseries(
        self,
        block: list[list[acb_mat]],
        k_taylor_order: int,
    ) -> list[list[acb_mat]]:
        """利用截断 ``(epsilon,k)`` 卷积结构递推矩阵逆。

        输入采用 ``[epsilon power][k Taylor degree]`` 布局。固定数值频率后，
        本路线的 ``000`` block 与 epsilon 无关且对 k 只含线性项；保留通用卷积
        写法是为了和完整矩阵输入接口交叉检查。
        """

        inverse = self._zero_biseries(k_taylor_order)
        identity8 = acb_mat(8, 8)
        for index in range(8):
            identity8[index, index] = acb(1)
        inverse000 = block[0][0].solve(identity8)
        inverse[0][0] = inverse000
        epsilon_linear = block[1][0] if self.epsilon_order >= 1 else None
        epsilon_quadratic = block[2][0] if self.epsilon_order >= 2 else None

        for k_degree in range(k_taylor_order):
            for epsilon_power in range(self.epsilon_order + 1):
                if k_degree == 0 and epsilon_power == 0:
                    continue
                right_hand_side = acb_mat(8, 8)
                if k_degree >= 1:
                    right_hand_side += inverse[epsilon_power][k_degree - 1]
                if epsilon_power >= 1:
                    right_hand_side -= (
                        epsilon_linear * inverse[epsilon_power - 1][k_degree]
                    )
                if epsilon_power >= 2:
                    right_hand_side -= (
                        epsilon_quadratic * inverse[epsilon_power - 2][k_degree]
                    )
                inverse[epsilon_power][k_degree] = inverse000 * right_hand_side
        return inverse

    def _multiply_biseries(
        self,
        left: list[list[acb_mat]],
        right: list[list[acb_mat]],
        k_taylor_order: int,
    ) -> list[list[acb_mat]]:
        """卷积两个同截断范围的 8×8 双级数。"""

        product = self._zero_biseries(k_taylor_order)
        for k_degree in range(k_taylor_order):
            for epsilon_power in range(self.epsilon_order + 1):
                total = acb_mat(8, 8)
                for left_k in range(k_degree + 1):
                    for left_epsilon in range(epsilon_power + 1):
                        total += (
                            left[left_epsilon][left_k]
                            * right[epsilon_power - left_epsilon][k_degree - left_k]
                        )
                product[epsilon_power][k_degree] = total
        return product

    def _assemble_de_matrix(
        self,
        left: acb_mat,
        right: acb_mat,
        add_identity: bool,
        z_block: acb_mat | None = None,
    ) -> acb_mat:
        """按 16/24 维物理 basis 顺序组装 k-DE 矩阵。"""

        if self.basis_mode == "subsector16":
            return acb_mat(
                [
                    [
                        (
                            left[row % 8, column]
                            if column < 8
                            else right[row % 8, column - 8]
                            + (
                                acb(1)
                                if add_identity
                                and row >= 8
                                and row - 8 == column - 8
                                else acb(0)
                            )
                        )
                        for column in range(16)
                    ]
                    for row in range(16)
                ]
            )

        if z_block is None:
            raise ValueError("top24 DE assembly requires the active-Z source block")
        inverse_two_z = acb(1) / (acb(2) * self.active_z_coefficient)
        z_diagonal = acb(1) + acb(3) * inverse_two_z
        rows: list[list[acb]] = []
        for row in range(24):
            source_row = row % 8
            values: list[acb] = []
            for column in range(24):
                block = column // 8
                source_column = column % 8
                if row < 8:
                    base = (z_block, left, right)[block][source_row, source_column]
                    value = base * inverse_two_z
                    if add_identity and block == 0 and source_row == source_column:
                        value += z_diagonal
                else:
                    value = (z_block, left, right)[block][source_row, source_column]
                    if (
                        add_identity
                        and row >= 16
                        and block == 2
                        and source_row == source_column
                    ):
                        value += acb(1)
                values.append(value)
            rows.append(values)
        return acb_mat(rows)

    def taylor_matrix_coefficients(
        self,
        center: acb,
        solution_order: int,
    ) -> list[list[acb_mat]]:
        """返回 ``A[epsilon_power][k_taylor_degree]`` 供 regulator-jet 后端调用。"""

        blocks = self._preinverse_block_series(center, solution_order)
        inverse000 = self._inverse_preinverse000_biseries(
            blocks["000"],
            solution_order,
        )
        de_left = self._multiply_biseries(
            inverse000,
            blocks["010"],
            solution_order,
        )
        de_right = self._multiply_biseries(
            inverse000,
            blocks["100"],
            solution_order,
        )
        de_z = self._multiply_biseries(
            inverse000,
            blocks["001"],
            solution_order,
        )
        return [
            [
                self._assemble_de_matrix(
                    de_left[epsilon_power][k_degree],
                    de_right[epsilon_power][k_degree],
                    epsilon_power == 0 and k_degree == 0,
                    de_z[epsilon_power][k_degree],
                )
                for k_degree in range(solution_order)
            ]
            for epsilon_power in range(self.epsilon_order + 1)
        ]


def _constant_epsilon_series(value: acb, epsilon_order: int) -> acb_series:
    """把 Acb 常数提升为指定截断阶数的 epsilon 级数。"""

    return acb_series([value], epsilon_order + 1)


def watson_epsilon_boundary_coefficients(
    system_payload: dict[str, Any],
    *,
    k_start: acb,
    taylor_order: int,
    epsilon_order: int,
    it0: int,
    it1: int,
) -> tuple[list[acb_mat], float]:
    """构造 fixed-EC QNM Watson 边界的解析 epsilon 系数。

    单腿 QNM 解与频率均不随 epsilon 改变；horizon regulator 只把每个积分分量的
    Laplace/Gamma 指数从 ``sigma`` 平移为 ``sigma+epsilon``。
    """

    required = {"physical_basis_keys", "numeric_point"}
    missing = sorted(required.difference(system_payload))
    if missing:
        raise ValueError(f"epsilon Watson payload is missing fields: {missing}")
    if taylor_order < 0 or epsilon_order < 0:
        raise ValueError("Watson Taylor/epsilon orders must be nonnegative")
    clock = time.perf_counter()
    point = system_payload["numeric_point"]
    iws_scalar = [
        exact_inputform_to_acb(point["iw1"]),
        exact_inputform_to_acb(point["iw2"]),
        exact_inputform_to_acb(point["iw3"]),
    ]
    if not abs(iws_scalar[2] - iws_scalar[1] - iws_scalar[0]).contains(0):
        raise ValueError("epsilon Watson boundary requires exact EC frequencies")
    iws = [_constant_epsilon_series(value, epsilon_order) for value in iws_scalar]
    cls = [exact_inputform_to_acb(point[f"cl2{leg}"]) for leg in range(1, 4)]
    leg_u = [
        _cu_qnm_bc_epsilon_series(
            iw,
            lam,
            taylor_order + 2,
            epsilon_order,
            it0=it0,
            it1=it1,
        )
        for iw, lam in zip(iws, cls)
    ]
    leg_du = [_derivative_series(coefficients) for coefficients in leg_u]
    z_value = acb(1) / k_start
    iw_sum = sum(iws, _constant_epsilon_series(acb(0), epsilon_order))
    epsilon_series = (
        acb_series([acb(0), acb(1)], epsilon_order + 1)
        if epsilon_order >= 1
        else _constant_epsilon_series(acb(0), 0)
    )
    component_series: list[acb_series] = []
    for key in system_payload["physical_basis_keys"]:
        n_values = key[:3]
        a_x, a_x_minus_one = key[3:5]
        z_powers = key[5:]
        active_legs = [index for index, power in enumerate(z_powers) if power != 0]
        if len(active_legs) > 1:
            raise ValueError(f"multiple active Z factors are unsupported: {key}")
        x_factor = [
            _constant_epsilon_series(value, epsilon_order)
            for value in _integer_binomial_series(-a_x, taylor_order)
        ]
        if not active_legs:
            z_factor = [_constant_epsilon_series(acb(1), epsilon_order)] + [
                _constant_epsilon_series(acb(0), epsilon_order)
                for _ in range(taylor_order)
            ]
        else:
            active_leg = active_legs[0]
            z_power = z_powers[active_leg]
            z_horizon = 3 + 2 * cls[active_leg]
            ratio = 2 * cls[active_leg] / z_horizon
            z_factor = [
                _constant_epsilon_series(
                    coefficient * ratio**degree * z_horizon ** (-z_power),
                    epsilon_order,
                )
                for degree, coefficient in enumerate(
                    _integer_binomial_series(-z_power, taylor_order)
                )
            ]
        modes = [
            leg_u[leg] if n_values[leg] == 0 else leg_du[leg]
            for leg in range(3)
        ]
        endpoint = _multiply_taylor_series_over_epsilon(
            [x_factor, z_factor, *modes],
            taylor_order,
            epsilon_order,
        )
        sigma = it1 * iw_sum - a_x_minus_one + epsilon_series
        component = _constant_epsilon_series(acb(0), epsilon_order)
        for degree in range(taylor_order + 1):
            gamma_argument = sigma + degree + 1
            if abs(gamma_argument[0]).contains(0):
                raise ArithmeticError(
                    f"epsilon Watson Gamma argument contains zero for {key}, degree {degree}"
                )
            laplace_power = (gamma_argument * z_value.log()).exp()
            component += endpoint[degree] * gamma_argument.gamma() * laplace_power
        component_series.append(component)
    coefficients = [
        acb_mat([[component[power]] for component in component_series])
        for power in range(epsilon_order + 1)
    ]
    return coefficients, time.perf_counter() - clock


def solve_epsilon_regular_path(
    metadata: dict[str, Any],
    path: list[acb],
    *,
    epsilon_order: int,
    k_taylor_order: int,
    watson_taylor_order: int,
    it0: int,
    it1: int,
    a_x_zero: int,
) -> EpsilonRegularTransportResult:
    """从解析 Watson 边界沿 ``path`` 输运全部 epsilon 系数。"""

    if len(path) < 2:
        raise ValueError("path must contain at least two k points")
    system = FlintLocalEpsilonPreinverseSystem(
        metadata,
        epsilon_order=epsilon_order,
        it0=it0,
        it1=it1,
        a_x_zero=a_x_zero,
    )
    boundary, boundary_seconds = watson_epsilon_boundary_coefficients(
        metadata,
        k_start=path[0],
        taylor_order=watson_taylor_order,
        epsilon_order=epsilon_order,
        it0=it0,
        it1=it1,
    )
    snapshots, reports, transport_seconds = transport_epsilon_path(
        system,
        boundary,
        path,
        k_taylor_order,
    )
    return EpsilonRegularTransportResult(
        boundary_vectors=boundary,
        snapshots=snapshots,
        segment_reports=reports,
        boundary_seconds=boundary_seconds,
        transport_seconds=transport_seconds,
    )


def _stack_epsilon_vectors(vectors: list[acb_mat]) -> acb_mat:
    """按 epsilon 幂次由低到高堆叠等长列向量。"""

    if not vectors:
        raise ValueError("epsilon vector list must not be empty")
    dimension = vectors[0].nrows()
    if any(vector.ncols() != 1 or vector.nrows() != dimension for vector in vectors):
        raise ValueError("epsilon vectors must be equal-sized columns")
    return acb_mat(
        [[vector[row, 0]] for vector in vectors for row in range(dimension)]
    )


def fit_epsilon_k0_local_basis(
    local_basis: EpsilonK0NilpotentLocalBasis,
    sample_points: Sequence[acb],
    sample_vectors: Sequence[list[acb_mat]],
) -> EpsilonK0FitResult:
    """用三个或更多普通点拟合、验证同一个 ``k=0`` power-log 局部基。

    第一个点确定结果常数，第二个点独立确定稳定性常数，其余点只验证第一组常数。
    每个普通点必须提供从 ``epsilon^0`` 到 ``epsilon^P`` 的完整 16 维列向量。
    """

    if len(sample_points) < 3 or len(sample_points) != len(sample_vectors):
        raise ValueError("k0 local fit requires equally sized lists with at least three samples")
    expected_count = local_basis.regulator_order + 1
    if any(len(vectors) != expected_count for vectors in sample_vectors):
        raise ValueError("sample vectors do not cover the local-basis epsilon order")

    import time

    local_matrices = [local_basis.evaluate(point) for point in sample_points]
    stacked_vectors = [_stack_epsilon_vectors(vectors) for vectors in sample_vectors]
    fit_clock = time.perf_counter()
    constants = local_matrices[0].solve(stacked_vectors[0])
    stability_constants = local_matrices[1].solve(stacked_vectors[1])
    validation_errors = [
        relative_difference_inf(matrix * constants, vector)
        for matrix, vector in zip(local_matrices[1:], stacked_vectors[1:])
    ]
    finite_parts = [
        local_basis.finite_part_matrix(power) * constants
        for power in range(expected_count)
    ]
    stability_parts = [
        local_basis.finite_part_matrix(power) * stability_constants
        for power in range(expected_count)
    ]
    stability_errors = [
        relative_difference_inf(primary, reference)
        for primary, reference in zip(finite_parts, stability_parts)
    ]
    return EpsilonK0FitResult(
        constants=constants,
        stability_constants=stability_constants,
        finite_part_vectors=finite_parts,
        validation_errors=validation_errors,
        stability_errors=stability_errors,
        fit_seconds=time.perf_counter() - fit_clock,
    )


def _epsilon_payload_for_jet_backend(payload: dict[str, Any]) -> dict[str, Any]:
    """把公开 epsilon schema 映射到共享 regulator-jet 后端的内部字段。"""

    required = {"epsilon_order", "epsilon_coefficients"}
    missing = sorted(required.difference(payload))
    if missing:
        raise ValueError(f"epsilon rational payload is missing fields: {missing}")
    records = []
    for record in payload["epsilon_coefficients"]:
        converted = dict(record)
        converted["regulator_power"] = int(record["epsilon_power"])
        records.append(converted)
    converted_payload = dict(payload)
    converted_payload["regulator_order"] = int(payload["epsilon_order"])
    converted_payload["regulator_coefficients"] = records
    return converted_payload


def solve_epsilon_k0_finite_parts(
    metadata: dict[str, Any],
    rational_payload: dict[str, Any],
    path: list[acb],
    sample_positions: Sequence[int],
    *,
    epsilon_order: int,
    k_taylor_order: int,
    watson_taylor_order: int,
    k0_series_order: int,
    it0: int,
    it1: int,
    a_x_zero: int,
) -> EpsilonK0MatchResult:
    """输运 ``I_0..I_P``，并在三个近零常点匹配 ``k=0`` power-log 基。

    第一个 sample 用于求齐次常数，第二个独立重求常数并检查有限部稳定性，其余
    sample 只验证局部级数。payload 必须至少覆盖请求的 ``epsilon_order``；任何身份、
    维数或阶数不一致都直接阻断。
    """

    if len(sample_positions) < 3:
        raise ValueError("k0 matching requires at least three sample positions")
    if any(position < 0 or position >= len(path) for position in sample_positions):
        raise ValueError("k0 sample position is outside the transport path")
    for key, expected in (
        ("dimension", 16),
        ("it0", it0),
        ("it1", it1),
        ("a_x_zero", a_x_zero),
    ):
        if int(rational_payload[key]) != expected:
            raise ValueError(
                f"rational payload {key}={rational_payload[key]!r}, expected {expected!r}"
            )
    if int(rational_payload["epsilon_order"]) < epsilon_order:
        raise ValueError("rational payload does not cover the requested epsilon order")

    regular = solve_epsilon_regular_path(
        metadata,
        path,
        epsilon_order=epsilon_order,
        k_taylor_order=k_taylor_order,
        watson_taylor_order=watson_taylor_order,
        it0=it0,
        it1=it1,
        a_x_zero=a_x_zero,
    )
    rational_system = FlintEpsilonRationalSystem.from_payload(
        _epsilon_payload_for_jet_backend(rational_payload)
    )
    local_basis, local_seconds = build_epsilon_k0_nilpotent_local_basis(
        rational_system,
        solution_regulator_order=epsilon_order,
        series_order=k0_series_order,
    )
    sample_vectors = [
        _stack_epsilon_vectors(regular.snapshots[position])
        for position in sample_positions
    ]
    # 局部拟合统一走预计算向量接口；完整 scalar 的有限部仍由上层统一提取。
    fit = fit_epsilon_k0_local_basis(
        local_basis,
        [path[position] for position in sample_positions],
        [regular.snapshots[position] for position in sample_positions],
    )
    maximum_error = max(
        [
            *(arb_midpoint_float(error) for error in local_basis.compatibility_residuals),
            *(arb_midpoint_float(error) for error in fit.validation_errors),
            *(arb_midpoint_float(error) for error in fit.stability_errors),
        ],
        default=0.0,
    )
    if not maximum_error < 1.0:
        raise ArithmeticError(f"epsilon k0 matching is numerically unstable: {maximum_error}")
    return EpsilonK0MatchResult(
        regular_transport=regular,
        local_basis=local_basis,
        constants=fit.constants,
        stability_constants=fit.stability_constants,
        finite_part_vectors=fit.finite_part_vectors,
        validation_errors=fit.validation_errors,
        stability_errors=fit.stability_errors,
        local_basis_seconds=local_seconds,
        fit_seconds=fit.fit_seconds,
    )
