#!/usr/bin/env python3
"""用 FlintNDE 从字面量 ``inf`` 和 horizon 奇点验证 Schwarzschild QNM。

示例把统一 ``u`` 二阶方程写成 exact ``RationalMatrixSystem``。两个端点都使用
``{a,b,C}`` 边界；无穷远的指数根由二阶 Laurent 矩阵和 ``C`` 自动识别。独立保留的
二阶标量递推只用于核对局部无穷远初始化，不参与 FlintNDE 主输运。
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
import warnings
from pathlib import Path
from typing import Any

from flint import acb, acb_mat


# 独立复制 example 时只需让该路径指向已安装源码或删除本段后使用正式安装。
EXAMPLE_DIR = Path(__file__).resolve().parent
# 发布布局直接取上级目录；发展仓布局取 versions/ 下版本号最大的当前版本目录。
PACKAGE_ROOT = (EXAMPLE_DIR / "..").resolve()
if not (PACKAGE_ROOT / "flintnde").is_dir():
    _version_dirs = sorted((PACKAGE_ROOT / "versions").glob("FlintNDE-*"))
    if not _version_dirs:
        raise RuntimeError("no versions/FlintNDE-* directory beside examples/")

    def _version_key(path: Path) -> tuple[int, ...]:
        parts = path.name.removeprefix("FlintNDE-").lstrip("v").split(".")
        try:
            return tuple(int(part) for part in parts)
        except ValueError:
            raise RuntimeError(f"cannot parse version directory name: {path}") from None

    PACKAGE_ROOT = max(_version_dirs, key=_version_key).resolve()
DEFAULT_CONFIG_PATH = EXAMPLE_DIR / ".." / "config" / "qnm_u_unified_it0_3_it1_minus1.json"
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from flintnde import (  # noqa: E402
    GaussianRational,
    NamedPoint,
    RationalMatrixSystem,
    build_adaptive_path,
    build_local_solution_basis,
    column_vector,
    configure_working_precision,
    exponential_boundary,
    frobenius_boundary,
    gaussian_rational,
    initialize_output_layout,
    plan_transport_path,
    rational_function,
    relative_difference_inf,
    transport_path_refined,
    transport_planned_path_refined,
)


configure_working_precision(70, 32)

ELL = 2
QNM_IW_EXACT = {
    "real": "0.17792463137787093",
    "imag": "0.74734336883608354",
}
CONTROL_RELATIVE_SHIFT = gaussian_rational("1/1000")
X_MATCH = acb(4)
INFINITY_MATCH_RATIO = 1.0 / 30.0
HORIZON_MATCH_RATIO = 0.12
FORMAL_MINIMUM_ORDER_FACTOR = 3
CAUCHY_RADIUS_FRACTION = 0.40
LOCAL_SERIES_ORDER = 55
PRIMARY_ORDER = 38
REFERENCE_ORDER = 55
TARGET_RELATIVE_ERROR = "1e-10"
PRIMARY_SAMPLES = 88
REFERENCE_SAMPLES = 120
CURVE_RSTAR_MIN = -6.0
CURVE_RSTAR_MAX = 160.0
CURVE_POINT_COUNT = 720
# QNM 解在无穷远按 e^{Re(w) r*} 增长，长输运链上每步截断误差会被后续增长
# 放大；规划步长收紧到半径的 6% 后，38/55 阶即可全程保持 1e-10 相对精度。
CURVE_PLANNING_RADIUS_FRACTION = 0.06


def load_convention_config(path: str | Path = DEFAULT_CONFIG_PATH) -> dict[str, Any]:
    """读取并验证 package-local QNM convention，返回含实际绝对路径的记录。"""

    config_path = Path(path).expanduser().resolve()
    data = json.loads(config_path.read_text(encoding="utf-8"))
    required = {"schema", "ConventionName", "it0", "it1", "EnergyCondition"}
    missing = sorted(required - set(data))
    if missing:
        raise ValueError(f"QNM convention config is missing fields: {missing}")
    if isinstance(data["it0"], bool) or not isinstance(data["it0"], int):
        raise TypeError("QNM convention it0 must be an integer")
    if isinstance(data["it1"], bool) or not isinstance(data["it1"], int):
        raise TypeError("QNM convention it1 must be an integer")
    return {**data, "ConfigPath": str(config_path)}


def _unified_parameters(
    iw: GaussianRational,
    convention: dict[str, Any],
    ell: int = ELL,
) -> tuple[GaussianRational, ...]:
    """返回统一方程中依赖 exact 频率、角动量和 convention 的五个系数。"""

    it0 = gaussian_rational(convention["it0"])
    it1 = gaussian_rational(convention["it1"])
    lam = gaussian_rational((ell - 1) * (ell + 2)) / 2
    a_x = -1 + 2 * it0
    a_h = 1 + 2 * it1 * iw
    b_h = 1 - it1 * iw - 2 * iw * iw - 2 * lam + it0 + 2 * it1 * iw * it0
    b_x = -1 + it1 * iw + 2 * lam - it0 - 2 * it1 * iw * it0
    return lam, a_x, a_h, b_h, b_x


def unified_system(
    iw: Any,
    convention: dict[str, Any],
    ell: int = ELL,
) -> RationalMatrixSystem:
    """构造 ``d(u,u')/dx=DEmatrix(x)(u,u')`` 的 exact 有理矩阵。"""

    exact_iw = gaussian_rational(iw)
    _lam, a_x, a_h, b_h, b_x = _unified_parameters(exact_iw, convention, ell)
    x = rational_function((0, 1))
    a_value = rational_function(a_x) / x + rational_function(a_h) / (x - 1)
    b_value = (
        rational_function(-(exact_iw**2))
        + rational_function(b_h) / (x - 1)
        + rational_function(b_x) / x
    )
    return RationalMatrixSystem(
        ((0, 1), (-b_value, -a_value)),
        variable_name="x",
        name=f"qnm-u-l{ell}",
    )


def horizon_coefficients(
    iw: GaussianRational,
    root: acb,
    order: int,
    convention: dict[str, Any],
    ell: int = ELL,
) -> list[acb]:
    """独立递推 horizon 标量 Frobenius 分支，固定 ``u`` 首系数为一。"""

    _lam, a_x_exact, a_h_exact, b_h_exact, b_x_exact = _unified_parameters(
        iw, convention, ell
    )
    a_x, a_h = a_x_exact.to_acb(), a_h_exact.to_acb()
    b_h, b_x = b_h_exact.to_acb(), b_x_exact.to_acb()
    iw_acb = iw.to_acb()
    a_regular = [a_x * acb((-1) ** degree) for degree in range(order + 1)]
    b_regular = [b_x * acb((-1) ** degree) for degree in range(order + 1)]
    b_regular[0] -= iw_acb * iw_acb
    coefficients = [acb(1)]
    for degree in range(1, order + 1):
        known = b_h * coefficients[degree - 1]
        for regular_degree in range(degree):
            previous = degree - 1 - regular_degree
            known += a_regular[regular_degree] * (root + previous) * coefficients[previous]
        for regular_degree in range(degree - 1):
            known += b_regular[regular_degree] * coefficients[degree - 2 - regular_degree]
        denominator = (root + degree) * (root + degree - 1 + a_h)
        if abs(denominator).contains(0):
            raise ValueError(f"horizon recurrence is resonant at degree {degree}")
        coefficients.append(-known / denominator)
    return coefficients


def horizon_vector(
    iw: GaussianRational,
    root: acb,
    point: acb,
    convention: dict[str, Any],
    ell: int = ELL,
) -> acb_mat:
    """独立计算 horizon 邻域的 ``(u,u')``，只用于局部交叉验证。"""

    t_value = point - 1
    coefficients = horizon_coefficients(
        iw, root, LOCAL_SERIES_ORDER, convention, ell
    )
    series = acb(0)
    derivative_series = acb(0)
    for degree, coefficient in enumerate(coefficients):
        series += coefficient * t_value**degree
        if degree > 0:
            derivative_series += degree * coefficient * t_value ** (degree - 1)
    common = t_value**root
    return column_vector(
        [common * series, common * (derivative_series + root * series / t_value)]
    )


def infinity_coefficients(
    iw: GaussianRational,
    branch: str,
    order: int,
    convention: dict[str, Any],
    ell: int = ELL,
) -> tuple[acb, acb, list[acb]]:
    """独立递推二阶标量方程的 infinity outgoing/incoming 广义级数。"""

    if branch not in {"outgoing", "incoming"}:
        raise ValueError(f"unsupported infinity branch: {branch}")
    _lam, a_x_exact, a_h_exact, b_h_exact, b_x_exact = _unified_parameters(
        iw, convention, ell
    )
    iw_acb = iw.to_acb()
    a_x, a_h = a_x_exact.to_acb(), a_h_exact.to_acb()
    b_h, b_x = b_h_exact.to_acb(), b_x_exact.to_acb()
    exponential = iw_acb if branch == "outgoing" else -iw_acb
    power = (
        acb(convention["it0"]) - 2 * iw_acb
        if branch == "outgoing"
        else acb(convention["it0"])
    )
    c_zero = a_x + a_h
    d_zero = b_x + b_h

    def p_coefficient(index: int) -> acb:
        if index == 2:
            return -2 * exponential
        if index == 3:
            return 2 * power + 2 - c_zero
        return -a_h if index >= 4 else acb(0)

    def q_coefficient(index: int) -> acb:
        if index < 1:
            return acb(0)
        first = (
            -2 * exponential * power + exponential * c_zero + d_zero
            if index == 1
            else exponential * a_h + b_h
        )
        second = power * (power + 1) - power * c_zero if index == 2 else -power * a_h
        return first + (second if index >= 2 else 0)

    if not abs(q_coefficient(1)).contains(0):
        raise ValueError("infinity exponent failed the independent indicial gate")
    coefficients = [acb(1)]
    for degree in range(1, order + 1):
        target_power = degree + 1
        residual = acb(0)
        lower = target_power - 2
        if 0 <= lower < degree:
            residual += lower * (lower - 1) * coefficients[lower]
        for index in range(2, target_power + 1):
            coefficient_index = target_power - index + 1
            if 0 <= coefficient_index < degree:
                residual += p_coefficient(index) * coefficient_index * coefficients[coefficient_index]
        for index in range(1, target_power + 1):
            coefficient_index = target_power - index
            if 0 <= coefficient_index < degree:
                residual += q_coefficient(index) * coefficients[coefficient_index]
        denominator = p_coefficient(2) * degree + q_coefficient(1)
        if abs(denominator).contains(0):
            raise ValueError(f"infinity recurrence is resonant at degree {degree}")
        coefficients.append(-residual / denominator)
    return exponential, power, coefficients


def infinity_vector(
    iw: GaussianRational,
    branch: str,
    point: acb,
    convention: dict[str, Any],
    ell: int = ELL,
) -> acb_mat:
    """独立求值 infinity 标量递推；不调用 FlintNDE 的矩阵形式递推。"""

    z_value = 1 / point
    exponential, power, coefficients = infinity_coefficients(
        iw, branch, LOCAL_SERIES_ORDER, convention, ell
    )
    series = acb(0)
    derivative_z = acb(0)
    for degree, coefficient in enumerate(coefficients):
        series += coefficient * z_value**degree
        if degree > 0:
            derivative_z += degree * coefficient * z_value ** (degree - 1)
    common = (exponential / z_value).exp() * z_value**power
    value = common * series
    derivative = common * (
        exponential * series - power * z_value * series - z_value * z_value * derivative_z
    )
    return column_vector([value, derivative])


def _matrix_from_columns(first: acb_mat, second: acb_mat) -> acb_mat:
    """把两个二维列向量拼成局部基矩阵。"""

    return acb_mat([[first[row, 0], second[row, 0]] for row in range(2)])


def _complex_record(value: acb, digits: int = 35) -> dict[str, str]:
    """保存 Acb 中点和 ball，避免 JSON 中退回 binary64。"""

    return {
        "re": value.real.str(digits, radius=False),
        "im": value.imag.str(digits, radius=False),
        "re_ball": value.real.str(digits),
        "im_ball": value.imag.str(digits),
    }


def _ratio(numerator: acb, denominator: acb) -> acb:
    """返回禁戒/允许系数的绝对值比。"""

    if abs(denominator).contains(0):
        raise ZeroDivisionError("allowed boundary coefficient contains zero")
    return acb(abs(numerator) / abs(denominator))


def _transport_from_endpoint(
    system: RationalMatrixSystem,
    boundary: Any,
    start: NamedPoint,
    max_step_over_radius: float,
) -> tuple[dict[str, Any], Any]:
    """从一个字面奇点起步，执行主/参考两条阶数链并返回自描述路径。"""

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        path = build_adaptive_path(
            system,
            start,
            NamedPoint("match", 4),
            max_step_over_radius=max_step_over_radius,
            formal_asymptotic_order=PRIMARY_ORDER,
            formal_minimum_order_factor=FORMAL_MINIMUM_ORDER_FACTOR,
        )
    result = transport_path_refined(
        system,
        boundary,
        path,
        primary_order=PRIMARY_ORDER,
        reference_order=REFERENCE_ORDER,
        primary_sample_count=PRIMARY_SAMPLES,
        reference_sample_count=REFERENCE_SAMPLES,
        # 更远的无穷远匹配点会放大 DFT 混叠误差，缩小采样圆比提高尾项抑制。
        radius_fraction=CAUCHY_RADIUS_FRACTION,
        target_relative_error=TARGET_RELATIVE_ERROR,
    )
    return result, path


def _selected_asymptotic_diagnostic(transport: dict[str, Any]) -> dict[str, Any]:
    """从奇点初始化报告提取用户实际选择的指数 sector 截断诊断。"""

    initialization = transport["primary_segments"][0]
    evaluation = initialization["local_evaluation"]
    diagnostic = dict(evaluation["selected_branch_diagnostics"][0])
    diagnostic["formal_accuracy_checks_passed"] = evaluation[
        "formal_accuracy_checks_passed"
    ]
    diagnostic["formal_accuracy_issues"] = evaluation["formal_accuracy_issues"]
    return diagnostic


def run_frequency(
    label: str,
    iw: Any,
    ell: int = ELL,
    *,
    config_path: str | Path = DEFAULT_CONFIG_PATH,
) -> dict[str, Any]:
    """从 literal horizon/inf 分别出发，并计算对端禁戒分量。"""

    clock = time.perf_counter()
    convention = load_convention_config(config_path)
    exact_iw = gaussian_rational(iw)
    system = unified_system(exact_iw, convention, ell)
    _lam, _a_x, a_h, b_h, _b_x = _unified_parameters(exact_iw, convention, ell)

    horizon_allowed_boundary = frobenius_boundary(
        [{"a": 0, "b": 0, "C": [1, -b_h / a_h]}]
    )
    horizon_root = 1 - a_h
    horizon_forbidden_boundary = frobenius_boundary(
        [{"a": -a_h, "b": 0, "C": [0, horizon_root]}]
    )
    infinity_outgoing_boundary = exponential_boundary(
        [
            {
                "phi": [{"power": -1, "coefficient": exact_iw}],
                "a": gaussian_rational(convention["it0"]) - 2 * exact_iw,
                "b": 0,
                "C": [1, exact_iw],
            }
        ]
    )
    infinity_incoming_boundary = exponential_boundary(
        [
            {
                "phi": [{"power": -1, "coefficient": -exact_iw}],
                "a": convention["it0"],
                "b": 0,
                "C": [1, -exact_iw],
            }
        ]
    )

    horizon_allowed, horizon_path = _transport_from_endpoint(
        system,
        horizon_allowed_boundary,
        NamedPoint("horizon", 1),
        HORIZON_MATCH_RATIO,
    )
    horizon_forbidden, _ = _transport_from_endpoint(
        system,
        horizon_forbidden_boundary,
        NamedPoint("horizon", 1),
        HORIZON_MATCH_RATIO,
    )
    infinity_allowed, infinity_path = _transport_from_endpoint(
        system,
        infinity_outgoing_boundary,
        NamedPoint("infinity", "inf"),
        INFINITY_MATCH_RATIO,
    )
    infinity_forbidden, _ = _transport_from_endpoint(
        system,
        infinity_incoming_boundary,
        NamedPoint("infinity", "inf"),
        INFINITY_MATCH_RATIO,
    )

    h_allowed = horizon_allowed["primary_snapshots"][-1]
    h_forbidden = horizon_forbidden["primary_snapshots"][-1]
    i_allowed = infinity_allowed["primary_snapshots"][-1]
    i_forbidden = infinity_forbidden["primary_snapshots"][-1]
    infinity_coefficients_at_match = _matrix_from_columns(i_allowed, i_forbidden).solve(h_allowed)
    horizon_coefficients_at_match = _matrix_from_columns(h_allowed, h_forbidden).solve(i_allowed)
    infinity_contamination = _ratio(
        infinity_coefficients_at_match[1, 0], infinity_coefficients_at_match[0, 0]
    )
    horizon_contamination = _ratio(
        horizon_coefficients_at_match[1, 0], horizon_coefficients_at_match[0, 0]
    )

    # 在 package 实际首匹配点，用独立二阶标量递推核对两个 infinity 列。
    first_z = infinity_path[0]
    manual_outgoing = infinity_vector(exact_iw, "outgoing", 1 / first_z, convention, ell)
    manual_incoming = infinity_vector(exact_iw, "incoming", 1 / first_z, convention, ell)
    outgoing_local_error = relative_difference_inf(
        infinity_allowed["primary_snapshots"][0], manual_outgoing
    )
    incoming_local_error = relative_difference_inf(
        infinity_forbidden["primary_snapshots"][0], manual_incoming
    )

    inverted_basis = build_local_solution_basis(system.inverted(), 0, REFERENCE_ORDER)
    nearby_z = acb(1) / 28
    nearby_out_constants, _ = inverted_basis.resolve_boundary(infinity_outgoing_boundary)
    nearby_in_constants, _ = inverted_basis.resolve_boundary(infinity_incoming_boundary)
    nearby_out = inverted_basis.evaluate(nearby_z) * nearby_out_constants.to_acb()
    nearby_in = inverted_basis.evaluate(nearby_z) * nearby_in_constants.to_acb()
    nearby_out_error = relative_difference_inf(
        nearby_out,
        infinity_vector(exact_iw, "outgoing", 1 / nearby_z, convention, ell),
    )
    nearby_in_error = relative_difference_inf(
        nearby_in,
        infinity_vector(exact_iw, "incoming", 1 / nearby_z, convention, ell),
    )

    return {
        "label": label,
        "ell": ell,
        "iw": _complex_record(exact_iw.to_acb()),
        "ConfigPath": convention["ConfigPath"],
        "ConventionName": convention["ConventionName"],
        "it0": convention["it0"],
        "it1": convention["it1"],
        "EnergyCondition": convention["EnergyCondition"],
        "horizon_start_to_infinity_basis": {
            "allowed_outgoing_coefficient": _complex_record(infinity_coefficients_at_match[0, 0]),
            "forbidden_incoming_coefficient": _complex_record(infinity_coefficients_at_match[1, 0]),
            "forbidden_over_allowed": _complex_record(infinity_contamination),
        },
        "infinity_start_to_horizon_basis": {
            "allowed_horizon_coefficient": _complex_record(horizon_coefficients_at_match[0, 0]),
            "forbidden_horizon_coefficient": _complex_record(horizon_coefficients_at_match[1, 0]),
            "forbidden_over_allowed": _complex_record(horizon_contamination),
        },
        "literal_singular_starts": {
            "horizon": "1",
            "infinity": "inf",
            "infinity_working_variable": "sinv=1/x",
            "infinity_first_ordinary_match_sinv": first_z.str(40),
            "infinity_first_ordinary_match_x": (1 / first_z).str(40),
            "infinity_method": "formal_exponential_asymptotic",
            "formal_asymptotic_match_estimate": (
                infinity_path.formal_asymptotic_match_estimate
            ),
        },
        "infinity_asymptotic_diagnostics": {
            "outgoing": _selected_asymptotic_diagnostic(infinity_allowed),
            "incoming": _selected_asymptotic_diagnostic(infinity_forbidden),
            "independent_scalar_recurrence_relative_difference_at_first_match": {
                "outgoing": outgoing_local_error.str(20),
                "incoming": incoming_local_error.str(20),
            },
            "nearby_sinv_1_over_28_relative_difference": {
                "outgoing": nearby_out_error.str(20),
                "incoming": nearby_in_error.str(20),
            },
            "classification": "genuinely formal/Gevrey asymptotic; no convergence radius is claimed",
        },
        "path_point_counts": {
            "horizon_to_match": len(horizon_path),
            "infinity_to_match": len(infinity_path),
        },
        "horizon_allowed_u_at_x_match": _complex_record(h_allowed[0, 0]),
        "transport_refinement_maximum": max(
            float(record["relative_difference_midpoint"])
            for record in (
                horizon_allowed,
                horizon_forbidden,
                infinity_allowed,
                infinity_forbidden,
            )
        ),
        "transport_target_relative_error": TARGET_RELATIVE_ERROR,
        "transport_target_relative_error_met": all(
            record["target_relative_error_met"]
            for record in (
                horizon_allowed,
                horizon_forbidden,
                infinity_allowed,
                infinity_forbidden,
            )
        ),
        "elapsed_seconds": time.perf_counter() - clock,
    }


def _x_from_tortoise(rstar: float) -> float:
    """Newton 求解 ``x + ln(x-1) = rstar``（``x > 1``），返回 float 坐标。"""

    target = rstar - 1.0
    t = target - math.log(target) if target > 2.0 else math.exp(target)
    for _ in range(60):
        step = (t + math.log(t) - target) / (1.0 + 1.0 / t)
        t -= step
        if abs(step) < 1.0e-15 * max(1.0, abs(t)):
            break
    return 1.0 + t


def run_curve(
    iw: Any,
    reference_u_at_match: acb,
    ell: int = ELL,
    *,
    config_path: str | Path = DEFAULT_CONFIG_PATH,
) -> dict[str, Any]:
    """从视界单端输运，输出 QNM 解全曲线并做三层交叉验证。

    曲线按 tortoise 坐标 ``rstar = x + ln(x-1)`` 均匀取点（近视界振荡在 ``x``
    坐标下无限加密，均匀 ``rstar`` 才能均匀采样振荡），再 Newton 反解 ``x``。
    除均匀网格外强制包含精确 ``x = 4``：该点的曲线值与 ``run_frequency`` 中
    独立规划的 horizon→match 链端点比较，验证不同节点划分的路径一致性。
    """

    clock = time.perf_counter()
    convention = load_convention_config(config_path)
    exact_iw = gaussian_rational(iw)
    iw_acb = exact_iw.to_acb()
    system = unified_system(exact_iw, convention, ell)
    _lam, _a_x, a_h, b_h, _b_x = _unified_parameters(exact_iw, convention, ell)
    horizon_allowed_boundary = frobenius_boundary(
        [{"a": 0, "b": 0, "C": [1, -b_h / a_h]}]
    )

    # tortoise 均匀网格 + 强制 x=4（rstar=4+ln3），按 rstar 排序去重。
    rstar_grid = [
        CURVE_RSTAR_MIN
        + (CURVE_RSTAR_MAX - CURVE_RSTAR_MIN) * index / (CURVE_POINT_COUNT - 1)
        for index in range(CURVE_POINT_COUNT)
    ]
    rstar_match = 4.0 + math.log(3.0)
    rstar_grid = sorted(set(rstar_grid + [rstar_match]))
    x_grid = [_x_from_tortoise(rstar) for rstar in rstar_grid]
    x_grid = [4.0 if abs(value - 4.0) < 1.0e-12 else value for value in x_grid]

    # 第一步：horizon 奇点起步到首个曲线点，得到普通初值并做独立递推种子核对。
    x_start = x_grid[0]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        seed_path = build_adaptive_path(
            system,
            NamedPoint("horizon", 1),
            NamedPoint("curve_start", repr(x_start)),
            max_step_over_radius=HORIZON_MATCH_RATIO,
        )
        seed_transport = transport_path_refined(
            system,
            horizon_allowed_boundary,
            seed_path,
            primary_order=PRIMARY_ORDER,
            reference_order=REFERENCE_ORDER,
            primary_sample_count=PRIMARY_SAMPLES,
            reference_sample_count=REFERENCE_SAMPLES,
            radius_fraction=CAUCHY_RADIUS_FRACTION,
            target_relative_error=TARGET_RELATIVE_ERROR,
        )
    seed_vector = seed_transport["primary_snapshots"][-1]
    seed_error = relative_difference_inf(
        seed_vector, horizon_vector(exact_iw, acb(0), acb(repr(x_start)), convention, ell)
    )

    # 第二步：从普通起点整体规划剩余曲线点；前瞻规则让一个节点覆盖收敛盘内
    # 的多个用户点，被跨过的点共享该段 Taylor 系数（桶式多点求值）。
    # 规划入口消费认证过的 P(x)+Sum R_j/(x-p_j) 公式系统；本例矩阵只有
    # 一阶有限极点加常数 -w^2，specialize 认证必然通过。
    planning_system, specialization_report = system.specialize_polynomial_simple_poles()
    if planning_system is None:
        raise RuntimeError(
            "QNM unified system failed the simple-pole specialization gate: "
            f"{specialization_report.get('reason')}"
        )
    user_points = [repr(value) for value in x_grid[1:]]
    plan = plan_transport_path(
        planning_system,
        acb(repr(x_start)),
        user_points,
        radius_fraction=CURVE_PLANNING_RADIUS_FRACTION,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        curve = transport_planned_path_refined(
            planning_system,
            seed_vector,
            plan,
            primary_order=PRIMARY_ORDER,
            reference_order=REFERENCE_ORDER,
            primary_sample_count=PRIMARY_SAMPLES,
            reference_sample_count=REFERENCE_SAMPLES,
            radius_fraction=CAUCHY_RADIUS_FRACTION,
            target_relative_error=TARGET_RELATIVE_ERROR,
        )

    # 第三步：按坐标收集采样值。被前瞻规则选为节点的用户点不在
    # sample_results 里，其值取主链该节点快照；两来源以 40 位坐标串对齐。
    values_by_key: dict[str, Any] = {}
    for record in curve["sample_results"]:
        values_by_key[record["coordinate"]] = record
    for node, snapshot in zip(plan.nodes, curve["primary_snapshots"]):
        values_by_key.setdefault(
            node.str(40),
            {"value": snapshot, "evaluation_algorithm": "node_snapshot"},
        )
    samples = []
    max_relative_radius = 0.0
    match_relative_difference = None
    for rstar, x_float in zip(rstar_grid[1:], x_grid[1:]):
        record = values_by_key[acb(repr(x_float)).str(40)]
        u_value = record["value"][0, 0]
        x_value = acb(repr(x_float))
        t_value = x_value - 1
        prefactor = (x_value**3) * (t_value ** (-iw_acb))
        u_physical = prefactor * u_value
        radius = max(u_value.real.rad(), u_value.imag.rad())
        magnitude = abs(u_value)
        relative_radius = float(radius / magnitude)
        max_relative_radius = max(max_relative_radius, relative_radius)
        if x_float == 4.0:
            match_relative_difference = float(
                (abs(u_value - reference_u_at_match) / abs(reference_u_at_match)).mid()
            )
        samples.append(
            {
                "rstar": rstar,
                "x": x_float,
                "re_u_physical": float(u_physical.real),
                "im_u_physical": float(u_physical.imag),
                "abs_u_physical": float(abs(u_physical)),
                "relative_radius": relative_radius,
                "evaluation_algorithm": record["evaluation_algorithm"],
            }
        )

    dense_buckets = [
        report["dense_evaluation"] for report in curve["primary_segments"]
        if "dense_evaluation" in report
    ]
    bucket_sizes = [bucket["point_count"] for bucket in dense_buckets]
    fast_buckets = sum(1 for bucket in dense_buckets if bucket["algorithm"] == "fast")
    curve_checks = {
        "seed_independent_recurrence": bool(float(seed_error.mid()) < 1.0e-20),
        "x_match_path_consistency": bool(
            match_relative_difference is not None
            and match_relative_difference < float(TARGET_RELATIVE_ERROR)
        ),
        "target_relative_error_met": bool(curve["target_relative_error_met"]),
    }
    return {
        "rstar_grid": {
            "minimum": CURVE_RSTAR_MIN,
            "maximum": CURVE_RSTAR_MAX,
            "uniform_points": CURVE_POINT_COUNT,
            "forced_x_match": 4.0,
            "total_points": len(x_grid),
        },
        "planning": {
            "node_count": len(plan.report["nodes"]),
            "segment_count": plan.report["segment_count"],
            "covered_sample_count": plan.report["covered_sample_count"],
            "effective_radius_fraction": plan.report["effective_radius_fraction"],
            "bucket_count": len(bucket_sizes),
            "fast_multipoint_buckets": fast_buckets,
            "largest_bucket": max(bucket_sizes),
            "mean_bucket": sum(bucket_sizes) / len(bucket_sizes),
        },
        "precision": {
            "target_relative_error": TARGET_RELATIVE_ERROR,
            "chain_relative_difference": curve["relative_difference_midpoint"],
            "seed_error_vs_independent_recurrence": seed_error.str(20),
            "x_match_relative_difference": match_relative_difference,
            "sample_value_note": (
                "dense evaluation returns midpoint matrices; per-sample radii are "
                "recorded in the samples file"
            ),
        },
        "timing_seconds": {
            "seed_transport": seed_transport["primary_seconds"]
            + seed_transport["reference_seconds"],
            "curve_primary": curve["primary_seconds"],
            "curve_reference": curve["reference_seconds"],
            "total": time.perf_counter() - clock,
        },
        "checks": curve_checks,
        "samples": samples,
    }


def main() -> None:
    """运行 QNM 与 exact 频率偏移反事实，并保存可复核摘要。"""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG_PATH)
    args = parser.parse_args()
    convention = load_convention_config(args.config)
    output_layout = initialize_output_layout(__file__)
    qnm_iw = gaussian_rational(QNM_IW_EXACT)
    qnm = run_frequency("qnm-frequency", qnm_iw, config_path=args.config)
    match_record = qnm["horizon_allowed_u_at_x_match"]
    reference_u_at_match = (
        acb(match_record["re"]) + acb(match_record["im"]) * acb(0, 1)
    )
    curve = run_curve(
        qnm_iw,
        reference_u_at_match,
        config_path=args.config,
    )
    control = run_frequency(
        "shifted-frequency-control",
        qnm_iw * (1 + CONTROL_RELATIVE_SHIFT),
        config_path=args.config,
    )
    qnm_inf = float(qnm["horizon_start_to_infinity_basis"]["forbidden_over_allowed"]["re"])
    qnm_hor = float(qnm["infinity_start_to_horizon_basis"]["forbidden_over_allowed"]["re"])
    control_inf = float(control["horizon_start_to_infinity_basis"]["forbidden_over_allowed"]["re"])
    control_hor = float(control["infinity_start_to_horizon_basis"]["forbidden_over_allowed"]["re"])
    status = (
        "passed"
        if max(qnm_inf, qnm_hor) < 1.0e-10
        and min(control_inf, control_hor) > 1.0e-8
        and all(curve["checks"].values())
        else "failed"
    )
    summary = {
        "schema": "flintnde_qnm_2x2_example_v3",
        "status": status,
        "package_root": str(PACKAGE_ROOT),
        "output_directory": str(output_layout.run_root),
        "ConfigPath": convention["ConfigPath"],
        "ConventionName": convention["ConventionName"],
        "it0": convention["it0"],
        "it1": convention["it1"],
        "EnergyCondition": convention["EnergyCondition"],
        "equation": "unified-u exact rational 2x2 first-order system in x",
        "ell": ELL,
        "boundary_test": "allowed literal singular endpoint branch transports into the allowed subspace at the other endpoint",
        "integration": {
            "backend": "FlintNDE exact local structure plus Acb Cauchy-DFT Taylor transport",
            "primary_order": PRIMARY_ORDER,
            "reference_order": REFERENCE_ORDER,
            "primary_samples": PRIMARY_SAMPLES,
            "reference_samples": REFERENCE_SAMPLES,
            "target_relative_error": TARGET_RELATIVE_ERROR,
            "formal_minimum_order_factor": FORMAL_MINIMUM_ORDER_FACTOR,
            "cauchy_radius_fraction": CAUCHY_RADIUS_FRACTION,
            "x_match": str(X_MATCH),
            "infinity_start": "inf",
            "infinity_match_step_over_radius": INFINITY_MATCH_RATIO,
            "horizon_start": "1",
            "horizon_match_step_over_radius": HORIZON_MATCH_RATIO,
        },
        "qnm": qnm,
        "control": control,
        "curve": curve,
    }
    output_path = output_layout.write_json("summary", "qnm_2x2_summary.json", summary)
    curve_path = output_layout.write_json(
        "transport", "qnm_2x2_curve.json", {key: curve[key] for key in curve if key != "samples"}
    )
    samples_path = output_layout.write_json(
        "transport", "qnm_2x2_curve_samples.json", curve["samples"]
    )
    print(
        json.dumps(
            {
                "status": status,
                "qnm_inf": qnm_inf,
                "qnm_hor": qnm_hor,
                "control_inf": control_inf,
                "control_hor": control_hor,
                "curve": {
                    "points": curve["rstar_grid"]["total_points"],
                    "nodes": curve["planning"]["node_count"],
                    "largest_bucket": curve["planning"]["largest_bucket"],
                    "chain_relative_difference": curve["precision"][
                        "chain_relative_difference"
                    ],
                    "x_match_relative_difference": curve["precision"][
                        "x_match_relative_difference"
                    ],
                    "checks": curve["checks"],
                },
            },
            indent=2,
        )
    )
    if status != "passed":
        raise RuntimeError(
            f"QNM 2x2 example failed: {output_path} / {curve_path} / {samples_path}"
        )


if __name__ == "__main__":
    main()
