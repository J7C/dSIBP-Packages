"""演示带 pole 的 ``ep`` Laurent 级数自适应重构。

本例的零矩阵 DE 保持边界值不变；符号边界
``1/ep + 2 + 3 ep + 4 ep^2 + 5 ep^3 + 6 ep^4`` 明确认证最低阶为 -1。
FlintNDE 只消费该证书，不从数值样本猜最低阶。用户只请求有限项，本例提供冗余的
生产候选池和独立验证点，并把首轮内部最高幂设为 2。首轮只计算所需候选前缀；独立
验证失败时再追加两点且复用已有生产点和验证点，绝不取用户池外的 regulator 值。
随后用 ``sample_angle_range=(-pi/3,pi/3)`` 再跑自动网格：程序只在该开角域内部
均匀选择三条射线，模长仍由目标精度自动决定且没有用户模长上限参数。

并行：``parallel_task_count`` 缺省为 12；实际进程数自动取任务数与 12 的较小者。
Windows 多进程要求 DE、边界和路径工厂都定义在模块顶层。
"""

from __future__ import annotations

import sys
import math
from pathlib import Path


EXAMPLE_ROOT = Path(__file__).resolve().parent
PACKAGE_ROOT = EXAMPLE_ROOT.parent
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from flint import acb, acb_mat, arb, fmpq  # noqa: E402
from flintnde import AnalyticMatrixSystem, reconstruct_series_solution  # noqa: E402


def zero_system(_ep: fmpq) -> AnalyticMatrixSystem:
    """返回一维零连接；远端哑奇点只用于限定普通点 Taylor 圆。"""

    return AnalyticMatrixSystem(
        lambda _z: acb_mat(1, 1),
        1,
        (acb(-10), acb(10)),
        "ep-series-zero-system",
    )


def laurent_boundary(ep: fmpq) -> list[acb]:
    """返回符号上已知最低阶为 ``ep^-1`` 的普通点边界。"""

    regulator = acb(ep)
    return [
        1 / regulator
        + 2
        + 3 * regulator
        + 4 * regulator**2
        + 5 * regulator**3
        + 6 * regulator**4
    ]


def fixed_path(_ep: fmpq, _system: AnalyticMatrixSystem) -> list[acb]:
    """返回每个正规化样本共用的普通点路径。"""

    return [acb(0), acb(1)]


def main() -> None:
    """重构 pole 与有限项，并检查增量扩阶和缓存复用诊断。"""

    production_candidates = (
        "1/10", "9/100", "2/25", "7/100", "3/50", "1/20", "1/25"
    )
    validation_points = ("3/100", "1/50")
    result = reconstruct_series_solution(
        DEmatrix=zero_system,
        boundary=laurent_boundary,
        path=fixed_path,
        maximum_power=0,
        leading_power=-1,
        leading_power_certificate={
            "status": "certified",
            "leading_power": -1,
            "method": "analytic-example-boundary-and-zero-DE",
        },
        goal_digits=12,
        sample_points=production_candidates,
        validation_points=validation_points,
        initial_internal_maximum_power=2,
        validation_tolerance="1e-20",
        transport_order=4,
        transport_extra_order=2,
        parallel_task_count=12,
    )
    angle_result = reconstruct_series_solution(
        DEmatrix=zero_system,
        boundary=laurent_boundary,
        path=fixed_path,
        maximum_power=0,
        leading_power=-1,
        leading_power_certificate={
            "status": "certified",
            "leading_power": -1,
            "method": "analytic-example-boundary-and-zero-DE",
        },
        goal_digits=12,
        sample_angle_range=(str(-math.pi / 3), str(math.pi / 3)),
        transport_order=4,
        transport_extra_order=2,
        parallel_task_count=12,
    )

    pole = result.coefficient(-1)[0, 0]
    finite_part = result.coefficient(0)[0, 0]
    history = result.diagnostics["fit_expansion_history"]
    checks = {
        "certified_leading_power": result.leading_power == -1,
        "pole_coefficient": bool(abs(pole - 1) < arb("1e-30")),
        "finite_part": bool(abs(finite_part - 2) < arb("1e-30")),
        "extra_fit_powers": result.internal_maximum_power >= 2,
        "independent_validation": history[-1]["passed"],
        "incremental_reuse": any(
            item["reused_production_sample_count"] > 0
            and item["reused_validation_sample_count"] > 0
            for item in history[1:]
        ),
        "candidate_prefix_only":
            result.effective_parameters["sample_candidate_count"] == 7
            and result.effective_parameters["unused_sample_candidate_count"] == 1,
        "user_range_respected": all(
            any(abs(point - acb(candidate)).contains(0) for candidate in production_candidates)
            for point in result.sample_points
        ),
        "default_parallel_limit":
            result.effective_parameters["parallel_task_count_requested"] == 12,
        "angle_range_precision_certified":
            angle_result.effective_parameters["precision_target_met"],
        "angle_range_open_and_three_rays": (
            all(
                -math.pi / 3 < math.atan2(
                    float(point.imag.mid()), float(point.real.mid())
                ) < math.pi / 3
                for point in angle_result.sample_points
            )
            and len({
                round(math.atan2(float(point.imag.mid()), float(point.real.mid())), 12)
                for point in angle_result.sample_points
            }) == 3
        ),
    }

    print(f"certified powers returned: {result.leading_power}..{result.maximum_power}")
    print(f"internal maximum power: {result.internal_maximum_power}")
    print(f"pole coefficient: {pole.str(30)}")
    print(f"finite part: {finite_part.str(30)}")
    print(f"fit expansion history: {history}")
    print(f"production candidates: {production_candidates}")
    print(f"used production points: {tuple(point.str(20) for point in result.sample_points)}")
    print(
        "automatic angle-range points: "
        f"{tuple(point.str(20) for point in angle_result.sample_points)}"
    )
    print(f"checks: {checks}")
    if not all(checks.values()):
        raise SystemExit(1)
    print("ep_series_reconstruction.py: PASSED")


if __name__ == "__main__":
    main()
