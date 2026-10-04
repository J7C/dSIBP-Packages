"""演示 FlintNDE 对不同固定 ep 取值进行有界多进程 NDE 计算。

方程为 ``y'(x)=ep/(1+x) y(x)``、``y(0)=1``，故 ``y(1)=2^ep``。
阶数与目标相对误差在下面的模块常量处声明一次；``passed`` 检查后端自己的
``target_relative_error_met``，并用同一个目标与闭式解 ``2^ep`` 比相对差。
``parallel_task_count`` 缺省为 12；实际 worker 数由程序自动取该值与 ep 数量的较小者，
任务更多时完成一个自动续交一个。任务函数必须定义在模块顶层以支持 Windows spawn。
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any


EXAMPLE_ROOT = Path(__file__).resolve().parent
PACKAGE_ROOT = EXAMPLE_ROOT.parent
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from flint import acb, arb  # noqa: E402
from flintnde import (  # noqa: E402
    RationalMatrixSystem,
    NamedPoint,
    build_adaptive_path,
    column_vector,
    configure_working_precision,
    rational_function,
    run_ep_tasks,
    transport_path_refined,
)


# 阶数与目标相对误差只在这里声明一次，闭式对照用的就是同一个目标；再另写一个
# 更松或更紧的阈值，会在用户改目标后检查一个与该示例无关的量。
PRIMARY_ORDER = 112
REFERENCE_ORDER = 128
TARGET_RELATIVE_ERROR = "1e-40"


def solve_ep(ep_value: str) -> dict[str, Any]:
    """在 worker 内构造 FLINT 对象并返回可进程传输的字符串摘要。"""

    configure_working_precision(60)
    system = RationalMatrixSystem(
        ((rational_function((ep_value,), (1, 1)),),),
        variable_name="x",
        name=f"ep-parallel-{ep_value}",
    )
    path = build_adaptive_path(
        system,
        NamedPoint("start", acb(0)),
        NamedPoint("target", acb(1)),
    )
    result = transport_path_refined(
        system,
        column_vector([1]),
        path,
        primary_order=PRIMARY_ORDER,
        reference_order=REFERENCE_ORDER,
        target_relative_error=TARGET_RELATIVE_ERROR,
    )
    value = result["primary_snapshots"][-1][0, 0]
    expected = acb(2) ** acb(ep_value)
    relative_error = abs(value - expected) / abs(expected)
    return {
        "ep": ep_value,
        "value": value.str(45),
        "expected": expected.str(45),
        "absolute_difference": abs(value - expected).str(20),
        "relative_difference": relative_error.str(20),
        # 对照的就是上面声明的那个目标：后端自认达标，且与闭式解的相对差同量级。
        "passed": bool(
            result["target_relative_error_met"]
            and relative_error < arb(TARGET_RELATIVE_ERROR)
        ),
    }


def main() -> None:
    """运行输入同序的 ep 任务；修改 parallel_task_count 即可设置并行上限。"""

    ep_values = ("1/5", "1/6", "1/7")
    parallel_task_count = 12  # 缺省值；例如改为 4 即最多同时运行四个 ep。
    batch = run_ep_tasks(
        ep_values,
        solve_ep,
        parallel_task_count=parallel_task_count,
    )
    print(
        "parallel requested/effective: "
        f"{batch.parallel_task_count_requested}/"
        f"{batch.parallel_task_count_effective}"
    )
    for result in batch.results:
        print(result)
    if not all(result["passed"] for result in batch.results):
        raise SystemExit(1)
    print("ep_parallel.py: PASSED")


if __name__ == "__main__":
    main()
