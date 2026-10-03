"""演示 FlintNDE 对不同固定 ep 取值进行有界多进程 NDE 计算。

方程为 ``y'(x)=ep/(1+x) y(x)``、``y(0)=1``，故 ``y(1)=2^ep``。
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
        primary_order=96,
        reference_order=112,
        target_relative_error="1e-40",
    )
    value = result["primary_snapshots"][-1][0, 0]
    expected = acb(2) ** acb(ep_value)
    return {
        "ep": ep_value,
        "value": value.str(45),
        "expected": expected.str(45),
        "absolute_difference": abs(value - expected).str(20),
        "passed": bool(abs(value - expected) < arb("1e-30")),
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
