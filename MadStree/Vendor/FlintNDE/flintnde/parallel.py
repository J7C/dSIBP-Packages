"""为互相独立的固定 regulator 任务提供有界多进程调度。

本模块只并行不同 regulator/ep 取值，不并行一条路径上前后依赖的输运节点。
缺省最多同时运行 12 个任务；任务更少时只启动任务数对应的 worker，任务更多时
由执行器在任一 worker 完成后自动续交下一个。独立进程隔离 python-flint 的全局上下文。
"""

from __future__ import annotations

from dataclasses import dataclass
from multiprocessing import get_context
from typing import Any, Callable, Sequence


DEFAULT_PARALLEL_TASK_COUNT = 12


@dataclass(frozen=True)
class EpTaskBatchResult:
    """保存输入同序结果和实际使用的并行调度参数。"""

    ep_values: tuple[Any, ...]
    results: tuple[Any, ...]
    parallel_task_count_requested: int
    parallel_task_count_effective: int


def _run_indexed_ep_task(
    task: Callable[[Any], Any], index: int, ep_value: Any
) -> tuple[int, Any]:
    """在 worker 中运行一个 ep 任务，并携带输入序号返回。"""

    return index, task(ep_value)


def run_ep_tasks(
    ep_values: Sequence[Any],
    task: Callable[[Any], Any],
    *,
    parallel_task_count: int = DEFAULT_PARALLEL_TASK_COUNT,
) -> EpTaskBatchResult:
    """按输入顺序返回独立 ep 任务结果，缺省并行数为 12。

    ``task`` 必须是可由 Python 多进程导入的模块顶层函数，输入 ep 和返回值也必须可
    pickle。NDE 任务应在 worker 内构造并使用 ``acb/acb_mat``，返回时改成字符串、数字、
    list/dict 或其它可传输摘要。若只需调试不可序列化的局部函数，可显式设为 1 串行运行。
    """

    values = tuple(ep_values)
    if not values:
        raise ValueError("ep_values must be a nonempty sequence")
    if not callable(task):
        raise TypeError("task must be callable")
    if (
        isinstance(parallel_task_count, bool)
        or not isinstance(parallel_task_count, int)
        or parallel_task_count < 1
    ):
        raise ValueError("parallel_task_count must be a positive integer")
    effective_count = min(parallel_task_count, len(values))
    print(
        "FlintNDE ep task pool: "
        f"requested={parallel_task_count}, effective={effective_count}; "
        "default parallel_task_count=12."
    )
    payloads = [(task, index, value) for index, value in enumerate(values)]
    with get_context("spawn").Pool(
        processes=effective_count, maxtasksperchild=1
    ) as pool:
        completed = pool.starmap(_run_indexed_ep_task, payloads, chunksize=1)
    ordered = tuple(result for _index, result in completed)
    return EpTaskBatchResult(
        ep_values=values,
        results=ordered,
        parallel_task_count_requested=parallel_task_count,
        parallel_task_count_effective=effective_count,
    )
