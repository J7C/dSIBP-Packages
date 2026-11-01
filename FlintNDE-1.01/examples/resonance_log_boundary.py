#!/usr/bin/env python3
"""演示共振对数：无 Jordan 块时 log 由 indicial 根的整数差触发。

系统 ``y1' = 0``、``y2' = y1 + y2/x``，即 ``A = [[0, 0], [1, 1/x]]``，
``x = 0`` 是正则奇点。留数矩阵 ``R = diag(0, 1)`` 无 Jordan 块，indicial 根
``{0, 1}`` 相差整数 1：``rho = 0`` 支在 Taylor 阶 ``n = 1`` 撞上根 1，
且 ``A_0 v_0 = [0, 1]`` 在根 1 本征方向上的投影非零（谱投影相容性不满足），
log 因此首现于 ``x^1 log x`` 项。闭式解为
``y1 = d``、``y2 = c x + d x log x``。示例分别输运两条 ``{a,b,C}`` 边界
记录——``{a=1, b=0, C=[0,1]}`` 选中纯幂支 ``[0, x]``、
``{a=0, b=0, C=[1,0]}`` 选中共振支 ``[1, x log x]``——到普通点 ``x = 2``
并与闭式比较。
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

from flint import acb, arb


# 仓库内示例集中设置源码路径；安装 FlintNDE 后可删除本段并直接 import。
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
OUTPUT_DIR = EXAMPLE_DIR / "results" / "resonance_log_boundary"
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from flintnde import (  # noqa: E402
    RationalMatrixSystem,
    build_adaptive_path,
    configure_working_precision,
    frobenius_boundary,
    rational_function,
    transport_path_refined,
)


def main() -> None:
    """分别输运纯幂支与共振支到 x=2，并与闭式 ``[0,2]``、``[1,2 log 2]`` 比较。"""

    configure_working_precision(50, 24)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    zero = rational_function(0)
    one = rational_function(1)
    one_over_x = rational_function(1, (0, 1))
    system = RationalMatrixSystem(
        ((zero, zero), (one, one_over_x)),
        variable_name="x",
        name="resonance-log-boundary-example",
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        path = build_adaptive_path(
            system,
            (0, "save"),
            (2, "save"),
            max_step_over_radius=0.30,
        )

    branches = (
        ("power_branch", [{"a": 1, "b": 0, "C": [0, 1]}]),
        ("resonance_branch", [{"a": 0, "b": 0, "C": [1, 0]}]),
    )
    log2 = acb(2).log()
    expected = {
        "power_branch": (acb(0), acb(2)),
        "resonance_branch": (acb(1), acb(2) * log2),
    }
    checks = {}
    for label, records in branches:
        result = transport_path_refined(
            system,
            frobenius_boundary(records),
            path,
            primary_order=48,
            reference_order=64,
            target_relative_error="1e-14",
            save_output_directory=OUTPUT_DIR / label,
        )
        final = result["primary_snapshots"][-1]
        print(f"{label} final value:", final.str(30))
        for row, value in enumerate(expected[label]):
            # abs(acb) 返回 arb 模长；零分量的期望值改用绝对差，避免 0/0。
            difference = abs(final[row, 0] - value)
            scale = abs(value)
            relative = difference / scale if not scale.is_zero() else difference
            print(f"{label} component {row} vs closed form:", relative.str(5))
            checks[f"{label}_component_{row}"] = bool(relative < arb("1e-20"))
        checks[f"{label}_target_met"] = bool(result["target_relative_error_met"])

    print("checks:", checks)
    print("saved files:", OUTPUT_DIR)
    if not all(checks.values()):
        raise SystemExit(1)
    print("resonance_log_boundary.py: PASSED")


if __name__ == "__main__":
    main()
