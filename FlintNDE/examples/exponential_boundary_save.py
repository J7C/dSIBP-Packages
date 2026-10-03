#!/usr/bin/env python3
"""展示已认证指数型奇点 ``{phi,a,b,C}`` 边界的保存与复用。

标量系统 ``y'=y/x^2`` 的解为 ``C exp(-1/x)``。示例显式输入
``phi=-1/x``，保存奇点边界和普通终点；保存的 terms 可直接再次交给
``exponential_boundary``。
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path


# 仓库内示例集中设置源码路径；安装 FlintNDE 后可删除本段并直接 import。
EXAMPLE_DIR = Path(__file__).resolve().parent
PACKAGE_ROOT = (EXAMPLE_DIR / "..").resolve()
OUTPUT_DIR = EXAMPLE_DIR / "results" / "exponential_boundary_save"
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from flintnde import (  # noqa: E402
    RationalMatrixSystem,
    build_adaptive_path,
    configure_working_precision,
    exponential_boundary,
    rational_function,
    transport_path,
)


def main() -> None:
    """输运指数型边界，并验证保存 JSON 能重建同一边界对象。"""

    configure_working_precision(50, 24)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    system = RationalMatrixSystem(
        ((rational_function(1, [0, 0, 1]),),),
        variable_name="x",
        name="exponential-boundary-save-example",
    )
    boundary = exponential_boundary(
        [{"phi": [{"power": -1, "coefficient": -1}], "a": 0, "b": 0, "C": [3]}]
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        path = build_adaptive_path(
            system,
            (0, "save"),
            (1, "save"),
            max_step_over_radius=0.25,
        )

    snapshots, _reports, _elapsed = transport_path(
        system,
        boundary,
        path,
        order=32,
        sample_count=96,
        save_output_directory=OUTPUT_DIR,
    )

    first_record = json.loads(
        sorted(OUTPUT_DIR.glob("flintnde_save_[0-9][0-9][0-9].json"))[0].read_text(
            encoding="utf-8"
        )
    )
    reused = exponential_boundary(first_record["result"]["terms"])

    print("final value:", snapshots[-1][0, 0])
    print("reused boundary:", reused.to_json())
    print("saved files:", OUTPUT_DIR)


if __name__ == "__main__":
    main()
