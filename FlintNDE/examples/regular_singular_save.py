#!/usr/bin/env python3
"""展示正则奇点 ``{a,b,C}`` 边界、refinement 和保存点输出。

系统 ``y'=2 y/x`` 在 ``x=0`` 有正则奇点。输入边界表示
``y(x)=3 x^2`` 的领头项；路径起点和普通终点都用无名 ``save`` 标签保存。
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path


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
OUTPUT_DIR = EXAMPLE_DIR / "results" / "regular_singular_save"
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
    """从正则奇点输运到普通点，并在调用者目录保存两端数据。"""

    configure_working_precision(50, 24)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    system = RationalMatrixSystem(
        ((rational_function(2, [0, 1]),),),
        variable_name="x",
        name="regular-singular-save-example",
    )
    boundary = frobenius_boundary([{"a": 2, "b": 0, "C": [3]}])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        path = build_adaptive_path(
            system,
            (0, "save"),
            (1, "save"),
            max_step_over_radius=0.30,
        )

    result = transport_path_refined(
        system,
        boundary,
        path,
        primary_order=20,
        reference_order=28,
        target_relative_error="1e-14",
        save_output_directory=OUTPUT_DIR,
    )

    print("final value:", result["primary_snapshots"][-1][0, 0])
    print("relative difference:", result["relative_difference_inf"])
    print("saved files:", OUTPUT_DIR)


if __name__ == "__main__":
    main()
