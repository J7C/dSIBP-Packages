#!/usr/bin/env python3
"""演示多项式指数相位 ``{phi,a,b,C}`` 边界：两个 sector 各含两项负幂。

系统为对角矩阵 ``A = diag(A11, A22)``，其中

    A11 = 3/(2 z^3) + 1/z^2 + (1/2)/z + 2,
    A22 = 1/z^3 + (1/3)/z^2 + (7/2)/z - 1.

``z = 0`` 是三阶极点（本性奇点）。两个对角元的高阶 Laurent 系数不同，因此拆成
两个指数 sector，指数多项式分别为 ``Phi_1 = -3/(4 z^2) - 1/z``（幂 ``a = 1/2``）与
``Phi_2 = -1/(2 z^2) - 1/(3 z)``（幂 ``a = 7/2``）；闭式解在 ``z = 1`` 为
``(e^{1/4}, e^{-11/6})``。示例用两条 ``{phi,a,b,C}`` 记录分别选中两个 sector，
输运到普通点并与闭式比较。
"""

from __future__ import annotations

import sys
from pathlib import Path

from flint import acb, arb


# 仓库内示例集中设置源码路径；安装 FlintNDE 后可删除本段并直接 import。
EXAMPLE_DIR = Path(__file__).resolve().parent
# 发布布局直接取上级目录；发展仓布局取 versions/ 下唯一的当前版本目录。
PACKAGE_ROOT = (EXAMPLE_DIR / "..").resolve()
if not (PACKAGE_ROOT / "flintnde").is_dir():
    _version_dirs = sorted((PACKAGE_ROOT / "versions").glob("FlintNDE-*"))
    if len(_version_dirs) != 1:
        raise RuntimeError(
            "expected exactly one versions/FlintNDE-* directory beside examples/"
        )
    PACKAGE_ROOT = _version_dirs[0].resolve()
OUTPUT_DIR = EXAMPLE_DIR / "results" / "exponential_polynomial_boundary"
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from flintnde import (  # noqa: E402
    RationalMatrixSystem,
    build_adaptive_path,
    configure_working_precision,
    exponential_boundary,
    rational_function,
    transport_path_refined,
)


def main() -> None:
    """从三阶极点输运双 sector 指数边界，并与闭式比较。"""

    configure_working_precision(50, 24)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # 分母为 z^3；分子按幂次升序，因此第 k 项对应 Laurent 幂 z^(k-3)：
    # 想要 Laurent 展开 (3/2)/z^3 + 1/z^2 + (1/2)/z + 2，分子须写
    # (3/2) + 1*z + (1/2)*z^2 + 2*z^3。
    a11 = rational_function(("3/2", 1, "1/2", 2), (0, 0, 0, 1))
    a22 = rational_function((1, "1/3", "7/2", "-1"), (0, 0, 0, 1))
    system = RationalMatrixSystem(
        ((a11, 0), (0, a22)),
        variable_name="z",
        name="exponential-polynomial-boundary-example",
    )

    # phi 系数 = -c_k/(k-1)：sector 1 由 c_3=3/2、c_2=1 给出 -3/4 与 -1；
    # sector 2 由 c_3=1、c_2=1/3 给出 -1/2 与 -1/3。C 分别选中各自 sector。
    boundary = exponential_boundary(
        [
            {
                "phi": [
                    {"power": -2, "coefficient": "-3/4"},
                    {"power": -1, "coefficient": "-1"},
                ],
                "a": "1/2",
                "b": 0,
                "C": [1, 0],
            },
            {
                "phi": [
                    {"power": -2, "coefficient": "-1/2"},
                    {"power": -1, "coefficient": "-1/3"},
                ],
                "a": "7/2",
                "b": 0,
                "C": [0, 1],
            },
        ]
    )

    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        path = build_adaptive_path(
            system,
            (0, "save"),
            (1, "save"),
            max_step_over_radius=0.25,
        )

    result = transport_path_refined(
        system,
        boundary,
        path,
        primary_order=48,
        reference_order=64,
        radius_fraction=0.45,
        target_relative_error="1e-14",
        save_output_directory=OUTPUT_DIR,
    )

    final = result["primary_snapshots"][-1]
    expected1 = (acb(1) / 4).exp()
    expected2 = (acb(-11) / 6).exp()
    relative1 = abs(final[0, 0] - expected1) / abs(expected1)
    relative2 = abs(final[1, 0] - expected2) / abs(expected2)
    checks = {
        "closed_form_sector_1": bool(relative1 < arb("1e-20")),
        "closed_form_sector_2": bool(relative2 < arb("1e-20")),
        "target_relative_error_met": bool(result["target_relative_error_met"]),
    }
    print("final value:", final.str(50))
    print("relative differences vs closed form:", relative1.str(10), relative2.str(10))
    print("checks:", checks)
    print("saved files:", OUTPUT_DIR)
    if not all(checks.values()):
        raise SystemExit(1)
    print("exponential_polynomial_boundary.py: PASSED")


if __name__ == "__main__":
    main()
