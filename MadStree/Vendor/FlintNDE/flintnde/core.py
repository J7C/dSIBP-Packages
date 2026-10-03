"""FlintNDE 的 FLINT 数值、精度和矩阵公共接口。

本模块集中管理工作精度、列向量构造、中点重启和范数比较。高层输运与局部级数
模块只通过这些接口处理数值约定，避免各自隐式降到 binary64。
"""

from __future__ import annotations

import math
import re
from fractions import Fraction
from typing import Any

from flint import acb, acb_mat, arb, ctx, fmpq


DEFAULT_WORKING_PRECISION_DIGITS = 200


def require_exact_keys(
    record: Any,
    expected_keys: set[str] | frozenset[str],
    field_name: str,
) -> dict[str, Any]:
    """要求 JSON 对象严格采用当前 schema，不接受缺字段或额外字段。"""

    if not isinstance(record, dict):
        raise ValueError(f"{field_name} must be an object")
    observed = set(record)
    missing = sorted(expected_keys - observed)
    unexpected = sorted(observed - expected_keys)
    if missing or unexpected:
        details = []
        if missing:
            details.append(f"missing={missing}")
        if unexpected:
            details.append(f"unexpected={unexpected}")
        raise ValueError(
            f"{field_name} does not match the current schema: " + "; ".join(details)
        )
    return record


def configure_working_precision(
    decimal_digits: int = DEFAULT_WORKING_PRECISION_DIGITS,
    guard_bits: int = 32,
) -> int:
    """设置 Acb 工作精度并返回实际二进制位数；缺省使用 200 位十进制精度。"""

    if decimal_digits <= 0 or guard_bits < 0:
        raise ValueError("decimal_digits must be positive and guard_bits nonnegative")
    ctx.prec = math.ceil(decimal_digits * math.log2(10)) + guard_bits
    return ctx.prec


# 导入 FlintNDE 时建立确定的 200 位缺省；公开配置函数仍允许用户显式覆盖。
configure_working_precision()


def arb_ball_to_json(value: arb, decimal_digits: int) -> dict[str, Any]:
    """把 Arb 球写成保证覆盖原区间的十进制中点、半径和指数。"""

    if decimal_digits <= 0:
        raise ValueError("decimal_digits must be positive")
    midpoint, radius, exponent = value.mid_rad_10exp(decimal_digits)
    return {
        "midpoint": str(midpoint),
        "radius": str(radius),
        "exponent": int(exponent),
    }


def arb_ball_from_json(record: Any, field_name: str) -> arb:
    """从中点--半径记录恢复 Arb 球，并拒绝含混或负半径输入。"""

    require_exact_keys(record, {"midpoint", "radius", "exponent"}, field_name)
    try:
        midpoint = int(record["midpoint"])
        radius = int(record["radius"])
        exponent = int(record["exponent"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(
            f"{field_name} must contain integer midpoint, radius and exponent"
        ) from error
    if radius < 0:
        raise ValueError(f"{field_name} radius must be nonnegative")
    return arb(f"[{midpoint}e{exponent} +/- {radius}e{exponent}]")


def exact_rational(value: Any) -> fmpq:
    """把整数、分数字符串或十进制字符串精确转换为 FLINT ``fmpq``。

    Python ``float`` 只按其十进制打印值转为有理数；需要保留输入误差时应使用 Acb
    ball，而不应把 float 传给 exact Frobenius gate。
    """

    if isinstance(value, fmpq):
        return value
    text = re.sub(r"`(?:\d+(?:\.\d*)?|\.\d+)?", "", str(value)).replace("*^", "e")
    try:
        rational = Fraction(text)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError(f"expected an exact rational scalar, got {value}") from error
    return fmpq(rational.numerator, rational.denominator)


def exact_to_acb(value: Any) -> acb:
    """把 exact 有理实数转换为 Acb，不经过 machine number。"""

    return acb(exact_rational(value))


def column_vector(values: list[acb | int | str]) -> acb_mat:
    """由标量列表构造 Acb 列向量。"""

    return acb_mat([[value if isinstance(value, acb) else exact_to_acb(value)] for value in values])


def identity_matrix(dimension: int) -> acb_mat:
    """构造指定维数的 Acb 单位矩阵。"""

    return acb_mat(
        [[acb(1 if row == column else 0) for column in range(dimension)] for row in range(dimension)]
    )


def acb_midpoint_matrix(matrix: acb_mat) -> acb_mat:
    """以当前高精度中点重建矩阵，控制多段输运的区间包裹膨胀。"""

    return acb_mat(
        [
            [acb(matrix[row, column].real.mid(), matrix[row, column].imag.mid()) for column in range(matrix.ncols())]
            for row in range(matrix.nrows())
        ]
    )


def _arb_midpoint(value: arb) -> float:
    """只为排序或摘要读取 Arb 中点。"""

    return float(value.mid())


def vector_norm_inf(vector: acb_mat) -> arb:
    """计算 Acb 列向量的无穷范数。"""

    if vector.ncols() != 1:
        raise ValueError("vector_norm_inf requires a column vector")
    values = [abs(vector[row, 0]) for row in range(vector.nrows())]
    return arb(0) if not values else max(values, key=_arb_midpoint)


def matrix_norm_inf(matrix: acb_mat) -> arb:
    """计算 Acb 矩阵的最大行和范数。"""

    rows = [
        sum((abs(matrix[row, column]) for column in range(matrix.ncols())), arb(0))
        for row in range(matrix.nrows())
    ]
    return arb(0) if not rows else max(rows, key=_arb_midpoint)


def relative_difference_inf(primary: acb_mat, reference: acb_mat) -> arb:
    """计算同维列向量或矩阵的无穷范数相对差。

    对列向量（``ncols==1``）与 :func:`vector_norm_inf` 行为一致；对多列矩阵
    使用 :func:`matrix_norm_inf` 的最大行和范数，使 batch Frobenius 输运的
    多列 snapshots 也能通过嵌入式截断认证路径。
    """

    if primary.nrows() != reference.nrows() or primary.ncols() != reference.ncols():
        raise ValueError("relative-difference dimensions do not match")
    denominator = matrix_norm_inf(reference)
    if denominator.contains(0):
        raise ZeroDivisionError("reference matrix norm contains zero")
    return matrix_norm_inf(reference - primary) / denominator
