"""形式渐近级数的公开尾项诊断。

本模块只分析 caller 已生成的标量级数项，不推断奇点类型或选择截断阶。返回对象保留
分子、分母、比值、阈值和严格门禁结果，供下游在写出物理量前 fail closed。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

from flint import acb, arb


@dataclass(frozen=True)
class FiveTermTailDiagnostic:
    """保存固定截断阶前后各五项的标量尾项诊断。"""

    evaluation_degree: int
    previous_five_absolute_sum: arb
    next_five_sum_absolute: arb
    next_over_previous: arb
    threshold: arb
    passed: bool

    def as_dict(self, digits: int = 30) -> dict[str, Any]:
        """返回可直接写入 JSON 的稳定字段。"""

        return {
            "evaluation_degree": self.evaluation_degree,
            "lookahead_order_count": 5,
            "numerator": self.next_five_sum_absolute.str(digits),
            "denominator": self.previous_five_absolute_sum.str(digits),
            "ratio": self.next_over_previous.str(digits),
            "threshold": self.threshold.str(digits),
            "passed": self.passed,
            "comparison": "strict_less_than",
            "definition": "abs(sum(T[N+1:N+6]))/sum(abs(T[N-4:N+1]))",
        }


def five_term_tail_diagnostic(
    terms: Sequence[acb | int | str],
    evaluation_degree: int,
    *,
    threshold: arb | int | float | str = "0.1",
) -> FiveTermTailDiagnostic:
    """计算 ``abs(sum(next five))/sum(abs(previous five))`` 严格门禁。

    ``terms[n]`` 必须是级数第 ``n`` 项在实际求值点的完整标量贡献。分母包含零时
    无法认证相对尾项，函数直接失败；区间比较不能严格证明小于阈值时也返回失败。
    """

    if evaluation_degree < 4:
        raise ValueError("five-term diagnostic requires evaluation_degree >= 4")
    required_count = evaluation_degree + 6
    if len(terms) < required_count:
        raise ValueError(
            "five-term diagnostic requires terms through degree "
            f"{evaluation_degree + 5}, got {len(terms) - 1}"
        )
    threshold_value = threshold if isinstance(threshold, arb) else arb(str(threshold))
    if not threshold_value > 0:
        raise ValueError("five-term diagnostic threshold must be positive")

    converted = [value if isinstance(value, acb) else acb(value) for value in terms]
    previous = converted[evaluation_degree - 4 : evaluation_degree + 1]
    following = converted[evaluation_degree + 1 : evaluation_degree + 6]
    denominator = sum((abs(value) for value in previous), arb(0))
    if denominator.contains(0):
        raise ZeroDivisionError("previous five absolute-term sum contains zero")
    numerator = abs(sum(following, acb(0)))
    ratio = numerator / denominator
    return FiveTermTailDiagnostic(
        evaluation_degree=evaluation_degree,
        previous_five_absolute_sum=denominator,
        next_five_sum_absolute=numerator,
        next_over_previous=ratio,
        threshold=threshold_value,
        passed=bool(ratio < threshold_value),
    )
