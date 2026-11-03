"""letter list 输入的认证根隔离、部分分式装配与 letter 级无穷远反演。

letter 记录 ``{"numerator": [...], "denominator": [...], "residue": [[...]]}`` 描述
``d Log[N(x)/D(x)]·R`` 的 dlog 贡献。实现思路：先在 letter 层完成根隔离（代数根
一律取认证小数包络，不进符号层）与整块留数装配，输出可作输运快速系统的
``PartialFractionSystem``；``x→1/t`` 反演同样只对 letter 系数做倒数多项式替换并
补一条 ``dlog t`` letter 吸收幂次差留数，不对整矩阵做变量代换。
"""

from __future__ import annotations

from typing import Any, Sequence

from flint import acb, acb_mat, acb_poly, arb, ctx

from .core import require_exact_keys
from .exact_gaussian import GaussianRational, gaussian_rational
from .systems import PartialFractionSystem


def _trimmed(coefficients: list[GaussianRational]) -> list[GaussianRational]:
    """去掉系数表末尾的零系数，使最高次系数非零（度数按实际值计）。"""

    while len(coefficients) > 1 and coefficients[-1].is_zero:
        coefficients.pop()
    return coefficients


def isolated_roots(coefficients: Sequence[GaussianRational]) -> list[tuple[acb, int]] | None:
    """把多项式系数表隔离成 (root, multiplicity)；任一环节不满足返回 None。

    先剪尾部零系数，零根再按首部零系数个数精确剥离；一次因子用精确 Q(i)
    除法；高次因子用 acb 球根隔离并验证 ``N(r)`` 含 0、根数等于次数
    （squarefree）。代数根一律以认证小数包络表示，不进入符号层。
    """

    core = _trimmed(list(coefficients))
    zero_multiplicity = 0
    while zero_multiplicity < len(core) - 1 and core[zero_multiplicity].is_zero:
        zero_multiplicity += 1
    core = core[zero_multiplicity:]
    degree = len(core) - 1
    if degree <= 0:
        roots = []
    elif degree == 1:
        roots = [(-core[0] / core[1]).to_acb()]
    else:
        polynomial = acb_poly([value.to_acb() for value in core])
        # tol 达到 prec 一半量级时，roots 内部默认 maxprec 不够，需显式放大重试。
        tolerance = arb(10) ** (-(ctx.prec // 2))
        isolated = None
        for maxprec in (4 * ctx.prec, 8 * ctx.prec):
            try:
                isolated = polynomial.roots(tolerance, maxprec)
                break
            except (ValueError, ZeroDivisionError):
                continue
        if isolated is None or len(isolated) != degree:
            return None
        for root in isolated:
            if not polynomial(root).contains(0):
                return None
        roots = isolated
    result = []
    if zero_multiplicity:
        result.append((acb(0), zero_multiplicity))
    result.extend((root, 1) for root in roots)
    return result


def _letter_coefficients(record: dict[str, Any], label: str) -> tuple[list[GaussianRational], list[GaussianRational]]:
    """解析并修剪一条 letter 的分子/分母精确系数表。"""

    numerator = [gaussian_rational(item) for item in record["numerator"]]
    denominator = [gaussian_rational(item) for item in record["denominator"]]
    if not numerator or not denominator:
        raise ValueError(f"{label} numerator and denominator must be nonempty")
    return _trimmed(numerator), _trimmed(denominator)


def _letter_matrix(record: dict[str, Any], dimension: int, label: str) -> list[list[GaussianRational]]:
    """解析 ``dimension×dimension`` 的精确留数矩阵。"""

    rows = record["residue"]
    if not isinstance(rows, list) or len(rows) != dimension:
        raise ValueError(f"{label}.residue must be a list of {dimension} rows")
    matrix = []
    for row_index, row in enumerate(rows):
        if not isinstance(row, list) or len(row) != dimension:
            raise ValueError(f"{label}.residue row {row_index} must have {dimension} entries")
        matrix.append([gaussian_rational(item) for item in row])
    return matrix


def _zero_matrix(dimension: int) -> list[list[GaussianRational]]:
    return [[GaussianRational() for _ in range(dimension)] for _ in range(dimension)]


def letter_partial_fraction_system(
    letters: Sequence[dict[str, Any]],
    dimension: int,
    *,
    name: str = "letter-partial-fraction-system",
    max_degree: int = 3,
) -> PartialFractionSystem | None:
    """按 ``d Log[N/D]`` 分解装配 letter 快速系统；前置条件不满足返回 None。

    ``N`` 的根带 ``+R`` 留数、``D`` 的根带 ``−R`` 留数，重数按因子幂次计入；
    不同 letter 的同名极点按球包含合并。留数在 Q(i) 上精确累加，只在装配完成
    时转成 acb 球。超过 ``max_degree`` 的 letter、根隔离失败都返回 None，由
    调用方回退一般路线。
    """

    if dimension <= 0:
        raise ValueError("dimension must be positive")
    if max_degree < 1:
        raise ValueError("max_degree must be at least one")
    poles: list[acb] = []
    buckets: list[list[list[GaussianRational]]] = []
    for index, record in enumerate(letters):
        label = f"letter {index}"
        require_exact_keys(record, {"numerator", "denominator", "residue"}, label)
        numerator, denominator = _letter_coefficients(record, label)
        if max(len(numerator), len(denominator)) - 1 > max_degree:
            return None
        residue = _letter_matrix(record, dimension, label)
        if all(value.is_zero for row in residue for value in row):
            continue
        for coefficients, sign in ((numerator, 1), (denominator, -1)):
            if len(coefficients) == 1:
                continue
            roots = isolated_roots(coefficients)
            if roots is None:
                return None
            for root, multiplicity in roots:
                position = next(
                    (existing for existing, pole in enumerate(poles)
                     if (root - pole).contains(0)),
                    None,
                )
                if position is None:
                    poles.append(root)
                    buckets.append(_zero_matrix(dimension))
                    position = len(poles) - 1
                target = buckets[position]
                for _ in range(multiplicity):
                    for row in range(dimension):
                        for column in range(dimension):
                            value = residue[row][column]
                            target[row][column] = (
                                target[row][column] + value if sign == 1
                                else target[row][column] - value
                            )
    residues = tuple(
        acb_mat([[value.to_acb() for value in row] for row in bucket])
        for bucket in buckets
    )
    return PartialFractionSystem(
        constant=acb_mat(dimension, dimension),
        residues=residues,
        poles=tuple(poles),
        name=name,
    )


def letter_system_matches_reference(
    fast_system: PartialFractionSystem,
    reference_system: Any,
    probes: Sequence[acb],
) -> bool:
    """在探针点比对 letter 快速系统与参照系统的矩阵球，全部包含才通过。

    两条求值路径的舍入差为工作精度量级（约 1 ulp），探针容差取
    ``10^-(prec/6)``：远高于舍入噪声，又比常规误差目标严数十个量级，只用来
    拦截装配性错误（错极点/漏 letter 的偏差是 O(1) 量级）。
    """

    tolerance = arb(10) ** (-(ctx.prec // 6))
    dimension = reference_system.dimension
    for probe in probes:
        fast = fast_system.evaluate(probe)
        reference = reference_system.evaluate(probe)
        for row in range(dimension):
            for column in range(dimension):
                difference = abs(fast[row, column] - reference[row, column])
                scale = abs(reference[row, column]) + arb(1)
                if not difference < tolerance * scale:
                    return False
    return True


def inverted_letter_system(
    letters: Sequence[dict[str, Any]],
    dimension: int,
    *,
    name: str = "letter-partial-fraction-system-at-infinity",
    max_degree: int = 3,
) -> PartialFractionSystem | None:
    """在 letter 层完成 ``x=1/t`` 拉回并装配 t 平面快速系统。

    恒等式：``-t^-2 dLog[N/D](1/t) dt = dLog[N*/D*] - (deg N - deg D) dLog t``，
    其中 ``N*(t)=t^deg N N(1/t)`` 即系数反转。逐 letter 只反转分子/分母系数，
    再把全部幂次差留数 ``-Σ (deg N - deg D) R`` 合成一条 ``dlog t`` letter 一并
    装配；不做任何整矩阵变换。前置条件不满足返回 None。
    """

    if dimension <= 0:
        raise ValueError("dimension must be positive")
    inverted: list[dict[str, Any]] = []
    shift = _zero_matrix(dimension)
    for index, record in enumerate(letters):
        label = f"letter {index}"
        require_exact_keys(record, {"numerator", "denominator", "residue"}, label)
        numerator, denominator = _letter_coefficients(record, label)
        power_difference = (len(numerator) - 1) - (len(denominator) - 1)
        residue = _letter_matrix(record, dimension, label)
        if not all(value.is_zero for row in residue for value in row):
            for row in range(dimension):
                for column in range(dimension):
                    shift[row][column] = (
                        shift[row][column] - residue[row][column] * power_difference
                    )
        inverted.append(
            {
                "numerator": list(reversed(numerator)),
                "denominator": list(reversed(denominator)),
                "residue": record["residue"],
            }
        )
    if any(not value.is_zero for row in shift for value in row):
        inverted.append(
            {
                "numerator": [GaussianRational(), GaussianRational(1)],
                "denominator": [GaussianRational(1)],
                "residue": shift,
            }
        )
    return letter_partial_fraction_system(
        inverted, dimension, name=name, max_degree=max_degree
    )
