"""用户点命中有限奇点时的局部解判定与可序列化输出。"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from flint import acb, acb_mat, arb

from .local_solutions import build_local_solution_basis
from .singularities import RationalMatrixSystem


@dataclass(frozen=True)
class SingularTargetResult:
    """保存奇点目标的分量值、分类和局部诊断。"""

    values: tuple[acb | str, ...]
    component_classifications: tuple[str, ...]
    classification: str
    report: dict[str, Any]


@dataclass(frozen=True)
class SingularTargetMatchPlan:
    """保存末端奇点的普通匹配点和局部收敛域证据。"""

    match_point: acb
    inserted: bool
    convergence_radius: arb | None
    report: dict[str, Any]


def plan_singular_target_match(
    system: RationalMatrixSystem,
    initial_point: Any,
    target_point: Any,
    *,
    radius_fraction: Any = "1/4",
) -> SingularTargetMatchPlan:
    """按目标奇点的局部收敛圆选择末端求值的普通匹配点。

    若当前普通点已落在保守半径内则直接复用；否则沿当前点到目标奇点的来向，
    在目标奇点收敛圆的 ``radius_fraction`` 处插入隐藏普通点。函数只规划点，
    不执行输运；调用方必须先把主链和参考链分别输运到同一个匹配点。
    """

    fraction = arb(radius_fraction)
    if not arb(0) < fraction < arb(1):
        raise ValueError("radius_fraction must lie strictly between zero and one")
    pole = target_point.to_acb() if hasattr(target_point, "to_acb") else acb(target_point)
    start = initial_point.to_acb() if hasattr(initial_point, "to_acb") else acb(initial_point)
    start_distance = abs(start - pole)
    if start_distance.contains(0):
        raise ValueError("initial_point must be ordinary and distinct from target_point")
    inventory = system.singularity_inventory()
    matching_records = [
        record
        for record in inventory.finite
        if record.location is not None and abs(record.location - pole).contains(0)
    ]
    if len(matching_records) != 1:
        raise ValueError("target_point must match exactly one finite system singularity")
    other_distances = [
        abs(record.location - pole)
        for record in inventory.finite
        if record.location is not None and not abs(record.location - pole).contains(0)
    ]
    convergence_radius = min(other_distances) if other_distances else None
    allowed_distance = (
        None if convergence_radius is None else convergence_radius * fraction
    )
    inserted = allowed_distance is not None and not start_distance < allowed_distance
    if inserted:
        direction = (start - pole) / start_distance
        match_point = acb(pole + direction * acb(allowed_distance))
    else:
        match_point = acb(start)
    match_distance = abs(match_point - pole)
    if convergence_radius is not None and not match_distance < convergence_radius:
        raise RuntimeError("planned singular-target match point lies outside its convergence disk")
    return SingularTargetMatchPlan(
        match_point,
        inserted,
        convergence_radius,
        {
            "target": pole.str(40),
            "sourcePoint": start.str(40),
            "matchPoint": match_point.str(40),
            "matchPointInserted": inserted,
            "convergenceRadius": (
                None if convergence_radius is None else convergence_radius.str(40)
            ),
            "matchDistance": match_distance.str(40),
            "radiusFraction": fraction.str(40),
            "planningBasis": "target singularity nearest-other-singularity distance",
        },
    )


def _component_classification(magnitudes: list[float], target_relative_error: float) -> str:
    """用局部幂对数解的缩小序列区分有限、pole 和 log 发散。"""

    if not magnitudes or max(magnitudes) == 0.0:
        return "removable_singularity"
    tail = magnitudes[-3:]
    finite_scale = max(1.0, max(tail))
    tail_spread = max(tail) - min(tail)
    if tail_spread / finite_scale <= max(target_relative_error, 1.0e-12) * 10:
        return "removable_singularity"
    first = max(magnitudes[0], 1.0e-300)
    last = max(magnitudes[-1], 1.0e-300)
    growth = last / first
    if growth > 1.0e4:
        # 幂次 pole 的增长随 1/|z|^p 明显快于零阶 log。
        return "true_pole"
    return "log_divergent_singularity"


def classify_singular_sample_vectors(
    vectors: list[acb_mat],
    *,
    target_point: Any,
    local_basis_method: str,
    local_basis_order: int,
    sample_exponents: tuple[int, ...],
    target_relative_error: Any,
    sample_scale: arb | None = None,
    convergence_radius: arb | None = None,
) -> SingularTargetResult:
    """把同一局部解沿趋近奇点的样本序列转换为用户可见分类。

    本函数不负责构造或匹配局部基；普通奇点目标和路径奇点 bucket 共用它，
    从而保证有限值、幂次 pole 与零阶 log 发散使用同一判据和输出 schema。
    """

    if not vectors:
        raise ValueError("singular target classification needs at least one sample vector")
    dimension = vectors[0].nrows()
    if any(vector.nrows() != dimension or vector.ncols() != 1 for vector in vectors):
        raise ValueError("singular target sample vectors must share one column-vector shape")
    tolerance = float(target_relative_error)
    component_classes: list[str] = []
    values: list[acb | str] = []
    magnitude_rows: list[list[float]] = []
    for row in range(dimension):
        magnitudes = [float(abs(vector[row, 0]).mid()) for vector in vectors]
        classification = _component_classification(magnitudes, tolerance)
        component_classes.append(classification)
        magnitude_rows.append(magnitudes)
        values.append(
            "Infinity"
            if classification != "removable_singularity"
            else vectors[-1][row, 0]
        )
    if "true_pole" in component_classes:
        overall = "true_pole"
    elif "log_divergent_singularity" in component_classes:
        overall = "log_divergent_singularity"
    else:
        overall = "removable_singularity"
    return SingularTargetResult(
        tuple(values),
        tuple(component_classes),
        overall,
        {
            "target": str(target_point),
            "classification": overall,
            "componentClassifications": component_classes,
            "sampleExponents": list(sample_exponents),
            "sampleScale": None if sample_scale is None else sample_scale.str(40),
            "convergenceRadius": (
                None if convergence_radius is None else convergence_radius.str(40)
            ),
            "sampleMagnitudes": magnitude_rows,
            "localBasisMethod": local_basis_method,
            "localBasisOrder": local_basis_order,
            "targetRelativeError": str(target_relative_error),
        },
    )


def evaluate_singular_target(
    system: RationalMatrixSystem,
    initial_point: Any,
    initial_vector: acb_mat,
    target_point: Any,
    *,
    order: int = 48,
    target_relative_error: Any = 1.0e-25,
    sample_exponents: tuple[int, ...] = (8, 12, 16, 20, 24),
) -> SingularTargetResult:
    """在 exact 有限奇点处沿局部幂对数解求值，不把 pole 当普通数值点。

    `initial_point` 必须是同一局部支上的普通点。有限分量返回最后一个局部样本的
    Acb 近似；真实 pole 或零阶 log 发散分量返回文本 ``Infinity``。该函数只依赖
    已认证 local basis，不调用数值积分或绕行路径。
    """

    if initial_vector.ncols() != 1 or initial_vector.nrows() != system.dimension:
        raise ValueError("initial_vector dimension does not match the system")
    pole = target_point.to_acb() if hasattr(target_point, "to_acb") else acb(target_point)
    start = initial_point.to_acb() if hasattr(initial_point, "to_acb") else acb(initial_point)
    if abs(start - pole).contains(0):
        raise ValueError("initial_point must be ordinary and distinct from target_point")
    inventory = system.singularity_inventory()
    other_distances = [
        abs(record.location - pole)
        for record in inventory.finite
        if record.location is not None and not abs(record.location - pole).contains(0)
    ]
    convergence_radius = min(other_distances) if other_distances else None
    start_distance = abs(start - pole)
    if convergence_radius is not None and not start_distance < convergence_radius:
        raise ValueError(
            "initial_point lies outside the singular local convergence disk; "
            "transport to an ordinary match point inside the disk first"
        )
    basis = build_local_solution_basis(system, target_point, order)
    constants = basis.evaluate(start - pole).solve(initial_vector)
    direction = (start - pole) / abs(start - pole)
    sample_scale = start_distance / arb(2)
    if convergence_radius is not None:
        sample_scale = min(sample_scale, convergence_radius / arb(4))
    samples = [
        direction * acb(sample_scale) * acb(10) ** (-exponent)
        for exponent in sample_exponents
    ]
    vectors = [basis.evaluate(local) * constants for local in samples]
    return classify_singular_sample_vectors(
        vectors,
        target_point=target_point,
        local_basis_method=basis.method,
        local_basis_order=order,
        sample_exponents=sample_exponents,
        target_relative_error=target_relative_error,
        sample_scale=sample_scale,
        convergence_radius=convergence_radius,
    )
