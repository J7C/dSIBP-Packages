"""一般矩阵微分方程的普通点 Taylor 输运与统一奇点局部桥接。

实现按 ``Y'=A(z)Y`` 的列向量约定递推解系数。路径可以由用户完整指定，也可以在
已声明奇点的前提下沿直线自动分段；自描述 adaptive path 遇到奇点时，先自动选择
Frobenius、shearing 后 Frobenius 或认证的指数广义幂级数，在入射侧反解局部常数，
再在出射侧求值。普通段的 Cauchy 圆和步长严格留在最近奇点以内。

快速路径包括：``PartialFractionSystem`` 极点--留数递推；``sample_points``
段内 dense output（保存点不再需要成为路径节点）；嵌入式单链截断认证
（``transport_path_refined(..., certification_mode="embedded")``）。
"""

from __future__ import annotations

import time
import warnings
from pathlib import Path
from typing import Any

from flint import acb, acb_mat, acb_poly, arb

from .boundary import ExponentialBoundary, FrobeniusBoundary
from .core import acb_midpoint_matrix, relative_difference_inf
from .local_solutions import build_local_solution_basis
from .savepoints import (
    SavePointWriter,
    singular_boundary_record_from_constants,
    vector_record,
)
from .singularities import RationalMatrixSystem
from .systems import AnalyticMatrixSystem, PartialFractionSystem


def vector_taylor_coefficients(
    matrix_coefficients: list[acb_mat], initial_vector: acb_mat
) -> list[acb_mat]:
    """递推 ``Y'=A Y`` 的一个或多个解列的 Taylor 系数。"""

    if not matrix_coefficients:
        raise ValueError("matrix coefficient list must not be empty")
    dimension = matrix_coefficients[0].nrows()
    if initial_vector.nrows() != dimension or initial_vector.ncols() < 1:
        raise ValueError("initial vector dimension mismatch")
    solution = [acb_mat(initial_vector)]
    for degree in range(len(matrix_coefficients)):
        total = acb_mat(dimension, initial_vector.ncols())
        for matrix_degree in range(degree + 1):
            total += matrix_coefficients[matrix_degree] * solution[degree - matrix_degree]
        solution.append(total / acb(degree + 1))
    return solution


def _matrix_column(matrix: acb_mat, column: int) -> acb_mat:
    """复制一个矩阵列，供逐物理解报告使用。"""

    result = acb_mat(matrix.nrows(), 1)
    for row in range(matrix.nrows()):
        result[row, 0] = matrix[row, column]
    return result


def evaluate_vector_series(coefficients: list[acb_mat], delta: acb) -> acb_mat:
    """用 Horner 法计算列向量 Taylor 多项式。"""

    if not coefficients:
        raise ValueError("solution coefficient list must not be empty")
    value = acb_mat(coefficients[-1])
    for coefficient in reversed(coefficients[:-1]):
        value = value * delta + coefficient
    return value


def evaluate_vector_series_many(
    coefficients: list[acb_mat],
    deltas: list[acb],
    *,
    fast_threshold: int = 8,
) -> tuple[list[acb_mat], str]:
    """批量计算同一 Taylor 向量多项式。

    同一输运节点覆盖的局部坐标按分量交给 FLINT ``acb_poly.evaluate``；大桶使用
    子积树/余数树快速多点求值，小桶使用 Horner，避免树构造固定开销。
    """

    if not coefficients:
        raise ValueError("solution coefficient list must not be empty")
    if not deltas:
        return [], "none"
    if fast_threshold < 2:
        raise ValueError("fast multipoint threshold must be at least 2")
    dimension = coefficients[0].nrows()
    column_count = coefficients[0].ncols()
    if any(
        coefficient.nrows() != dimension or coefficient.ncols() != column_count
        for coefficient in coefficients
    ):
        raise ValueError("solution coefficient matrices must have one common shape")
    algorithm = "fast" if len(deltas) >= fast_threshold else "iter"
    values = [acb_mat(dimension, column_count) for _ in deltas]
    for row in range(dimension):
        for column in range(column_count):
            polynomial = acb_poly(
                [coefficient[row, column] for coefficient in coefficients]
            )
            evaluated = polynomial.evaluate(deltas, algorithm=algorithm)
            for point_index, value in enumerate(evaluated):
                values[point_index][row, column] = value
    return [acb_midpoint_matrix(value) for value in values], algorithm


def _assign_dense_sample_points(
    sample_points: list[acb],
    path: list[acb],
    segment_radii: list[arb | None],
) -> list[tuple[int | None, acb | None]]:
    """把段内采样点路由到覆盖它的输运段，返回 ``(段号, 局部坐标)``。

    采样点必须严格位于某普通段的 Cauchy 圆盘内并落在线段范围内；恰为路径节点
    的点由节点 snapshot/保存机制覆盖，这里拒绝重复路由。采样点不改变路径，也
    不产生新的递推起点；绕行点则属于路径节点并建立新的递推起点。
    """

    assignments: list[tuple[int | None, acb | None]] = []
    for point in sample_points:
        chosen: tuple[int | None, acb | None] = (None, None)
        for index, (start, target) in enumerate(zip(path[:-1], path[1:])):
            local = point - start
            if abs(local).contains(0):
                raise ValueError(
                    f"sample point {point.str(40)} coincides with a path node; "
                    "use the node snapshot instead"
                )
            radius = segment_radii[index]
            if radius is None:
                continue
            step = target - start
            if not (abs(local) < radius and abs(local) <= abs(step)):
                continue
            ratio = local / step
            # 端点允许球半径级容差；与节点重合的点仍由 contains(0) 分支拒绝
            if not (
                abs(ratio.imag) <= arb("1e-12")
                and ratio.real >= arb("-1e-12")
                and ratio.real <= arb("1") + arb("1e-12")
            ):
                continue
            if chosen[0] is not None:
                raise ValueError(
                    f"sample point {point.str(40)} is covered by more than one segment"
                )
            chosen = (index, local)
        if chosen[0] is None:
            raise ValueError(
                f"sample point {point.str(40)} lies outside every segment Cauchy disk"
            )
        assignments.append(chosen)
    return assignments


def _sample_request(point: acb) -> dict[str, Any]:
    """为段内采样点构造保存记录请求。"""

    return {
        "coordinate": point.str(40),
        "working_coordinate": point.str(40),
        "classification": "ordinary",
        "singularity_identifier": None,
        "role": "sample",
    }


def build_straight_path(
    system: AnalyticMatrixSystem | PartialFractionSystem,
    start: acb,
    target: acb,
    *,
    step_fraction: float = 0.20,
) -> list[acb]:
    """按最近奇点距离的一定比例构造直线路径。"""

    if not 0 < step_fraction < 1:
        raise ValueError("step_fraction must lie between zero and one")
    points = [acb(start)]
    current = acb(start)
    direction = target - start
    total_distance = abs(direction)
    if total_distance.contains(0):
        return points
    if isinstance(system, PartialFractionSystem) and not system.poles:
        return [acb(start), acb(target)]
    unit = direction / total_distance
    for _ in range(100000):
        remaining = abs(target - current)
        if remaining.contains(0):
            return points
        radius = system.nearest_singularity_distance(current)
        proposed = radius * arb(str(step_fraction))
        step = remaining if float(remaining.mid()) <= float(proposed.mid()) else proposed
        current = acb(target) if step is remaining else current + unit * step
        points.append(acb(current))
        if abs(target - current).contains(0):
            return points
    raise RuntimeError("automatic path construction exceeded the segment limit")


def _resolve_transport_boundary(
    system: AnalyticMatrixSystem | RationalMatrixSystem,
    boundary: Any,
    path: list[acb],
    order: int,
    target_relative_error: arb | None,
) -> tuple[acb_mat, dict[str, Any] | None]:
    """把普通列向量或奇点 ``{a,b,C}`` 边界转换为首个执行点的列向量。"""

    start = getattr(path, "start_classification", None)
    if start is None or start.kind == "ordinary":
        if isinstance(boundary, (FrobeniusBoundary, ExponentialBoundary)):
            raise TypeError("singular boundary data require a singular path start")
        if not isinstance(boundary, acb_mat):
            raise TypeError("ordinary path start requires an acb_mat initial column vector")
        return acb_mat(boundary), None
    if not isinstance(system, RationalMatrixSystem):
        raise TypeError("singular {a,b,C} boundary requires its RationalMatrixSystem")
    working_system = getattr(path, "_working_system", None)
    if not isinstance(working_system, RationalMatrixSystem):
        raise TypeError("singular-start path lacks its exact working RationalMatrixSystem")
    inventory = working_system.singularity_inventory()
    record = inventory.find_finite(start.singularity_identifier)
    if record is None or record.location is None or record.location_exact is None:
        raise ValueError("singular boundary requires an exact Q(i) start point")
    basis = build_local_solution_basis(working_system, record.location_exact, order)
    constants, boundary_report = basis.resolve_boundary(boundary)
    local_point = path[0] - record.location
    vector = acb_midpoint_matrix(basis.evaluate(local_point) * constants.to_acb())
    method = (
        "regular_singular_boundary_initialization"
        if basis.method == "regular_singular_power_log"
        else f"{basis.method}_boundary_initialization"
    )
    report = _singular_initialization_report(
        basis,
        record,
        path,
        working_system,
        local_point,
        order,
        boundary_report,
        method,
        target_relative_error,
    )
    return vector, report


def _accuracy_target(value: Any | None) -> arb | None:
    """规范化可选相对精度目标，并拒绝无意义的阈值。"""

    if value is None:
        return None
    target = arb(str(value))
    if not 0 < target < 1:
        raise ValueError("target_relative_error must lie strictly between zero and one")
    return target


def _assess_formal_evaluation(
    evaluation: dict[str, Any],
    target_relative_error: arb | None,
) -> None:
    """评估实际选中分支的五阶诊断；失败只警告并保留计算结果。"""

    diagnostics = evaluation.get(
        "selected_branch_diagnostics", evaluation.get("branch_diagnostics", ())
    )
    issues: list[str] = []
    for index, diagnostic in enumerate(diagnostics, 1):
        sector_index = diagnostic.get("sector_index", index)
        ratio_text = diagnostic.get("next_five_over_previous_five")
        decreasing = diagnostic.get("next_five_over_previous_five_below_one")
        if decreasing is not True:
            ratio_label = "undefined" if ratio_text is None else ratio_text
            issues.append(
                f"sector {sector_index}: five-order block ratio={ratio_label} is not < 1"
            )
        refinement_text = diagnostic.get("five_order_relative_refinement")
        if target_relative_error is not None:
            meets_target = (
                refinement_text is not None
                and arb(refinement_text) < target_relative_error
            )
            diagnostic["target_relative_error"] = target_relative_error.str(30)
            diagnostic["five_order_relative_refinement_meets_target"] = bool(meets_target)
            if not meets_target:
                refinement_label = "undefined" if refinement_text is None else refinement_text
                issues.append(
                    f"sector {sector_index}: five-order relative refinement={refinement_label} "
                    f"does not meet {target_relative_error.str(12)}"
                )
    evaluation["formal_accuracy_checks_passed"] = not issues
    evaluation["formal_accuracy_issues"] = issues
    if issues:
        warnings.warn(
            "formal-asymptotic convergence warning; result retained: " + "; ".join(issues),
            UserWarning,
            stacklevel=3,
        )


def _transport_points_and_save_requests(path: list[Any]) -> tuple[list[acb], tuple[dict[str, Any], ...]]:
    """解析普通 ``list`` 中的 ``(coordinate,"save")``，并保留 AdaptivePath 元数据。"""

    adaptive_requests = getattr(path, "save_requests", None)
    if adaptive_requests is not None:
        return path, tuple(adaptive_requests)
    points: list[acb] = []
    requests: list[dict[str, Any]] = []
    for index, item in enumerate(path):
        tagged = isinstance(item, (list, tuple)) and len(item) == 2 and item[1] == "save"
        value = item[0] if tagged else item
        point = acb(value)
        points.append(point)
        if tagged:
            requests.append(
                {
                    "coordinate": point.str(40),
                    "working_coordinate": point.str(40),
                    "classification": "ordinary",
                    "singularity_identifier": None,
                    "role": "path",
                    "execution_index": index,
                }
            )
    return points, tuple(requests)


def _writer_for_requests(
    requests: tuple[dict[str, Any], ...],
    output_directory: str | Path | None,
    summary_filename: str,
    digits: int,
    enabled: bool,
) -> SavePointWriter | None:
    """只在路径确有 save 标签时建立输出根。"""

    if not enabled or not requests:
        return None
    return SavePointWriter(
        Path.cwd() if output_directory is None else Path(output_directory),
        summary_filename=summary_filename,
        digits=digits,
    )


def _write_ordinary_requests(
    writer: SavePointWriter | None,
    requests: tuple[dict[str, Any], ...],
    execution_index: int,
    vector: acb_mat,
) -> None:
    """在指定执行点完成后立即保存全部对应 ordinary 标签。"""

    if writer is None:
        return
    for request in requests:
        if request["classification"] == "ordinary" and request["execution_index"] == execution_index:
            writer.write(
                request,
                {
                    "resultType": "ordinary_vector",
                    "result": vector_record(vector, writer.digits),
                },
            )


def _write_singular_start_request(
    writer: SavePointWriter | None,
    requests: tuple[dict[str, Any], ...],
    path: list[acb],
    boundary: Any,
    order: int,
) -> None:
    """把已验证的奇点起点恢复为对应局部基的可复用边界。"""

    if writer is None:
        return
    start_requests = [
        request
        for request in requests
        if request["classification"] != "ordinary" and request["role"] == "start"
    ]
    if not start_requests:
        return
    working_system = getattr(path, "_working_system", None)
    if not isinstance(working_system, RationalMatrixSystem):
        raise TypeError("saved singular start requires its exact RationalMatrixSystem")
    inventory = working_system.singularity_inventory()
    for request in start_requests:
        if not isinstance(boundary, (FrobeniusBoundary, ExponentialBoundary)):
            raise TypeError("saved singular start requires a certified singular boundary")
        record = inventory.find_finite(request["singularity_identifier"])
        if record is None or record.location_exact is None:
            raise ValueError("saved singular start requires an exact Q(i) finite center")
        basis = build_local_solution_basis(working_system, record.location_exact, order)
        constants, _report = basis.resolve_boundary(boundary)
        if basis.method == "regular_singular_power_log" and isinstance(
            boundary, FrobeniusBoundary
        ):
            result_type = "frobenius_boundary"
            result = {
                "schema": "flintnde_frobenius_boundary_v1",
                "terms": boundary.to_json(),
                "localBasisMethod": basis.method,
            }
        else:
            result_type, result = singular_boundary_record_from_constants(
                basis, constants.to_acb(), writer.digits
            )
        writer.write(
            request,
            {
                "resultType": result_type,
                "result": result,
            },
        )


def _write_singular_bridge_requests(
    writer: SavePointWriter | None,
    requests: tuple[dict[str, Any], ...],
    singularity_identifier: str,
    basis: Any,
    constants: acb_mat,
) -> None:
    """在 continuation-ready bridge 的入射侧反解后即时保存中间奇点边界。"""

    if writer is None:
        return
    for request in requests:
        if (
            request["classification"] != "ordinary"
            and request["role"] == "detour"
            and request["singularity_identifier"] == singularity_identifier
        ):
            result_type, result = singular_boundary_record_from_constants(
                basis, constants, writer.digits
            )
            writer.write(request, {"resultType": result_type, "result": result})


def _write_singular_target_request(
    writer: SavePointWriter | None,
    requests: tuple[dict[str, Any], ...],
    path: list[acb],
    vector: acb_mat,
    order: int,
) -> None:
    """由终点入射匹配值反解局部常数，再输出可复用的 Frobenius 领头数据。"""

    if writer is None:
        return
    target_requests = [
        request
        for request in requests
        if request["classification"] != "ordinary" and request["role"] == "target"
    ]
    if not target_requests:
        return
    working_system = getattr(path, "_working_system", None)
    if not isinstance(working_system, RationalMatrixSystem):
        raise TypeError("saved singular target requires its exact RationalMatrixSystem")
    inventory = working_system.singularity_inventory()
    for request in target_requests:
        record = inventory.find_finite(request["singularity_identifier"])
        if record is None or record.location is None or record.location_exact is None:
            raise ValueError("saved singular target requires an exact Q(i) finite center")
        basis = build_local_solution_basis(working_system, record.location_exact, order)
        if not basis.continuation_ready:
            raise NotImplementedError(
                f"saved singular target requires continuation-ready local data, got {basis.method}"
            )
        constants = basis.evaluate(path[-1] - record.location).solve(vector)
        result_type, result = singular_boundary_record_from_constants(
            basis, constants, writer.digits
        )
        writer.write(
            request,
            {
                "resultType": result_type,
                "result": result,
            },
        )


def _singular_initialization_report(
    basis: Any,
    record: Any,
    path: list[acb],
    working_system: RationalMatrixSystem,
    local_point: acb,
    order: int,
    boundary_report: dict[str, Any],
    method: str,
    target_relative_error: arb | None,
) -> dict[str, Any]:
    """构造奇点起点的统一初始化报告，并附加局部求值诊断。

    单列 ``transport_path`` 与批量 ``_resolve_frobenius_boundary_batch`` 共用；
    两者只在 ``method`` 标签上不同。
    """

    report = {
        "segment_index": 0,
        "method": method,
        "singularity_identifier": record.identifier,
        "singularity": record.location_exact,
        "working_variable": getattr(path, "working_variable", working_system.variable_name),
        "ordinary_match_point": path[0].str(40),
        "local_coordinate": local_point.str(40),
        "order": order,
        "maximum_log_degree": basis.maximum_log_degree,
        "local_basis": basis.manifest,
        "branch_convention": "python-flint principal acb log/power branch",
        "boundary": boundary_report,
    }
    evaluation_report = basis.evaluation_report(local_point)
    if evaluation_report is not None:
        selected_sectors = sorted(
            {
                item["verified_sector"]
                for item in boundary_report.get("term_resolutions", ())
                if "verified_sector" in item
            }
        )
        if selected_sectors:
            diagnostics = evaluation_report.get("branch_diagnostics", ())
            evaluation_report["selected_sector_indices"] = selected_sectors
            evaluation_report["selected_branch_diagnostics"] = [
                {**diagnostics[index - 1], "sector_index": index}
                for index in selected_sectors
                if 0 < index <= len(diagnostics)
            ]
        _assess_formal_evaluation(evaluation_report, target_relative_error)
        report["local_evaluation"] = evaluation_report
    return report


def _transport_ordinary_segment(
    analytic_system: AnalyticMatrixSystem | PartialFractionSystem,
    vector_or_vectors: acb_mat,
    start: acb,
    target: acb,
    index: int,
    *,
    order: int,
    radius_fraction: float,
    sample_count: int | None,
    method_label: str,
    extra_report: dict[str, Any] | None = None,
) -> tuple[acb_mat, dict[str, Any], list[acb_mat]]:
    """执行一个普通点 Taylor 段并返回推进后的矢量、段报告与解系数。

    单列 ``transport_path`` 与多列 ``_transport_frobenius_boundary_batch`` 共用；
    保存点写入、snapshot 收集和逐列报告由调用方负责。返回的解系数供段内 dense
    output 求值与嵌入式截断认证的前缀求值复用；``PartialFractionSystem`` 走极点
    状态递推直接得到解系数，不做 Cauchy 采样。
    """

    unbounded_polynomial = (
        isinstance(analytic_system, PartialFractionSystem) and not analytic_system.poles
    )
    nearest = (
        None
        if unbounded_polynomial
        else analytic_system.nearest_singularity_distance(start)
    )
    radius = None if nearest is None else acb(nearest * arb(str(radius_fraction)))
    delta = target - start
    if radius is not None and float(abs(delta).mid()) >= float(abs(radius).mid()):
        raise ValueError(
            f"{analytic_system.name}: segment {index} leaves its Cauchy disk; refine the path"
        )
    segment_clock = time.perf_counter()
    if isinstance(analytic_system, PartialFractionSystem):
        solution_coefficients = analytic_system.solution_taylor_coefficients(
            start, order, vector_or_vectors
        )
        method = f"{method_label}_pole_recurrence"
        reported_sample_count = None
    else:
        matrix_coefficients = analytic_system.taylor_matrix_coefficients(
            start,
            order,
            radius=radius,
            sample_count=sample_count,
        )
        solution_coefficients = vector_taylor_coefficients(
            matrix_coefficients, vector_or_vectors
        )
        method = method_label
        reported_sample_count = sample_count or max(32, 2 * order)
    result = acb_midpoint_matrix(evaluate_vector_series(solution_coefficients, delta))
    elapsed = time.perf_counter() - segment_clock
    report: dict[str, Any] = {
        "segment_index": index,
        "method": method,
        "order": order,
        "sample_count": reported_sample_count,
        "step_over_convergence_radius": (
            0.0 if nearest is None else float((abs(delta) / nearest).mid())
        ),
        "cauchy_radius_over_nearest_singularity": (
            None if nearest is None else radius_fraction
        ),
        "entire_polynomial_connection": nearest is None,
        "elapsed_seconds": elapsed,
    }
    if extra_report:
        report.update(extra_report)
    return result, report, solution_coefficients


def transport_path(
    system: AnalyticMatrixSystem | PartialFractionSystem | RationalMatrixSystem,
    initial_vector: Any,
    path: list[acb],
    *,
    order: int,
    sample_count: int | None = None,
    radius_fraction: float = 0.60,
    target_relative_error: Any | None = None,
    save_output_directory: str | Path | None = None,
    save_summary_filename: str = "flintnde_save_points.json",
    save_digits: int = 40,
    write_save_points: bool = True,
    sample_points: list[Any] | None = None,
    collect_patches: bool = False,
) -> tuple[list[acb_mat], list[dict[str, Any]], float, list[dict[str, Any]] | None]:
    """沿普通路径或含统一奇点 bridge 的 adaptive path 输运列向量。

    普通 ``list[acb]`` 保持原接口。由 ``build_adaptive_path`` 返回的列表子类可和
    ``RationalMatrixSystem`` 直接配合；主函数在开始计算前读取全部 transition，确认
    哪些线段需要从奇点局部基反向求常数。

    ``sample_points`` 是段内 dense output 求值点：不改变路径，逐段复用局部解
    Taylor 系数 Horner 求值；启用保存时其结果以 ``role="sample"`` 记录写入。
    ``collect_patches=True`` 时额外返回逐段解系数补丁，供嵌入式截断认证复用。
    启用 ``sample_points`` 或 ``collect_patches`` 时返回四元组，第四个元素为
    ``{"sample_results": [...], "patches": [...] | None}``；否则保持既有三返回値。
    """

    path, save_requests = _transport_points_and_save_requests(path)
    if len(path) < 2:
        raise ValueError("transport path needs at least two points")
    if not 0 < radius_fraction < 1:
        raise ValueError("radius_fraction must lie between zero and one")
    accuracy_target = _accuracy_target(target_relative_error)
    writer = _writer_for_requests(
        save_requests,
        save_output_directory,
        save_summary_filename,
        save_digits,
        write_save_points,
    )
    transitions = getattr(path, "_transitions", None)
    specialization_report: dict[str, Any] | None = None
    if isinstance(system, RationalMatrixSystem):
        if transitions is None:
            raise TypeError(
                "RationalMatrixSystem transport requires a path returned by build_adaptive_path"
            )
        working_system = getattr(path, "_working_system")
        inventory = working_system.singularity_inventory()
        specialized, specialization_report = (
            working_system.specialize_polynomial_simple_poles(inventory)
        )
        analytic_system = (
            specialized
            if specialized is not None
            else working_system.to_analytic_system(inventory)
        )
    else:
        if transitions is not None and any(
            item.get("method") == "regular_singular_bridge" for item in transitions
        ):
            raise TypeError(
                "a path with singular local bridges must receive its RationalMatrixSystem"
            )
        working_system = None
        inventory = None
        analytic_system = system

    total_clock = time.perf_counter()
    dense_points = (
        None if sample_points is None else [acb(value) for value in sample_points]
    )
    dense_assignments: list[tuple[int | None, acb | None]] = []
    if dense_points:
        segment_radii: list[arb | None] = []
        for segment_index, segment_start in enumerate(path[:-1]):
            transition = (
                transitions[segment_index]
                if transitions is not None
                else {"method": "ordinary_taylor"}
            )
            if transition["method"] == "regular_singular_bridge":
                segment_radii.append(None)
            else:
                segment_radii.append(
                    (
                        abs(path[segment_index + 1] - segment_start) * arb(2)
                        if isinstance(analytic_system, PartialFractionSystem)
                        and not analytic_system.poles
                        else analytic_system.nearest_singularity_distance(segment_start)
                        * arb(str(radius_fraction))
                    )
                )
        dense_assignments = _assign_dense_sample_points(
            dense_points, path, segment_radii
        )
    patches: list[dict[str, Any]] | None = [] if collect_patches else None
    sample_results: list[dict[str, Any]] = []
    vector, boundary_report = _resolve_transport_boundary(
        system, initial_vector, path, order, accuracy_target
    )
    snapshots = [acb_mat(vector)]
    reports: list[dict[str, Any]] = (
        [] if boundary_report is None else [boundary_report]
    )
    _write_singular_start_request(writer, save_requests, path, initial_vector, order)
    _write_ordinary_requests(writer, save_requests, 0, vector)
    for index, (start, target) in enumerate(zip(path[:-1], path[1:]), start=1):
        transition = (
            transitions[index - 1]
            if transitions is not None
            else {"method": "ordinary_taylor"}
        )
        step_reports = [
            report
            for report in getattr(path, "step_reports", ())
            if report["execution_segment_index"] == index
        ]
        if transition["method"] == "regular_singular_bridge":
            if working_system is None or inventory is None:
                raise TypeError("singular local bridge requires an exact Q(i) working system")
            singularity = transition["singularity"]
            record = inventory.find_finite(singularity.singularity_identifier)
            if record is None or record.location is None:
                raise ValueError("local bridge singularity is absent from the working inventory")
            if record.location_exact is None:
                raise ValueError(
                    f"{record.identifier}: exact local bridge currently requires a Q(i) root"
                )
            segment_clock = time.perf_counter()
            basis = build_local_solution_basis(working_system, record.location_exact, order)
            if not basis.continuation_ready:
                raise NotImplementedError(
                    f"{basis.method} cannot bridge an internal singularity without Stokes data"
                )
            incoming_basis = basis.evaluate(start - record.location)
            constants = incoming_basis.solve(vector)
            _write_singular_bridge_requests(
                writer,
                save_requests,
                record.identifier,
                basis,
                constants,
            )
            vector = acb_midpoint_matrix(
                basis.evaluate(target - record.location) * constants
            )
            elapsed = time.perf_counter() - segment_clock
            snapshots.append(acb_mat(vector))
            _write_ordinary_requests(writer, save_requests, index, vector)
            ratios = [
                report["step_over_convergence_radius"] for report in step_reports
            ]
            bridge_report = {
                    "segment_index": index,
                    "method": f"{basis.method}_bridge",
                    "singularity_identifier": record.identifier,
                    "singularity": record.location_exact,
                    "order": order,
                    "maximum_log_degree": basis.maximum_log_degree,
                    "local_basis": basis.manifest,
                    "branch_convention": "python-flint principal acb log/power branch",
                    "step_over_convergence_radius": ratios,
                    "maximum_step_over_convergence_radius": max(
                        (value for value in ratios if value is not None), default=None
                    ),
                    "elapsed_seconds": elapsed,
                }
            incoming_evaluation = basis.evaluation_report(start - record.location)
            outgoing_evaluation = basis.evaluation_report(target - record.location)
            if incoming_evaluation is not None or outgoing_evaluation is not None:
                if incoming_evaluation is not None:
                    _assess_formal_evaluation(incoming_evaluation, accuracy_target)
                if outgoing_evaluation is not None:
                    _assess_formal_evaluation(outgoing_evaluation, accuracy_target)
                bridge_report["local_evaluation"] = {
                    "incoming": incoming_evaluation,
                    "outgoing": outgoing_evaluation,
                }
            reports.append(bridge_report)
            if patches is not None:
                patches.append({"segment_index": index, "method": "singular_bridge"})
            continue

        vector, segment_report, solution_coefficients = _transport_ordinary_segment(
            analytic_system,
            vector,
            start,
            target,
            index,
            order=order,
            radius_fraction=radius_fraction,
            sample_count=sample_count,
            method_label="ordinary_taylor",
            extra_report=(
                {"system_specialization": specialization_report}
                if specialization_report is not None
                else None
            ),
        )
        snapshots.append(acb_mat(vector))
        _write_ordinary_requests(writer, save_requests, index, vector)
        reports.append(segment_report)
        if patches is not None:
            patches.append(
                {
                    "segment_index": index,
                    "method": segment_report["method"],
                    "center": start,
                    "delta": target - start,
                    "order": order,
                    "solution_coefficients": solution_coefficients,
                }
            )
        segment_dense = [
            (dense_index, local)
            for dense_index, (segment_index, local) in enumerate(dense_assignments)
            if segment_index == index - 1 and local is not None
        ]
        dense_values, dense_algorithm = evaluate_vector_series_many(
            solution_coefficients,
            [local for _dense_index, local in segment_dense],
        )
        segment_report["dense_evaluation"] = {
            "algorithm": dense_algorithm,
            "point_count": len(segment_dense),
        }
        for (dense_index, _local), dense_value in zip(segment_dense, dense_values):
            if writer is not None:
                writer.write(
                    _sample_request(dense_points[dense_index]),
                    {
                        "resultType": "ordinary_vector",
                        "result": vector_record(dense_value, writer.digits),
                    },
                )
            sample_results.append(
                {
                    "coordinate": dense_points[dense_index].str(40),
                    "segment_index": index,
                    "evaluation_algorithm": dense_algorithm,
                    "value": acb_mat(dense_value),
                }
            )
    _write_singular_target_request(writer, save_requests, path, vector, order)
    if writer is not None:
        writer.finalize()
    if dense_points or patches is not None:
        return (
            snapshots,
            reports,
            time.perf_counter() - total_clock,
            {"sample_results": sample_results, "patches": patches},
        )
    return snapshots, reports, time.perf_counter() - total_clock


def _resolve_frobenius_boundary_batch(
    system: RationalMatrixSystem,
    boundaries: list[FrobeniusBoundary],
    path: list[acb],
    order: int,
    target_relative_error: arb | None,
) -> tuple[acb_mat, list[dict[str, Any]]]:
    """一次构造局部基，并把多个 Frobenius 边界变成同一匹配点的解列。"""

    start = getattr(path, "start_classification", None)
    if start is None or start.kind != "regular_singular":
        raise TypeError("Frobenius batch requires a regular-singular path start")
    working_system = getattr(path, "_working_system", None)
    if not isinstance(working_system, RationalMatrixSystem):
        raise TypeError("Frobenius batch path lacks its exact working RationalMatrixSystem")
    inventory = working_system.singularity_inventory()
    record = inventory.find_finite(start.singularity_identifier)
    if record is None or record.location is None or record.location_exact is None:
        raise ValueError("Frobenius batch requires an exact Q(i) start point")
    basis = build_local_solution_basis(working_system, record.location_exact, order)
    local_point = path[0] - record.location
    evaluated_basis = basis.evaluate(local_point)
    vectors = acb_mat(system.dimension, len(boundaries))
    reports: list[dict[str, Any]] = []
    for column, boundary in enumerate(boundaries):
        constants, boundary_report = basis.resolve_boundary(boundary)
        vector = acb_midpoint_matrix(evaluated_basis * constants.to_acb())
        for row in range(system.dimension):
            vectors[row, column] = vector[row, 0]
        report = _singular_initialization_report(
            basis,
            record,
            path,
            working_system,
            local_point,
            order,
            boundary_report,
            "regular_singular_boundary_initialization_batch",
            target_relative_error,
        )
        reports.append(report)
    return vectors, reports


def _transport_frobenius_boundary_batch(
    system: RationalMatrixSystem,
    boundaries: list[FrobeniusBoundary],
    path: list[acb],
    *,
    order: int,
    radius_fraction: float,
    target_relative_error: arb | None,
    collect_patches: bool = False,
) -> tuple[
    list[acb_mat], list[dict[str, Any]], list[dict[str, Any]], float,
    list[dict[str, Any]] | None,
]:
    """共享局部基和普通段矩阵系数，批量输运多个 Frobenius 解列。

    ``collect_patches=True`` 时额外返回逐段解系数补丁，供嵌入式截断认证的
    前缀求值复用。
    """

    path, save_requests = _transport_points_and_save_requests(path)
    if save_requests:
        raise NotImplementedError("Frobenius batch does not yet support save points")
    if len(path) < 2:
        raise ValueError("transport path needs at least two points")
    if not 0 < radius_fraction < 1:
        raise ValueError("radius_fraction must lie between zero and one")
    transitions = getattr(path, "_transitions", None)
    if transitions is None or any(
        item.get("method") != "ordinary_taylor" for item in transitions
    ):
        raise NotImplementedError(
            "Frobenius batch currently requires an ordinary path after the singular launch"
        )
    working_system = getattr(path, "_working_system", None)
    if not isinstance(working_system, RationalMatrixSystem):
        raise TypeError("Frobenius batch path lacks its exact working RationalMatrixSystem")
    inventory = working_system.singularity_inventory()
    analytic_system = working_system.to_analytic_system(inventory)
    total_clock = time.perf_counter()
    vectors, boundary_reports = _resolve_frobenius_boundary_batch(
        system, boundaries, path, order, target_relative_error
    )
    snapshots = [acb_mat(vectors)]
    reports: list[dict[str, Any]] = []
    patches: list[dict[str, Any]] | None = [] if collect_patches else None
    for index, (start, target) in enumerate(zip(path[:-1], path[1:]), start=1):
        vectors, segment_report, solution_coefficients = _transport_ordinary_segment(
            analytic_system,
            vectors,
            start,
            target,
            index,
            order=order,
            radius_fraction=radius_fraction,
            sample_count=None,
            method_label="ordinary_taylor_batch",
            extra_report={"column_count": len(boundaries)},
        )
        snapshots.append(acb_mat(vectors))
        reports.append(segment_report)
        if patches is not None:
            patches.append(
                {
                    "segment_index": index,
                    "method": segment_report["method"],
                    "center": start,
                    "delta": target - start,
                    "order": order,
                    "solution_coefficients": solution_coefficients,
                }
            )
    return (
        snapshots,
        reports,
        boundary_reports,
        time.perf_counter() - total_clock,
        patches,
    )


def _embedded_prefix_snapshots(
    snapshots: list[acb_mat],
    patches: list[dict[str, Any]],
    primary_order: int,
) -> tuple[list[acb_mat], list[arb]]:
    """从单链高阶运行的解系数补丁提取主阶前缀 snapshot 与逐段截断差。

    递推的前缀性质保证前 ``primary_order+1`` 个解系数与独立低阶运行逐位相同，
    所以主链结果只需对前缀多项式重新 Horner 求值；奇点 bridge 段没有 Taylor
    补丁，其主链値直接复用参考链値。
    """

    primary_snapshots = [acb_mat(snapshots[0])]
    truncations: list[arb] = []
    for node_index in range(1, len(snapshots)):
        patch = patches[node_index - 1]
        if patch.get("method") == "singular_bridge":
            raise ValueError(
                "embedded prefix extraction cannot certify a singular bridge"
            )
        coefficients = patch["solution_coefficients"][: primary_order + 1]
        prefix_value = acb_midpoint_matrix(
            evaluate_vector_series(coefficients, patch["delta"])
        )
        primary_snapshots.append(acb_mat(prefix_value))
        truncations.append(
            relative_difference_inf(prefix_value, snapshots[node_index])
        )
    return primary_snapshots, truncations


def _embedded_prefix_sample_results(
    reference_results: list[dict[str, Any]],
    patches: list[dict[str, Any]],
    primary_order: int,
    local_points: list[acb],
) -> list[dict[str, Any]]:
    """从单段参考链系数前缀重算用户可见的主阶 dense 值。"""

    if len(reference_results) != len(local_points):
        raise ValueError("embedded dense result count does not match sample points")
    if not reference_results:
        return []
    if len(patches) != 1 or patches[0].get("method") == "singular_bridge":
        raise ValueError("embedded dense extraction requires one ordinary Taylor patch")
    coefficients = patches[0]["solution_coefficients"][: primary_order + 1]
    values, algorithm = evaluate_vector_series_many(coefficients, local_points)
    return [
        {
            **record,
            "evaluation_algorithm": algorithm,
            "value": acb_midpoint_matrix(value),
        }
        for record, value in zip(reference_results, values)
    ]


def transport_frobenius_boundaries_refined(
    system: RationalMatrixSystem,
    boundaries: list[FrobeniusBoundary],
    path: list[acb],
    *,
    primary_order: int,
    reference_order: int,
    radius_fraction: float = 0.60,
    target_relative_error: Any | None = None,
    certification_mode: str = "certified",
) -> dict[str, Any]:
    """对同一系统的多个 Frobenius 初值共享构造，并逐列报告 refinement。

    ``certification_mode="embedded"`` 只跑 reference_order 一条链，主链
    结果取解系数前缀；``"certified"``（缺省）保留独立双链完整重算用于论文级认证。
    """

    if not boundaries:
        raise ValueError("Frobenius boundary batch must not be empty")
    if reference_order <= primary_order:
        raise ValueError("reference_order must exceed primary_order")
    if certification_mode not in ("embedded", "certified"):
        raise ValueError('certification_mode must be "embedded" or "certified"')
    requested_certification_mode = certification_mode
    singular_boundary_forced_certified = certification_mode == "embedded"
    if singular_boundary_forced_certified:
        warnings.warn(
            "embedded certification cannot reuse a Frobenius boundary initialization; "
            "running independent primary/reference chains instead",
            UserWarning,
            stacklevel=2,
        )
        certification_mode = "certified"
    accuracy_target = _accuracy_target(target_relative_error)
    if certification_mode == "embedded":
        chain = _transport_frobenius_boundary_batch(
            system,
            boundaries,
            path,
            order=reference_order,
            radius_fraction=radius_fraction,
            target_relative_error=accuracy_target,
            collect_patches=True,
        )
        primary_snapshots, truncations = _embedded_prefix_snapshots(
            chain[0], chain[4], primary_order
        )
        primary_segments = [
            {**report, "method": f"{report['method']}_embedded_prefix"}
            if "method" in report
            else report
            for report in chain[1]
        ]
        differences = [
            relative_difference_inf(
                _matrix_column(primary_snapshots[-1], column),
                _matrix_column(chain[0][-1], column),
            )
            for column in range(len(boundaries))
        ]
        result = {
            "certification_mode": "embedded",
            "primary_snapshots": primary_snapshots,
            "reference_snapshots": chain[0],
            "primary_segments": primary_segments,
            "reference_segments": chain[1],
            "primary_boundary_reports": chain[2],
            "reference_boundary_reports": chain[2],
            "primary_seconds": chain[3],
            "reference_seconds": chain[3],
            "segment_truncation_differences_inf": truncations,
            "segment_truncation_differences_midpoint": [
                float(value.mid()) for value in truncations
            ],
        }
    else:
        primary = _transport_frobenius_boundary_batch(
            system,
            boundaries,
            path,
            order=primary_order,
            radius_fraction=radius_fraction,
            target_relative_error=accuracy_target,
        )
        reference = _transport_frobenius_boundary_batch(
            system,
            boundaries,
            path,
            order=reference_order,
            radius_fraction=radius_fraction,
            target_relative_error=accuracy_target,
        )
        differences = [
            relative_difference_inf(
                _matrix_column(primary[0][-1], column),
                _matrix_column(reference[0][-1], column),
            )
            for column in range(len(boundaries))
        ]
        result = {
            "certification_mode": "certified",
            "primary_snapshots": primary[0],
            "reference_snapshots": reference[0],
            "primary_segments": primary[1],
            "reference_segments": reference[1],
            "primary_boundary_reports": primary[2],
            "reference_boundary_reports": reference[2],
            "primary_seconds": primary[3],
            "reference_seconds": reference[3],
        }
    meets_target = [
        None if accuracy_target is None else bool(difference < accuracy_target)
        for difference in differences
    ]
    if any(value is False for value in meets_target):
        warnings.warn(
            "NDE batch refinement accuracy warning; results retained",
            UserWarning,
            stacklevel=2,
        )
    result.update(
        {
            "relative_differences_inf": differences,
            "certification_mode_requested": requested_certification_mode,
            "singular_boundary_forced_certified": singular_boundary_forced_certified,
            "relative_differences_midpoint": [float(value.mid()) for value in differences],
            "target_relative_error": None if accuracy_target is None else accuracy_target.str(20),
            "target_relative_error_met": meets_target,
        }
    )
    return result


def transport_path_refined(
    system: AnalyticMatrixSystem | PartialFractionSystem | RationalMatrixSystem,
    initial_vector: Any,
    path: list[acb],
    *,
    primary_order: int,
    reference_order: int,
    primary_sample_count: int | None = None,
    reference_sample_count: int | None = None,
    radius_fraction: float = 0.60,
    target_relative_error: Any | None = None,
    save_output_directory: str | Path | None = None,
    save_summary_filename: str = "flintnde_save_points.json",
    save_digits: int = 40,
    certification_mode: str = "certified",
    sample_points: list[Any] | None = None,
) -> dict[str, Any]:
    """输运并认证末点精度。

    ``certification_mode="certified"``（缺省）运行独立双链。用户结果、dense sample
    和保存点均来自 ``primary_order`` 主链；``reference_order`` 链只用于精度核验。
    ``"embedded"`` 仅对单段普通 Taylor 路径复用 reference_order 解系数前缀；
    多段路径无法据此传播低阶初值误差，含奇点局部段也不能复用，因此两者均升级为
    独立 ``"certified"`` 双链。含保存点的请求也升级为双链，使保存文件严格来自主链。
    ``sample_points`` 为段内 dense output 求值点，不改变路径。
    """

    if reference_order <= primary_order:
        raise ValueError("reference_order must exceed primary_order")
    if certification_mode not in ("embedded", "certified"):
        raise ValueError('certification_mode must be "embedded" or "certified"')
    requested_certification_mode = certification_mode
    transport_nodes, save_requests = _transport_points_and_save_requests(path)
    transitions = getattr(path, "_transitions", ())
    start_classification = getattr(path, "start_classification", None)
    target_classification = getattr(path, "target_classification", None)
    has_singular_local_step = (
        any(item.get("method") == "regular_singular_bridge" for item in transitions)
        or (
            start_classification is not None
            and start_classification.kind != "ordinary"
        )
        or (
            target_classification is not None
            and target_classification.kind != "ordinary"
        )
    )
    singular_local_forced_certified = (
        certification_mode == "embedded" and has_singular_local_step
    )
    multi_segment_forced_certified = (
        certification_mode == "embedded" and len(transport_nodes) > 2
    )
    save_points_forced_certified = (
        certification_mode == "embedded" and bool(save_requests)
    )
    if (
        singular_local_forced_certified
        or multi_segment_forced_certified
        or save_points_forced_certified
    ):
        reasons = []
        if singular_local_forced_certified:
            reasons.append("singular local data cannot be reused")
        if multi_segment_forced_certified:
            reasons.append("low-order state errors cannot be propagated by segment prefixes")
        if save_points_forced_certified:
            reasons.append("saved values must be produced by the primary chain")
        warnings.warn(
            "embedded certification upgraded to independent primary/reference chains: "
            + "; ".join(reasons),
            UserWarning,
            stacklevel=2,
        )
        certification_mode = "certified"
    accuracy_target = _accuracy_target(target_relative_error)
    if certification_mode == "embedded":
        snapshots, reports, elapsed, extra = transport_path(
            system,
            initial_vector,
            path,
            order=reference_order,
            sample_count=reference_sample_count,
            radius_fraction=radius_fraction,
            target_relative_error=accuracy_target,
            save_output_directory=save_output_directory,
            save_summary_filename=save_summary_filename,
            save_digits=save_digits,
            sample_points=sample_points,
            collect_patches=True,
        )
        primary_snapshots, truncations = _embedded_prefix_snapshots(
            snapshots, extra["patches"], primary_order
        )
        primary_segments = [
            {**report, "method": f"{report['method']}_embedded_prefix"}
            if "method" in report
            else report
            for report in reports
        ]
        difference = relative_difference_inf(primary_snapshots[-1], snapshots[-1])
        result = {
            "certification_mode": "embedded",
            "certification_mode_requested": requested_certification_mode,
            "singular_local_forced_certified": singular_local_forced_certified,
            "multi_segment_forced_certified": multi_segment_forced_certified,
            "save_points_forced_certified": save_points_forced_certified,
            "primary_snapshots": primary_snapshots,
            "reference_snapshots": snapshots,
            "primary_segments": primary_segments,
            "reference_segments": reports,
            "primary_seconds": elapsed,
            "reference_seconds": elapsed,
            "segment_truncation_differences_inf": truncations,
            "segment_truncation_differences_midpoint": [
                float(value.mid()) for value in truncations
            ],
        }
        if sample_points is not None:
            result["sample_results"] = _embedded_prefix_sample_results(
                extra["sample_results"],
                extra["patches"],
                primary_order,
                [acb(point) - transport_nodes[0] for point in sample_points],
            )
        meets_target = None if accuracy_target is None else bool(difference < accuracy_target)
        if meets_target is False:
            warnings.warn(
                "NDE embedded truncation accuracy warning; result retained: "
                f"truncation difference {difference.str(20)} does not meet "
                f"{accuracy_target.str(20)}",
                UserWarning,
                stacklevel=2,
            )
        result.update(
            {
                "relative_difference_inf": difference,
                "relative_difference_midpoint": float(difference.mid()),
                "target_relative_error": (
                    None if accuracy_target is None else accuracy_target.str(30)
                ),
                "target_relative_error_met": meets_target,
            }
        )
        return result

    primary = transport_path(
        system,
        initial_vector,
        path,
        order=primary_order,
        sample_count=primary_sample_count,
        radius_fraction=radius_fraction,
        target_relative_error=accuracy_target,
        save_output_directory=save_output_directory,
        save_summary_filename=save_summary_filename,
        save_digits=save_digits,
        sample_points=sample_points,
    )
    reference = transport_path(
        system,
        initial_vector,
        path,
        order=reference_order,
        sample_count=reference_sample_count,
        radius_fraction=radius_fraction,
        target_relative_error=accuracy_target,
        write_save_points=False,
    )
    difference = relative_difference_inf(primary[0][-1], reference[0][-1])
    meets_target = None if accuracy_target is None else bool(difference < accuracy_target)
    if meets_target is False:
        warnings.warn(
            "NDE refinement accuracy warning; result retained: relative difference "
            f"{difference.str(20)} does not meet {accuracy_target.str(20)}",
            UserWarning,
            stacklevel=2,
        )
    result = {
        "certification_mode": "certified",
        "certification_mode_requested": requested_certification_mode,
        "singular_local_forced_certified": singular_local_forced_certified,
        "multi_segment_forced_certified": multi_segment_forced_certified,
        "save_points_forced_certified": save_points_forced_certified,
        "primary_snapshots": primary[0],
        "reference_snapshots": reference[0],
        "primary_segments": primary[1],
        "reference_segments": reference[1],
        "primary_seconds": primary[2],
        "reference_seconds": reference[2],
        "relative_difference_inf": difference,
        "relative_difference_midpoint": float(difference.mid()),
        "target_relative_error": (
            None if accuracy_target is None else accuracy_target.str(30)
        ),
        "target_relative_error_met": meets_target,
    }
    if sample_points is not None:
        result["sample_results"] = primary[3]["sample_results"]
    return result
