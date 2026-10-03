"""FlintNDE 的 Wolfram/Mathematica 嵌入桥接。

命令行用法：``python -m flintnde.mathematica_bridge request.json output.m``。

桥接只负责"最终输出"：读取 Wolfram 端写出的 JSON 请求，运行一次
``transport_path_refined``，把末点结果、可选段内采样点与认证摘要写成 Wolfram
``Get`` 可直接加载的 ``.m`` 文件。中间泰勒系数补丁不写盘、不回传；是否把结果
落盘为 MMA 文件由用户在 Wolfram 端用自带输出命令决定。
"""

from __future__ import annotations

import argparse
import json
from multiprocessing import get_context
import re
import sys
import warnings
from pathlib import Path
from typing import Any

from flint import acb, acb_mat, arb

from .core import configure_working_precision, exact_rational, require_exact_keys
from .singularity_jump import (
    SingularPathError,
    plan_transport_path,
    planned_path_from_json,
    transport_planned_path_refined,
)
from .routing import (
    AdaptivePathSingularityError,
    NamedPoint,
    adaptive_path_from_json,
    adaptive_path_to_json,
    build_adaptive_path,
)
from .singularities import RationalMatrixSystem, rational_function
from .systems import PartialFractionSystem
from .transport import transport_path_refined

REQUEST_SCHEMA = "flintnde_mathematica_request_v1"
RESULT_SCHEMA = "flintnde_mathematica_bridge_v1"
DEFAULT_PARALLEL_TASK_COUNT = 12


def _mma_scalar(text: str) -> Any:
    """把 MMA 风格数值字符串规整为 exact 有理转换可接受的十进制/分数字符串。

    支持 Wolfram ``*^`` 指数记号和反引号精度后缀；``1/3`` 等分数按 exact 有理
    处理，不经过 machine float。
    """

    cleaned = re.sub(r"`(?:\d+(?:\.\d*)?|\.\d+)?", "", str(text).strip())
    if "*^" in cleaned:
        mantissa, exponent = cleaned.split("*^", 1)
        cleaned = f"{mantissa}e{exponent}"
    return cleaned


def _entry_to_acb(value: Any) -> acb:
    """把当前 bridge 标量记录转换为 Acb，不接受缺分量的复数对象。"""

    if isinstance(value, dict):
        require_exact_keys(value, {"re", "im"}, "complex scalar")
        return acb(
            acb(exact_rational(_mma_scalar(value["re"]))).real,
            acb(exact_rational(_mma_scalar(value["im"]))).real,
        )
    if isinstance(value, (list, tuple)):
        if len(value) != 2:
            raise ValueError(f"complex entry must have two components, got {value}")
        return acb(
            acb(exact_rational(_mma_scalar(value[0]))).real,
            acb(exact_rational(_mma_scalar(value[1]))).real,
        )
    return acb(exact_rational(_mma_scalar(value)))

def _column_vector(entries: list[Any], name: str) -> acb_mat:
    if not isinstance(entries, list) or not entries:
        raise ValueError(f"{name} must be a nonempty list")
    return acb_mat([[_entry_to_acb(entry)] for entry in entries])


def _matrix_from_records(records: list[list[Any]], name: str) -> acb_mat:
    if not isinstance(records, list) or not records:
        raise ValueError(f"{name} must be a nonempty matrix record list")
    return acb_mat([[_entry_to_acb(entry) for entry in row] for row in records])


def _exact_coefficient(value: Any, field_name: str) -> Any:
    """把当前 Wolfram exact 标量记录改成高斯有理数字段。"""

    if isinstance(value, dict):
        require_exact_keys(value, {"re", "im"}, field_name)
        return {
            "real": _mma_scalar(value["re"]),
            "imag": _mma_scalar(value["im"]),
        }
    return _mma_scalar(value)


def _build_partial_fraction_system(record: dict[str, Any]) -> PartialFractionSystem:
    """构造唯一的显式多项式加简单极点系统记录。"""

    require_exact_keys(
        record,
        {"type", "polynomialCoefficients", "residues", "poles"},
        "partialFraction system",
    )
    if record["type"] != "partialFraction":
        raise ValueError("partialFraction system has an inconsistent type")
    polynomial_records = record["polynomialCoefficients"]
    residue_records = record["residues"]
    pole_records = record["poles"]
    if not isinstance(polynomial_records, list) or not polynomial_records:
        raise ValueError("partialFraction system needs nonempty polynomialCoefficients")
    if not isinstance(residue_records, list) or not isinstance(pole_records, list):
        raise ValueError("partialFraction residues and poles must be lists")
    if len(residue_records) != len(pole_records):
        raise ValueError("partialFraction residue and pole counts must agree")
    polynomial = tuple(
        _matrix_from_records(item, f"polynomial coefficient {degree}")
        for degree, item in enumerate(polynomial_records)
    )
    residues = tuple(
        _matrix_from_records(item, f"residue {index}")
        for index, item in enumerate(residue_records)
    )
    poles = tuple(_entry_to_acb(item) for item in pole_records)
    return PartialFractionSystem(
        constant=polynomial[0],
        polynomial_coefficients=polynomial[1:],
        residues=residues,
        poles=poles,
        name="Mathematica-partial-fraction-system",
    )


def _build_rational_matrix_system(record: dict[str, Any]) -> RationalMatrixSystem:
    """构造 exact Q(i)(x) 有理矩阵；奇点由 FlintNDE 内部发现。"""

    require_exact_keys(
        record,
        {"type", "variable", "name", "matrix"},
        "rationalMatrix system",
    )
    if record["type"] != "rationalMatrix":
        raise ValueError("rationalMatrix system has an inconsistent type")
    if not isinstance(record["variable"], str) or not isinstance(record["name"], str):
        raise ValueError("rationalMatrix variable and name must be strings")
    matrix_records = record["matrix"]
    if not isinstance(matrix_records, list) or not matrix_records:
        raise ValueError("rationalMatrix system needs a nonempty matrix")
    entries = []
    for row_index, row in enumerate(matrix_records):
        if not isinstance(row, list) or not row:
            raise ValueError(f"rationalMatrix row {row_index} must be nonempty")
        output_row = []
        for column_index, entry in enumerate(row):
            field_name = f"matrix[{row_index}][{column_index}]"
            require_exact_keys(entry, {"numerator", "denominator"}, field_name)
            numerator = entry["numerator"]
            denominator = entry["denominator"]
            if not isinstance(numerator, list) or not numerator:
                raise ValueError(f"{field_name}.numerator must be a nonempty list")
            if not isinstance(denominator, list) or not denominator:
                raise ValueError(f"{field_name}.denominator must be a nonempty list")
            output_row.append(
                rational_function(
                    tuple(
                        _exact_coefficient(item, f"{field_name}.numerator")
                        for item in numerator
                    ),
                    tuple(
                        _exact_coefficient(item, f"{field_name}.denominator")
                        for item in denominator
                    ),
                )
            )
        entries.append(tuple(output_row))
    return RationalMatrixSystem(
        tuple(entries),
        variable_name=record["variable"],
        name=record["name"],
    )


def _build_system(
    record: dict[str, Any],
) -> PartialFractionSystem | RationalMatrixSystem:
    """按唯一公开 system type 构造单变量矩阵微分方程。"""

    if not isinstance(record, dict) or "type" not in record:
        raise ValueError("bridge system must be an object containing type")
    system_type = record["type"]
    if system_type == "partialFraction":
        return _build_partial_fraction_system(record)
    if system_type == "rationalMatrix":
        return _build_rational_matrix_system(record)
    raise ValueError(f"unsupported Mathematica bridge system type: {system_type}")

def _mma_float_literal(value: float) -> str:
    """把 machine float 写成 Wolfram 可解析的 ``*^`` 科学记号。"""

    mantissa, exponent = f"{value:.16e}".split("e", 1)
    return f"{mantissa}*^{int(exponent)}"


def _mma_string(text: str) -> str:
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def _acb_record(value: acb, digits: int) -> dict[str, str]:
    """把 Acb 球的中点与半径保存为十进制字符串。"""

    return {
        "real": value.real.mid().str(digits, radius=False, more=True),
        "imag": value.imag.mid().str(digits, radius=False, more=True),
        "realRadius": value.real.rad().str(digits),
        "imagRadius": value.imag.rad().str(digits),
    }


def _vector_records(vector: acb_mat, digits: int) -> list[dict[str, str]]:
    """把列向量写成 Wolfram 端可恢复的高精度复数记录。"""

    if vector.ncols() != 1:
        raise ValueError("bridge output vector must be a column")
    return [_acb_record(vector[row, 0], digits) for row in range(vector.nrows())]


def _singular_value_records(values: tuple[Any, ...], digits: int) -> list[dict[str, str]]:
    """序列化奇点分量；发散态保留文本，有限分量保留 Acb 球。"""

    return [
        {"text": value} if isinstance(value, str) else _acb_record(value, digits)
        for value in values
    ]


def _resolve_singularity_mode(request: dict[str, Any]) -> str:
    """读取公开奇点模式值并严格验证允许集合。"""

    try:
        mode = request["singularityMode"]
    except KeyError as error:
        raise ValueError("plan request must provide singularityMode") from error
    if mode not in {"avoid", "singularity_jump"}:
        raise ValueError(
            'singularityMode must be exactly "avoid" or "singularity_jump"'
        )
    return mode


def _planning_route(
    system: PartialFractionSystem | RationalMatrixSystem,
) -> tuple[
    str,
    PartialFractionSystem | RationalMatrixSystem,
    dict[str, Any] | None,
    dict[str, Any],
]:
    """自动选择公式快速路线或通用有理矩阵路线。"""

    if isinstance(system, PartialFractionSystem):
        return (
            "explicitPolynomialSimplePoles",
            system,
            None,
            {
                "schema": "flintnde_polynomial_simple_pole_specialization_v1",
                "eligible": True,
                "reason": "explicit partialFraction representation",
                "polynomial_degree": system.polynomial_degree,
            },
        )
    inventory = system.singularity_inventory()
    specialized, report = system.specialize_polynomial_simple_poles(inventory)
    if specialized is not None:
        return "certifiedPolynomialSimplePoles", specialized, inventory.to_json(), report
    return "generalRationalMatrix", system, inventory.to_json(), report


def _run_plan_request(request: dict[str, Any]) -> dict[str, Any]:
    """从当前 plan 请求的原始用户点生成一次可序列化计划。"""

    digits = request["outputDigits"]
    working_digits = request["workingPrecisionDigits"]
    if any(
        isinstance(value, bool) or not isinstance(value, int) or value <= 0
        for value in (digits, working_digits)
    ):
        raise ValueError("plan precision fields must be positive integers")
    configure_working_precision(working_digits)
    system = _build_system(request["system"])
    solver_kind, planning_system, inventory, specialization = _planning_route(system)
    start = _entry_to_acb(request["start"])
    point_records = request["points"]
    if not isinstance(point_records, list) or not point_records:
        raise ValueError("FlintNDE plan action needs a nonempty point list")
    points = [_entry_to_acb(item) for item in point_records]
    mode = _resolve_singularity_mode(request)
    language = request["messageLanguage"]
    if language not in {"EN", "CN"}:
        raise ValueError('messageLanguage must be exactly "EN" or "CN"')
    caught_messages: list[str] = []
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", UserWarning)
            if isinstance(planning_system, PartialFractionSystem):
                plan = plan_transport_path(
                    planning_system,
                    start,
                    points,
                    radius_fraction=float(request["radiusFraction"]),
                    singularity_mode=mode,
                    singularity_jump_threshold=float(
                        request["singularityJumpThreshold"]
                    ),
                    match_fraction=float(request["matchFraction"]),
                    max_singularity_jumps=int(request["maxSingularityJumps"]),
                    message_language=language,
                )
                serialized_path = plan.to_json(digits=max(20, working_digits))
                path_kind = "plannedPath"
                plan_report = plan.report
            else:
                named_detours = tuple(
                    NamedPoint(f"user_point_{index:03d}", point)
                    for index, point in enumerate(points[:-1], start=1)
                )
                path = build_adaptive_path(
                    planning_system,
                    NamedPoint("start", start),
                    NamedPoint("target", points[-1]),
                    detour_points=named_detours,
                    max_step_over_radius=float(request["maxStepOverRadius"]),
                    singularity_mode=mode,
                    message_language=language,
                )
                serialized_path = adaptive_path_to_json(
                    path, digits=max(20, working_digits)
                )
                path_kind = "adaptivePath"
                plan_report = {
                    "schema": "flintnde_adaptive_path_plan_v1",
                    "planning_action": path.planning_action,
                    "singularity_mode": path.singularity_mode,
                    "message_language": path.message_language,
                    "point_count": len(path),
                    "internal_singularity_count": len(path.internal_singularities),
                    "messages": [],
                }
        caught_messages = [str(item.message) for item in caught]
    except SingularPathError as error:
        return {
            "schema": RESULT_SCHEMA,
            "status": "singularPathRefused",
            "operation": "plan",
            "message": str(error),
            "singularityMode": mode,
            "messageLanguage": language,
            "singularPathPairs": [
                [_acb_record(first, digits), _acb_record(second, digits)]
                for first, second in error.singular_path_pairs
            ],
            "singularityInventory": inventory,
            "specialization": specialization,
        }
    except AdaptivePathSingularityError as error:
        return {
            "schema": RESULT_SCHEMA,
            "status": "singularPathRefused",
            "operation": "plan",
            "message": str(error),
            "singularityMode": mode,
            "messageLanguage": language,
            "singularPathReport": error.report,
            "singularityInventory": inventory,
            "specialization": specialization,
        }
    return {
        "schema": RESULT_SCHEMA,
        "status": "complete",
        "operation": "plan",
        "planningAction": "raw_points_automatic_plan",
        "solverKind": solver_kind,
        "pathKind": path_kind,
        "singularityMode": mode,
        "messageLanguage": language,
        "singularityInventory": inventory,
        "specialization": specialization,
        "planReport": plan_report,
        "warnings": caught_messages,
        "plan": {
            "schema": "flintnde_bridge_execution_plan_v1",
            "planningPrecisionDigits": working_digits,
            "solverKind": solver_kind,
            "pathKind": path_kind,
            "path": serialized_path,
            "specialization": specialization,
        },
    }

def _execution_plan_precision(plan_record: dict[str, Any]) -> int:
    """读取并交叉核对当前执行计划精度；冲突时要求重新规划。"""

    require_exact_keys(
        plan_record,
        {
            "schema",
            "planningPrecisionDigits",
            "solverKind",
            "pathKind",
            "path",
            "specialization",
        },
        "execution plan",
    )
    if plan_record["schema"] != "flintnde_bridge_execution_plan_v1":
        raise ValueError("unsupported FlintNDE bridge execution-plan schema")
    planning_digits = plan_record["planningPrecisionDigits"]
    if (
        isinstance(planning_digits, bool)
        or not isinstance(planning_digits, int)
        or planning_digits <= 0
    ):
        raise ValueError(
            "execution plan lacks a positive integer planning precision; replan"
        )
    path_record = plan_record["path"]
    if not isinstance(path_record, dict):
        raise ValueError("execution plan path must be an object; replan")
    serialized_digits = path_record.get("planning_precision_digits")
    if (
        isinstance(serialized_digits, bool)
        or not isinstance(serialized_digits, int)
        or serialized_digits <= 0
    ):
        raise ValueError(
            "serialized path lacks a positive integer planning precision; replan"
        )
    if serialized_digits != planning_digits:
        raise ValueError(
            "execution-plan and serialized-path precisions disagree; replan"
        )
    return planning_digits

def _restore_execution_path(
    system: PartialFractionSystem | RationalMatrixSystem,
    plan_record: dict[str, Any],
) -> tuple[
    PartialFractionSystem | RationalMatrixSystem,
    Any,
]:
    """按当前计划声明的 solver 恢复路径，并核对系统结构。"""

    solver_kind = plan_record["solverKind"]
    path_kind = plan_record["pathKind"]
    if solver_kind == "explicitPolynomialSimplePoles":
        if not isinstance(system, PartialFractionSystem):
            raise ValueError("execution plan expects an explicit partialFraction system")
        execution_system: PartialFractionSystem | RationalMatrixSystem = system
    elif solver_kind == "certifiedPolynomialSimplePoles":
        if not isinstance(system, RationalMatrixSystem):
            raise ValueError("execution plan expects a rationalMatrix source system")
        inventory = system.singularity_inventory()
        specialized, report = system.specialize_polynomial_simple_poles(inventory)
        if specialized is None or report["eligible"] is not True:
            raise ValueError("the current rational system no longer passes specialization")
        execution_system = specialized
    elif solver_kind == "generalRationalMatrix":
        if not isinstance(system, RationalMatrixSystem):
            raise ValueError("execution plan expects a rationalMatrix source system")
        execution_system = system
    else:
        raise ValueError(f"unsupported execution-plan solver kind: {solver_kind}")
    if path_kind == "plannedPath":
        if not isinstance(execution_system, PartialFractionSystem):
            raise ValueError("plannedPath execution needs a polynomial-simple-pole system")
        return execution_system, planned_path_from_json(
            plan_record["path"], system=execution_system
        )
    if path_kind == "adaptivePath":
        if not isinstance(execution_system, RationalMatrixSystem):
            raise ValueError("adaptivePath execution needs a RationalMatrixSystem")
        return execution_system, adaptive_path_from_json(
            execution_system, plan_record["path"]
        )
    raise ValueError(f"unsupported execution-plan path kind: {path_kind}")

def _run_execute_request(request: dict[str, Any]) -> dict[str, Any]:
    """直接消费 plan 动作的完整成功结果；本函数不调用路径规划器。"""

    planned_result = request["plannedResult"]
    require_exact_keys(
        planned_result,
        {
            "schema",
            "status",
            "operation",
            "planningAction",
            "solverKind",
            "pathKind",
            "singularityMode",
            "messageLanguage",
            "singularityInventory",
            "specialization",
            "planReport",
            "warnings",
            "plan",
        },
        "plannedResult",
    )
    if (
        planned_result["schema"] != RESULT_SCHEMA
        or planned_result["status"] != "complete"
        or planned_result["operation"] != "plan"
    ):
        raise ValueError("plannedResult must be a complete FlintNDE plan result")
    plan_record = planned_result["plan"]
    planning_digits = _execution_plan_precision(plan_record)
    if (
        planned_result["solverKind"] != plan_record["solverKind"]
        or planned_result["pathKind"] != plan_record["pathKind"]
    ):
        raise ValueError("plannedResult and its execution plan disagree")
    digits = request["outputDigits"]
    working_digits = request["workingPrecisionDigits"]
    if any(
        isinstance(value, bool) or not isinstance(value, int) or value <= 0
        for value in (digits, working_digits)
    ):
        raise ValueError("execute precision fields must be positive integers")
    if working_digits > planning_digits:
        raise ValueError(
            "execution working precision "
            f"{working_digits} exceeds planned precision {planning_digits}; "
            "replan at the requested working precision"
        )
    configure_working_precision(working_digits)
    source_system = _build_system(request["system"])
    execution_system, path = _restore_execution_path(source_system, plan_record)
    initial_vector = _column_vector(request["initialVector"], "initialVector")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", UserWarning)
        if plan_record["pathKind"] == "plannedPath":
            result = transport_planned_path_refined(
                execution_system,
                initial_vector,
                path,
                primary_order=int(request["primaryOrder"]),
                reference_order=int(request["referenceOrder"]),
                radius_fraction=float(request["radiusFraction"]),
                target_relative_error=request["targetRelativeError"],
                certification_mode=request["certificationMode"],
            )
        else:
            result = transport_path_refined(
                execution_system,
                initial_vector,
                path,
                primary_order=int(request["primaryOrder"]),
                reference_order=int(request["referenceOrder"]),
                radius_fraction=float(request["radiusFraction"]),
                target_relative_error=request["targetRelativeError"],
                certification_mode=request["certificationMode"],
            )
    output = {
        "schema": RESULT_SCHEMA,
        "status": "complete",
        "operation": "execute",
        "executionAction": "execute_existing_plan_without_replanning",
        "planningPrecisionDigits": planning_digits,
        "workingPrecisionDigits": working_digits,
        "solverKind": plan_record["solverKind"],
        "pathKind": plan_record["pathKind"],
        "certificationMode": result["certification_mode"],
        "certificationModeRequested": result["certification_mode_requested"],
        "primaryFinalVector": _vector_records(
            result["primary_snapshots"][-1], digits
        ),
        "referenceFinalVector": _vector_records(
            result["reference_snapshots"][-1], digits
        ),
        "relativeDifferenceInf": result["relative_difference_inf"].mid().str(digits),
        "targetRelativeErrorMet": result["target_relative_error_met"],
        "primarySeconds": float(result["primary_seconds"]),
        "referenceSeconds": float(result["reference_seconds"]),
        "warnings": [str(item.message) for item in caught],
    }
    if "segment_truncation_differences_midpoint" in result:
        output["segmentTruncationDifferences"] = result[
            "segment_truncation_differences_midpoint"
        ]
    if "sample_results" in result:
        output["samplePoints"] = [
            {
                "coordinate": item["coordinate"],
                "userPointIndex": item.get("user_point_index"),
                "source": item.get("source"),
                "value": _vector_records(item["value"], digits),
            }
            for item in result["sample_results"]
        ]
    if "singular_target_results" in result:
        output["singularTargets"] = [
            {
                "coordinate": item["coordinate"],
                "userPointIndex": item.get("user_point_index"),
                "source": item.get("source"),
                "classification": item["classification"],
                "componentClassifications": list(
                    item["component_classifications"]
                ),
                "values": _singular_value_records(item["values"], digits),
                "report": item["report"],
            }
            for item in result["singular_target_results"]
        ]
    return output


def _run_evaluate_request(request: dict[str, Any]) -> dict[str, Any]:
    """对一个固定 ep job 连续执行规划与输运，供有界任务池调用。"""

    plan_request = {
        "schema": request["schema"],
        "action": "plan",
        "system": request["system"],
        "start": request["start"],
        "points": request["points"],
        "workingPrecisionDigits": request["workingPrecisionDigits"],
        "outputDigits": request["outputDigits"],
        "messageLanguage": request["messageLanguage"],
        "singularityMode": request["singularityMode"],
        "radiusFraction": request["radiusFraction"],
        "maxStepOverRadius": request["maxStepOverRadius"],
        "singularityJumpThreshold": request["singularityJumpThreshold"],
        "matchFraction": request["matchFraction"],
        "maxSingularityJumps": request["maxSingularityJumps"],
    }
    planned = _run_plan_request(plan_request)
    if planned.get("status") != "complete":
        return {"ep": request["ep"], "status": "planFailed", "plan": planned}
    execute_request = {
        "schema": request["schema"],
        "action": "execute",
        "system": request["system"],
        "initialVector": request["initialVector"],
        "plannedResult": planned,
        "workingPrecisionDigits": request["workingPrecisionDigits"],
        "outputDigits": request["outputDigits"],
        "primaryOrder": request["primaryOrder"],
        "referenceOrder": request["referenceOrder"],
        "targetRelativeError": request["targetRelativeError"],
        "certificationMode": request["certificationMode"],
        "radiusFraction": request["radiusFraction"],
        "messageLanguage": request["messageLanguage"],
    }
    executed = _run_execute_request(execute_request)
    return {
        "ep": request["ep"],
        "status": executed.get("status", "error"),
        "plan": planned,
        "execution": executed,
    }

def run_request(request: dict[str, Any]) -> dict[str, Any]:
    """只分发当前 plan 或 execute 请求，不接受其它字段或 action。"""

    if not isinstance(request, dict):
        raise ValueError("FlintNDE bridge request must be an object")
    if request.get("schema") != REQUEST_SCHEMA:
        raise ValueError(f"unsupported bridge request schema: {request.get('schema')}")
    action = request.get("action")
    if action == "ep_batch":
        require_exact_keys(
            request,
            {"schema", "action", "requests", "parallelTaskCount"},
            "ep batch request",
        )
        return run_ep_requests(
            request["requests"],
            parallel_task_count=request["parallelTaskCount"],
        )
    if action == "plan":
        require_exact_keys(
            request,
            {
                "schema",
                "action",
                "system",
                "start",
                "points",
                "workingPrecisionDigits",
                "outputDigits",
                "messageLanguage",
                "singularityMode",
                "radiusFraction",
                "maxStepOverRadius",
                "singularityJumpThreshold",
                "matchFraction",
                "maxSingularityJumps",
            },
            "plan request",
        )
        return _run_plan_request(request)
    if action == "execute":
        require_exact_keys(
            request,
            {
                "schema",
                "action",
                "system",
                "initialVector",
                "plannedResult",
                "workingPrecisionDigits",
                "outputDigits",
                "primaryOrder",
                "referenceOrder",
                "targetRelativeError",
                "certificationMode",
                "radiusFraction",
                "messageLanguage",
            },
            "execute request",
        )
        if request["messageLanguage"] not in {"EN", "CN"}:
            raise ValueError('messageLanguage must be exactly "EN" or "CN"')
        return _run_execute_request(request)
    if action == "evaluate":
        require_exact_keys(
            request,
            {
                "schema",
                "action",
                "ep",
                "system",
                "start",
                "points",
                "initialVector",
                "workingPrecisionDigits",
                "outputDigits",
                "primaryOrder",
                "referenceOrder",
                "targetRelativeError",
                "certificationMode",
                "messageLanguage",
                "singularityMode",
                "radiusFraction",
                "maxStepOverRadius",
                "singularityJumpThreshold",
                "matchFraction",
                "maxSingularityJumps",
            },
            "evaluate request",
        )
        return _run_evaluate_request(request)
    raise ValueError(
        'FlintNDE bridge request action must be exactly "plan", "execute", '
        '"evaluate", or "ep_batch"'
    )


def _run_indexed_request(item: tuple[int, dict[str, Any]]) -> tuple[int, dict[str, Any]]:
    """在独立 worker 进程执行一个请求，并保留原始输入序号。"""

    index, request = item
    return index, run_request(request)


def run_ep_requests(
    requests: list[dict[str, Any]],
    *,
    parallel_task_count: int = DEFAULT_PARALLEL_TASK_COUNT,
) -> dict[str, Any]:
    """用有界进程池执行互相独立的固定 ep NDE 请求。

    ``parallel_task_count`` 缺省为 12；实际 worker 数严格取请求数与该值的较小者。
    任务数超过 worker 数时，``ProcessPoolExecutor`` 会在任一任务结束后自动提交队列中的
    后续任务。独立进程隔离 python-flint 的全局精度和线程上下文，返回结果按输入顺序重排。
    """

    if not isinstance(requests, list) or not requests:
        raise ValueError("requests must be a nonempty list")
    if (
        isinstance(parallel_task_count, bool)
        or not isinstance(parallel_task_count, int)
        or parallel_task_count < 1
    ):
        raise ValueError("parallel_task_count must be a positive integer")
    effective_count = min(parallel_task_count, len(requests))
    print(
        "FlintNDE ep task pool: "
        f"requested={parallel_task_count}, effective={effective_count}; "
        "default parallel_task_count=12."
    )
    payloads = [(index, request) for index, request in enumerate(requests)]
    with get_context("spawn").Pool(
        processes=effective_count, maxtasksperchild=1
    ) as pool:
        completed = pool.map(_run_indexed_request, payloads, chunksize=1)
    results = [result for _index, result in completed]
    return {
        "schema": "flintnde_ep_batch_v1",
        "status": "complete",
        "parallelTaskCountRequested": parallel_task_count,
        "parallelTaskCountEffective": effective_count,
        "taskCount": len(results),
        "results": results,
    }


def _mma_expression(value: Any) -> str:
    """递归把 JSON 型结果写成 Wolfram Association/List 表达式。"""

    if value is None:
        return "Null"
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, str):
        return _mma_string(value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return _mma_float_literal(value)
    if isinstance(value, dict):
        rules = ", ".join(
            f"{_mma_string(str(key))} -> {_mma_expression(item)}"
            for key, item in value.items()
        )
        return "<|" + rules + "|>"
    if isinstance(value, (list, tuple)):
        return "{" + ", ".join(_mma_expression(item) for item in value) + "}"
    raise TypeError(f"cannot serialize bridge result value of type {type(value).__name__}")


def _serialize_generic_result(output: dict[str, Any]) -> str:
    """序列化两阶段 bridge 的任意嵌套 Association 结果。"""

    return "FlintNDEBridgeResult = " + _mma_expression(output) + ";\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="FlintNDE Wolfram bridge")
    parser.add_argument("request", help="JSON request file written by the Wolfram loader")
    parser.add_argument("output", help="target .m file loadable by Wolfram Get")
    args = parser.parse_args(argv)
    try:
        request = json.loads(Path(args.request).read_text(encoding="utf-8"))
        output = run_request(request)
        serialized = _serialize_generic_result(output)
        Path(args.output).write_text(serialized, encoding="utf-8")
    except Exception as error:  # noqa: BLE001 - 桥接把任何失败回传给 Wolfram 端
        failure = f'FlintNDEBridgeResult = <| "schema" -> "{RESULT_SCHEMA}", "status" -> "error", "message" -> {_mma_string(str(error))} |>;\n'
        Path(args.output).write_text(failure, encoding="utf-8")
        print(f"flintnde.mathematica_bridge: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
