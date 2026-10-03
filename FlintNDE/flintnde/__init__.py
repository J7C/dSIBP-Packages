"""FlintNDE 的稳定公开接口。

程序包只在本文件汇总面向用户的对象；内部实现可在不改变这些导入名的前提下继续
拆分。所有数值矩阵均使用 python-flint 的 ``acb``/``acb_mat``。
"""

import os
import sys


FLINTNDE_CITATION_NOTICE = "FlintNDE package paper, arXiv identifier pending."

# 普通 import 每个进程只提醒一次；文件协议适配器可显式静默并自行保存引用 metadata。
if (
    not globals().get("_FLINTNDE_CITATION_REMINDER_SHOWN", False)
    and os.environ.get("FLINTNDE_SUPPRESS_CITATION_NOTICE") != "1"
):
    print(
        "To support FlintNDE, please acknowledge our work when you use FlintNDE "
        "in your research. Please cite: " + FLINTNDE_CITATION_NOTICE,
        file=sys.stderr,
    )
    _FLINTNDE_CITATION_REMINDER_SHOWN = True

from .core import (
    DEFAULT_WORKING_PRECISION_DIGITS,
    acb_midpoint_matrix,
    column_vector,
    configure_working_precision,
    relative_difference_inf,
    vector_norm_inf,
)
from .asymptotics import FiveTermTailDiagnostic, five_term_tail_diagnostic
from .boundary import (
    ExponentialBoundary,
    ExponentialBoundaryTerm,
    FrobeniusBoundary,
    FrobeniusBoundaryTerm,
    exponential_boundary,
    frobenius_boundary,
)
from .exact_gaussian import GaussianRational, gaussian_rational
from .frobenius import (
    PowerLogBasis,
    RegularSingularSystem,
    build_exact_frobenius_manifest,
    build_frobenius_manifest,
    build_power_log_basis,
)
from .fuchsian import (
    FuchsianReductionResult,
    MoserBalanceTransformation,
    attempt_fuchsian_reduction,
)
from .local_solutions import (
    LocalReductionError,
    LocalSolutionBasis,
    build_local_solution_basis,
)
from .regularization import (
    SeriesReconstructionResult,
    SeriesValidationError,
    fit_sampled_series,
    reconstruct_series_solution,
)
from .parallel import (
    DEFAULT_PARALLEL_TASK_COUNT,
    EpTaskBatchResult,
    run_ep_tasks,
)
from .output_layout import OutputLayout, initialize_output_layout
from .numeric_structure import (
    NumericalFrobeniusOptions,
    NumericalRegularSingularSystem,
    build_numerical_frobenius_manifest,
)
from .systems import AnalyticMatrixSystem, PartialFractionSystem
from .singularity_jump import (
    PlannedPath,
    SingularityJumpBasis,
    SingularityJumpError,
    SingularityJumpSpec,
    SingularPathError,
    direct_user_point_path,
    plan_transport_path,
    planned_path_from_json,
    planned_path_to_json,
    transport_planned_path,
    transport_planned_path_refined,
)
from .singularities import (
    RationalFunction,
    RationalMatrixSystem,
    SingularityInventory,
    SingularityRecord,
    analyze_singularities,
    rational_function,
)
from .singular_targets import (
    SingularTargetMatchPlan,
    SingularTargetResult,
    evaluate_singular_target,
    plan_singular_target_match,
)
from .routing import (
    AdaptivePath,
    AdaptivePathSingularityError,
    adaptive_path_from_json,
    adaptive_path_to_json,
    LocalExpansion,
    NamedPoint,
    PathPlan,
    PointClassification,
    build_adaptive_path,
    build_adaptive_path_plan,
    classify_point,
    prepare_local_expansion,
)
from .transport import (
    build_straight_path,
    transport_frobenius_boundaries_refined,
    transport_path,
    transport_path_refined,
)

__all__ = [
    "AdaptivePath",
    "AdaptivePathSingularityError",
    "adaptive_path_from_json",
    "adaptive_path_to_json",
    "AnalyticMatrixSystem",
    "ExponentialBoundary",
    "ExponentialBoundaryTerm",
    "FuchsianReductionResult",
    "FrobeniusBoundary",
    "FrobeniusBoundaryTerm",
    "FiveTermTailDiagnostic",
    "FLINTNDE_CITATION_NOTICE",
    "GaussianRational",
    "DEFAULT_PARALLEL_TASK_COUNT",
    "DEFAULT_WORKING_PRECISION_DIGITS",
    "EpTaskBatchResult",
    "run_ep_tasks",
    "LocalExpansion",
    "LocalReductionError",
    "LocalSolutionBasis",
    "MoserBalanceTransformation",
    "NamedPoint",
    "NumericalFrobeniusOptions",
    "NumericalRegularSingularSystem",
    "OutputLayout",
    "PartialFractionSystem",
    "PlannedPath",
    "SingularityJumpBasis",
    "SingularityJumpError",
    "SingularityJumpSpec",
    "PowerLogBasis",
    "SingularPathError",
    "PathPlan",
    "PointClassification",
    "RationalFunction",
    "RationalMatrixSystem",
    "RegularSingularSystem",
    "SeriesReconstructionResult",
    "SeriesValidationError",
    "SingularityInventory",
    "SingularityRecord",
    "acb_midpoint_matrix",
    "attempt_fuchsian_reduction",
    "build_exact_frobenius_manifest",
    "planned_path_from_json",
    "planned_path_to_json",
    "build_frobenius_manifest",
    "build_local_solution_basis",
    "build_adaptive_path",
    "build_adaptive_path_plan",
    "build_numerical_frobenius_manifest",
    "build_power_log_basis",
    "build_straight_path",
    "classify_point",
    "column_vector",
    "configure_working_precision",
    "direct_user_point_path",
    "frobenius_boundary",
    "five_term_tail_diagnostic",
    "fit_sampled_series",
    "exponential_boundary",
    "gaussian_rational",
    "analyze_singularities",
    "initialize_output_layout",
    "plan_transport_path",
    "relative_difference_inf",
    "rational_function",
    "SingularTargetResult",
    "SingularTargetMatchPlan",
    "evaluate_singular_target",
    "plan_singular_target_match",
    "reconstruct_series_solution",
    "transport_path",
    "transport_path_refined",
    "transport_planned_path",
    "transport_planned_path_refined",
    "transport_frobenius_boundaries_refined",
    "prepare_local_expansion",
    "vector_norm_inf",
]
