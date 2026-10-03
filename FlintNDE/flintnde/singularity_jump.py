"""路径预规划、前瞻节点选择与单极点数值奇点折跃。

本模块实现路径规划路线：

- ``plan_transport_path``：从起点与用户点序列整体规划输运节点。规划分两步：
  先做奇点节点化——线段穿过（或以用户阈值接近）某极点时，把该极点选为奇点折跃
  支点，在线段上放入射/出射匹配点；再沿折线做前瞻节点选择——从当前节点向后
  按用户顺序检查调用方已放入同一个单变量复参数平面的点，到哪一个为止都落在
  本节点一步可达范围，取最后一个做下一节点，步内覆盖的用户点成为该步的求值点；
  不要求这些复参数点额外落在同一条实直线上。
- ``SingularityJumpBasis``：``PartialFractionSystem`` 简单极点的数值幂对数局部基。
  留数可对角化且无共振时构造 ``Y_j(t)=t^{λ_j}Σ a_n t^n`` 型基本解组；结构
  不确定（缺陷留数、共振、特征向量矩阵奇异、残差超限）时一律 fail closed，
  不猜测 log 结构。
- ``transport_planned_path`` / ``transport_planned_path_refined``：执行规划路径，
  奇点折跃段走局部基反解，普通段复用极点--留数递推，步内求值点用该步缓存的解
  多项式系数按节点分桶；大桶使用 FLINT 子积树/余数树快速多点求值，小桶使用
  Horner，避免树构造的固定开销。

约定：不预缓存任何"每极点收敛半径"——极点位置来自系统本身，奇点折跃走廊半径在
规划时由极点到最近其它极点实时求出；步长控制用实时的最近极点距离。每极点至多
奇点折跃一次，奇点折跃总数有硬上限，超限 fail closed。
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass
import warnings
from typing import Any

from flint import acb, acb_mat, arb, ctx

from .core import (
    acb_midpoint_matrix,
    arb_ball_from_json,
    arb_ball_to_json,
    relative_difference_inf,
    require_exact_keys,
)
from .singular_targets import classify_singular_sample_vectors
from .systems import AnalyticMatrixSystem, PartialFractionSystem
from .transport import (
    _accuracy_target,
    _transport_ordinary_segment,
    _embedded_prefix_sample_results,
    _embedded_prefix_snapshots,
    evaluate_vector_series,
    evaluate_vector_series_many,
)


# 注意：规划安全系数不在模块级构造 arb——import 时 ctx.prec 还是默认
# 53 bit，模块级 arb("0.9") 会带 ~1e-17 球半径，后续所有球运算会被
# FLINT 自适应精度封顶在 ~90 bit；改在 _plan_node_walk 内按当时工作精度构造
_PLANNING_REACH_MARGIN_TEXT = "0.9"
_COLLINEARITY_TOLERANCE = 1.0e-12
# 匹配点方向相对线段方向旋转 pi/4：主支幂的割线是极点的负实方向射线。
# 旋转因子必须在调用时按当前 Arb 精度构造，不能在 import 阶段冻结为 binary64。
_MATCH_ROTATION_DIVISOR = 4
# 步长/最近奇点距离比的规划上限：每步截断尾项约为该比的 ``order`` 次幂，
# 奇点邻近段的半径就是到极点的距离，比值过大时高阶截断也压不下去。远离
# 奇点时可达范围不构成约束，此上限不增加节点数
_NEAR_SINGULARITY_RADIUS_FRACTION = 0.25
_BRANCH_CONVENTION = (
    "principal acb log/power branch; a singularity jump detours on the left of the travel "
    "direction, i.e. through the upper half-plane for a real-axis path "
    "traversed left to right"
)

def _normalize_message_language(value: str) -> str:
    """规范化运行时提示语言；只接受公开合同中的 EN/CN。"""

    if not isinstance(value, str):
        raise TypeError("message_language must be a string")
    language = value
    if language not in {"EN", "CN"}:
        raise ValueError('message_language must be "EN" or "CN"')
    return language



def _normalize_singularity_mode(value: str) -> str:
    """校验公开奇点模式是否属于允许集合。"""

    if not isinstance(value, str):
        raise TypeError("singularity_mode must be a string")
    if value not in {"avoid", "singularity_jump"}:
        raise ValueError(
            'singularity_mode must be exactly "avoid" or "singularity_jump"'
        )
    return value


def _match_rotation(sign: int) -> acb:
    """按当前 FLINT 精度构造匹配点旋转因子，避免 binary64 精度封顶。"""

    if sign not in {-1, 1}:
        raise ValueError("match rotation sign must be -1 or 1")
    angle = arb.pi() / arb(_MATCH_ROTATION_DIVISOR)
    return acb(0, arb(sign) * angle).exp()


class SingularityJumpError(NotImplementedError):
    """奇点折跃局部基的结构门禁失败；调用方应改用精确有理系统或细化路径。"""


class SingularPathError(SingularityJumpError):
    """无奇点模式下极点落在用户折线上，拒绝执行。

    ``singular_path_pairs`` 携带经过极点的相邻端点对列表，供调用方直接作为
    结构化结果输出（用户不需要复制提示文本）。
    """

    def __init__(self, message: str, singular_path_pairs: list[tuple[acb, acb]]) -> None:
        super().__init__(message)
        self.singular_path_pairs = singular_path_pairs


def _input_precision_digits(matrix: acb_mat) -> int:
    """由当前工作精度和矩阵球精度估计可靠十进制位数。"""

    digits = max(1, math.floor(int(ctx.prec) / math.log2(10)))
    for row in range(matrix.nrows()):
        for column in range(matrix.ncols()):
            accuracy = int(matrix[row, column].rel_accuracy_bits())
            if 0 < accuracy < 1_000_000_000:
                digits = min(digits, max(1, math.floor(accuracy / math.log2(10))))
    return digits


def _snap_midpoint(value: acb) -> acb:
    """以当前工作精度的中点重建标量。

    规划几何用球算术做包含性判断；但节点若带着双精度来源的球半径进入输运，
    FLINT 球算术会把后续所有运算的中点精度压到输入球精度，误差地板与阶数
    无关。规划器产出的节点与匹配点因此一律取中点。
    """

    return acb(value.real.mid(), value.imag.mid())




def _argument_sweep(start: acb, target: acb, pole: acb) -> arb:
    """用 Arb 计算不穿过极点的直线段连续辐角变化。

    输入段必须已由规划器证明不含极点；端点比值的主值辐角正是该直线段
    在 ``(-pi, pi]`` 内的连续扫掠，整个计算保留当前 FLINT 工作精度。
    """

    incoming = start - pole
    outgoing = target - pole
    if incoming.contains(0) or outgoing.contains(0):
        raise SingularityJumpError("argument sweep endpoints must not coincide with the pole")
    return (outgoing / incoming).arg()


def _nearest_integer(value: arb, label: str) -> int:
    """从 Arb 小区间提取唯一最近整数；不确定时拒绝猜测 winding。"""

    rounded = (value + arb("0.5")).floor()
    try:
        integer = rounded.unique_fmpz()
    except ValueError as error:
        raise SingularityJumpError(
            f"{label} is not a uniquely certified winding number: {value.str(30)}"
        ) from error
    return int(integer)


class SingularityJumpBasis:
    """``PartialFractionSystem`` 一个简单极点的数值幂对数局部基。

    系统写成 ``A(s)=R/(s-p)+B(s)``；``B`` 的泰勒系数由其余极点闭式给出，
    收敛半径即极点到最近其它极点的距离（单极点时 ``B`` 为常矩阵，级数全域
    收敛）。基本解列取 ``Y_j(t)=t^{λ_j}Σ_{n≤N} a_n^{(j)} t^n``（``t=s-p``，
    主支幂），其中 ``a_0^{(j)}`` 是留数 ``R`` 的特征向量，高阶系数由
    ``((n+λ_j)I-R)a_n=Σ_{m<n}B_m a_{n-1-m}`` 递推。任一结构门禁（特征向量
    残差超限、特征向量矩阵奇异、递推算子奇异、求值残差超限）失败即抛出
    ``SingularityJumpError``。
    """

    def __init__(
        self,
        system: PartialFractionSystem,
        pole_index: int,
        order: int,
    ) -> None:
        if order <= 0:
            raise ValueError("singularity-jump basis order must be positive")
        self.system = system
        self.pole_index = pole_index
        self.order = order
        self.pole = acb(system.poles[pole_index])
        self.residue = acb_mat(system.residues[pole_index])
        dimension = self.residue.nrows()
        self.dimension = dimension
        digits = _input_precision_digits(self.residue)
        gate = arb(10) ** (-arb(max(3, digits // 2)))
        self.gate_tolerance = gate

        eigenvalues, eigenvector_matrix = self.residue.eig(
            right=True, algorithm="approx"
        )
        if any("nan" in str(value).lower() for value in eigenvalues):
            raise SingularityJumpError(
                f"{system.name}: pole {self.pole.str(30)} residue spectrum is uncertain"
            )
        if (
            eigenvector_matrix.nrows() != dimension
            or eigenvector_matrix.ncols() != dimension
        ):
            raise SingularityJumpError(
                f"{system.name}: pole {self.pole.str(30)} residue did not return "
                "a full eigenvector matrix"
            )
        residue_scale = arb(0)
        for row in range(dimension):
            for column in range(dimension):
                residue_scale = max(residue_scale, abs(self.residue[row, column]).mid())
        for column in range(dimension):
            vector = acb_mat(dimension, 1)
            for row in range(dimension):
                vector[row, 0] = eigenvector_matrix[row, column]
            residual = self.residue * vector - eigenvalues[column] * vector
            residual_size = arb(0)
            vector_size = arb(0)
            for row in range(dimension):
                residual_size = max(residual_size, abs(residual[row, 0]).mid())
                vector_size = max(vector_size, abs(vector[row, 0]).mid())
            reference = residue_scale * vector_size
            if reference > 0 and arb(residual_size) > gate * reference:
                raise SingularityJumpError(
                    f"{system.name}: pole {self.pole.str(30)} eigenvector residual "
                    f"{residual_size.str(10)} exceeds the numeric gate; the residue "
                    "structure needs the exact rational route"
                )
        determinant = eigenvector_matrix.det()
        if determinant.contains(0):
            raise SingularityJumpError(
                f"{system.name}: pole {self.pole.str(30)} residue eigenvector matrix "
                "is numerically singular; the local structure needs the exact route"
            )
        self.exponents = [acb(value.mid()) for value in eigenvalues]
        self.eigenvector_matrix = eigenvector_matrix

        # B(s)=C+Σ_{j≠*} R_j/(s-p_j) 在极点 p 处的泰勒系数，闭式生成
        other_terms = [
            (acb_mat(system.residues[index]), acb(system.poles[index]))
            for index in range(len(system.poles))
            if index != pole_index
        ]
        corridor: arb | None = None
        for _residue, other_pole in other_terms:
            distance = abs(other_pole - self.pole)
            if distance.contains(0):
                raise SingularityJumpError(
                    f"{system.name}: pole {self.pole.str(30)} collides with another pole"
                )
            corridor = distance if corridor is None else min(corridor, distance)
        self.corridor_radius = corridor
        self.regular_coefficients = system._shifted_polynomial_coefficients(
            self.pole, order
        )
        for degree, polynomial_coefficient in enumerate(self.regular_coefficients):
            coefficient = acb_mat(polynomial_coefficient)
            for residue, other_pole in other_terms:
                offset = other_pole - self.pole
                power = acb(1) / offset
                for _repeat in range(degree):
                    power /= offset
                coefficient -= residue * power
            self.regular_coefficients[degree] = coefficient

        self.series: list[list[acb_mat]] = []
        for column in range(dimension):
            seed = acb_mat(dimension, 1)
            for row in range(dimension):
                seed[row, 0] = eigenvector_matrix[row, column]
            coefficients = [seed]
            for degree in range(1, order + 1):
                right_hand = acb_mat(dimension, 1)
                for regular_degree in range(degree):
                    right_hand += (
                        self.regular_coefficients[regular_degree]
                        * coefficients[degree - 1 - regular_degree]
                    )
                operator = acb_mat(dimension, dimension)
                for row in range(dimension):
                    operator[row, row] = acb(degree) + self.exponents[column]
                operator -= self.residue
                candidate = operator.solve(right_hand)
                # 球算术下奇异算子可能不抛异常：用回代残差做共振门禁
                back_residual = operator * candidate - right_hand
                back_size = arb(0)
                right_size = arb(0)
                for row in range(dimension):
                    back_size = max(back_size, abs(back_residual[row, 0]).mid())
                    right_size = max(right_size, abs(right_hand[row, 0]).mid())
                if right_size > 0 and back_size > gate * right_size:
                    raise SingularityJumpError(
                        f"{system.name}: pole {self.pole.str(30)} recurrence operator "
                        f"is singular at power {degree} for exponent "
                        f"{self.exponents[column].str(20)}; resonance needs the exact route"
                    )
                coefficients.append(candidate)
            self.series.append(coefficients)

    def _column_value(self, column: int, local: acb) -> tuple[acb_mat, acb_mat]:
        """返回一列基本解在局部点的值与其导数，供残差诊断使用。"""

        coefficients = self.series[column]
        phi = acb_mat(coefficients[-1])
        phi_derivative = acb_mat(self.dimension, 1)
        for degree in range(len(coefficients) - 1, 0, -1):
            phi_derivative = phi_derivative * local + coefficients[degree] * acb(degree)
            phi = phi * local + coefficients[degree - 1]
        power = local ** self.exponents[column]
        value = power * phi
        derivative = power * (phi_derivative + phi * (self.exponents[column] / local))
        return value, derivative

    def evaluate(self, local: acb) -> acb_mat:
        """在局部点 ``t=s-p`` 计算基本解矩阵（主支幂）。"""

        return self.lifted_evaluate(local, 0)

    def lifted_evaluate(self, local: acb, winding: int) -> acb_mat:
        """在局部点计算基本解矩阵，幂因子用 ``Log(t)+2πi·winding`` 提升。

        主支幂的割线是极点的负实方向射线；输运路径绕极点若干圈后，同一
        局部点上的解属于不同的分支，用整数提升表示。
        """

        if local.contains(0):
            raise ZeroDivisionError("a singularity-jump basis cannot be evaluated at the pole")
        matrix = acb_mat(self.dimension, self.dimension)
        for column in range(self.dimension):
            value, _derivative = self._column_value(column, local)
            if winding:
                factor = (acb(0, 2) * arb.pi() * winding * self.exponents[column]).exp()
                value = value * factor
            for row in range(self.dimension):
                matrix[row, column] = value[row, 0]
        return matrix

    def solve(self, local: acb, vector: acb_mat) -> acb_mat:
        """由入射侧匹配值反解局部常数列。"""

        return self.evaluate(local).solve(vector)

    def residual_report(self, local: acb) -> dict[str, Any]:
        """比较 ``F'(t)`` 与 ``A(p+t)F(t)``，返回相对残差诊断。"""

        if local.contains(0):
            raise ZeroDivisionError("a singularity-jump basis cannot be checked at the pole")
        value_matrix = acb_mat(self.dimension, self.dimension)
        derivative_matrix = acb_mat(self.dimension, self.dimension)
        for column in range(self.dimension):
            value, derivative = self._column_value(column, local)
            for row in range(self.dimension):
                value_matrix[row, column] = value[row, 0]
                derivative_matrix[row, column] = derivative[row, 0]
        connection = self.system.evaluate(self.pole + local)
        defect = derivative_matrix - connection * value_matrix
        defect_size = arb(0)
        value_size = arb(0)
        for row in range(self.dimension):
            for column in range(self.dimension):
                defect_size = max(defect_size, abs(defect[row, column]).mid())
                value_size = max(value_size, abs(value_matrix[row, column]).mid())
        relative = defect_size / value_size if value_size > 0 else defect_size
        return {
            "local_coordinate": local.str(40),
            "relative_residual": relative.str(30),
            "relative_residual_midpoint": float(relative.mid()),
        }

    @property
    def manifest(self) -> dict[str, Any]:
        return {
            "schema": "flintnde_singularity_jump_basis_v1",
            "pole": self.pole.str(40),
            "exponents": [exponent.str(40) for exponent in self.exponents],
            "order": self.order,
            "corridor_radius": (
                None if self.corridor_radius is None else self.corridor_radius.str(40)
            ),
            "gate_tolerance": self.gate_tolerance.str(30),
            "branch_convention": _BRANCH_CONVENTION,
        }


@dataclass(frozen=True)
class SingularityJumpSpec:
    """一条奇点折跃记录：极点、两侧匹配点与走廊几何。"""

    pole_index: int
    pole: acb
    incoming: acb
    outgoing: acb
    match_distance: arb
    corridor_radius: arb | None
    outgoing_route: tuple[acb, ...] = ()


@dataclass(frozen=True)
class PlannedPath:
    """规划结果：执行节点、奇点折跃段、步内求值点路由与规划报告。"""

    nodes: list[acb]
    singularity_jump_segments: dict[int, SingularityJumpSpec]
    sample_assignments: list[dict[str, Any]]
    report: dict[str, Any]


    def to_json(self, digits: int = 80) -> dict[str, Any]:
        """序列化完整计划，供后续进程直接执行而不重新规划。"""

        return planned_path_to_json(self, digits=digits)


def _planned_acb_record(value: acb, digits: int) -> dict[str, Any]:
    """保存 Acb 展示中点及可严格恢复的实部、虚部 Arb 球。"""

    return {
        "real": value.real.mid().str(digits, radius=False, more=True),
        "imag": value.imag.mid().str(digits, radius=False, more=True),
        "real_ball": arb_ball_to_json(value.real, digits),
        "imag_ball": arb_ball_to_json(value.imag, digits),
    }


def _planned_acb_from_record(record: Any, field_name: str) -> acb:
    """从当前计划 JSON 的严格 Arb 球字段恢复 Acb。"""

    require_exact_keys(
        record,
        {"real", "imag", "real_ball", "imag_ball"},
        field_name,
    )
    return acb(
        arb_ball_from_json(record["real_ball"], f"{field_name}.real_ball"),
        arb_ball_from_json(record["imag_ball"], f"{field_name}.imag_ball"),
    )

def planned_path_to_json(plan: PlannedPath, *, digits: int = 80) -> dict[str, Any]:
    """把 ``PlannedPath`` 转为可移植的执行计划。

    计划保存所有普通节点、奇点折跃几何和步内求值路由。执行端只需验证计划与当前
    系统的极点一致，再调用 ``transport_planned_path[_refined]``；不得再次调用规划器。
    """

    if digits < 20:
        raise ValueError("planned-path serialization needs at least 20 decimal digits")
    return {
        "schema": "flintnde_planned_path_serialized_v1",
        "planning_precision_digits": digits,
        "nodes": [_planned_acb_record(point, digits) for point in plan.nodes],
        "singularity_jump_segments": [
            {
                "segment_index": segment_index,
                "pole_index": spec.pole_index,
                "pole": _planned_acb_record(spec.pole, digits),
                "incoming": _planned_acb_record(spec.incoming, digits),
                "outgoing": _planned_acb_record(spec.outgoing, digits),
                "match_distance": spec.match_distance.mid().str(
                    digits, radius=False, more=True
                ),
                "corridor_radius": (
                    None
                    if spec.corridor_radius is None
                    else spec.corridor_radius.mid().str(
                        digits, radius=False, more=True
                    )
                ),
                "outgoing_route": [
                    _planned_acb_record(point, digits) for point in spec.outgoing_route
                ],
            }
            for segment_index, spec in sorted(plan.singularity_jump_segments.items())
        ],
        "sample_assignments": [
            {
                "segment_index": int(record["segment_index"]),
                "point": _planned_acb_record(record["point"], digits),
                "local": _planned_acb_record(record["local"], digits),
                "user_point_index": record["user_point_index"],
                "source": str(record["source"]),
            }
            for record in plan.sample_assignments
        ],
        "report": plan.report,
    }


def planned_path_from_json(
    record: dict[str, Any],
    *,
    system: AnalyticMatrixSystem | PartialFractionSystem | None = None,
) -> PlannedPath:
    """恢复并校验当前计划；只反序列化，不执行路径规划。

    若给出 ``system``，每个奇点折跃的极点编号和坐标必须与当前系统一致，避免把
    一套微分方程的计划用于另一套系统。当前 schema 的字段必须完整且不得夹带其它字段。
    """

    require_exact_keys(
        record,
        {
            "schema",
            "planning_precision_digits",
            "nodes",
            "singularity_jump_segments",
            "sample_assignments",
            "report",
        },
        "planned path",
    )
    if record["schema"] != "flintnde_planned_path_serialized_v1":
        raise ValueError("unsupported planned-path schema")
    planning_digits = record["planning_precision_digits"]
    if isinstance(planning_digits, bool) or not isinstance(planning_digits, int):
        raise ValueError("planned-path planning_precision_digits must be an integer")
    if planning_digits < 20:
        raise ValueError("planned-path planning_precision_digits must be at least 20")
    node_records = record["nodes"]
    if not isinstance(node_records, list):
        raise ValueError("planned-path nodes must be a list")
    nodes = [
        _planned_acb_from_record(item, f"nodes[{index}]")
        for index, item in enumerate(node_records)
    ]
    if len(nodes) < 2:
        raise ValueError("a serialized planned path needs at least two nodes")
    jump_records = record["singularity_jump_segments"]
    if not isinstance(jump_records, list):
        raise ValueError("planned-path singularity_jump_segments must be a list")
    singularity_jump_segments: dict[int, SingularityJumpSpec] = {}
    system_poles = None if system is None else tuple(system.singularities)
    for position, item in enumerate(jump_records):
        require_exact_keys(
            item,
            {
                "segment_index",
                "pole_index",
                "pole",
                "incoming",
                "outgoing",
                "match_distance",
                "corridor_radius",
                "outgoing_route",
            },
            f"singularity_jump_segments[{position}]",
        )
        segment_index = int(item["segment_index"])
        pole_index = int(item["pole_index"])
        if segment_index < 0 or segment_index >= len(nodes) - 1:
            raise ValueError("serialized singularity-jump segment index is outside the node chain")
        if segment_index in singularity_jump_segments:
            raise ValueError("serialized planned path contains a duplicate singularity-jump segment")
        pole = _planned_acb_from_record(item["pole"], "singularity-jump pole")
        incoming = _planned_acb_from_record(item["incoming"], "singularity-jump incoming")
        outgoing = _planned_acb_from_record(item["outgoing"], "singularity-jump outgoing")
        if not abs(nodes[segment_index] - incoming).contains(0):
            raise ValueError("serialized singularity-jump incoming point does not match its path node")
        if not abs(nodes[segment_index + 1] - outgoing).contains(0):
            raise ValueError("serialized singularity-jump outgoing point does not match its path node")
        if system_poles is not None:
            if pole_index < 0 or pole_index >= len(system_poles):
                raise ValueError("serialized singularity-jump pole index is absent from the current system")
            if not abs(system_poles[pole_index] - pole).contains(0):
                raise ValueError("serialized singularity-jump pole does not match the current system")
        route_records = item["outgoing_route"]
        if not isinstance(route_records, list):
            raise ValueError("serialized singularity-jump outgoing_route must be a list")
        singularity_jump_segments[segment_index] = SingularityJumpSpec(
            pole_index,
            pole,
            incoming,
            outgoing,
            arb(str(item["match_distance"])),
            (
                None
                if item["corridor_radius"] is None
                else arb(str(item["corridor_radius"]))
            ),
            tuple(
                _planned_acb_from_record(point, "singularity-jump outgoing route point")
                for point in route_records
            ),
        )
    assignment_records = record["sample_assignments"]
    if not isinstance(assignment_records, list):
        raise ValueError("planned-path sample_assignments must be a list")
    sample_assignments: list[dict[str, Any]] = []
    for position, item in enumerate(assignment_records):
        require_exact_keys(
            item,
            {
                "segment_index",
                "point",
                "local",
                "user_point_index",
                "source",
            },
            f"sample_assignments[{position}]",
        )
        segment_index = int(item["segment_index"])
        if segment_index < 0 or segment_index >= len(nodes) - 1:
            raise ValueError("serialized sample segment index is outside the node chain")
        user_point_index = item["user_point_index"]
        if user_point_index is not None and (
            isinstance(user_point_index, bool) or not isinstance(user_point_index, int)
        ):
            raise ValueError("serialized sample user_point_index must be an integer or null")
        if not isinstance(item["source"], str):
            raise ValueError("serialized sample source must be a string")
        sample_assignments.append(
            {
                "segment_index": segment_index,
                "point": _planned_acb_from_record(item["point"], "sample point"),
                "local": _planned_acb_from_record(item["local"], "sample local coordinate"),
                "user_point_index": user_point_index,
                "source": item["source"],
            }
        )
    report = record["report"]
    if not isinstance(report, dict):
        raise ValueError("serialized planned-path report must be an object")
    return PlannedPath(nodes, singularity_jump_segments, sample_assignments, report)

def _segment_projection(start: acb, target: acb, pole: acb) -> tuple[arb, arb]:
    """用 Arb 求极点的实投影参数和球距离，并把参数截断到 [0,1]。"""

    direction = target - start
    length_squared = (
        direction.real * direction.real + direction.imag * direction.imag
    )
    if length_squared.contains(0):
        return arb(0), abs(pole - start)
    offset = pole - start
    parameter = (
        offset.real * direction.real + offset.imag * direction.imag
    ) / length_squared
    parameter_midpoint = parameter.mid()
    if parameter_midpoint <= 0:
        clamped = arb(0)
    elif parameter_midpoint >= 1:
        clamped = arb(1)
    else:
        clamped = parameter
    projected = start + direction * acb(clamped)
    return clamped, abs(pole - projected)

def _plan_singularity_jump_segments(
    system: AnalyticMatrixSystem | PartialFractionSystem,
    chain: list[acb],
    *,
    singularity_mode: str,
    singularity_jump_threshold: float,
    match_fraction: float,
    max_singularity_jumps: int,
    message_language: str = "EN",
) -> tuple[list[dict[str, Any]], list[str], arb | None]:
    """沿折线扫描极点，返回奇点折跃段记录、规划消息与极点到路径最小距离。

    ``singularity_mode="avoid"`` 不建立奇点节点：极点落在路径上时
    报错拒绝执行，并携带经过极点的相邻端点对。只有显式设置
    ``singularity_mode="singularity_jump"`` 时，穿过极点或进入阈值范围的线段才
    建立局部基桥。每个极点至多使用一次；匹配点必须落在合理位置，否则 fail closed。
    未命中奇点时仍报告极点到路径的最小距离。
    """

    language = _normalize_message_language(message_language)
    mode = _normalize_singularity_mode(singularity_mode)
    if not 0 < match_fraction < 1:
        raise ValueError("match_fraction must lie between zero and one")
    if singularity_jump_threshold <= 0:
        raise ValueError("singularity_jump_threshold must be positive")
    rotation_left = _match_rotation(1)
    rotation_right = _match_rotation(-1)
    jumps: list[dict[str, Any]] = []
    messages: list[str] = []
    used_poles: set[int] = set()
    refused_segment_indices: set[int] = set()
    singularities = list(system.singularities)
    minimum_path_distance: arb | None = None
    singular_pairs: list[tuple[acb, acb]] = []
    for segment_index, (start, target) in enumerate(zip(chain[:-1], chain[1:])):
        length = abs(target - start)
        if length.contains(0):
            continue
        unit = (target - start) / length
        candidates: list[tuple[arb, arb, bool, int, acb, arb, arb | None]] = []
        for pole_index, pole in enumerate(singularities):
            parameter, distance = _segment_projection(start, target, pole)
            parameter_midpoint = parameter.mid()
            corridor: arb | None = None
            for other_index, other in enumerate(singularities):
                if other_index == pole_index:
                    continue
                separation = abs(other - pole)
                corridor = (
                    separation if corridor is None else min(corridor, separation)
                )
            on_segment = (
                0 < parameter_midpoint < 1
                and (
                    distance.contains(0)
                    or distance < arb(str(_COLLINEARITY_TOLERANCE)) * length
                )
            )
            if minimum_path_distance is None or distance < minimum_path_distance:
                minimum_path_distance = distance
            if on_segment and mode == "avoid":
                if segment_index not in refused_segment_indices:
                    singular_pairs.append((acb(start), acb(target)))
                    refused_segment_indices.add(segment_index)
                continue
            if mode == "avoid":
                continue
            if pole_index in used_poles:
                if on_segment:
                    raise SingularityJumpError(
                        f"singularity jump at pole {pole.str(30)} would be required "
                        "more than once; split or reroute the input polyline"
                    )
                continue
            # 阈值尺度取 min(走廊半径, 段长)：远离路径的孤立极点走廊很大，
            # 若直接以走廊半径为尺度，任何有限线段都会被判成 near miss
            threshold_scale = (
                min(corridor, length) if corridor is not None else length
            )
            near_miss = (
                0 < parameter_midpoint < 1
                and not on_segment
                and distance < arb(str(singularity_jump_threshold)) * threshold_scale
            )
            if not (on_segment or near_miss):
                continue
            if len(used_poles) + len(candidates) + 1 > max_singularity_jumps:
                raise SingularityJumpError(
                    f"singularity-jump budget exceeded: more than {max_singularity_jumps} "
                    "singularity jumps required; refine the input point sequence "
                    "or raise max_singularity_jumps"
                )
            candidates.append(
                (
                    parameter_midpoint,
                    parameter,
                    on_segment,
                    pole_index,
                    pole,
                    distance,
                    corridor,
                )
            )
        candidates.sort(key=lambda item: item[0])
        for (
            parameter_midpoint,
            parameter,
            on_segment,
            pole_index,
            pole,
            distance,
            corridor,
        ) in candidates:
            projected = start + (target - start) * acb(parameter)
            endpoint_clearance = min(abs(projected - start), abs(target - projected))
            if corridor is not None:
                match_distance = arb(str(match_fraction)) * corridor
            else:
                match_distance = arb(str(match_fraction)) * endpoint_clearance
            match_distance = min(
                match_distance, arb("0.5") * endpoint_clearance
            )
            if match_distance <= 0 or match_distance.contains(0):
                raise SingularityJumpError(
                    f"pole {pole.str(30)} leaves no room for match points inside its segment"
                )
            # 匹配点分居极点前后、整体放在绕行侧：穿过极点的线段用行进方向
            # 左侧约定（实轴左到右路径的上半平面）；near miss 线段必须与原
            # 路径保持同侧，否则奇点折跃会改变解的分支。绕行侧是拓扑选择，
            # 绝对不可翻转；分支割线的跨越由 winding 解析追踪处理
            direction = unit * rotation_left
            pole_offset = pole - projected
            pole_side = (
                unit.real * pole_offset.imag - unit.imag * pole_offset.real
            )
            if not on_segment and pole_side.contains(0):
                raise SingularityJumpError(
                    f"pole {pole.str(30)} has an uncertain detour side; "
                    "increase the working precision"
                )
            detour_left = on_segment or pole_side < 0
            if detour_left:
                incoming = projected - match_distance * unit * rotation_right
                outgoing = projected + match_distance * direction
            else:
                incoming = projected - match_distance * direction
                outgoing = projected + match_distance * unit * rotation_right
            used_poles.add(pole_index)
            jumps.append(
                {
                    "segment_index": segment_index,
                    "spec": SingularityJumpSpec(
                        pole_index,
                        acb(pole),
                        _snap_midpoint(incoming),
                        _snap_midpoint(outgoing),
                        match_distance,
                        corridor,
                    ),
                    "on_segment": on_segment,
                    "closest_distance": distance,
                }
            )
            if language == "CN":
                messages.append(
                    f"已在第 {segment_index} 段的奇点 {pole.str(30)} 规划奇点折跃"
                    f"（匹配距离 {match_distance.str(20)}）"
                )
            else:
                messages.append(
                    f"planned singularity jump at {pole.str(30)} on segment {segment_index} "
                    f"(match distance {match_distance.str(20)})"
                )
    if singular_pairs:
        pair_text = "; ".join(
            f"({start_point.str(30)} -> {end_point.str(30)})"
            for start_point, end_point in singular_pairs
        )
        message = (
            "缺省避开奇点模式拒绝执行：以下用户折线路段经过奇点："
            f"{pair_text}。请修改用户路径使各段避开奇点；也可显式设置 "
            "singularity_mode='singularity_jump' 启用奇点折跃。奇点折跃选择的"
            "多值分支等价于某一绕行路径，必须由用户确认。"
            if language == "CN"
            else
            "Default avoid-singularity mode refuses to run: poles lie on user "
            f"polyline segments {pair_text}. Reroute the user points, or explicitly "
            "set singularity_mode='singularity_jump'. A singularity jump selects a "
            "multivalued branch equivalent to a detour path; the user must confirm it."
        )
        raise SingularPathError(message, singular_pairs)
    return jumps, messages, minimum_path_distance


def _plan_node_walk(
    system: AnalyticMatrixSystem | PartialFractionSystem,
    chain: list[acb],
    jumps: list[dict[str, Any]],
    radius_fraction: float,
    user_point_indices: list[int] | None = None,
) -> tuple[list[acb], dict[int, SingularityJumpSpec], list[dict[str, Any]]]:
    """按奇点折跃边界与前瞻规则走折线，返回节点、奇点折跃段号与步内求值点路由。

    前瞻规则：从当前节点向后按用户顺序检查途经点，到哪一个为止都在本节点一步
    可达范围（``radius_fraction`` 倍最近奇点距离再乘规划安全系数）且与当前节点
    复参数平面内可共享系数，取最后一个做下一节点；被跨过的用户点成为该步的求值点。奇点折跃入射/出射
    匹配点是不可跨越的必经点。一步够不到最近必经点时沿其方向插入中间节点。
    """

    must_pass: list[tuple[acb, str, SingularityJumpSpec | None, int | None]] = []
    if user_point_indices is None:
        user_point_indices = list(range(len(chain) - 1))
    if len(user_point_indices) != len(chain) - 1:
        raise ValueError("user_point_indices must align with the chain targets")
    reach_scale = arb(_PLANNING_REACH_MARGIN_TEXT) * arb(str(radius_fraction))
    incoming_bucket_by_spec: dict[int, list[tuple[acb, int]]] = {}
    incoming_bucket_target_indices: set[int] = set()
    for jump in jumps:
        spec = jump["spec"]
        if spec.corridor_radius is None:
            continue
        singular_reach = reach_scale * spec.corridor_radius
        covered: list[tuple[acb, int]] = []
        # jump 位于 chain[segment_index] -> chain[segment_index+1]；入射侧
        # 用户点最多到前一段终点，即 target_index=segment_index-1。
        for target_index in range(jump["segment_index"] - 1, -1, -1):
            target = chain[target_index + 1]
            if abs(target - spec.pole) > singular_reach:
                break
            covered.append((target, user_point_indices[target_index]))
            incoming_bucket_target_indices.add(target_index)
        if covered:
            incoming_bucket_by_spec[id(spec)] = list(reversed(covered))
    # user_point_index 保留原始用户序号；移除精确奇点规划点后也不能重编号。
    for segment_index, (_start, target) in enumerate(zip(chain[:-1], chain[1:])):
        segment_jumps = [jump for jump in jumps if jump["segment_index"] == segment_index]
        if len(segment_jumps) > 1:
            segment_jumps.sort(
                key=lambda item: _segment_projection(
                    chain[segment_index], target, item["spec"].pole
                )[0].mid()
            )
        for jump in segment_jumps:
            must_pass.append((jump["spec"].incoming, "jump_in", jump["spec"], None))
            must_pass.append((jump["spec"].outgoing, "jump_out", jump["spec"], None))
        if segment_index not in incoming_bucket_target_indices:
            must_pass.append((target, "user", None, user_point_indices[segment_index]))
    nodes = [acb(chain[0])]
    singularity_jump_segments: dict[int, SingularityJumpSpec] = {}
    samples: list[dict[str, Any]] = []
    has_singularities = bool(system.singularities)
    position = 0
    subdivision_guard = 0
    while position < len(must_pass):
        current = nodes[-1]
        # 先消耗与当前节点重合的用户点（值即节点值，挂到下一步局部坐标 0）
        while position < len(must_pass):
            point, kind, _spec, user_index = must_pass[position]
            if kind != "user" or not abs(point - current).contains(0):
                break
            samples.append(
                {
                    "segment_index": len(nodes) - 1,
                    "point": point,
                    "local": acb(0),
                    "user_point_index": user_index,
                    "source": "coincident_with_node",
                }
            )
            position += 1
        if position >= len(must_pass):
            break
        radius = (
            system.nearest_singularity_distance(current)
            if has_singularities
            else None
        )
        reach = None if radius is None else reach_scale * radius
        best: int | None = None
        scan = position
        while scan < len(must_pass):
            point, kind, _spec, _user_index = must_pass[scan]
            delta = point - current
            within_reach = reach is None or abs(delta).contains(0) or abs(delta) <= reach
            if not within_reach:
                break
            if kind in ("jump_in", "jump_out"):
                if abs(delta).contains(0):
                    raise SingularityJumpError(
                        "a singularity-jump match point coincides with the current path node"
                    )
                if all(must_pass[item][1] == "user" for item in range(position, scan)):
                    best = scan
                break
            # points 已属于调用方声明的同一个复参数平面。只要候选点与
            # 其前面的用户点都落在当前 Taylor 收敛圆盘内，就可共享本节点
            # 系数；复平面内不要求这些求值点再落在同一条实线段上。奇点折跃
            # 匹配点仍是不可跨越的必经点。
            if any(
                must_pass[item][1] != "user"
                for item in range(position, scan)
            ):
                break
            best = scan
            scan += 1
        if best is None:
            point, _kind, _spec, _user_index = must_pass[position]
            direction = point - current
            nodes.append(
                _snap_midpoint(acb(current + direction / abs(direction) * reach))
            )
            subdivision_guard += 1
            if subdivision_guard > 100000:
                raise RuntimeError(
                    "lookahead planning exceeded the subdivision limit; a pole may lie "
                    "on the walk without a planned singularity jump"
                )
            continue
        subdivision_guard = 0
        next_segment = len(nodes) - 1
        node_point, kind, spec, _user_index = must_pass[best]
        incoming_covered_positions: set[int] = set()
        if kind == "jump_in" and spec is not None and spec.corridor_radius is not None:
            singular_reach = reach_scale * spec.corridor_radius
            for candidate_position in range(best - 1, position - 1, -1):
                candidate, candidate_kind, _candidate_spec, _candidate_user_index = (
                    must_pass[candidate_position]
                )
                if candidate_kind != "user" or abs(candidate - spec.pole) > singular_reach:
                    break
                incoming_covered_positions.add(candidate_position)
        for item in range(position, best):
            if item in incoming_covered_positions:
                continue
            interior_point, _interior_kind, _interior_spec, user_index = must_pass[item]
            local = _snap_midpoint(interior_point - current)
            samples.append(
                {
                    "segment_index": next_segment,
                    "point": interior_point,
                    "local": local,
                    "user_point_index": user_index,
                    "source": (
                        "coincident_with_node"
                        if local.contains(0)
                        else "covered_by_step"
                    ),
                }
            )
        if abs(node_point - current).contains(0):
            # 选中的下一节点就是当前节点：不造零长段，重合点直接记录
            samples.append(
                {
                    "segment_index": next_segment,
                    "point": node_point,
                    "local": acb(0),
                    "user_point_index": _user_index,
                    "source": "coincident_with_node",
                }
            )
            position = best + 1
            continue
        nodes.append(acb(node_point))
        if kind == "jump_in":
            if best + 1 >= len(must_pass):
                raise SingularityJumpError("a planned singularity jump lacks its outgoing match point")
            out_point, out_kind, out_spec, _index = must_pass[best + 1]
            if out_kind != "jump_out" or out_spec is not spec:
                raise SingularityJumpError("planned singularity-jump match points are not paired")
            jump_segment_index = len(nodes) - 1
            for covered_point, user_index in incoming_bucket_by_spec.get(id(spec), []):
                samples.append(
                    {
                        "segment_index": jump_segment_index,
                        "point": covered_point,
                        "local": _snap_midpoint(covered_point - spec.pole),
                        "user_point_index": user_index,
                        "source": "covered_by_singularity_jump_incoming",
                    }
                )
            for covered_position in sorted(incoming_covered_positions):
                covered_point, _kind, _covered_spec, user_index = must_pass[
                    covered_position
                ]
                samples.append(
                    {
                        "segment_index": jump_segment_index,
                        "point": covered_point,
                        "local": _snap_midpoint(covered_point - spec.pole),
                        "user_point_index": user_index,
                        "source": "covered_by_singularity_jump_incoming",
                    }
                )
            next_position = best + 2
            covered_positions: list[int] = []
            outgoing_route: list[acb] = []
            scan_after_pole = next_position
            singular_reach = (
                None
                if spec.corridor_radius is None
                else reach_scale * spec.corridor_radius
            )
            while scan_after_pole < len(must_pass):
                candidate, candidate_kind, _candidate_spec, _candidate_user_index = (
                    must_pass[scan_after_pole]
                )
                if candidate_kind != "user":
                    break
                if singular_reach is not None and abs(candidate - spec.pole) > singular_reach:
                    break
                covered_positions.append(scan_after_pole)
                outgoing_route.append(acb(candidate))
                scan_after_pole += 1

            if covered_positions:
                out_point = outgoing_route[-1]
                for covered_position in covered_positions[:-1]:
                    covered_point, _kind, _covered_spec, user_index = must_pass[
                        covered_position
                    ]
                    samples.append(
                        {
                            "segment_index": jump_segment_index,
                            "point": covered_point,
                            "local": _snap_midpoint(covered_point - spec.pole),
                            "user_point_index": user_index,
                            "source": "covered_by_singularity_jump",
                        }
                    )
                position = covered_positions[-1] + 1
            elif next_position < len(must_pass) and must_pass[next_position][1] == "user":
                next_user = must_pass[next_position][0]
                direction_from_pole = next_user - spec.pole
                out_point = _snap_midpoint(
                    spec.pole
                    + direction_from_pole
                    / abs(direction_from_pole)
                    * singular_reach
                )
                outgoing_route = [acb(out_point)]
                position = next_position
            else:
                outgoing_route = [acb(out_point)]
                position = next_position

            routed_spec = SingularityJumpSpec(
                spec.pole_index,
                spec.pole,
                spec.incoming,
                acb(out_point),
                spec.match_distance,
                spec.corridor_radius,
                tuple(outgoing_route),
            )
            singularity_jump_segments[jump_segment_index] = routed_spec
            nodes.append(acb(out_point))
        else:
            position = best + 1
    return nodes, singularity_jump_segments, samples


def plan_transport_path(
    system: AnalyticMatrixSystem | PartialFractionSystem,
    start: Any,
    points: list[Any],
    *,
    radius_fraction: float = 0.60,
    singularity_mode: str = "avoid",
    singularity_jump_threshold: float = 0.5,
    match_fraction: float = 0.6,
    max_singularity_jumps: int = 16,
    message_language: str = "EN",
) -> PlannedPath:
    """整体规划从 ``start`` 依次经过 ``points`` 的折线输运路径。

    返回的 ``PlannedPath`` 只含可直接执行的普通节点；奇点折跃极点作为段元数据存在，
    步内覆盖的用户点带段路由。显式 ``singularity_jump`` 模式允许一个中间用户点
    精确等于简单极点：该点从普通节点链中移除，保留为局部基目标，并要求其前后
    普通点给出的穿越段共线。末端奇点由 ``evaluate_singular_target`` 单独处理。

    ``singularity_mode="avoid"`` 是缺省模式：不采取奇点节点；极点落在路径上时抛
    ``SingularPathError`` 拒绝执行（携带经过极点的相邻端点对）；显式选择
    ``"singularity_jump"`` 才允许局部基穿越奇点。规划器同时
    并在报告中给出奇点到路径与到节点距离的最小值——该值过小会导致所需节点数
    大幅增加，但规划不阻拦。

    ``radius_fraction`` 会被收敛上限 ``_NEAR_SINGULARITY_RADIUS_FRACTION`` 封顶：
    奇点邻近段的步长/半径比决定截断收敛速度，比值过大时高阶截断也压不下去；
    封顶只增加奇点邻近的节点数，实际生效值写入报告 ``effective_radius_fraction``。
    """

    language = _normalize_message_language(message_language)
    mode = _normalize_singularity_mode(singularity_mode)
    if not 0 < radius_fraction < 1:
        raise ValueError("radius_fraction must lie between zero and one")
    if not points:
        raise ValueError("a planned path needs at least one target point")
    start_point = acb(start)
    waypoints = [acb(point) for point in points]
    singularities = list(system.singularities)
    for pole in singularities:
        if (start_point - pole).contains(0):
            raise ValueError(
                f"the path start coincides with pole {pole.str(30)}; singular data required"
            )
    regular_waypoints: list[acb] = []
    regular_user_indices: list[int] = []
    singular_user_targets: list[dict[str, Any]] = []
    for waypoint_index, waypoint in enumerate(waypoints):
        matching_poles = [
            (pole_index, pole)
            for pole_index, pole in enumerate(singularities)
            if (waypoint - pole).contains(0)
        ]
        if not matching_poles:
            regular_waypoints.append(waypoint)
            regular_user_indices.append(waypoint_index)
            continue
        if len(matching_poles) != 1:
            raise SingularityJumpError(
                f"user point {waypoint_index} matches more than one unresolved pole"
            )
        pole_index, pole = matching_poles[0]
        if waypoint_index == len(waypoints) - 1:
            raise ValueError(
                f"terminal user point {waypoint_index} coincides with pole {pole.str(30)}; "
                "use evaluate_singular_target for a terminal singular value"
            )
        if mode != "singularity_jump":
            previous = start_point if not regular_waypoints else regular_waypoints[-1]
            raise SingularPathError(
                f"user point {waypoint_index} coincides with pole {pole.str(30)}; "
                "select singularity_mode='singularity_jump' explicitly",
                [(previous, waypoint)],
            )
        next_waypoint = waypoints[waypoint_index + 1]
        if any((next_waypoint - other).contains(0) for other in singularities):
            raise SingularityJumpError("consecutive singular user points are not continuation-ready")
        previous = start_point if not regular_waypoints else regular_waypoints[-1]
        parameter, distance = _segment_projection(previous, next_waypoint, pole)
        chord_length = abs(next_waypoint - previous)
        parameter_midpoint = parameter.mid()
        collinear_crossing = (
            0 < parameter_midpoint < 1
            and (
                distance.contains(0)
                or distance < arb(str(_COLLINEARITY_TOLERANCE)) * chord_length
            )
        )
        if not collinear_crossing:
            raise SingularityJumpError(
                f"singular user point {waypoint_index} changes direction at pole "
                f"{pole.str(30)}; the continuation branch is not fixed by that input"
            )
        singular_user_targets.append(
            {
                "user_point_index": waypoint_index,
                "pole_index": pole_index,
                "pole": acb(pole),
            }
        )
    chain = [start_point, *regular_waypoints]
    if len(chain) < 2:
        raise ValueError("a planned path needs at least one ordinary target point")
    jumps, messages, minimum_path_distance = _plan_singularity_jump_segments(
        system,
        chain,
        singularity_mode=mode,
        singularity_jump_threshold=singularity_jump_threshold,
        match_fraction=match_fraction,
        max_singularity_jumps=max_singularity_jumps,
        message_language=language,
    )
    effective_fraction = min(radius_fraction, _NEAR_SINGULARITY_RADIUS_FRACTION)
    if effective_fraction < radius_fraction:
        if language == "CN":
            messages.append(
                f"为控制奇点附近截断尾项，radius_fraction 已从 {radius_fraction} "
                f"收紧为 {effective_fraction}"
            )
        else:
            messages.append(
                f"radius_fraction tightened from {radius_fraction} to "
                f"{effective_fraction} to keep the truncation tail small near singularities"
            )
    nodes, singularity_jump_segments, samples = _plan_node_walk(
        system,
        chain,
        jumps,
        effective_fraction,
        user_point_indices=regular_user_indices,
    )
    for target_record in singular_user_targets:
        matching_segments = [
            segment_index
            for segment_index, spec in singularity_jump_segments.items()
            if spec.pole_index == target_record["pole_index"]
            and abs(spec.pole - target_record["pole"]).contains(0)
        ]
        if len(matching_segments) != 1:
            raise SingularityJumpError(
                f"singular user point {target_record['user_point_index']} did not map to "
                "exactly one singularity-jump segment"
            )
        samples.append(
            {
                "segment_index": matching_segments[0],
                "point": target_record["pole"],
                "local": acb(0),
                "user_point_index": target_record["user_point_index"],
                "source": "singular_target",
            }
        )
    for index in range(len(nodes) - 1):
        if index in singularity_jump_segments:
            continue
        if not system.singularities:
            continue
        step = abs(nodes[index + 1] - nodes[index])
        radius = system.nearest_singularity_distance(nodes[index])
        if not step < radius * arb(str(effective_fraction)):
            raise RuntimeError(
                f"planned segment {index} exceeds its Cauchy disk; this indicates a planner bug"
            )
    minimum_node_distance: arb | None = None
    for node in nodes:
        if not system.singularities:
            break
        node_distance = system.nearest_singularity_distance(node)
        if minimum_node_distance is None or node_distance < minimum_node_distance:
            minimum_node_distance = node_distance
    if mode == "avoid":
        mode_message = (
            "已从原始用户点生成路径计划。当前使用避开奇点模式（缺省）；"
            "如需奇点折跃，请显式设置 singularity_mode='singularity_jump'。"
            if language == "CN"
            else
            "The path was planned from raw user points in avoid-singularity mode "
            "(default). Set singularity_mode='singularity_jump' explicitly to request "
            "a singularity jump."
        )
    else:
        mode_message = (
            "已从原始用户点生成路径计划。当前显式使用奇点折跃模式；多值分支"
            "等价于某一绕行路径，必须由用户确认。"
            if language == "CN"
            else
            "The path was planned from raw user points in explicit singularity-jump "
            "mode. The selected multivalued branch is equivalent to a detour path "
            "and must be confirmed by the user."
        )
    messages.insert(0, mode_message)
    report = {
        "schema": "flintnde_planned_path_v1",
        "planning_action": "raw_points_automatic_plan",
        "singularity_mode": mode,
        "message_language": language,
        "working_precision_bits": int(ctx.prec),
        "nodes": [node.str(40) for node in nodes],
        "segment_count": len(nodes) - 1,
        "singularity_jump_count": len(singularity_jump_segments),
        "singular_target_count": len(singular_user_targets),
        "singular_target_user_indices": [
            record["user_point_index"] for record in singular_user_targets
        ],
        "covered_sample_count": len(samples),
        "effective_radius_fraction": str(effective_fraction),
        "minimum_pole_path_distance": (
            minimum_path_distance.str(30) if minimum_path_distance is not None else None
        ),
        "minimum_pole_node_distance": (
            minimum_node_distance.str(30) if minimum_node_distance is not None else None
        ),
        "messages": messages,
    }
    return PlannedPath(nodes, singularity_jump_segments, samples, report)


def direct_user_point_path(
    system: AnalyticMatrixSystem | PartialFractionSystem,
    start: Any,
    points: list[Any],
    *,
    message_language: str = "EN",
) -> PlannedPath:
    """严格把用户点作为输运节点，不插点、不删点，也不调用路径规划器。

    本入口只执行必要的奇点安全扫描；若用户折线穿过奇点或端点落在奇点上，
    立即 fail closed。节点间距是否落在局部级数收敛圆内由输运层继续校验。
    """

    language = _normalize_message_language(message_language)
    if not points:
        raise ValueError("a direct user-point path needs at least one target point")
    start_point = acb(start)
    waypoints = [acb(point) for point in points]
    for pole in system.singularities:
        if (start_point - pole).contains(0):
            raise ValueError(
                f"the path start coincides with pole {pole.str(30)}; singular data required"
            )
        for waypoint_index, waypoint in enumerate(waypoints):
            if (waypoint - pole).contains(0):
                raise ValueError(
                    f"user point {waypoint_index} coincides with pole {pole.str(30)}; "
                    "values diverge on the singular locus"
                )
    nodes = [start_point, *waypoints]
    _, _, minimum_path_distance = _plan_singularity_jump_segments(
        system,
        nodes,
        singularity_mode="avoid",
        singularity_jump_threshold=0.5,
        match_fraction=0.6,
        max_singularity_jumps=0,
        message_language=language,
    )
    message = (
        "已关闭路径规划：FlintNDE 将按输入顺序严格使用用户点作为输运节点。"
        if language == "CN"
        else
        "Path planning is disabled: FlintNDE uses the supplied user points as "
        "transport nodes in their input order."
    )
    report = {
        "schema": "flintnde_direct_user_nodes_v1",
        "planning_action": "disabled_use_user_points_as_nodes",
        "singularity_mode": "avoid",
        "message_language": language,
        "working_precision_bits": int(ctx.prec),
        "nodes": [node.str(40) for node in nodes],
        "segment_count": len(nodes) - 1,
        "singularity_jump_count": 0,
        "covered_sample_count": 0,
        "minimum_pole_path_distance": (
            minimum_path_distance.str(30) if minimum_path_distance is not None else None
        ),
        "messages": [message],
    }
    return PlannedPath(nodes, {}, [], report)


def transport_planned_path(
    system: AnalyticMatrixSystem | PartialFractionSystem,
    initial_vector: acb_mat,
    plan: PlannedPath,
    *,
    order: int,
    sample_count: int | None = None,
    radius_fraction: float = 0.60,
    target_relative_error: Any | None = None,
    collect_patches: bool = False,
) -> tuple[list[acb_mat], list[dict[str, Any]], float, dict[str, Any]]:
    """执行规划路径：奇点折跃段走局部基，普通段走递推，步内求值点做 Horner 代入。

    第四个返回值汇总 ``sample_results``、``patches``（可选）。奇点折跃基按极点缓存，
    同一极点只构造一次；每个奇点折跃段在入射和出射匹配点各做一次相对残差诊断并
    写入报告。
    """

    if not isinstance(initial_vector, acb_mat):
        raise TypeError("a planned path requires an acb_mat initial column vector")
    _accuracy_target(target_relative_error)
    nodes = plan.nodes
    if len(nodes) < 2:
        raise ValueError("a planned path needs at least two nodes")
    total_clock = time.perf_counter()
    bases: dict[int, SingularityJumpBasis] = {}
    snapshots = [acb_mat(initial_vector)]
    reports: list[dict[str, Any]] = []
    patches: list[dict[str, Any]] | None = [] if collect_patches else None
    sample_results: list[dict[str, Any]] = []
    singular_target_results: list[dict[str, Any]] = []
    vector = acb_mat(initial_vector)
    for index, (start, target) in enumerate(zip(nodes[:-1], nodes[1:]), start=1):
        jump = plan.singularity_jump_segments.get(index - 1)
        if jump is not None:
            if not isinstance(system, PartialFractionSystem):
                raise TypeError("singularity jumps currently require a PartialFractionSystem")
            segment_clock = time.perf_counter()
            basis = bases.get(jump.pole_index)
            if basis is None:
                basis = SingularityJumpBasis(system, jump.pole_index, order)
                bases[jump.pole_index] = basis
            incoming_local = start - jump.pole
            outgoing_local = target - jump.pole
            # 分支提升：主支幂的割线是极点的负实方向射线。入射侧提升用已输运
            # 折线对极点的解析辐角扫掠逐段累加（含穿过割线的情形）；出射侧
            # 提升 = 入射侧提升加上匹配点主支辐角差——对实轴上从左到右穿过
            # 实极点的路径，这等于上半平面绕行约定
            two_pi = arb(2) * arb.pi()
            lift = (nodes[0] - jump.pole).arg()
            for previous_node, node in zip(nodes[: index - 1], nodes[1:index]):
                lift += _argument_sweep(previous_node, node, jump.pole)
            incoming_argument = incoming_local.arg()
            incoming_winding = _nearest_integer(
                (lift - incoming_argument) / two_pi, "incoming branch lift"
            )
            constants = basis.lifted_evaluate(
                incoming_local, incoming_winding
            ).solve(vector)
            incoming_assignments = [
                assignment
                for assignment in plan.sample_assignments
                if assignment["segment_index"] == index - 1
                and assignment["source"] == "covered_by_singularity_jump_incoming"
            ]
            for assignment in incoming_assignments:
                local_point = assignment["point"] - jump.pole
                point_lift = lift + _argument_sweep(start, assignment["point"], jump.pole)
                point_winding = _nearest_integer(
                    (point_lift - local_point.arg()) / two_pi,
                    "incoming singularity-bucket branch lift",
                )
                dense_value = acb_midpoint_matrix(
                    basis.lifted_evaluate(local_point, point_winding) * constants
                )
                sample_results.append(
                    {
                        "coordinate": assignment["point"].str(40),
                        "segment_index": index,
                        "user_point_index": assignment["user_point_index"],
                        "source": assignment["source"],
                        "value": acb_mat(dense_value),
                    }
                )
            singular_assignments = [
                assignment
                for assignment in plan.sample_assignments
                if assignment["segment_index"] == index - 1
                and assignment["source"] == "singular_target"
            ]
            if singular_assignments:
                sample_exponents = (8, 12, 16, 20, 24)
                incoming_direction = incoming_local / abs(incoming_local)
                sample_scale = abs(incoming_local) / arb(2)
                if basis.corridor_radius is not None:
                    sample_scale = min(
                        sample_scale,
                        basis.corridor_radius / arb(4),
                    )
                singular_vectors = [
                    basis.lifted_evaluate(
                        incoming_direction
                        * acb(sample_scale)
                        * acb(10) ** (-exponent),
                        incoming_winding,
                    )
                    * constants
                    for exponent in sample_exponents
                ]
                singular_result = classify_singular_sample_vectors(
                    singular_vectors,
                    target_point=jump.pole.str(40),
                    local_basis_method=basis.manifest["method"]
                    if "method" in basis.manifest
                    else basis.manifest["schema"],
                    local_basis_order=order,
                    sample_exponents=sample_exponents,
                    target_relative_error=(
                        target_relative_error
                        if target_relative_error is not None
                        else "1e-25"
                    ),
                    sample_scale=sample_scale,
                    convergence_radius=basis.corridor_radius,
                )
                for assignment in singular_assignments:
                    singular_target_results.append(
                        {
                            "coordinate": assignment["point"].str(40),
                            "segment_index": index,
                            "user_point_index": assignment["user_point_index"],
                            "source": assignment["source"],
                            "values": singular_result.values,
                            "component_classifications": (
                                singular_result.component_classifications
                            ),
                            "classification": singular_result.classification,
                            "report": singular_result.report,
                        }
                    )
            route_points = jump.outgoing_route or (acb(target),)
            route_start = acb(start)
            route_lift = lift
            outgoing_winding = incoming_winding
            jump_assignments = [
                assignment
                for assignment in plan.sample_assignments
                if assignment["segment_index"] == index - 1
                and assignment["source"] == "covered_by_singularity_jump"
            ]
            for route_point in route_points:
                route_lift += _argument_sweep(route_start, route_point, jump.pole)
                local_point = route_point - jump.pole
                point_winding = _nearest_integer(
                    (route_lift - local_point.arg()) / two_pi,
                    "singularity-jump route branch lift",
                )
                for assignment in jump_assignments:
                    if not abs(assignment["point"] - route_point).contains(0):
                        continue
                    dense_value = acb_midpoint_matrix(
                        basis.lifted_evaluate(local_point, point_winding) * constants
                    )
                    sample_results.append(
                        {
                            "coordinate": assignment["point"].str(40),
                            "segment_index": index,
                            "user_point_index": assignment["user_point_index"],
                            "source": assignment["source"],
                            "value": acb_mat(dense_value),
                        }
                    )
                route_start = acb(route_point)
                outgoing_winding = point_winding
            if not abs(route_start - target).contains(0):
                raise SingularityJumpError("singularity-jump outgoing route does not end at its target node")
            outgoing_lift = route_lift
            vector = acb_midpoint_matrix(
                basis.lifted_evaluate(outgoing_local, outgoing_winding) * constants
            )
            elapsed = time.perf_counter() - segment_clock
            reports.append(
                {
                    "segment_index": index,
                    "method": "singularity_jump_bridge",
                    "pole": jump.pole.str(40),
                    "pole_index": jump.pole_index,
                    "order": order,
                    "local_basis": basis.manifest,
                    "branch_convention": _BRANCH_CONVENTION,
                    "incoming_winding": incoming_winding,
                    "outgoing_winding": outgoing_winding,
                    "incoming_residual": basis.residual_report(incoming_local),
                    "outgoing_residual": basis.residual_report(outgoing_local),
                    "elapsed_seconds": elapsed,
                }
            )
            snapshots.append(acb_mat(vector))
            if patches is not None:
                # 与既有 bridge 一致：奇点折跃段无 Taylor 补丁，前缀求值直接复用参考链值
                patches.append({"segment_index": index, "method": "singular_bridge"})
            continue
        vector, segment_report, solution_coefficients = _transport_ordinary_segment(
            system,
            vector,
            start,
            target,
            index,
            order=order,
            radius_fraction=radius_fraction,
            sample_count=sample_count,
            method_label="ordinary_taylor",
        )
        snapshots.append(acb_mat(vector))
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
        segment_assignments = [
            assignment
            for assignment in plan.sample_assignments
            if assignment["segment_index"] == index - 1
        ]
        dense_values, dense_algorithm = evaluate_vector_series_many(
            solution_coefficients,
            [assignment["local"] for assignment in segment_assignments],
        )
        segment_report["dense_evaluation"] = {
            "algorithm": dense_algorithm,
            "point_count": len(segment_assignments),
        }
        for assignment, dense_value in zip(segment_assignments, dense_values):
            sample_results.append(
                {
                    "coordinate": assignment["point"].str(40),
                    "segment_index": index,
                    "user_point_index": assignment["user_point_index"],
                    "source": assignment["source"],
                    "evaluation_algorithm": dense_algorithm,
                    "value": acb_mat(dense_value),
                }
            )
    elapsed_total = time.perf_counter() - total_clock
    return (
        snapshots,
        reports,
        elapsed_total,
        {
            "sample_results": sample_results,
            "singular_target_results": singular_target_results,
            "patches": patches,
        },
    )


def transport_planned_path_refined(
    system: AnalyticMatrixSystem | PartialFractionSystem,
    initial_vector: acb_mat,
    plan: PlannedPath,
    *,
    reference_initial_vector: acb_mat | None = None,
    primary_order: int,
    reference_order: int,
    primary_sample_count: int | None = None,
    reference_sample_count: int | None = None,
    radius_fraction: float = 0.60,
    target_relative_error: Any | None = None,
    certification_mode: str = "certified",
) -> dict[str, Any]:
    """对规划路径做双链或单段嵌入式截断认证。

    多段计划无法仅靠各段参考初值的前缀传播主阶误差，含奇点折跃段也不能复用局部基；
    两种情形请求 ``embedded`` 时均升级为独立 ``certified`` 双链。用户可见节点与
    dense 值来自主链；参考链只用于精度核验。跨调用续算时可显式传入独立的
    ``reference_initial_vector``，此时也必须运行双链。
    """

    if reference_order <= primary_order:
        raise ValueError("reference_order must exceed primary_order")
    if certification_mode not in ("embedded", "certified"):
        raise ValueError('certification_mode must be "embedded" or "certified"')
    requested_certification_mode = certification_mode
    singular_bridge_forced_certified = (
        certification_mode == "embedded" and bool(plan.singularity_jump_segments)
    )
    multi_segment_forced_certified = (
        certification_mode == "embedded" and len(plan.nodes) > 2
    )
    separate_initials_forced_certified = (
        certification_mode == "embedded" and reference_initial_vector is not None
    )
    if (
        singular_bridge_forced_certified
        or multi_segment_forced_certified
        or separate_initials_forced_certified
    ):
        reasons = []
        if singular_bridge_forced_certified:
            reasons.append("singular bridges cannot be reused")
        if multi_segment_forced_certified:
            reasons.append("low-order state errors cannot be propagated by segment prefixes")
        if separate_initials_forced_certified:
            reasons.append("primary and reference chains have independent initial states")
        warnings.warn(
            "embedded certification upgraded to independent primary/reference chains: "
            + "; ".join(reasons),
            UserWarning,
            stacklevel=2,
        )
        certification_mode = "certified"
    accuracy_target = _accuracy_target(target_relative_error)
    if certification_mode == "embedded":
        snapshots, reports, elapsed, extra = transport_planned_path(
            system,
            initial_vector,
            plan,
            order=reference_order,
            sample_count=reference_sample_count,
            radius_fraction=radius_fraction,
            target_relative_error=accuracy_target,
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
        # 前缀链各段与参考链共用段初值，端点差只反映最后一段的截断；总误差
        # 以逐段截断差之和为上界（不计传播放大），认证量取两者更保守者
        accumulated_truncation = arb(0)
        for value in truncations:
            accumulated_truncation += value
        difference = max(
            relative_difference_inf(primary_snapshots[-1], snapshots[-1]),
            accumulated_truncation,
        )
        result = {
            "certification_mode": "embedded",
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
            "sample_results": _embedded_prefix_sample_results(
                extra["sample_results"],
                extra["patches"],
                primary_order,
                [assignment["local"] for assignment in plan.sample_assignments],
            ),
            "singular_target_results": extra["singular_target_results"],
        }
    else:
        primary = transport_planned_path(
            system,
            initial_vector,
            plan,
            order=primary_order,
            sample_count=primary_sample_count,
            radius_fraction=radius_fraction,
            target_relative_error=accuracy_target,
        )
        reference = transport_planned_path(
            system,
            initial_vector if reference_initial_vector is None else reference_initial_vector,
            plan,
            order=reference_order,
            sample_count=reference_sample_count,
            radius_fraction=radius_fraction,
            target_relative_error=accuracy_target,
        )
        difference = relative_difference_inf(primary[0][-1], reference[0][-1])
        primary_singular = primary[3]["singular_target_results"]
        reference_singular = reference[3]["singular_target_results"]
        if len(primary_singular) != len(reference_singular):
            raise SingularityJumpError(
                "primary/reference singular-target result counts do not agree"
            )
        singular_difference = arb(0)
        for primary_item, reference_item in zip(primary_singular, reference_singular):
            if (
                primary_item["user_point_index"] != reference_item["user_point_index"]
                or primary_item["component_classifications"]
                != reference_item["component_classifications"]
            ):
                raise SingularityJumpError(
                    "primary/reference singular-target classifications do not agree"
                )
            for primary_value, reference_value in zip(
                primary_item["values"], reference_item["values"]
            ):
                if isinstance(primary_value, str) or isinstance(reference_value, str):
                    if primary_value != reference_value:
                        raise SingularityJumpError(
                            "primary/reference singular-target divergence states do not agree"
                        )
                    continue
                denominator = max(arb(1), abs(reference_value), key=lambda item: float(item.mid()))
                singular_difference = max(
                    singular_difference,
                    abs(primary_value - reference_value) / denominator,
                )
        difference = max(difference, singular_difference)
        result = {
            "certification_mode": "certified",
            "primary_snapshots": primary[0],
            "reference_snapshots": reference[0],
            "primary_segments": primary[1],
            "reference_segments": reference[1],
            "primary_seconds": primary[2],
            "reference_seconds": reference[2],
            "sample_results": primary[3]["sample_results"],
            "singular_target_results": primary[3]["singular_target_results"],
            "reference_singular_target_results": reference[3][
                "singular_target_results"
            ],
            "singular_target_relative_difference_inf": singular_difference,
            "singular_target_relative_difference_midpoint": float(
                singular_difference.mid()
            ),
        }
    meets_target = None if accuracy_target is None else bool(difference < accuracy_target)
    result.update(
        {
            "relative_difference_inf": difference,
            "relative_difference_midpoint": float(difference.mid()),
            "target_relative_error": (
                None if accuracy_target is None else accuracy_target.str(30)
            ),
            "target_relative_error_met": meets_target,
            "plan": plan.report,
            "certification_mode_requested": requested_certification_mode,
            "singular_bridge_forced_certified": singular_bridge_forced_certified,
            "multi_segment_forced_certified": multi_segment_forced_certified,
            "separate_initials_forced_certified": separate_initials_forced_certified,
            "execution_action": "execute_existing_plan_without_replanning",
        }
    )
    return result
