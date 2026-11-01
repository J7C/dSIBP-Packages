"""一般解析矩阵系统与单极点部分分式系统的解系数生成。

本模块提供两条解泰勒系数路线：

- ``AnalyticMatrixSystem``：黑盒矩阵函数的 Cauchy--DFT 圆周采样重建，只要求采样圆
  内解析；单位根缓存与稀疏采样掩码避免逐段重复构造 roots 与遍历
  恒零矩阵元。
- ``PartialFractionSystem``：``A(z)=P(z)+Σ R_j/(z-p_j)`` 型多项式加有限简单
  极点系统的公式递推；``P`` 可为任意有限次数。递推全程 acb 球算术，不做
  Cauchy 采样，也不丢弃误差球。
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

from flint import acb, acb_mat, arb


MatrixFunction = Callable[[acb], acb_mat]

_ROOTS_CACHE: dict[tuple[int, int], list[acb]] = {}


def _unit_roots(count: int) -> list[acb]:
    """按 ``(count, ctx.prec)`` 缓存单位根，避免逐段重复指数求值。"""

    from flint import ctx

    key = (count, int(ctx.prec))
    cached = _ROOTS_CACHE.get(key)
    if cached is None:
        cached = [
            acb(0, 2 * arb.pi() * index / count).exp() for index in range(count)
        ]
        _ROOTS_CACHE[key] = cached
    return cached


@dataclass(frozen=True)
class AnalyticMatrixSystem:
    """保存一般解析矩阵函数、维数和用于收敛圆控制的奇点。

    ``sparsity_mask`` 列出需要求值的 ``(row, column)`` 矩阵元；其余矩阵元按恒零
    处理，采样时直接跳过。缺省 ``None`` 表示稠密求值。
    """

    matrix_function: MatrixFunction
    dimension: int
    singularities: tuple[acb, ...] = ()
    name: str = "analytic-matrix-system"
    sparsity_mask: tuple[tuple[int, int], ...] | None = None

    def evaluate(self, point: acb) -> acb_mat:
        """计算矩阵并验证输出维数。"""

        matrix = self.matrix_function(point)
        if matrix.nrows() != self.dimension or matrix.ncols() != self.dimension:
            raise ValueError(f"{self.name}: matrix dimension mismatch at {point}")
        return matrix

    def _sparse_samples(self, point: acb) -> acb_mat:
        """按稀疏掩码只计算声明的非零矩阵元，其余置零球。"""

        full = self.evaluate(point)
        sparse = acb_mat(self.dimension, self.dimension)
        for row, column in self.sparsity_mask:
            sparse[row, column] = full[row, column]
        return sparse

    def nearest_singularity_distance(self, point: acb) -> arb:
        """返回到最近已知奇点的距离；未声明奇点时要求调用者显式给半径。"""

        if not self.singularities:
            raise ValueError(f"{self.name}: no singularities declared for automatic radius control")
        return min((abs(point - singularity) for singularity in self.singularities), key=lambda x: float(x.mid()))

    def taylor_matrix_coefficients(
        self,
        center: acb,
        order: int,
        *,
        radius: acb,
        sample_count: int | None = None,
    ) -> list[acb_mat]:
        """用 Cauchy 圆周样本和离散 Fourier 投影生成 ``A_0..A_(order-1)``。

        `sample_count` 必须大于目标阶数；提高它可以压低高阶解析项的离散混叠。
        收敛误差应由另一组阶数和样本数独立复算，而不是把 Acb 采样半径当成完整
        截断误差证明。
        """

        if order <= 0:
            raise ValueError("Taylor order must be positive")
        count = sample_count or max(32, 2 * order)
        if count <= order:
            raise ValueError("Cauchy sample_count must exceed the Taylor order")
        if abs(radius).contains(0):
            raise ValueError("Cauchy radius must be nonzero")
        roots = _unit_roots(count)
        sample_points = [center + radius * root for root in roots]
        if self.sparsity_mask is not None:
            samples = [self._sparse_samples(point) for point in sample_points]
        else:
            samples = [self.evaluate(point) for point in sample_points]
        coefficients: list[acb_mat] = []
        denominator = acb(count)
        radius_power = acb(1)
        for degree in range(order):
            total = acb_mat(self.dimension, self.dimension)
            for root, sample in zip(roots, samples):
                total += sample * (root ** (-degree))
            coefficients.append(total / denominator / radius_power)
            radius_power *= radius
        return coefficients


@dataclass(frozen=True)
class PartialFractionSystem:
    """``A(z)=P(z)+Σ_j R_j/(z-p_j)`` 型系统的公式解系数递推。

    ``constant`` 是 ``P_0``；``polynomial_coefficients`` 依次保存
    ``P_1,P_2,...``。所有有限极点
    必须是一阶且互不相同；结构是否满足这一条件由 exact rational 输入的自动认证器
    决定，显式构造者则承担该前置条件。
    """

    constant: acb_mat
    residues: tuple[acb_mat, ...]
    poles: tuple[acb, ...]
    polynomial_coefficients: tuple[acb_mat, ...] = ()
    name: str = "partial-fraction-system"

    def __post_init__(self) -> None:
        dimension = self.constant.nrows()
        if self.constant.ncols() != dimension:
            raise ValueError(f"{self.name}: constant matrix must be square")
        if len(self.residues) != len(self.poles):
            raise ValueError(f"{self.name}: residues and poles must have equal length")
        for index, residue in enumerate(self.residues):
            if residue.nrows() != dimension or residue.ncols() != dimension:
                raise ValueError(f"{self.name}: residue {index} dimension mismatch")
        for degree, coefficient in enumerate(self.polynomial_coefficients, start=1):
            if coefficient.nrows() != dimension or coefficient.ncols() != dimension:
                raise ValueError(
                    f"{self.name}: polynomial coefficient {degree} dimension mismatch"
                )
        seen: list[str] = []
        for pole in self.poles:
            key = pole.str(40)
            if key in seen:
                raise ValueError(f"{self.name}: duplicate pole at {key}; merge residues first")
            seen.append(key)

    @property
    def dimension(self) -> int:
        return self.constant.nrows()

    @property
    def singularities(self) -> tuple[acb, ...]:
        return self.poles
    @property
    def polynomial_degree(self) -> int:
        """返回显式保存的多项式最高次数。"""

        return len(self.polynomial_coefficients)

    def _shifted_polynomial_coefficients(
        self, center: acb, order: int
    ) -> list[acb_mat]:
        """用二项式公式生成 ``P(center+s)`` 的前 ``order`` 个系数。"""

        coefficients = (self.constant, *self.polynomial_coefficients)
        shifted = [acb_mat(self.dimension, self.dimension) for _ in range(order)]
        for source_degree, matrix in enumerate(coefficients):
            for local_degree in range(min(source_degree, order - 1) + 1):
                scalar = (
                    acb(math.comb(source_degree, local_degree))
                    * center ** (source_degree - local_degree)
                )
                shifted[local_degree] += matrix * scalar
        return shifted


    def evaluate(self, point: acb) -> acb_mat:
        """直接求值 ``P(z)+Σ R_j/(z-p_j)``，供对照与调试使用。"""

        polynomial = (self.constant, *self.polynomial_coefficients)
        matrix = acb_mat(polynomial[-1])
        for coefficient in reversed(polynomial[:-1]):
            matrix = matrix * point + coefficient
        for residue, pole in zip(self.residues, self.poles):
            matrix += residue / (point - pole)
        return matrix

    def nearest_singularity_distance(self, point: acb) -> arb:
        if not self.poles:
            raise ValueError(f"{self.name}: no poles declared for automatic radius control")
        return min((abs(point - pole) for pole in self.poles), key=lambda x: float(x.mid()))

    def solution_taylor_coefficients(
        self,
        center: acb,
        order: int,
        initial_vector: acb_mat,
    ) -> list[acb_mat]:
        """用公式递推直接生成 ``Y'=A Y`` 的解系数 ``c_0..c_order``。

        每步状态更新 ``u_j[n]=(c[n]-u_j[n-1])/(center-p_j)``、
        ``c[n+1]=(Σ_m P_m(center)c[n-m]+Σ_j R_j·u_j[n])/(n+1)``。
        多项式局部系数由二项式公式一次生成；递推前缀性质保证高阶运行给出的低阶
        系数与低阶运行逐位相同。
        """

        if order <= 0:
            raise ValueError("Taylor order must be positive")
        if initial_vector.nrows() != self.dimension or initial_vector.ncols() < 1:
            raise ValueError(f"{self.name}: initial vector dimension mismatch")
        solution = [acb_mat(initial_vector)]
        states = [
            acb_mat(self.dimension, initial_vector.ncols()) for _ in self.poles
        ]
        # 与几何级数展开一致取 offset=p_j-center；u_j[n]=(c[n]-u_j[n-1])/(center-p_j)
        # 的负号由 u_j=-s_j 吸收，递推贡献仍为 +R_j·u_j[n]
        offsets = [pole - center for pole in self.poles]
        for offset in offsets:
            if offset.contains(0):
                raise ValueError(
                    f"{self.name}: expansion center coincides with a pole; use local singular data"
                )
        polynomial = self._shifted_polynomial_coefficients(center, order)
        for degree in range(order):
            current = solution[degree]
            total = acb_mat(self.dimension, initial_vector.ncols())
            for polynomial_degree in range(
                min(degree, len(polynomial) - 1) + 1
            ):
                total += polynomial[polynomial_degree] * solution[
                    degree - polynomial_degree
                ]
            for index, offset in enumerate(offsets):
                states[index] = (states[index] - current) / offset
                total += self.residues[index] * states[index]
            solution.append(total / acb(degree + 1))
        return solution

    def taylor_matrix_coefficients(
        self,
        center: acb,
        order: int,
        *,
        radius: acb | None = None,
        sample_count: int | None = None,
    ) -> list[acb_mat]:
        """生成 ``A_0..A_(order-1)``；多项式与简单极点均按公式展开。

        ``radius``/``sample_count`` 仅为与 ``AnalyticMatrixSystem`` 接口一致而保留。
        """

        if order <= 0:
            raise ValueError("Taylor order must be positive")
        coefficients = self._shifted_polynomial_coefficients(center, order)
        for residue, pole in zip(self.residues, self.poles):
            offset = pole - center
            if offset.contains(0):
                raise ValueError(
                    f"{self.name}: expansion center coincides with a pole; use local singular data"
                )
            power = acb(1) / offset
            for degree in range(order):
                # R/(z-p) = -R/((p-center)-(z-center))，逐项 -R/(p-center)^{n+1}
                coefficients[degree] -= residue * power
                power /= offset
        return coefficients
