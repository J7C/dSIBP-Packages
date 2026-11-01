"""exact Gaussian-rational 矩阵系统及可审计的奇点发现。

公开输入由 :class:`RationalFunction` 矩阵组成，多项式系数按幂次升序保存为 Q(i) 精确数。
每个完整矩阵元先在 Q(i)[x] 中约分，再做 square-free 分解；这既能严格验证相消，也把
运行依赖限制为 python-flint。有限根用 Acb ball 隔离，无穷远则由变换后的系数矩阵
``-sinv**(-2) A(1/sinv)`` 分类。高阶极点只认证当前输入基底非 Fuchsian；在没有
Moser/shearing 分析时，不把它提升为基底不变量意义下的不规则奇点结论。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable

from flint import acb, acb_mat, acb_poly, arb

from .exact_gaussian import (
    GaussianMatrix,
    GaussianPolynomial,
    GaussianRational,
    gaussian_rational,
    series_quotient,
)


def _exact_scalar(value: Any, field_name: str) -> GaussianRational:
    """转换一个精确系数，并拒绝机器数或 ball 输入。"""

    try:
        return gaussian_rational(value)
    except (TypeError, ValueError) as error:
        raise type(error)(
            f"{field_name} received invalid exact input {value!r}; automatic singularity "
            "discovery requires Gaussian-rational coefficients"
        ) from error


def _exact_polynomial(coefficients: Iterable[Any], field_name: str) -> GaussianPolynomial:
    """由按幂次升序排列的系数构造精确多项式。"""

    records = tuple(coefficients)
    if not records:
        raise ValueError(f"{field_name} must contain at least one coefficient")
    return GaussianPolynomial.from_coefficients(
        [_exact_scalar(value, field_name) for value in records]
    )


def _reduced_pair(
    numerator_coefficients: Iterable[Any], denominator_coefficients: Iterable[Any]
) -> tuple[GaussianPolynomial, GaussianPolynomial]:
    """约分一个有理函数，并把分母归一化为首一多项式。"""

    numerator = _exact_polynomial(numerator_coefficients, "numerator")
    denominator = _exact_polynomial(denominator_coefficients, "denominator")
    if denominator.is_zero:
        raise ZeroDivisionError("rational-function denominator must not be zero")
    if numerator.is_zero:
        return (
            GaussianPolynomial.from_coefficients([0]),
            GaussianPolynomial.from_coefficients([1]),
        )
    common = numerator.gcd(denominator)
    numerator, numerator_remainder = divmod(numerator, common)
    denominator, denominator_remainder = divmod(denominator, common)
    if not numerator_remainder.is_zero or not denominator_remainder.is_zero:
        raise ArithmeticError("Gaussian polynomial gcd did not divide exactly")
    leading = denominator.coefficient(denominator.degree)
    numerator = numerator / leading
    denominator = denominator / leading
    return numerator, denominator


def _polynomial_records(polynomial: GaussianPolynomial) -> tuple[str, ...]:
    """按幂次升序序列化精确多项式。"""

    return polynomial.records()


def _horner(polynomial: GaussianPolynomial, point: acb) -> acb:
    """不经过 SymPy，在 Acb 点计算精确有理多项式。"""

    return polynomial.evaluate(point)


def _shift_polynomial(
    polynomial: GaussianPolynomial, center: GaussianRational
) -> GaussianPolynomial:
    """返回 ``p(center + z)`` 的精确局部系数。"""

    return polynomial.shifted(center)


def _series_quotient(
    numerator: GaussianPolynomial, denominator: GaussianPolynomial, order: int
) -> list[GaussianRational]:
    """对常数项非零的分母递推商的 ``order`` 个精确 Taylor 系数。"""

    return series_quotient(numerator, denominator, order)


@dataclass(frozen=True)
class RationalFunction:
    """分子和分母系数均按幂次升序保存的 Q(i) 精确有理函数。"""

    numerator: tuple[Any, ...]
    denominator: tuple[Any, ...] = (1,)
    _exact_pair: tuple[GaussianPolynomial, GaussianPolynomial] = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        """固定系数容器，并验证可精确约分的表示。"""

        object.__setattr__(self, "numerator", tuple(self.numerator))
        object.__setattr__(self, "denominator", tuple(self.denominator))
        object.__setattr__(
            self,
            "_exact_pair",
            _reduced_pair(self.numerator, self.denominator),
        )

    def exact_polynomials(self) -> tuple[GaussianPolynomial, GaussianPolynomial]:
        """返回约分后的精确分子和首一分母。"""

        return self._exact_pair

    def evaluate(self, point: acb) -> acb:
        """在 Acb 点计算约分后的有理函数。"""

        numerator, denominator = self.exact_polynomials()
        denominator_value = _horner(denominator, point)
        if denominator_value.contains(0):
            raise ZeroDivisionError(f"rational function is singular at {point}")
        return _horner(numerator, point) / denominator_value

    def __add__(self, other: RationalFunction | Any) -> RationalFunction:
        """相加两个有理项，并立即返回精确约分后的和。"""

        right = other if isinstance(other, RationalFunction) else RationalFunction((other,))
        left_numerator, left_denominator = self.exact_polynomials()
        right_numerator, right_denominator = right.exact_polynomials()
        numerator = left_numerator * right_denominator + right_numerator * left_denominator
        denominator = left_denominator * right_denominator
        return RationalFunction(_polynomial_records(numerator), _polynomial_records(denominator))

    def __radd__(self, other: RationalFunction | Any) -> RationalFunction:
        """支持 ``sum`` 以及标量加有理函数的输入写法。"""

        return self.__add__(other)

    @property
    def is_zero(self) -> bool:
        """返回完整约分后的分子是否严格为零。"""

        return self.exact_polynomials()[0].is_zero

    def __neg__(self) -> RationalFunction:
        """返回 exact 相反数。"""

        numerator, denominator = self.exact_polynomials()
        return RationalFunction(_polynomial_records(-numerator), _polynomial_records(denominator))

    def __sub__(self, other: RationalFunction | Any) -> RationalFunction:
        """相减两个有理函数并立即约分。"""

        return self + (-_coerce_entry(other))

    def __rsub__(self, other: RationalFunction | Any) -> RationalFunction:
        """支持标量或有理函数减去当前对象。"""

        return _coerce_entry(other) - self

    def __mul__(self, other: RationalFunction | Any) -> RationalFunction:
        """相乘两个有理函数并立即约分。"""

        right = _coerce_entry(other)
        left_numerator, left_denominator = self.exact_polynomials()
        right_numerator, right_denominator = right.exact_polynomials()
        return RationalFunction(
            _polynomial_records(left_numerator * right_numerator),
            _polynomial_records(left_denominator * right_denominator),
        )

    __rmul__ = __mul__

    def __truediv__(self, other: RationalFunction | Any) -> RationalFunction:
        """exact 相除并拒绝零除数。"""

        right = _coerce_entry(other)
        right_numerator, right_denominator = right.exact_polynomials()
        if right_numerator.is_zero:
            raise ZeroDivisionError("rational-function division by zero")
        left_numerator, left_denominator = self.exact_polynomials()
        return RationalFunction(
            _polynomial_records(left_numerator * right_denominator),
            _polynomial_records(left_denominator * right_numerator),
        )

    def __rtruediv__(self, other: RationalFunction | Any) -> RationalFunction:
        """支持标量或有理函数除以当前对象。"""

        return _coerce_entry(other) / self

    def __pow__(self, power: int) -> RationalFunction:
        """用平方乘算法计算整数次幂。"""

        if isinstance(power, bool) or not isinstance(power, int):
            raise TypeError("rational-function power must be an integer")
        if power < 0:
            return (RationalFunction((1,)) / self) ** (-power)
        result = RationalFunction((1,))
        base = self
        exponent = power
        while exponent:
            if exponent & 1:
                result *= base
            base *= base
            exponent //= 2
        return result

    def valuation_at(self, point: Any) -> int | None:
        """返回局部变量 ``z=x-point`` 的 exact valuation；零函数返回 ``None``。"""

        numerator, denominator = self.exact_polynomials()
        if numerator.is_zero:
            return None
        center = _exact_scalar(point, "local point")
        shifted_numerator = _shift_polynomial(numerator, center)
        shifted_denominator = _shift_polynomial(denominator, center)

        def zero_order(polynomial: GaussianPolynomial) -> int:
            return next(
                degree
                for degree in range(polynomial.degree + 1)
                if not polynomial.coefficient(degree).is_zero
            )

        return zero_order(shifted_numerator) - zero_order(shifted_denominator)

    def pole_order_at(self, point: Any) -> int:
        """返回指定 exact 点的约分后 pole 阶。"""

        valuation = self.valuation_at(point)
        return 0 if valuation is None else max(0, -valuation)

    def local_laurent_coefficients(
        self,
        point: Any,
        minimum_power: int,
        maximum_power: int,
    ) -> tuple[GaussianRational, ...]:
        """返回闭区间幂次上的 exact Laurent 系数。"""

        if minimum_power > maximum_power:
            raise ValueError("minimum_power must not exceed maximum_power")
        width = maximum_power - minimum_power + 1
        numerator, denominator = self.exact_polynomials()
        if numerator.is_zero:
            return tuple(GaussianRational() for _ in range(width))
        center = _exact_scalar(point, "local point")
        shifted_numerator = _shift_polynomial(numerator, center)
        shifted_denominator = _shift_polynomial(denominator, center)
        numerator_order = next(
            degree
            for degree in range(shifted_numerator.degree + 1)
            if not shifted_numerator.coefficient(degree).is_zero
        )
        denominator_order = next(
            degree
            for degree in range(shifted_denominator.degree + 1)
            if not shifted_denominator.coefficient(degree).is_zero
        )
        valuation = numerator_order - denominator_order
        if maximum_power < valuation:
            return tuple(GaussianRational() for _ in range(width))
        local_numerator = GaussianPolynomial.from_coefficients(
            shifted_numerator.coefficients()[numerator_order:]
        )
        local_denominator = GaussianPolynomial.from_coefficients(
            shifted_denominator.coefficients()[denominator_order:]
        )
        quotient = _series_quotient(
            local_numerator,
            local_denominator,
            maximum_power - valuation + 1,
        )
        return tuple(
            (
                quotient[power - valuation]
                if valuation <= power <= maximum_power
                else GaussianRational()
            )
            for power in range(minimum_power, maximum_power + 1)
        )

    def inverted_variable(self) -> RationalFunction:
        """返回含微分 Jacobian 的 ``-sinv**(-2) f(1/sinv)``。"""

        numerator, denominator = self.exact_polynomials()
        if numerator.is_zero:
            return RationalFunction((0,))
        numerator_reversed = list(reversed(numerator.coefficients()))
        denominator_reversed = list(reversed(denominator.coefficients()))
        power = denominator.degree - numerator.degree - 2
        if power >= 0:
            numerator_reversed = [GaussianRational()] * power + numerator_reversed
        else:
            denominator_reversed = [GaussianRational()] * (-power) + denominator_reversed
        numerator_reversed = [-value for value in numerator_reversed]
        return RationalFunction(
            tuple(numerator_reversed),
            tuple(denominator_reversed),
        )

    def pole_order_at_infinity(self) -> int:
        """返回 ``sinv=1/s`` 及微分算符变换后的 pole 阶。"""

        numerator, denominator = self.exact_polynomials()
        if numerator.is_zero:
            return 0
        return max(0, numerator.degree - denominator.degree + 2)


def rational_function(
    numerator: Iterable[Any] | Any,
    denominator: Iterable[Any] | Any = (1,),
) -> RationalFunction:
    """构造有理函数；标量参数自动提升为单项系数列表。"""

    numerator_records = (
        tuple(numerator) if isinstance(numerator, (list, tuple)) else (numerator,)
    )
    denominator_records = (
        tuple(denominator) if isinstance(denominator, (list, tuple)) else (denominator,)
    )
    return RationalFunction(numerator_records, denominator_records)


@dataclass(frozen=True)
class SingularityRecord:
    """保存有限或无穷远点及当前输入基底中已验证的矩阵 pole 阶。"""

    identifier: str
    location: acb | None
    location_exact: str | None
    pole_order: int
    kind: str
    factor_coefficients_exact: tuple[str, ...] | None
    witness_entries: tuple[tuple[int, int], ...]

    def to_json(self) -> dict[str, Any]:
        """返回可移植 JSON 记录，不把 Acb 位置降为 binary64。"""

        location = None
        if self.location is not None:
            location = {
                "exact": self.location_exact,
                "real_midpoint": self.location.real.str(40, radius=False),
                "imag_midpoint": self.location.imag.str(40, radius=False),
                "real_ball": self.location.real.str(40),
                "imag_ball": self.location.imag.str(40),
            }
        return {
            "identifier": self.identifier,
            "location": location,
            "pole_order": self.pole_order,
            "classification": self.kind,
            "reduced_denominator_factor_coefficients": (
                list(self.factor_coefficients_exact)
                if self.factor_coefficients_exact is not None
                else None
            ),
            "witness_matrix_entries_zero_based": [list(index) for index in self.witness_entries],
        }


@dataclass(frozen=True)
class SingularityInventory:
    """保存全部已验证有限奇点，以及单独分类的无穷远点。"""

    system_name: str
    variable_name: str
    finite: tuple[SingularityRecord, ...]
    infinity: SingularityRecord

    def to_json(self) -> dict[str, Any]:
        """序列化发现路线和全部奇点位置。"""

        return {
            "schema": "flintnde_singularity_inventory_v1",
            "system_name": self.system_name,
            "variable_name": self.variable_name,
            "discovery_route": "exact Q(i)[x] entry reduction and square-free decomposition, Acb root isolation",
            "cancellation_validation": "exact Gaussian-polynomial gcd before denominator collection",
            "infinity_validation": "pole order of -sinv^(-2) A(1/sinv)",
            "finite_singularities": [record.to_json() for record in self.finite],
            "infinity": self.infinity.to_json(),
        }

    def find_finite(self, identifier: str) -> SingularityRecord | None:
        """按 identifier 查找有限奇点记录；未命中返回 ``None``。"""

        return next(
            (item for item in self.finite if item.identifier == identifier),
            None,
        )


def _coerce_entry(value: RationalFunction | dict[str, Any] | Any) -> RationalFunction:
    """接收公开有理函数对象、系数字典或精确常数。"""

    if isinstance(value, RationalFunction):
        return value
    if isinstance(value, dict):
        if "numerator" not in value:
            raise ValueError("rational-function dictionary requires a numerator field")
        return rational_function(value["numerator"], value.get("denominator", (1,)))
    return rational_function(value)


@dataclass(frozen=True)
class RationalMatrixSystem:
    """支持自动奇点发现的方形 exact Gaussian-rational 系数矩阵。"""

    entries: tuple[tuple[RationalFunction | dict[str, Any] | Any, ...], ...]
    variable_name: str = "s"
    name: str = "rational-matrix-system"

    def __post_init__(self) -> None:
        """统一转换所有矩阵元，并要求矩阵非空且为方阵。"""

        rows = tuple(tuple(_coerce_entry(value) for value in row) for row in self.entries)
        if not rows or any(len(row) != len(rows) for row in rows):
            raise ValueError("rational matrix must be nonempty and square")
        if not self.variable_name:
            raise ValueError("variable_name must not be empty")
        object.__setattr__(self, "entries", rows)

    @property
    def dimension(self) -> int:
        """返回方阵维数。"""

        return len(self.entries)

    def evaluate(self, point: acb) -> acb_mat:
        """在 Acb 点计算每个已约分矩阵元。"""

        return acb_mat([[entry.evaluate(point) for entry in row] for row in self.entries])

    def pole_order_at(self, point: Any) -> int:
        """返回完整矩阵在指定 exact 点的最大 pole 阶。"""

        return max(entry.pole_order_at(point) for row in self.entries for entry in row)

    def local_laurent_matrices_at(
        self,
        point: Any,
        minimum_power: int,
        maximum_power: int,
    ) -> dict[int, GaussianMatrix]:
        """按幂次返回完整矩阵的 exact 局部 Laurent 系数。"""

        coefficient_rows = [
            [
                entry.local_laurent_coefficients(point, minimum_power, maximum_power)
                for entry in row
            ]
            for row in self.entries
        ]
        return {
            power: GaussianMatrix.from_records(
                [
                    [
                        coefficient_rows[row][column][power - minimum_power]
                        for column in range(self.dimension)
                    ]
                    for row in range(self.dimension)
                ]
            )
            for power in range(minimum_power, maximum_power + 1)
        }

    def constant_basis_transform(self, transformation: GaussianMatrix) -> RationalMatrixSystem:
        """对 ``I=S J`` 施加 exact 常数换基并返回 ``S^-1 A S``。"""

        if transformation.nrows != self.dimension or transformation.ncols != self.dimension:
            raise ValueError("constant basis transformation dimension mismatch")
        inverse = transformation.inverse()
        rows: list[list[RationalFunction]] = []
        for row in range(self.dimension):
            output_row: list[RationalFunction] = []
            for column in range(self.dimension):
                value = RationalFunction((0,))
                for left in range(self.dimension):
                    for right in range(self.dimension):
                        value += (
                            self.entries[left][right]
                            * inverse.scalar(row, left)
                            * transformation.scalar(right, column)
                        )
                output_row.append(value)
            rows.append(output_row)
        return RationalMatrixSystem(
            tuple(tuple(row) for row in rows),
            variable_name=self.variable_name,
            name=f"{self.name}-constant-basis",
        )

    def inverted(self, variable_name: str = "sinv") -> RationalMatrixSystem:
        """对矩阵系统施加 ``s=1/sinv`` 和完整微分 Jacobian。"""

        return RationalMatrixSystem(
            tuple(
                tuple(entry.inverted_variable() for entry in row)
                for row in self.entries
            ),
            variable_name=variable_name,
            name=f"{self.name}-at-infinity",
        )

    def singularity_inventory(self, *, root_tolerance: float = 1.0e-30) -> SingularityInventory:
        """逐元约分、square-free 分解分母、隔离根并分类无穷远。"""

        if root_tolerance <= 0:
            raise ValueError("root_tolerance must be positive")
        factors: dict[
            tuple[str, ...], dict[str, Any]
        ] = {}
        for row_index, row in enumerate(self.entries):
            for column_index, entry in enumerate(row):
                numerator, denominator = entry.exact_polynomials()
                if numerator.is_zero or denominator.degree == 0:
                    continue
                for factor, exponent in denominator.squarefree_factors():
                    factor = factor.monic()
                    key = _polynomial_records(factor)
                    item = factors.setdefault(
                        key,
                        {"factor": factor, "pole_order": 0, "witnesses": set()},
                    )
                    item["pole_order"] = max(item["pole_order"], int(exponent))
                    item["witnesses"].add((row_index, column_index))

        unresolved: list[
            tuple[acb, str | None, tuple[str, ...], int, tuple[tuple[int, int], ...]]
        ] = []
        for key, item in factors.items():
            factor = item["factor"]
            # square-free 分组不等于不可约分解；若同一因子还含代数根，先单独抽出
            # 精确零根，避免把边界起点 x=0 连同其它根一起降为无 exact 身份的 Acb 根。
            if factor.coefficient(0).is_zero:
                zero_factor = GaussianPolynomial.from_coefficients([0, 1])
                factor, remainder = divmod(factor, zero_factor)
                if not remainder.is_zero:
                    raise ArithmeticError("exact zero-root factor did not divide polynomial")
                unresolved.append(
                    (
                        gaussian_rational(0).to_acb(),
                        "0",
                        _polynomial_records(zero_factor),
                        item["pole_order"],
                        tuple(sorted(item["witnesses"])),
                    )
                )
                if factor.degree == 0:
                    continue
            try:
                exact_root_records = factor.gaussian_rational_roots()
                roots = [root.to_acb() for root, _multiplicity in exact_root_records]
                exact_roots = [str(root) for root, _multiplicity in exact_root_records]
            except ValueError:
                roots = acb_poly(
                    [coefficient.to_acb() for coefficient in factor.coefficients()]
                ).roots(tol=root_tolerance)
                exact_roots = [None for _ in roots]
            for root, exact_root in zip(roots, exact_roots):
                unresolved.append(
                    (
                        root,
                        exact_root,
                        key,
                        item["pole_order"],
                        tuple(sorted(item["witnesses"])),
                    )
                )
        merged: list[dict[str, Any]] = []
        for root, exact_root, factor_records, pole_order, witnesses in unresolved:
            match = next(
                (
                    item
                    for item in merged
                    if (
                        exact_root is not None
                        and item["exact_root"] is not None
                        and exact_root == item["exact_root"]
                    )
                    or abs(root - item["root"]).contains(0)
                    or float(abs(root - item["root"]).mid()) <= 10.0 * root_tolerance
                ),
                None,
            )
            if match is None:
                merged.append(
                    {
                        "root": root,
                        "exact_root": exact_root,
                        "factor_records": factor_records,
                        "pole_order": pole_order,
                        "witnesses": set(witnesses),
                    }
                )
                continue
            match["witnesses"].update(witnesses)
            if pole_order > match["pole_order"]:
                match["pole_order"] = pole_order
                match["factor_records"] = factor_records
            if match["exact_root"] is None and exact_root is not None:
                match["exact_root"] = exact_root
                match["root"] = root
        merged.sort(
            key=lambda item: (
                float(item["root"].real.mid()),
                float(item["root"].imag.mid()),
            )
        )
        finite = tuple(
            SingularityRecord(
                identifier=f"finite_{index:03d}",
                location=item["root"],
                location_exact=item["exact_root"],
                pole_order=item["pole_order"],
                kind=(
                    "regular_singular"
                    if item["pole_order"] == 1
                    else "non_fuchsian_input_basis"
                ),
                factor_coefficients_exact=item["factor_records"],
                witness_entries=tuple(sorted(item["witnesses"])),
            )
            for index, item in enumerate(merged, start=1)
        )
        infinity_order = max(
            entry.pole_order_at_infinity() for row in self.entries for entry in row
        )
        infinity = SingularityRecord(
            identifier="infinity",
            location=None,
            location_exact="inf",
            pole_order=infinity_order,
            kind=(
                "ordinary"
                if infinity_order == 0
                else "regular_singular"
                if infinity_order == 1
                else "non_fuchsian_input_basis"
            ),
            factor_coefficients_exact=None,
            witness_entries=tuple(
                (row_index, column_index)
                for row_index, row in enumerate(self.entries)
                for column_index, entry in enumerate(row)
                if entry.pole_order_at_infinity() == infinity_order and infinity_order > 0
            ),
        )
        return SingularityInventory(self.name, self.variable_name, finite, infinity)
    def specialize_polynomial_simple_poles(
        self,
        inventory: SingularityInventory | None = None,
    ) -> tuple[Any | None, dict[str, Any]]:
        """自动认证并构造 P(x)+Sum R_j/(x-p_j) 公式快速系统。

        认证先使用 exact Q(i)[x] 约分后的奇点阶数；只要任一有限极点高于一阶，
        就返回通用路线。通过结构门禁后，多项式商用 exact 除法得到，留数用
        N(p)/D'(p) 的 Acb 球值构造，并在固定普通探针上验证候选矩阵包含原矩阵。
        """

        from .systems import PartialFractionSystem

        resolved = inventory or self.singularity_inventory()
        report: dict[str, Any] = {
            "schema": "flintnde_polynomial_simple_pole_specialization_v1",
            "system_name": self.name,
            "eligible": False,
            "finite_pole_count": len(resolved.finite),
        }
        high_poles = [
            record.identifier for record in resolved.finite if record.pole_order != 1
        ]
        if high_poles:
            report["reason"] = (
                "finite poles are not all simple: " + ", ".join(high_poles)
            )
            return None, report

        quotient_rows: list[list[GaussianPolynomial]] = []
        exact_residue_rows: list[list[list[GaussianRational]]] = [
            [
                [GaussianRational() for _record in resolved.finite]
                for _column in range(self.dimension)
            ]
            for _row in range(self.dimension)
        ]
        maximum_degree = 0
        for row_index, row in enumerate(self.entries):
            quotient_row: list[GaussianPolynomial] = []
            for column_index, entry in enumerate(row):
                numerator, denominator = entry.exact_polynomials()
                quotient, _remainder = divmod(numerator, denominator)
                quotient_row.append(quotient)
                maximum_degree = max(maximum_degree, max(0, quotient.degree))
                reconstruction = RationalFunction(quotient.records())
                for pole_index, record in enumerate(resolved.finite):
                    if record.location_exact is None:
                        report["reason"] = (
                            f"{record.identifier} is not an exact Q(i) root; "
                            "formula specialization cannot be certified exactly"
                        )
                        return None, report
                    exact_pole = gaussian_rational(record.location_exact)
                    exact_residue = entry.local_laurent_coefficients(
                        exact_pole, -1, -1
                    )[0]
                    exact_residue_rows[row_index][column_index][pole_index] = (
                        exact_residue
                    )
                    reconstruction += RationalFunction(
                        (exact_residue,), (-exact_pole, GaussianRational(1))
                    )
                if not (entry - reconstruction).is_zero:
                    report["reason"] = (
                        f"exact formula reconstruction failed for entry "
                        f"({row_index}, {column_index})"
                    )
                    return None, report
            quotient_rows.append(quotient_row)
        polynomial = tuple(
            acb_mat(
                [
                    [
                        quotient_rows[row][column].coefficient(degree).to_acb()
                        for column in range(self.dimension)
                    ]
                    for row in range(self.dimension)
                ]
            )
            for degree in range(maximum_degree + 1)
        )

        poles: list[acb] = []
        residues: list[acb_mat] = []
        for record in resolved.finite:
            if record.location is None:
                report["reason"] = f"{record.identifier} has no finite root enclosure"
                return None, report
            pole = acb(record.location)
            residue = acb_mat(self.dimension, self.dimension)
            for row_index, row in enumerate(self.entries):
                for column_index, _entry in enumerate(row):
                    residue[row_index, column_index] = (
                        exact_residue_rows[row_index][column_index][len(poles)].to_acb()
                    )
            poles.append(pole)
            residues.append(residue)

        candidate = PartialFractionSystem(
            constant=polynomial[0],
            polynomial_coefficients=polynomial[1:],
            residues=tuple(residues),
            poles=tuple(poles),
            name=f"{self.name}-polynomial-simple-poles",
        )
        probe_records: list[str] = []
        # 上述 Q(i)(x) 恒等式负责认证；Acb 探针只作为独立实现诊断。
        for probe in (acb(0), acb(1), acb(0, 1), acb(2, 1)):
            if any(abs(probe - pole).contains(0) for pole in poles):
                continue
            difference = candidate.evaluate(probe) - self.evaluate(probe)
            if any(
                not difference[row, column].contains(0)
                for row in range(self.dimension)
                for column in range(self.dimension)
            ):
                report["reason"] = (
                    f"formula reconstruction failed at probe {probe.str(30)}"
                )
                report["verification_probes"] = probe_records
                return None, report
            probe_records.append(probe.str(30))

        report.update(
            {
                "eligible": True,
                "reason": "all finite poles are simple and the formula reconstruction passed",
                "polynomial_degree": maximum_degree,
                "verification_probes": probe_records,
                "ordinary_solver": "polynomial_simple_pole_recurrence",
            }
        )
        return candidate, report


    def to_analytic_system(self, inventory: SingularityInventory | None = None):
        """通过已有普通点 evaluator API 暴露 exact rational 系统。"""

        from .systems import AnalyticMatrixSystem

        resolved = inventory or self.singularity_inventory()
        return AnalyticMatrixSystem(
            self.evaluate,
            self.dimension,
            tuple(record.location for record in resolved.finite if record.location is not None),
            self.name,
        )

    def regular_singular_system_at(
        self,
        point: Any,
        series_order: int,
        *,
        allow_ordinary: bool = False,
    ):
        """在 Q(i) 点提取 residue 和解析部分；可显式允许零 residue。"""

        from .frobenius import RegularSingularSystem

        if series_order <= 0:
            raise ValueError("series_order must be positive")
        center = _exact_scalar(point, "singular point")
        residue_records = [["0" for _ in range(self.dimension)] for _ in range(self.dimension)]
        regular_records = [
            [["0" for _ in range(self.dimension)] for _ in range(self.dimension)]
            for _ in range(series_order)
        ]
        maximum_pole_order = 0
        for row_index, row in enumerate(self.entries):
            for column_index, entry in enumerate(row):
                numerator, denominator = entry.exact_polynomials()
                shifted_numerator = _shift_polynomial(numerator, center)
                shifted_denominator = _shift_polynomial(denominator, center)
                pole_order = 0
                while (
                    pole_order <= shifted_denominator.degree
                    and shifted_denominator.coefficient(pole_order).is_zero
                ):
                    pole_order += 1
                numerator_zero_order = 0
                while (
                    numerator_zero_order <= shifted_numerator.degree
                    and shifted_numerator.coefficient(numerator_zero_order).is_zero
                ):
                    numerator_zero_order += 1
                if shifted_numerator.is_zero:
                    effective_pole_order = 0
                else:
                    effective_pole_order = max(0, pole_order - numerator_zero_order)
                maximum_pole_order = max(maximum_pole_order, effective_pole_order)
                if effective_pole_order > 1:
                    continue
                if effective_pole_order == 1:
                    numerator_local = GaussianPolynomial.from_coefficients(
                        shifted_numerator.coefficients()[numerator_zero_order:]
                    )
                    denominator_local = GaussianPolynomial.from_coefficients(
                        shifted_denominator.coefficients()[pole_order:]
                    )
                    quotient = _series_quotient(
                        numerator_local,
                        denominator_local,
                        series_order + 1,
                    )
                    residue_records[row_index][column_index] = str(quotient[0])
                    for degree in range(series_order):
                        regular_records[degree][row_index][column_index] = str(quotient[degree + 1])
                else:
                    quotient = _series_quotient(
                        shifted_numerator,
                        shifted_denominator,
                        series_order,
                    )
                    for degree in range(series_order):
                        regular_records[degree][row_index][column_index] = str(quotient[degree])
        if maximum_pole_order > 1 or (maximum_pole_order != 1 and not allow_ordinary):
            raise ValueError(
                f"{self.name}: point {center} is not a regular singular point "
                f"(matrix pole order {maximum_pole_order})"
            )
        return RegularSingularSystem(
            residue_records,
            tuple(regular_records),
            name=f"{self.name}-at-{center}",
        )


def analyze_singularities(
    system: RationalMatrixSystem,
    *,
    root_tolerance: float = 1.0e-30,
    output_layout: Any | None = None,
    filename: str = "singularity_inventory.json",
) -> SingularityInventory:
    """构造完整奇点清单，并可保存到调用脚本本地的奇点目录。"""

    inventory = system.singularity_inventory(root_tolerance=root_tolerance)
    if output_layout is not None:
        output_layout.write_json("singularities", filename, inventory.to_json())
    return inventory
