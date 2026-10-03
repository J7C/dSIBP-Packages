"""FlintNDE 的精确高斯有理数、多项式和原维矩阵后端。

本模块只依赖 python-flint。复数 ``a+b I`` 由两个 ``fmpq`` 保存；复多项式由
实部/虚部两张 ``fmpq_poly`` 保存；复矩阵由实部/虚部两张同维 ``fmpq_mat`` 保存。
所有结构判定均保持原矩阵维数，不把 N 维复系统实块化为 2N 维系统。
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Iterable

from flint import acb, acb_mat, arb, fmpq, fmpq_mat, fmpq_poly


def _real_rational(value: Any, field_name: str) -> fmpq:
    """把一个无虚部精确标量转换为 ``fmpq``，并拒绝机器数。"""

    if isinstance(value, fmpq):
        return value
    if isinstance(value, (float, complex, acb, arb)):
        raise TypeError(f"{field_name} received floating input")
    if isinstance(value, Fraction):
        return fmpq(value.numerator, value.denominator)
    text = re.sub(r"`(?:\d+(?:\.\d*)?|\.\d+)?", "", str(value)).replace("*^", "e")
    try:
        rational = Fraction(text)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError(f"expected an exact rational scalar, got {value}") from error
    return fmpq(rational.numerator, rational.denominator)


@dataclass(frozen=True)
class GaussianRational:
    """保存一个精确高斯有理数 ``real + imag I``。"""

    real: fmpq = fmpq(0)
    imag: fmpq = fmpq(0)

    def __init__(self, real: Any = 0, imag: Any = 0) -> None:
        object.__setattr__(self, "real", _real_rational(real, "real part"))
        object.__setattr__(self, "imag", _real_rational(imag, "imaginary part"))

    @property
    def is_zero(self) -> bool:
        """返回实部和虚部是否严格为零。"""

        return self.real == 0 and self.imag == 0

    def to_acb(self) -> acb:
        """无损送入当前工作精度的 Acb ball。"""

        return acb(self.real, self.imag)

    def conjugate(self) -> GaussianRational:
        """返回复共轭。"""

        return GaussianRational(self.real, -self.imag)

    def sqrt_exact(self) -> GaussianRational:
        """若平方根仍在 ``Q(i)`` 中则返回，否则 fail closed。"""

        try:
            if self.imag == 0:
                if self.real >= 0:
                    return GaussianRational(self.real.sqrt(), 0)
                return GaussianRational(0, (-self.real).sqrt())
            norm = (self.real * self.real + self.imag * self.imag).sqrt()
            real_part = ((norm + self.real) / 2).sqrt()
            if real_part == 0:
                raise ValueError("nonreal Gaussian rational has zero real square-root part")
            imag_part = self.imag / (2 * real_part)
            candidate = GaussianRational(real_part, imag_part)
        except Exception as error:
            raise ValueError(f"{self} has no square root in Q(i)") from error
        if candidate * candidate != self:
            raise ValueError(f"{self} has no square root in Q(i)")
        return candidate

    def __add__(self, other: Any) -> GaussianRational:
        right = gaussian_rational(other)
        return GaussianRational(self.real + right.real, self.imag + right.imag)

    __radd__ = __add__

    def __neg__(self) -> GaussianRational:
        return GaussianRational(-self.real, -self.imag)

    def __sub__(self, other: Any) -> GaussianRational:
        return self + (-gaussian_rational(other))

    def __rsub__(self, other: Any) -> GaussianRational:
        return gaussian_rational(other) - self

    def __mul__(self, other: Any) -> GaussianRational:
        right = gaussian_rational(other)
        return GaussianRational(
            self.real * right.real - self.imag * right.imag,
            self.real * right.imag + self.imag * right.real,
        )

    __rmul__ = __mul__

    def __truediv__(self, other: Any) -> GaussianRational:
        right = gaussian_rational(other)
        norm = right.real * right.real + right.imag * right.imag
        if norm == 0:
            raise ZeroDivisionError("division by zero Gaussian rational")
        return GaussianRational(
            (self.real * right.real + self.imag * right.imag) / norm,
            (self.imag * right.real - self.real * right.imag) / norm,
        )

    def __rtruediv__(self, other: Any) -> GaussianRational:
        return gaussian_rational(other) / self

    def __pow__(self, power: int) -> GaussianRational:
        if not isinstance(power, int):
            return NotImplemented
        if power < 0:
            return GaussianRational(1) / (self ** (-power))
        result = GaussianRational(1)
        base = self
        exponent = power
        while exponent:
            if exponent & 1:
                result *= base
            base *= base
            exponent //= 2
        return result

    def __str__(self) -> str:
        if self.imag == 0:
            return str(self.real)
        if self.real == 0:
            return f"{self.imag}*I"
        sign = "+" if self.imag > 0 else "-"
        return f"{self.real} {sign} {abs(self.imag)}*I"


def gaussian_rational(value: Any = 0, imag: Any | None = None) -> GaussianRational:
    """解析精确高斯有理数；支持对象、二元组、字典和 ``a+b*I`` 字符串。"""

    if imag is not None:
        return GaussianRational(value, imag)
    if isinstance(value, GaussianRational):
        return value
    if isinstance(value, dict):
        unknown = set(value) - {"real", "imag"}
        if unknown:
            raise ValueError(f"unknown Gaussian-rational fields: {sorted(unknown)}")
        return GaussianRational(value.get("real", 0), value.get("imag", 0))
    if isinstance(value, tuple) and len(value) == 2:
        return GaussianRational(value[0], value[1])
    if isinstance(value, (float, complex, acb, arb)):
        raise TypeError("exact Gaussian-rational input cannot be floating")
    text = str(value).strip()
    if "I" not in text and "j" not in text:
        return GaussianRational(_real_rational(value, "exact scalar"), 0)
    compact = text.replace(" ", "").replace("(", "").replace(")", "").replace("j", "I")
    pieces = re.findall(r"[+-]?[^+-]+", compact)
    if not pieces or "".join(pieces) != compact:
        raise ValueError(f"unsupported exact Gaussian-rational syntax: {value}")
    real = fmpq(0)
    imaginary = fmpq(0)
    for piece in pieces:
        if "I" in piece:
            coefficient = piece.replace("*I", "").replace("I", "")
            if coefficient in {"", "+"}:
                coefficient = "1"
            elif coefficient == "-":
                coefficient = "-1"
            elif coefficient.startswith("/"):
                coefficient = "1" + coefficient
            elif coefficient.startswith("-/"):
                coefficient = "-1" + coefficient[1:]
            imaginary += _real_rational(coefficient, "imaginary coefficient")
        else:
            real += _real_rational(piece, "real coefficient")
    return GaussianRational(real, imaginary)


@dataclass(frozen=True)
class GaussianPolynomial:
    """保存一个 ``Q(i)[x]`` 多项式。"""

    real: fmpq_poly
    imag: fmpq_poly

    @classmethod
    def from_coefficients(cls, coefficients: Iterable[Any]) -> GaussianPolynomial:
        records = tuple(gaussian_rational(value) for value in coefficients)
        if not records:
            raise ValueError("polynomial must contain at least one coefficient")
        return cls(
            fmpq_poly([value.real for value in records]),
            fmpq_poly([value.imag for value in records]),
        )

    @classmethod
    def from_real_poly(cls, polynomial: fmpq_poly) -> GaussianPolynomial:
        return cls(polynomial, fmpq_poly([0]))

    @property
    def degree(self) -> int:
        return max(self.real.degree(), self.imag.degree())

    @property
    def is_zero(self) -> bool:
        return self.real == 0 and self.imag == 0

    @property
    def is_one(self) -> bool:
        return self.real == 1 and self.imag == 0

    def coefficient(self, degree: int) -> GaussianRational:
        return GaussianRational(self.real[degree], self.imag[degree])

    def coefficients(self) -> list[GaussianRational]:
        if self.is_zero:
            return [GaussianRational()]
        return [self.coefficient(index) for index in range(self.degree + 1)]

    def records(self) -> tuple[str, ...]:
        return tuple(str(value) for value in self.coefficients())

    def monic(self) -> GaussianPolynomial:
        if self.is_zero:
            return self
        return self / self.coefficient(self.degree)

    def derivative(self) -> GaussianPolynomial:
        if self.degree <= 0:
            return GaussianPolynomial.from_coefficients([0])
        return GaussianPolynomial.from_coefficients(
            [self.coefficient(degree) * degree for degree in range(1, self.degree + 1)]
        )

    def shifted(self, center: Any) -> GaussianPolynomial:
        point = gaussian_rational(center)
        shifted = [GaussianRational() for _ in range(max(1, self.degree + 1))]
        for source_degree, coefficient in enumerate(self.coefficients()):
            for local_degree in range(source_degree + 1):
                shifted[local_degree] += (
                    coefficient
                    * math.comb(source_degree, local_degree)
                    * point ** (source_degree - local_degree)
                )
        return GaussianPolynomial.from_coefficients(shifted)

    def evaluate(self, point: acb) -> acb:
        value = acb(0)
        for coefficient in reversed(self.coefficients()):
            value = value * point + coefficient.to_acb()
        return value

    def __add__(self, other: Any) -> GaussianPolynomial:
        right = as_gaussian_polynomial(other)
        return GaussianPolynomial(self.real + right.real, self.imag + right.imag)

    __radd__ = __add__

    def __neg__(self) -> GaussianPolynomial:
        return GaussianPolynomial(-self.real, -self.imag)

    def __sub__(self, other: Any) -> GaussianPolynomial:
        return self + (-as_gaussian_polynomial(other))

    def __mul__(self, other: Any) -> GaussianPolynomial:
        if isinstance(other, GaussianPolynomial):
            return GaussianPolynomial(
                self.real * other.real - self.imag * other.imag,
                self.real * other.imag + self.imag * other.real,
            )
        scalar = gaussian_rational(other)
        return GaussianPolynomial(
            self.real * scalar.real - self.imag * scalar.imag,
            self.real * scalar.imag + self.imag * scalar.real,
        )

    __rmul__ = __mul__

    def __truediv__(self, scalar: Any) -> GaussianPolynomial:
        value = gaussian_rational(scalar)
        inverse = GaussianRational(1) / value
        return self * inverse

    def __divmod__(self, divisor: GaussianPolynomial) -> tuple[GaussianPolynomial, GaussianPolynomial]:
        if divisor.is_zero:
            raise ZeroDivisionError("polynomial division by zero")
        if self.is_zero or self.degree < divisor.degree:
            return GaussianPolynomial.from_coefficients([0]), self
        remainder = self.coefficients()
        quotient = [GaussianRational() for _ in range(self.degree - divisor.degree + 1)]
        divisor_coefficients = divisor.coefficients()
        while len(remainder) - 1 >= divisor.degree and any(not value.is_zero for value in remainder):
            while len(remainder) > 1 and remainder[-1].is_zero:
                remainder.pop()
            if len(remainder) - 1 < divisor.degree:
                break
            shift = len(remainder) - 1 - divisor.degree
            factor = remainder[-1] / divisor_coefficients[-1]
            quotient[shift] += factor
            for degree, coefficient in enumerate(divisor_coefficients):
                remainder[shift + degree] -= factor * coefficient
        return (
            GaussianPolynomial.from_coefficients(quotient),
            GaussianPolynomial.from_coefficients(remainder),
        )

    def gcd(self, other: GaussianPolynomial) -> GaussianPolynomial:
        left, right = self, other
        while not right.is_zero:
            _quotient, remainder = divmod(left, right)
            left, right = right, remainder
        return left.monic() if not left.is_zero else GaussianPolynomial.from_coefficients([0])

    def squarefree_factors(self) -> list[tuple[GaussianPolynomial, int]]:
        """返回 characteristic-zero square-free factors 及其重数。"""

        if self.degree <= 0:
            return []
        polynomial = self.monic()
        repeated = polynomial.gcd(polynomial.derivative())
        squarefree, remainder = divmod(polynomial, repeated)
        if not remainder.is_zero:
            raise ArithmeticError("square-free decomposition division failed")
        factors: list[tuple[GaussianPolynomial, int]] = []
        multiplicity = 1
        while not squarefree.is_one:
            overlap = squarefree.gcd(repeated)
            factor, remainder = divmod(squarefree, overlap)
            if not remainder.is_zero:
                raise ArithmeticError("square-free factor division failed")
            if not factor.is_one:
                factors.append((factor.monic(), multiplicity))
            squarefree = overlap
            repeated, remainder = divmod(repeated, overlap)
            if not remainder.is_zero:
                raise ArithmeticError("repeated-factor division failed")
            multiplicity += 1
        if not repeated.is_one:
            raise ArithmeticError("square-free decomposition did not terminate")
        return factors

    def gaussian_rational_roots(self) -> list[tuple[GaussianRational, int]]:
        """精确分解完全落在 ``Q(i)`` 中的根，否则 fail closed。"""

        roots: list[tuple[GaussianRational, int]] = []
        for squarefree, multiplicity in self.squarefree_factors():
            remaining = squarefree
            norm = squarefree.real * squarefree.real + squarefree.imag * squarefree.imag
            _content, rational_factors = norm.factor()
            for rational_factor, _exponent in rational_factors:
                lifted = GaussianPolynomial.from_real_poly(rational_factor)
                factor = remaining.gcd(lifted)
                if factor.degree <= 0:
                    continue
                local_roots = factor._linear_roots()
                roots.extend((root, multiplicity) for root in local_roots)
                remaining, remainder = divmod(remaining, factor)
                if not remainder.is_zero:
                    raise ArithmeticError("Gaussian root factor did not divide exactly")
            if not remaining.is_one:
                raise ValueError("indicial polynomial does not split over Q(i)")
        if sum(multiplicity for _root, multiplicity in roots) != self.degree:
            raise ValueError("indicial polynomial does not split completely over Q(i)")
        return sorted(roots, key=lambda item: (float(item[0].real), float(item[0].imag)))

    def _linear_roots(self) -> list[GaussianRational]:
        if self.degree == 1:
            return [-self.coefficient(0) / self.coefficient(1)]
        if self.degree == 2:
            constant = self.coefficient(0)
            linear = self.coefficient(1)
            quadratic = self.coefficient(2)
            discriminant = linear * linear - 4 * quadratic * constant
            square_root = discriminant.sqrt_exact()
            roots = [
                (-linear - square_root) / (2 * quadratic),
                (-linear + square_root) / (2 * quadratic),
            ]
            product = GaussianPolynomial.from_coefficients([1])
            for root in roots:
                product *= GaussianPolynomial.from_coefficients([-root, 1])
            if product.monic().records() != self.monic().records():
                raise ValueError("quadratic factor did not split exactly over Q(i)")
            return roots
        raise ValueError("polynomial factor has non-Gaussian-rational roots")


def as_gaussian_polynomial(value: Any) -> GaussianPolynomial:
    if isinstance(value, GaussianPolynomial):
        return value
    return GaussianPolynomial.from_coefficients([value])


def series_quotient(
    numerator: GaussianPolynomial, denominator: GaussianPolynomial, order: int
) -> list[GaussianRational]:
    """递推常数项非零的 ``Q(i)`` Taylor 商。"""

    if order < 0:
        raise ValueError("series order must be nonnegative")
    denominator_zero = denominator.coefficient(0)
    if denominator_zero.is_zero:
        raise ZeroDivisionError("Taylor quotient denominator vanishes")
    result: list[GaussianRational] = []
    for degree in range(order):
        right = numerator.coefficient(degree)
        for denominator_degree in range(1, degree + 1):
            right -= denominator.coefficient(denominator_degree) * result[
                degree - denominator_degree
            ]
        result.append(right / denominator_zero)
    return result


@dataclass(frozen=True)
class GaussianMatrix:
    """用两张同维 ``fmpq_mat`` 保存原维精确复矩阵。"""

    real: fmpq_mat
    imag: fmpq_mat

    @classmethod
    def from_records(cls, records: list[list[Any]]) -> GaussianMatrix:
        values = [[gaussian_rational(value) for value in row] for row in records]
        if not values or any(len(row) != len(values[0]) for row in values):
            raise ValueError("exact Gaussian matrix must be nonempty and rectangular")
        return cls(
            fmpq_mat([[value.real for value in row] for row in values]),
            fmpq_mat([[value.imag for value in row] for row in values]),
        )

    @classmethod
    def zero(cls, nrows: int, ncols: int) -> GaussianMatrix:
        return cls(fmpq_mat(nrows, ncols), fmpq_mat(nrows, ncols))

    @classmethod
    def identity(cls, dimension: int) -> GaussianMatrix:
        return cls(
            fmpq_mat(
                [[1 if row == column else 0 for column in range(dimension)] for row in range(dimension)]
            ),
            fmpq_mat(dimension, dimension),
        )

    @property
    def nrows(self) -> int:
        return self.real.nrows()

    @property
    def ncols(self) -> int:
        return self.real.ncols()

    @property
    def is_zero(self) -> bool:
        return all(value == 0 for value in self.real.entries()) and all(
            value == 0 for value in self.imag.entries()
        )

    def scalar(self, row: int, column: int) -> GaussianRational:
        return GaussianRational(self.real[row, column], self.imag[row, column])

    def to_records(self) -> list[list[str]]:
        return [
            [str(self.scalar(row, column)) for column in range(self.ncols)]
            for row in range(self.nrows)
        ]

    def to_acb(self) -> acb_mat:
        return acb_mat(
            [[self.scalar(row, column).to_acb() for column in range(self.ncols)] for row in range(self.nrows)]
        )

    def column(self, column: int) -> GaussianMatrix:
        return GaussianMatrix(
            fmpq_mat([[self.real[row, column]] for row in range(self.nrows)]),
            fmpq_mat([[self.imag[row, column]] for row in range(self.nrows)]),
        )

    def __add__(self, other: GaussianMatrix) -> GaussianMatrix:
        return GaussianMatrix(self.real + other.real, self.imag + other.imag)

    def __neg__(self) -> GaussianMatrix:
        return GaussianMatrix(-self.real, -self.imag)

    def __sub__(self, other: GaussianMatrix) -> GaussianMatrix:
        return self + (-other)

    def __mul__(self, other: Any) -> GaussianMatrix:
        if isinstance(other, GaussianMatrix):
            return GaussianMatrix(
                self.real * other.real - self.imag * other.imag,
                self.real * other.imag + self.imag * other.real,
            )
        scalar = gaussian_rational(other)
        return GaussianMatrix(
            self.real * scalar.real - self.imag * scalar.imag,
            self.real * scalar.imag + self.imag * scalar.real,
        )

    def __rmul__(self, other: Any) -> GaussianMatrix:
        return self * other

    def __truediv__(self, scalar: Any) -> GaussianMatrix:
        return self * (GaussianRational(1) / gaussian_rational(scalar))

    def __pow__(self, power: int) -> GaussianMatrix:
        if power < 0 or self.nrows != self.ncols:
            raise ValueError("Gaussian matrix power requires a square matrix and nonnegative exponent")
        result = GaussianMatrix.identity(self.nrows)
        base = self
        exponent = power
        while exponent:
            if exponent & 1:
                result = result * base
            base = base * base
            exponent //= 2
        return result

    def trace(self) -> GaussianRational:
        return sum((self.scalar(index, index) for index in range(self.nrows)), GaussianRational())

    def charpoly(self) -> GaussianPolynomial:
        """用 Faddeev-LeVerrier 构造 exact characteristic polynomial。"""

        if self.nrows != self.ncols:
            raise ValueError("characteristic polynomial requires a square matrix")
        identity = GaussianMatrix.identity(self.nrows)
        auxiliary = identity
        coefficients = [GaussianRational(1)]
        for degree in range(1, self.nrows + 1):
            product = self * auxiliary
            coefficient = -product.trace() / degree
            coefficients.append(coefficient)
            auxiliary = product + identity * coefficient
        if not auxiliary.is_zero:
            raise ArithmeticError("exact characteristic recurrence failed Cayley-Hamilton check")
        return GaussianPolynomial.from_coefficients(reversed(coefficients))

    def rref(self) -> tuple[GaussianMatrix, int, list[int]]:
        """用 Q(i) Gauss-Jordan 消元返回 RREF、rank 和 pivot 列。"""

        values = [
            [self.scalar(row, column) for column in range(self.ncols)]
            for row in range(self.nrows)
        ]
        pivot_columns: list[int] = []
        pivot_row = 0
        for column in range(self.ncols):
            pivot = next(
                (row for row in range(pivot_row, self.nrows) if not values[row][column].is_zero),
                None,
            )
            if pivot is None:
                continue
            values[pivot_row], values[pivot] = values[pivot], values[pivot_row]
            pivot_value = values[pivot_row][column]
            values[pivot_row] = [value / pivot_value for value in values[pivot_row]]
            for row in range(self.nrows):
                if row == pivot_row or values[row][column].is_zero:
                    continue
                factor = values[row][column]
                values[row] = [
                    value - factor * pivot_entry
                    for value, pivot_entry in zip(values[row], values[pivot_row])
                ]
            pivot_columns.append(column)
            pivot_row += 1
            if pivot_row == self.nrows:
                break
        return GaussianMatrix.from_records(values), pivot_row, pivot_columns

    def nullspace(self) -> list[GaussianMatrix]:
        reduced, _rank, pivot_columns = self.rref()
        free_columns = [column for column in range(self.ncols) if column not in pivot_columns]
        basis: list[GaussianMatrix] = []
        for free_column in free_columns:
            values = [GaussianRational() for _ in range(self.ncols)]
            values[free_column] = GaussianRational(1)
            for row, column in enumerate(pivot_columns):
                values[column] = -reduced.scalar(row, free_column)
            basis.append(GaussianMatrix.from_records([[value] for value in values]))
        return basis

    def inverse(self) -> GaussianMatrix:
        """用 exact Q(i) Gauss-Jordan 消元返回方阵逆矩阵。"""

        if self.nrows != self.ncols:
            raise ValueError("exact matrix inverse requires a square matrix")
        dimension = self.nrows
        augmented = GaussianMatrix.from_records(
            [
                [self.scalar(row, column) for column in range(dimension)]
                + [1 if row == column else 0 for column in range(dimension)]
                for row in range(dimension)
            ]
        )
        reduced, rank, pivots = augmented.rref()
        if rank < dimension or pivots[:dimension] != list(range(dimension)):
            raise ZeroDivisionError("exact Gaussian matrix is singular")
        inverse = GaussianMatrix.from_records(
            [
                [reduced.scalar(row, dimension + column) for column in range(dimension)]
                for row in range(dimension)
            ]
        )
        if not (self * inverse - GaussianMatrix.identity(dimension)).is_zero:
            raise ArithmeticError("exact Gaussian matrix inverse failed substitution check")
        return inverse


def exact_column(values: Iterable[Any]) -> GaussianMatrix:
    """把 Q(i) 标量序列转换为 exact 列矩阵。"""

    return GaussianMatrix.from_records([[value] for value in values])


def exact_matrix_from_columns(columns: list[GaussianMatrix], dimension: int) -> GaussianMatrix:
    """按给定顺序把 exact 列向量拼成矩阵。

    ``dimension`` 指定行数；列数取 ``len(columns)``，允许少于行数的非方阵
    （用于本征子空间投影等场景）。需要完整方阵的调用方自行校验列数。
    """

    return GaussianMatrix.from_records(
        [
            [columns[column].scalar(row, 0) for column in range(len(columns))]
            for row in range(dimension)
        ]
    )


def canonical_exact_solve(
    left: GaussianMatrix,
    right: GaussianMatrix,
    *,
    incompatible_error: type[Exception] | None = None,
    incompatible_message: str = "exact linear system is incompatible",
    solve_error: type[Exception] = ArithmeticError,
    solve_message: str = "exact linear system failed substitution check",
    dimension_message: str = "exact linear solve dimension mismatch",
) -> GaussianMatrix:
    """对相容的奇异 exact 系统取自由变量为零的规范特解。

    ``incompatible_error``/``solve_error`` 允许调用方注入自己层次的异常类型；
    RREF 和代入验证逻辑为所有调用方共享，避免多份实现漂移。
    """

    if left.nrows != right.nrows or right.ncols != 1:
        raise ValueError(dimension_message)
    augmented = GaussianMatrix.from_records(
        [
            [left.scalar(row, column) for column in range(left.ncols)]
            + [right.scalar(row, 0)]
            for row in range(left.nrows)
        ]
    )
    reduced, _rank, pivots = augmented.rref()
    if left.ncols in pivots:
        error_type = incompatible_error or ValueError
        raise error_type(incompatible_message)
    solution = [GaussianRational() for _ in range(left.ncols)]
    for pivot_row, pivot_column in enumerate(pivots):
        if pivot_column < left.ncols:
            solution[pivot_column] = reduced.scalar(pivot_row, left.ncols)
    candidate = exact_column(solution)
    if not (left * candidate - right).is_zero:
        raise solve_error(solve_message)
    return candidate
