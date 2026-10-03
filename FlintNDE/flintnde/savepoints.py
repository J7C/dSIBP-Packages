"""把显式标记的路径点结果即时写入调用目录。

逐点文件是中断可恢复的计算记录；汇总文件只在整条输运链完成后写入。默认输出根固定为
调用 Python 进程的当前工作目录，FlintNDE 源码和安装目录不属于输出位置。
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from flint import acb, acb_mat

from .exact_gaussian import GaussianMatrix, gaussian_rational
from .output_layout import validate_single_path_name


def acb_record(value: acb, digits: int) -> dict[str, str]:
    """保留 Acb 中点和误差球，不经过 binary64。"""

    return {
        "real": value.real.mid().str(digits, radius=False, more=True),
        "imag": value.imag.mid().str(digits, radius=False, more=True),
        "realRadius": value.real.rad().str(digits),
        "imagRadius": value.imag.rad().str(digits),
    }


def acb_exact_input_record(value: acb, digits: int) -> dict[str, str]:
    """返回可直接作为下一次 ``C`` 分量使用的十进制 Q(i) 中点记录。"""

    return {
        "real": value.real.mid().str(digits, radius=False, more=True),
        "imag": value.imag.mid().str(digits, radius=False, more=True),
    }


def vector_record(vector: acb_mat, digits: int) -> list[dict[str, str]]:
    """序列化列向量并验证形状。"""

    if vector.ncols() != 1:
        raise ValueError("saved ordinary result must be a column vector")
    return [acb_record(vector[row, 0], digits) for row in range(vector.nrows())]


def _acb_column_zero(vector: acb_mat) -> bool:
    """按 Acb 包含零语义判断数值列是否可认证为零。"""

    return all(vector[row, 0].contains(0) for row in range(vector.nrows()))


def _canonical_range_preimage(matrix: GaussianMatrix, right: acb_mat) -> acb_mat:
    """对 ``matrix.x=right`` 使用与 exact 边界解析相同的自由变量为零特解。"""

    _reduced, _rank, pivot_columns = matrix.rref()
    if not pivot_columns:
        if _acb_column_zero(right):
            return acb_mat(matrix.ncols, 1)
        raise ValueError("numeric Frobenius coefficient is outside the nilpotent range")
    basis = GaussianMatrix.from_records(
        [
            [matrix.scalar(row, column) for column in pivot_columns]
            for row in range(matrix.nrows)
        ]
    )
    transposed = GaussianMatrix.from_records(
        [
            [basis.scalar(row, column) for row in range(basis.nrows)]
            for column in range(basis.ncols)
        ]
    )
    _row_reduced, _row_rank, pivot_rows = transposed.rref()
    selected_rows = tuple(pivot_rows[: len(pivot_columns)])
    square = GaussianMatrix.from_records(
        [
            [basis.scalar(row, column) for column in range(basis.ncols)]
            for row in selected_rows
        ]
    )
    selected_right = acb_mat([[right[row, 0]] for row in selected_rows])
    coefficients = square.inverse().to_acb() * selected_right
    result = acb_mat(matrix.ncols, 1)
    for local_index, column in enumerate(pivot_columns):
        result[column, 0] = coefficients[local_index, 0]
    return result


def _frobenius_terms_from_manifest(
    manifest: dict[str, Any],
    dimension: int,
    constants: acb_mat,
    digits: int,
) -> list[dict[str, Any]]:
    """从一个未作 meromorphic gauge 的 exact Frobenius manifest 恢复领头项。"""

    route = manifest.get("route")
    terms: list[dict[str, Any]] = []
    if route == "diagonalizable_roots_exact_gate":
        roots = [gaussian_rational(value) for value in manifest["solution_roots_exact"]]
        initial = [
            acb_mat([[gaussian_rational(value).to_acb()] for value in vector])
            for vector in manifest["initial_vectors_exact"]
        ]
        for root in dict.fromkeys(roots):
            coefficient = acb_mat(dimension, 1)
            for index, candidate_root in enumerate(roots):
                if candidate_root == root:
                    coefficient += initial[index] * constants[index, 0]
            if not _acb_column_zero(coefficient):
                terms.append(
                    {
                        "a": str(root),
                        "b": 0,
                        "C": [
                            acb_exact_input_record(coefficient[row, 0], digits)
                            for row in range(dimension)
                        ],
                        "CBalls": [
                            acb_record(coefficient[row, 0], digits)
                            for row in range(dimension)
                        ],
                    }
                )
        return terms
    if route != "single_root_jordan_exact_gate":
        raise NotImplementedError(f"unsupported saved Frobenius route: {route}")

    nilpotent = GaussianMatrix.from_records(manifest["nilpotent_exact"])
    remaining = acb_mat(constants)
    for log_degree in range(int(manifest["maximum_log_degree"]), -1, -1):
        power = nilpotent ** log_degree
        coefficient = power.to_acb() * remaining / acb(math.factorial(log_degree))
        if _acb_column_zero(coefficient):
            continue
        terms.append(
            {
                "a": str(gaussian_rational(manifest["root_exact"])),
                "b": log_degree,
                "C": [
                    acb_exact_input_record(coefficient[row, 0], digits)
                    for row in range(dimension)
                ],
                "CBalls": [
                    acb_record(coefficient[row, 0], digits)
                    for row in range(dimension)
                ],
            }
        )
        preimage = _canonical_range_preimage(
            power,
            coefficient * acb(math.factorial(log_degree)),
        )
        remaining -= preimage
    if not _acb_column_zero(remaining):
        raise ArithmeticError("saved Jordan boundary decomposition left nonzero constants")
    return terms


def frobenius_terms_from_constants(
    local_basis: Any,
    constants: acb_mat,
    digits: int,
) -> list[dict[str, Any]]:
    """把 direct regular-singular 常数恢复为可再次输入的 ``{a,b,C}`` 项。"""

    if local_basis.method != "regular_singular_power_log":
        raise NotImplementedError(
            "saved Frobenius output requires a direct regular-singular power-log basis"
        )
    return _frobenius_terms_from_manifest(
        local_basis.manifest["frobenius"], local_basis.dimension, constants, digits
    )


def _exponential_term(
    phi: list[dict[str, Any]],
    power_log_term: dict[str, Any],
    lifted: acb_mat,
    digits: int,
) -> dict[str, Any]:
    """把 transformed-sector 领头项换回原积分基并附上 exact 指数。"""

    return {
        "phi": phi,
        "a": power_log_term["a"],
        "b": power_log_term["b"],
        "C": [acb_exact_input_record(lifted[row, 0], digits) for row in range(lifted.nrows())],
        "CBalls": [acb_record(lifted[row, 0], digits) for row in range(lifted.nrows())],
    }


def _exponential_power_log_terms(
    local_basis: Any,
    constants: acb_mat,
    digits: int,
) -> list[dict[str, Any]]:
    """恢复严格解耦指数 sectors 的 ``{phi,a,b,C}`` 原基领头数据。"""

    manifest = local_basis.manifest
    constant_basis = GaussianMatrix.from_records(manifest["constant_basis_exact"]).to_acb()
    terms: list[dict[str, Any]] = []
    for sector in manifest["sectors"]:
        indices = tuple(int(index) for index in sector["component_indices_zero_based"])
        block_constants = acb_mat([[constants[index, 0]] for index in indices])
        if _acb_column_zero(block_constants):
            continue
        block_terms = _frobenius_terms_from_manifest(
            sector["frobenius"], len(indices), block_constants, digits
        )
        phi = [
            {
                "power": int(record["phi_power"]),
                "coefficient": record["phi_coefficient"],
            }
            for record in sector["exponential"]
        ]
        for block_term in block_terms:
            transformed = acb_mat(local_basis.dimension, 1)
            for local_index, global_index in enumerate(indices):
                value = block_term["C"][local_index]
                transformed[global_index, 0] = acb(value["real"], value["imag"])
            terms.append(
                _exponential_term(
                    phi,
                    block_term,
                    constant_basis * transformed,
                    digits,
                )
            )
    if not terms:
        raise ValueError("saved exponential boundary selects the zero solution")
    return terms


def _formal_exponential_terms(
    local_basis: Any,
    constants: acb_mat,
    digits: int,
) -> list[dict[str, Any]]:
    """恢复 start-only simple-spectrum 形式渐近分支的领头边界。"""

    terms: list[dict[str, Any]] = []
    for index, branch in enumerate(local_basis.manifest["branches"]):
        scalar = constants[index, 0]
        if scalar.contains(0):
            continue
        k = gaussian_rational(branch["k_exact"])
        leading = GaussianMatrix.from_records(branch["leading_vector_exact"]).to_acb() * scalar
        terms.append(
            _exponential_term(
                [] if k.is_zero else [{"power": -1, "coefficient": str(-k)}],
                {"a": branch["rho_exact"], "b": 0},
                leading,
                digits,
            )
        )
    if not terms:
        raise ValueError("saved formal exponential boundary selects the zero solution")
    return terms


def singular_boundary_record_from_constants(
    local_basis: Any,
    constants: acb_mat,
    digits: int,
) -> tuple[str, dict[str, Any]]:
    """按已认证局部基类型生成可复用的奇点边界记录。

    Lee--Moser 变换后的原基仍需逆局部 jet，未实现时直接拒绝；start-only formal
    分支可被保存，但调用方不得把它用于中间或终点 continuation。
    """

    if local_basis.method == "regular_singular_power_log":
        return "frobenius_boundary", {
            "schema": "flintnde_frobenius_boundary_v1",
            "terms": frobenius_terms_from_constants(local_basis, constants, digits),
            "canonicalBasisConstants": vector_record(constants, digits),
            "localBasisMethod": local_basis.method,
        }
    if local_basis.method == "exponential_power_log":
        terms = _exponential_power_log_terms(local_basis, constants, digits)
    elif local_basis.method == "formal_exponential_asymptotic":
        terms = _formal_exponential_terms(local_basis, constants, digits)
    else:
        raise NotImplementedError(
            f"saved singular output does not support local basis method {local_basis.method}"
        )
    return "exponential_boundary", {
        "schema": "flintnde_exponential_boundary_v1",
        "localVariable": "z=working_coordinate-singularity_center",
        "branchConvention": "python-flint principal acb log/power branch",
        "terms": terms,
        "canonicalBasisConstants": vector_record(constants, digits),
        "localBasisMethod": local_basis.method,
        "continuationReady": bool(local_basis.continuation_ready),
    }


@dataclass
class SavePointWriter:
    """管理一次输运的即时逐点文件和最终有序汇总。"""

    output_directory: Path
    summary_filename: str = "flintnde_save_points.json"
    digits: int = 40
    records: list[dict[str, Any]] = field(default_factory=list)

    def __post_init__(self) -> None:
        """冻结输出根并拒绝跨目录汇总文件名。"""

        self.output_directory = Path(self.output_directory).resolve()
        validate_single_path_name(self.summary_filename, "save summary filename")
        self.output_directory.mkdir(parents=True, exist_ok=True)

    def write(self, request: dict[str, Any], payload: dict[str, Any]) -> Path:
        """完成一个标记点后立即写文件，并把同一记录加入内存汇总。"""

        sequence = len(self.records) + 1
        record = {
            "schema": "flintnde_saved_point_v1",
            "sequence": sequence,
            "coordinate": request["coordinate"],
            "workingCoordinate": request.get("working_coordinate", request["coordinate"]),
            "classification": request["classification"],
            "singularityIdentifier": request.get("singularity_identifier"),
            "role": request.get("role", "path"),
            **payload,
        }
        filename = f"flintnde_save_{sequence:03d}.json"
        target = self.output_directory / filename
        target.write_text(
            json.dumps(record, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        record["file"] = filename
        self.records.append(record)
        return target

    def finalize(self) -> Path:
        """整条链完成后写一次汇总；逐点文件此前已经存在。"""

        target = self.output_directory / self.summary_filename
        target.write_text(
            json.dumps(
                {
                    "schema": "flintnde_saved_points_v1",
                    "status": "complete",
                    "pointCount": len(self.records),
                    "points": self.records,
                },
                ensure_ascii=False,
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        return target
