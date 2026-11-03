"""MadStree 到 FlintNDE 的数值适配器。

MadStree 只提交连续复仿射单变量段、exact dlog 拉回、边界和 master 顺序。
FlintNDE 在本进程内完成边界输运、各段自动规划或顺序直输、dense 多点求值和
奇点领头阶解析。正规化级数控制协议接收 MadStree 的符号最低阶证书，生成自适应
exact ep 网格，并对 MadStree 回传的 Acb 终点球作 Laurent 重构和独立验证。
满足引擎合同（有限边界＋avoid 直输运＋无极点接触）的请求整体优先交给
vendor 内置的 C++ 引擎子进程，其余请求与运行失败一律走/回 Python 车道。
"""

from __future__ import annotations

from multiprocessing import get_context
from decimal import Decimal
from fractions import Fraction
import json
import math
import os
import re
import subprocess
import sys
import tempfile
import time
import traceback
import warnings
from pathlib import Path
from typing import Any


EVALUATE_SCHEMA = "madstree_flintnde_evaluate_v2"
EVALUATE_BATCH_SCHEMA = "madstree_flintnde_ep_batch_v1"
SERIES_CONTROL_SCHEMA = "madstree_flintnde_ep_series_control_v1"
DEFAULT_PARALLEL_TASK_COUNT = 12
REQUEST_KEYS = {
    "schema", "backendPackagePath", "masterDigest", "dimension", "tensorSectors",
    "segments", "pathPlanning", "singularityMode", "boundary",
    "workingPrecisionDigits", "primaryOrder", "referenceOrder",
    "targetRelativeError", "certificationMode", "messageLanguage",
    "columnVectorConvention", "dlogStatus",
}
MADSTREE_CITATION_NOTICE = {
    "intro": (
        "To support MadStree, please acknowledge our work when you use MadStree "
        "in your research. You can do so by citing the following papers:"
    ),
    "entries": [
        {
            "authors": "Jiaqi Chen and Bo Feng",
            "title": (
                "Towards Systematic Evaluation of de Sitter Correlators via "
                "Generalized Integration-By-Parts Relations"
            ),
            "arXiv": "2401.00129",
            "url": "https://arxiv.org/abs/2401.00129",
        },
        {
            "authors": "Jiaqi Chen, Bo Feng and Yi-Xiao Tao",
            "title": (
                "Multivariate hypergeometric solutions of cosmological (dS) "
                "correlators by d log-form differential equations"
            ),
            "arXiv": "2411.03088",
            "url": "https://arxiv.org/abs/2411.03088",
        },
    ],
    "packagePaper": "MadStree package paper, arXiv identifier pending.",
}

def _citation_text() -> str:
    """生成与 Wolfram 入口同序的人工 CLI 引用文本。"""

    lines = [MADSTREE_CITATION_NOTICE["intro"]]
    for index, entry in enumerate(MADSTREE_CITATION_NOTICE["entries"], start=1):
        lines.append(
            f'{index}. {entry["authors"]}, "{entry["title"]}", '
            f'arXiv:{entry["arXiv"]}, {entry["url"]}.'
        )
    lines.append(f'3. {MADSTREE_CITATION_NOTICE["packagePaper"]}')
    return "\n".join(lines)


def _with_citation_metadata(result: dict[str, Any]) -> dict[str, Any]:
    """为每个文件协议响应附加稳定引用元数据。"""

    enriched = dict(result)
    enriched["metadata"] = {"citationNotice": MADSTREE_CITATION_NOTICE}
    return enriched


def _exact_keys(value: Any, expected: set[str], label: str) -> dict[str, Any]:
    """要求对象字段与当前唯一 schema 完全一致。"""

    if not isinstance(value, dict):
        raise TypeError(f"{label} must be a JSON object")
    actual = set(value)
    if actual != expected:
        raise ValueError(
            f"{label} fields do not match the current schema; "
            f"missing={sorted(expected - actual)}, unexpected={sorted(actual - expected)}"
        )
    return value


def _integer(value: Any, label: str, minimum: int = 1) -> int:
    """读取严格整数，拒绝 bool 冒充整数。"""

    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise TypeError(f"{label} must be an integer >= {minimum}")
    return value


def _string(value: Any, label: str) -> str:
    """读取非空字符串。"""

    if not isinstance(value, str) or not value:
        raise TypeError(f"{label} must be a nonempty string")
    return value


def _sector_key(value: Any, label: str) -> str:
    """读取 sectorKey；无根传播子的 sector 用空串作键。"""

    if not isinstance(value, str):
        raise TypeError(f"{label} must be a string")
    return value


def _complex_record(value: Any, label: str) -> None:
    """验证 Q(i) 或十进制复数记录。"""

    item = _exact_keys(value, {"real", "imag"}, label)
    _string(item["real"], f"{label}.real")
    _string(item["imag"], f"{label}.imag")


def _matrix(value: Any, dimension: int, label: str, validator: Any) -> None:
    """验证指定元素合同的方阵。"""

    if not isinstance(value, list) or len(value) != dimension:
        raise ValueError(f"{label} must be a {dimension} by {dimension} matrix")
    for row_index, row in enumerate(value):
        if not isinstance(row, list) or len(row) != dimension:
            raise ValueError(f"{label}[{row_index}] must have length {dimension}")
        for column_index, entry in enumerate(row):
            validator(entry, f"{label}[{row_index}][{column_index}]")


def _rational_function_record(value: Any, label: str) -> None:
    """验证一个 Q(i)(s) 有理函数系数记录。"""

    item = _exact_keys(value, {"numerator", "denominator"}, label)
    for key in ("numerator", "denominator"):
        coefficients = item[key]
        if not isinstance(coefficients, list) or not coefficients:
            raise ValueError(f"{label}.{key} must be a nonempty list")
        for index, coefficient in enumerate(coefficients):
            _complex_record(coefficient, f"{label}.{key}[{index}]")


TENSOR_SECTOR_KEYS = {
    "sectorKey", "dimension", "offset", "rawCount", "slotBases",
    "rawToMasterColumn", "masterRepresentativeRaw",
}
TENSOR_LETTER_KEYS = {"alpha", "beta", "contributions"}
BOUNDARY_LETTER_KEYS = {"numerator", "denominator", "contributions"}
TENSOR_MIDDLE_KEYS = {"rows", "columns", "entries"}
TENSOR_DIAGONAL_KEYS = {"kind", "scope", "sectorKey", "vector"}
TENSOR_CONJUGATED_KEYS = {
    "kind", "scope", "sectorKey", "preVector", "letterIndicator", "middle",
    "postFactor",
}
LOCAL_BASE_TAGS = {"massiveT", "masslessH", "identity"}


def _tensor_sectors(value: Any, dimension: int) -> dict[str, int]:
    """验证每个 sector 的张量因子记录，并返回 sectorKey 到局部维数的映射。"""

    if not isinstance(value, list) or not value:
        raise ValueError("tensorSectors must be a nonempty list")
    dimensions: dict[str, int] = {}
    offsets: list[tuple[int, int]] = []
    for index, record in enumerate(value):
        label = f"tensorSectors[{index}]"
        item = _exact_keys(record, TENSOR_SECTOR_KEYS, label)
        key = _sector_key(item["sectorKey"], f"{label}.sectorKey")
        if key in dimensions:
            raise ValueError(f"{label}: duplicate sectorKey {key}")
        local = _integer(item["dimension"], f"{label}.dimension")
        offset = _integer(item["offset"], f"{label}.offset")
        raw_count = _integer(item["rawCount"], f"{label}.rawCount")
        bases = item["slotBases"]
        if not isinstance(bases, list) or any(
            base not in LOCAL_BASE_TAGS for base in bases
        ):
            raise ValueError(f"{label}.slotBases must be known local base tags")
        if raw_count != 1 << len(bases):
            raise ValueError(f"{label}.rawCount must equal 2**len(slotBases)")
        columns = item["rawToMasterColumn"]
        if not isinstance(columns, list) or len(columns) != raw_count:
            raise ValueError(f"{label}.rawToMasterColumn length mismatch")
        for position, column in enumerate(columns):
            if isinstance(column, bool) or not isinstance(column, int) or not (
                0 <= column <= local
            ):
                raise ValueError(f"{label}.rawToMasterColumn[{position}] out of range")
        representatives = item["masterRepresentativeRaw"]
        if not isinstance(representatives, list) or len(representatives) != local:
            raise ValueError(f"{label}.masterRepresentativeRaw length mismatch")
        for position, raw in enumerate(representatives):
            if isinstance(raw, bool) or not isinstance(raw, int) or not (
                1 <= raw <= raw_count
            ):
                raise ValueError(f"{label}.masterRepresentativeRaw[{position}] out of range")
            if columns[raw - 1] != position + 1:
                raise ValueError(
                    f"{label}: raw representative does not map to its master column"
                )
        dimensions[key] = local
        offsets.append((offset, local))
    expected = 1
    for offset, local in sorted(offsets):
        if offset != expected:
            raise ValueError("tensorSectors must tile the global master order")
        expected += local
    if expected - 1 != dimension:
        raise ValueError("tensorSectors dimensions do not match dimension")
    return dimensions


def _tensor_sector_dimension(
    dimensions: dict[str, int], key: str, label: str
) -> int:
    """按 sectorKey 读取局部维数。"""

    if key not in dimensions:
        raise ValueError(f"{label} refers to unknown sector {key}")
    return dimensions[key]


def _tensor_vector(value: Any, label: str, length: int) -> None:
    """验证一个固定长度的 Q(i) 稠密向量。"""

    if not isinstance(value, list) or len(value) != length:
        raise ValueError(f"{label} must have length {length}")
    for index, entry in enumerate(value):
        _complex_record(entry, f"{label}[{index}]")


def _tensor_indicator(value: Any, label: str, length: int) -> None:
    """验证 master 空间的 0/1 指示符。"""

    if not isinstance(value, list) or len(value) != length:
        raise ValueError(f"{label} must have length {length}")
    for index, bit in enumerate(value):
        if isinstance(bit, bool) or bit not in (0, 1):
            raise ValueError(f"{label}[{index}] must be 0 or 1")


def _tensor_middle(value: Any, local: int, dimension: int, label: str) -> None:
    """验证 contact 块的稀疏 (sector x global) 记录，索引均为 1-based。"""

    item = _exact_keys(value, TENSOR_MIDDLE_KEYS, label)
    for key in ("rows", "columns", "entries"):
        if not isinstance(item[key], list):
            raise ValueError(f"{label}.{key} must be a list")
    if not (len(item["rows"]) == len(item["columns"]) == len(item["entries"])):
        raise ValueError(f"{label} sparse fields must have equal length")
    for position, row in enumerate(item["rows"]):
        if isinstance(row, bool) or not isinstance(row, int) or not (1 <= row <= local):
            raise ValueError(f"{label}.rows[{position}] out of range")
    for position, column in enumerate(item["columns"]):
        if isinstance(column, bool) or not isinstance(column, int) or not (
            1 <= column <= dimension
        ):
            raise ValueError(f"{label}.columns[{position}] out of range")
    for position, entry in enumerate(item["entries"]):
        _complex_record(entry, f"{label}.entries[{position}]")


def _tensor_contributions(
    value: Any, dimensions: dict[str, int], dimension: int, label: str
) -> None:
    """验证一条 letter 的结构化贡献列表。"""

    if not isinstance(value, list) or not value:
        raise ValueError(f"{label} must be a nonempty list")
    for index, record in enumerate(value):
        item_label = f"{label}[{index}]"
        if not isinstance(record, dict):
            raise TypeError(f"{item_label} must be a JSON object")
        kind = record.get("kind")
        if kind == "diagonal":
            item = _exact_keys(record, TENSOR_DIAGONAL_KEYS, item_label)
            scope = item["scope"]
            if scope == "sector":
                key = _sector_key(item["sectorKey"], f"{item_label}.sectorKey")
                length = _tensor_sector_dimension(dimensions, key, item_label)
            elif scope == "global":
                if item["sectorKey"] is not None:
                    raise ValueError(
                        f"{item_label}.sectorKey must be null for global scope"
                    )
                length = dimension
            else:
                raise ValueError(f'{item_label}.scope must be "sector" or "global"')
            _tensor_vector(item["vector"], f"{item_label}.vector", length)
            continue
        if kind != "conjugated":
            raise ValueError(f'{item_label}.kind must be "diagonal" or "conjugated"')
        item = _exact_keys(record, TENSOR_CONJUGATED_KEYS, item_label)
        key = _sector_key(item["sectorKey"], f"{item_label}.sectorKey")
        local = _tensor_sector_dimension(dimensions, key, item_label)
        _tensor_indicator(item["letterIndicator"], f"{item_label}.letterIndicator", local)
        _complex_record(item["postFactor"], f"{item_label}.postFactor")
        if item["scope"] == "sector":
            if item["preVector"] is None:
                raise ValueError(f"{item_label}.preVector is required for sector scope")
            _tensor_vector(item["preVector"], f"{item_label}.preVector", local)
            if item["middle"] is not None:
                raise ValueError(
                    f"{item_label}.middle must be null for sector scope"
                )
        elif item["scope"] == "global":
            if item["preVector"] is not None:
                raise ValueError(
                    f"{item_label}.preVector must be null for global scope"
                )
            if item["middle"] is None:
                raise ValueError(f"{item_label}.middle is required for global scope")
            _tensor_middle(item["middle"], local, dimension, f"{item_label}.middle")
        else:
            raise ValueError(f'{item_label}.scope must be "sector" or "global"')


def _letters(value: Any, dimensions: dict[str, int], dimension: int, label: str) -> None:
    """验证结构化仿射 dlog letters。"""

    if not isinstance(value, list):
        raise TypeError(f"{label} must be a list")
    for index, record in enumerate(value):
        item = _exact_keys(record, TENSOR_LETTER_KEYS, f"{label}[{index}]")
        _complex_record(item["alpha"], f"{label}[{index}].alpha")
        _complex_record(item["beta"], f"{label}[{index}].beta")
        _tensor_contributions(
            item["contributions"], dimensions, dimension, f"{label}[{index}].contributions"
        )


def _validate_segment(
    value: Any, dimensions: dict[str, int], dimension: int, label: str
) -> None:
    """验证一个连续复仿射段。"""

    item = _exact_keys(
        value, {"start", "points", "letters", "fromUserIndex", "userIndices"}, label
    )
    _string(item["start"], f"{label}.start")
    _integer(item["fromUserIndex"], f"{label}.fromUserIndex", 0)
    if not isinstance(item["points"], list) or not item["points"]:
        raise ValueError(f"{label}.points must be a nonempty list")
    if not isinstance(item["userIndices"], list) or len(item["userIndices"]) != len(item["points"]):
        raise ValueError(f"{label}.userIndices must match points")
    for index, point in enumerate(item["points"]):
        _complex_record(point, f"{label}.points[{index}]")
        _integer(item["userIndices"][index], f"{label}.userIndices[{index}]")
    if len(set(item["userIndices"])) != len(item["userIndices"]):
        raise ValueError(f"{label}.userIndices must be distinct")
    _letters(item["letters"], dimensions, dimension, f"{label}.letters")


def _boundary_letters(
    value: Any, dimensions: dict[str, int], dimension: int, label: str
) -> None:
    """验证边界曲线 letter 拉回记录：分子/分母系数与张量贡献。"""

    if not isinstance(value, list):
        raise TypeError(f"{label} must be a list")
    for index, record in enumerate(value):
        item = _exact_keys(record, BOUNDARY_LETTER_KEYS, f"{label}[{index}]")
        for key in ("numerator", "denominator"):
            coefficients = item[key]
            if not isinstance(coefficients, list) or not coefficients:
                raise ValueError(f"{label}[{index}].{key} must be a nonempty list")
            for position, coefficient in enumerate(coefficients):
                _complex_record(coefficient, f"{label}[{index}].{key}[{position}]")
        _tensor_contributions(
            item["contributions"], dimensions, dimension,
            f"{label}[{index}].contributions",
        )


def _validate_boundary(
    value: Any, dimensions: dict[str, int], dimension: int
) -> None:
    """验证有限边界或正则奇点边界。"""

    if not isinstance(value, dict):
        raise TypeError("boundary must be a JSON object")
    kind = value.get("kind")
    if kind == "finite":
        item = _exact_keys(value, {"kind", "values"}, "boundary")
        if not isinstance(item["values"], list) or len(item["values"]) != dimension:
            raise ValueError("finite boundary vector dimension mismatch")
        for index, record in enumerate(item["values"]):
            _complex_record(record, f"boundary.values[{index}]")
        return
    item = _exact_keys(
        value,
        {"kind", "variable", "matrix", "branches", "weights", "start", "target"}
        | ({"certificate"} if isinstance(value, dict) and "certificate" in value else set())
        | ({"letters"} if isinstance(value, dict) and "letters" in value else set()),
        "boundary",
    )
    if kind != "regular_singular":
        raise ValueError('boundary.kind must be "finite" or "regular_singular"')
    _string(item["variable"], "boundary.variable")
    _string(item["start"], "boundary.start")
    _string(item["target"], "boundary.target")
    _matrix(item["matrix"], dimension, "boundary.matrix", _rational_function_record)
    if "letters" in item:
        _boundary_letters(
            item["letters"], dimensions, dimension, "boundary.letters"
        )
    if not isinstance(item["branches"], list) or not item["branches"]:
        raise ValueError("boundary.branches must be a nonempty list")
    if not isinstance(item["weights"], list) or len(item["weights"]) != len(item["branches"]):
        raise ValueError("boundary.weights must match branches")
    for index, branch in enumerate(item["branches"]):
        record = _exact_keys(branch, {"a", "b", "C"}, f"boundary.branches[{index}]")
        _complex_record(record["a"], f"boundary.branches[{index}].a")
        _integer(record["b"], f"boundary.branches[{index}].b", 0)
        if not isinstance(record["C"], list) or len(record["C"]) != dimension:
            raise ValueError("boundary branch vector dimension mismatch")
        for column, coefficient in enumerate(record["C"]):
            _complex_record(coefficient, f"boundary.branches[{index}].C[{column}]")
        _complex_record(item["weights"][index], f"boundary.weights[{index}]")
    if "certificate" in item:
        certificate = _exact_keys(
            item["certificate"], {"roots", "multiplicities"}, "boundary.certificate"
        )
        if not isinstance(certificate["roots"], list) or not certificate["roots"]:
            raise ValueError("boundary.certificate.roots must be a nonempty list")
        for index, root in enumerate(certificate["roots"]):
            _complex_record(root, f"boundary.certificate.roots[{index}]")
        if (
            not isinstance(certificate["multiplicities"], list)
            or len(certificate["multiplicities"]) != len(certificate["roots"])
        ):
            raise ValueError("boundary.certificate.multiplicities must match roots")
        for index, multiplicity in enumerate(certificate["multiplicities"]):
            _integer(multiplicity, f"boundary.certificate.multiplicities[{index}]")
        if sum(certificate["multiplicities"]) != dimension:
            raise ValueError(
                "boundary.certificate multiplicities must sum to the master dimension"
            )


def _validate_request(value: Any) -> dict[str, Any]:
    """验证唯一请求格式；旧 schema 一律拒绝。"""

    data = _exact_keys(value, REQUEST_KEYS, "request")
    if data["schema"] != EVALUATE_SCHEMA:
        raise ValueError("unsupported MadStree-FlintNDE schema")
    _string(data["backendPackagePath"], "backendPackagePath")
    _string(data["masterDigest"], "masterDigest")
    dimension = _integer(data["dimension"], "dimension")
    if not isinstance(data["pathPlanning"], bool):
        raise TypeError("pathPlanning must be a JSON boolean")
    if data["singularityMode"] not in {"avoid", "singularity_jump"}:
        raise ValueError('singularityMode must be "avoid" or "singularity_jump"')
    if not data["pathPlanning"] and data["singularityMode"] != "avoid":
        raise ValueError("direct user-point transport supports avoid mode only")
    if data["messageLanguage"] not in {"EN", "CN"}:
        raise ValueError('messageLanguage must be "EN" or "CN"')
    _integer(data["workingPrecisionDigits"], "workingPrecisionDigits", 20)
    primary = _integer(data["primaryOrder"], "primaryOrder")
    reference = _integer(data["referenceOrder"], "referenceOrder")
    if reference <= primary:
        raise ValueError("referenceOrder must exceed primaryOrder")
    _string(data["targetRelativeError"], "targetRelativeError")
    if data["certificationMode"] not in {"embedded", "certified"}:
        raise ValueError('certificationMode must be "embedded" or "certified"')
    if data["columnVectorConvention"] != "Y'=A(s)Y":
        raise ValueError("columnVectorConvention mismatch")
    if data["dlogStatus"] != "certifiedByFormulaChecks":
        raise ValueError("uncertified dlog input")
    sector_dimensions = _tensor_sectors(data["tensorSectors"], dimension)
    if not isinstance(data["segments"], list) or not data["segments"]:
        raise ValueError("segments must be a nonempty list")
    for index, segment in enumerate(data["segments"]):
        _validate_segment(
            segment, sector_dimensions, dimension, f"segments[{index}]"
        )
    _validate_boundary(data["boundary"], sector_dimensions, dimension)
    return data


def _acb_record(value: Any, digits: int) -> dict[str, Any]:
    """同时写出用户可读中点和用于正规化认证的完整 Arb 球。"""

    from flintnde.core import arb_ball_to_json

    return {
        "real": value.real.mid().str(digits, radius=False, more=True),
        "imag": value.imag.mid().str(digits, radius=False, more=True),
        "real_ball": arb_ball_to_json(value.real, digits),
        "imag_ball": arb_ball_to_json(value.imag, digits),
    }


def _point_record(value: Any, digits: int) -> dict[str, Any]:
    """序列化实际输运节点坐标。"""

    return _acb_record(value, digits)


def _resolve_backend_path(data: dict[str, Any]) -> Path:
    """解析 FlintNDE 包并把它加入本进程导入路径。"""

    backend_path = Path(data["backendPackagePath"]).resolve()
    package_file = backend_path / "flintnde" / "__init__.py"
    if not package_file.is_file():
        raise FileNotFoundError(f"FlintNDE package not found: {package_file}")
    if str(backend_path) not in sys.path:
        sys.path.insert(0, str(backend_path))
    return backend_path


ENGINE_SCHEMA = "flintnde_engine_evaluate_v1"
ENGINE_OVERRIDE_ENV = "MADSTREE_FLINTNDE_ENGINE"
ENGINE_EXECUTABLE = "flintnde_engine.exe" if os.name == "nt" else "flintnde_engine"


def _engine_executable(backend_path: Path) -> Path | None:
    """探测内置 C++ 引擎可执行文件；返回 None 表示本请求走 Python 车道。

    显式环境变量优先（"python" 或空串＝强制关闭引擎车道，指向不存在文件
    则报错）；缺省探测 vendor engine/build 下的构建产物，未构建即回退。
    """

    override = os.environ.get(ENGINE_OVERRIDE_ENV)
    if override is not None:
        if override.lower() in {"", "python"}:
            return None
        candidate = Path(override)
        if not candidate.is_file():
            raise FileNotFoundError(
                f"{ENGINE_OVERRIDE_ENV} points to a missing engine: {override}"
            )
        return candidate
    candidate = backend_path / "engine" / "build" / ENGINE_EXECUTABLE
    return candidate if candidate.is_file() else None


def _engine_poles(letters: list[dict[str, Any]]) -> list[Any]:
    """按精确 Q(i) 收集段 letters 的有限极点（镜像 _letters_to_system）。"""

    from flintnde import gaussian_rational

    poles = []
    for record in letters:
        beta = gaussian_rational(record["beta"])
        if beta.is_zero:
            continue
        poles.append(-gaussian_rational(record["alpha"]) / beta)
    return poles


def _engine_environment() -> dict[str, str]:
    """为引擎子进程定位 libflint 运行库：优先 FLINTNDE_ENGINE_LIBPATH，否则
    探测常见 MSYS2 UCRT64 位置；都不命中时不改 PATH，依赖当前 PATH。

    与泛用版引擎 smoke 测试的 engine_environment 同一约定。
    """

    environment = dict(os.environ)
    override = os.environ.get("FLINTNDE_ENGINE_LIBPATH")
    candidates = [override] if override else [
        "F:/codinglib/msys64/ucrt64/bin",
        "C:/msys64/ucrt64/bin",
    ]
    for directory in candidates:
        if directory and Path(directory).is_dir():
            environment["PATH"] = directory + os.pathsep + environment.get("PATH", "")
            break
    return environment


def _engine_lane(data: dict[str, Any], backend_path: Path) -> dict[str, Any] | None:
    """把整个请求交给 C++ 引擎执行；不满足合同任一条即返回 None 回退 Python。

    引擎合同：pathPlanning=false＋avoid、有限边界、段起点与全部用户点不落
    极点（正则奇点边界与奇点目标求值是 Python 专属能力）。翻译只需换
    schema——结构化 letters 与 tensorSectors 由引擎原生精确解析。引擎运行期
    错误 fail-closed 上抛，不静默回退。
    """

    if data["pathPlanning"] or data["singularityMode"] != "avoid":
        return None
    if data["boundary"]["kind"] != "finite":
        return None
    engine = _engine_executable(backend_path)
    if engine is None:
        return None
    from flintnde import gaussian_rational

    for segment in data["segments"]:
        poles = _engine_poles(segment["letters"])
        values = [gaussian_rational(segment["start"])]
        values.extend(gaussian_rational(point) for point in segment["points"])
        if any((value - pole).is_zero for value in values for pole in poles):
            return None
    request = dict(data, schema=ENGINE_SCHEMA)
    with tempfile.TemporaryDirectory(prefix="madstree-flintnde-engine-") as scratch:
        input_file = Path(scratch) / "request.json"
        output_file = Path(scratch) / "response.json"
        input_file.write_text(json.dumps(request, ensure_ascii=False), encoding="utf-8")
        completed = subprocess.run(
            [str(engine), str(input_file), str(output_file)], capture_output=True,
            env=_engine_environment(),
        )
        if not output_file.is_file():
            raise RuntimeError(
                "FlintNDE engine produced no output file; "
                f"exit={completed.returncode}, "
                f"stderr={completed.stderr.decode('utf-8', 'replace')[-2000:]}"
            )
        result = json.loads(output_file.read_text(encoding="utf-8"))
    result["schema"] = data["schema"]
    return result


def _tensor_system(data: dict[str, Any]) -> Any:
    """按请求中的 sector 因子重建 exact 张量 dlog 系统。"""

    from flintnde import TensorDLogSystem, TensorSector

    sectors = tuple(
        TensorSector(
            sector_key=item["sectorKey"],
            dimension=item["dimension"],
            offset=item["offset"],
            raw_count=item["rawCount"],
            slot_bases=tuple(item["slotBases"]),
            raw_to_master_column=tuple(item["rawToMasterColumn"]),
            master_representative_raw=tuple(item["masterRepresentativeRaw"]),
        )
        for item in data["tensorSectors"]
    )
    return TensorDLogSystem(sectors, data["dimension"])


def _tensor_contributions_of(records: list[dict[str, Any]]) -> list[Any]:
    """把 JSON 贡献记录转成 exact 张量贡献对象；稀疏索引改为 0-based。"""

    from flintnde import TensorContribution, gaussian_rational

    contributions = []
    for record in records:
        if record["kind"] == "diagonal":
            contributions.append(TensorContribution(
                kind="diagonal",
                scope=record["scope"],
                sector_key=record["sectorKey"],
                vector=tuple(gaussian_rational(item) for item in record["vector"]),
            ))
            continue
        middle = record["middle"]
        contributions.append(TensorContribution(
            kind="conjugated",
            scope=record["scope"],
            sector_key=record["sectorKey"],
            pre_vector=(
                None if record["preVector"] is None
                else tuple(gaussian_rational(item) for item in record["preVector"])
            ),
            letter_indicator=tuple(record["letterIndicator"]),
            middle_rows=() if middle is None else tuple(row - 1 for row in middle["rows"]),
            middle_columns=(
                () if middle is None else tuple(column - 1 for column in middle["columns"])
            ),
            middle_values=(
                () if middle is None
                else tuple(gaussian_rational(item) for item in middle["entries"])
            ),
            post_factor=gaussian_rational(record["postFactor"]),
        ))
    return contributions


def _letters_to_system(
    records: list[dict[str, Any]], dimension: int, tensor_system: Any, name: str
) -> Any:
    """按精确极点合并整块直积算子，普通段不构造全局稠密留数。"""

    from flint import acb_mat
    from flintnde import PartialFractionSystem, gaussian_rational
    from flintnde.block_operators import BlockResidues, compile_block_residue

    merged: dict[Any, list[Any]] = {}
    for record in records:
        beta = gaussian_rational(record["beta"])
        if beta.is_zero:
            continue
        alpha = gaussian_rational(record["alpha"])
        if alpha.is_zero:
            raise ValueError("dlog letter vanishes at the segment anchor")
        pole = -alpha / beta
        terms = merged.setdefault(pole, [])
        # 先保留各项的 active 支撑，避免相加后掩盖不同的单位直积。
        for contribution in _tensor_contributions_of(record["contributions"]):
            terms.append(tensor_system.residue_entries((contribution,)))
    partition = tuple(
        (sector.first_global_index, sector.dimension) for sector in tensor_system.sectors
    )
    residues = BlockResidues(tuple(
        compile_block_residue(dimension, partition, terms) for terms in merged.values()
    ))
    return PartialFractionSystem(
        constant=acb_mat(dimension, dimension),
        residues=residues,
        poles=tuple(pole.to_acb() for pole in merged),
        name=name,
    )


def _rational_system(
    boundary: dict[str, Any], dimension: int, name: str,
    block_partition: tuple[tuple[int, int], ...] = (),
) -> Any:
    """重建正则奇点边界的一般 Q(i)(t) 系统。"""

    from flintnde import RationalMatrixSystem, rational_function

    matrix = tuple(
        tuple(rational_function(entry["numerator"], entry["denominator"]) for entry in row)
        for row in boundary["matrix"]
    )
    return RationalMatrixSystem(
        matrix, variable_name=boundary["variable"], name=name, block_partition=block_partition,
    )


def _isolated_roots(coefficients: list[Any]) -> list[tuple[Any, int]] | None:
    """把互素多项式的根隔离成 (root, multiplicity)；任一环节不满足返回 None。

    零根按系数表前导零个数精确剥离；一次因子用精确 Q(i) 除法；二次、三次因子
    用 acb 球根隔离并验证 |N(r)| 含 0、根数等于次数（squarefree）。代数根一律
    以认证小数包络表示，不进入符号层。
    """

    from flint import acb, acb_poly, arb, ctx

    zero_multiplicity = 0
    while zero_multiplicity < len(coefficients) - 1 and coefficients[zero_multiplicity].is_zero:
        zero_multiplicity += 1
    core = coefficients[zero_multiplicity:]
    degree = len(core) - 1
    if degree <= 0:
        roots = []
    elif degree == 1:
        roots = [(-core[0] / core[1]).to_acb()]
    else:
        values = [value.to_acb() for value in core]
        polynomial = acb_poly(values)
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


def _boundary_letters_to_system(
    records: list[dict[str, Any]], dimension: int, tensor_system: Any, name: str
) -> Any | None:
    """按 d Log[N/D] 分解装配边界 letter 快速系统；不满足前置条件返回 None。

    d Log[N(t)/D(t)] 的极点在 N 的根（留数 +R）与 D 的根（留数 -R），重数按
    因子幂次计入；不同 letter 的同名极点按球包含合并。矩阵结构保留整块直积
    张量形式，不构造全局稠密留数。
    """

    from flint import acb_mat
    from flintnde import PartialFractionSystem, gaussian_rational
    from flintnde.block_operators import BlockResidues, compile_block_residue

    max_degree = 3
    poles: list[Any] = []
    term_buckets: list[list[Any]] = []
    for record in records:
        numerator = [gaussian_rational(item) for item in record["numerator"]]
        denominator = [gaussian_rational(item) for item in record["denominator"]]
        if max(len(numerator), len(denominator)) - 1 > max_degree:
            return None
        contribution_terms = [
            tensor_system.residue_entries((contribution,))
            for contribution in _tensor_contributions_of(record["contributions"])
        ]
        if not contribution_terms:
            continue
        for coefficients, sign in ((numerator, 1), (denominator, -1)):
            if len(coefficients) == 1:
                continue
            roots = _isolated_roots(coefficients)
            if roots is None:
                return None
            for root, multiplicity in roots:
                index = next(
                    (position for position, pole in enumerate(poles)
                     if (root - pole).contains(0)),
                    None,
                )
                if index is None:
                    poles.append(root)
                    term_buckets.append([])
                    index = len(poles) - 1
                terms = (
                    contribution_terms if sign == 1
                    else [{position: -value for position, value in term.items()}
                          for term in contribution_terms]
                )
                for _ in range(multiplicity):
                    term_buckets[index].extend(terms)
    partition = tuple(
        (sector.first_global_index, sector.dimension) for sector in tensor_system.sectors
    )
    try:
        residues = BlockResidues(tuple(
            compile_block_residue(dimension, partition, terms) for terms in term_buckets
        ))
        return PartialFractionSystem(
            constant=acb_mat(dimension, dimension),
            residues=residues,
            poles=tuple(poles),
            name=name,
        )
    except ValueError:
        return None


def _boundary_letters_match(fast_system: Any, system: Any, probes: list[Any]) -> bool:
    """在探针点比对 letter 快速系统与有理系统的矩阵球，全部包含才通过。

    两条求值路径的舍入差为工作精度量级（约 1 ulp），探针容差取
    10^-(prec/6)：远高于舍入噪声、又比目标相对误差 1e-20 严数十个量级，
    只用来拦截装配性错误（错极点/漏 letter 的偏差是 O(1) 量级）。
    """

    from flint import arb, ctx

    tolerance = arb(10) ** (-(ctx.prec // 6))
    for probe in probes:
        fast = fast_system.evaluate(probe)
        reference = system.evaluate(probe)
        for row in range(system.dimension):
            for column in range(system.dimension):
                difference = abs(fast[row, column] - reference[row, column])
                scale = abs(reference[row, column]) + arb(1)
                if not difference < tolerance * scale:
                    return False
    return True


def _initial_vectors(
    data: dict[str, Any], digits: int, tensor_system: Any
) -> tuple[Any, Any | None, dict[str, Any]]:
    """分别构造主链与参考链初值；有限精确边界只需保存一份。

    边界请求携带 letter list 时，先装配极点递推快速系统并在探针点对表；
    任一环节失败即回退现行 Cauchy 路线，路径规划与 Frobenius 构造始终
    使用精确有理系统。
    """

    from flint import acb, acb_mat
    from flintnde import (
        IndicialSpectrumError, NamedPoint, build_adaptive_path, frobenius_boundary,
        gaussian_rational, transport_frobenius_boundaries_refined,
    )

    boundary = data["boundary"]
    dimension = data["dimension"]
    if boundary["kind"] == "finite":
        vector = acb_mat([[acb(item["real"], item["imag"])] for item in boundary["values"]])
        return vector, None, {"kind": "finite", "path": [], "seconds": 0.0}
    partition = tuple((item["offset"] - 1, item["dimension"]) for item in data["tensorSectors"])
    system = _rational_system(
        boundary, dimension, data["masterDigest"] + "-boundary", partition,
    )
    start = gaussian_rational(boundary["start"]).to_acb()
    target = gaussian_rational(boundary["target"]).to_acb()
    branches = [frobenius_boundary([record]) for record in boundary["branches"]]
    known_roots = None
    certificate = boundary.get("certificate")
    if certificate is not None:
        # 重数以内展开成含重复的 roots 列表；后端只做验证式消费，不再求谱。
        known_roots = [
            gaussian_rational(root)
            for root, multiplicity in zip(
                certificate["roots"], certificate["multiplicities"]
            )
            for _ in range(multiplicity)
        ]
    numerical_digits = max(1, digits // 2)

    def _launch(numerical: bool):
        """按 exact 或小数化数值路线建立边界路径、letter 快速系统并完成输运。"""

        launch_path = build_adaptive_path(
            system, NamedPoint("madstree_boundary", start), NamedPoint("madstree_anchor", target),
            path_name="madstree_boundary_to_anchor", max_step_over_radius=0.45,
            singularity_mode="avoid", message_language=data["messageLanguage"],
            numerical_roots=known_roots if numerical else None,
            numerical_precision_digits=numerical_digits if numerical else None,
        )
        launch_fast = None
        if boundary.get("letters"):
            launch_fast = _boundary_letters_to_system(
                boundary["letters"], dimension, tensor_system,
                data["masterDigest"] + "-boundary-letters",
            )
            if launch_fast is not None:
                probes = [(launch_path[index] + launch_path[index + 1]) / 2
                          for index in (0, len(launch_path) - 2)]
                if not _boundary_letters_match(launch_fast, system, probes):
                    launch_fast = None
        launch_result = transport_frobenius_boundaries_refined(
            system, branches, launch_path,
            primary_order=data["primaryOrder"], reference_order=data["referenceOrder"],
            radius_fraction=0.60, target_relative_error=data["targetRelativeError"],
            certification_mode=data["certificationMode"],
            known_roots=None if numerical else known_roots,
            fast_system=launch_fast,
            numerical_roots=known_roots if numerical else None,
            numerical_precision_digits=numerical_digits if numerical else None,
        )
        return launch_path, launch_fast, launch_result

    clock = time.perf_counter()
    frobenius_route = "exact"
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            path, fast_system, result = _launch(False)
        except IndicialSpectrumError:
            if known_roots is None:
                raise
            # exact 谱门禁失败（根式参数等非 Q(i) 谱）：按用户裁定自动改走小数化
            # 数值 Frobenius 路线；证书谱与后端独立特征值计算逐根交叉核对。
            # 声称谱的小数化扰动约 10^-digits，判定精度取一半保证斩杀线高于扰动。
            frobenius_route = "numerical"
            caught.clear()
            path, fast_system, result = _launch(True)
    weights = [acb(item["real"], item["imag"]) for item in boundary["weights"]]

    def weighted_vector(matrix: Any) -> Any:
        """按同一物理边界权重合并一条输运链的基解矩阵。"""

        vector = acb_mat(dimension, 1)
        for row in range(dimension):
            vector[row, 0] = sum(
                (weights[column] * matrix[row, column] for column in range(len(weights))),
                acb(0),
            )
        return vector

    primary_vector = weighted_vector(result["primary_snapshots"][-1])
    reference_vector = weighted_vector(result["reference_snapshots"][-1])
    return primary_vector, reference_vector, {
        "kind": "regular_singular",
        "frobeniusRoute": frobenius_route,
        "path": [_point_record(point, digits) for point in path],
        "pathPointCount": len(path),
        "seconds": time.perf_counter() - clock,
        "warnings": [str(item.message) for item in caught],
        "letterFastPath": fast_system is not None,
        "letterPoleCount": len(fast_system.poles) if fast_system is not None else 0,
    }


def _point_assignments(plan: Any, points: list[Any], user_indices: list[int]) -> list[dict[str, Any]]:
    """把每个用户点映射到节点快照或 dense sample。"""

    sample_by_index = {
        int(item["user_point_index"]): item for item in plan.sample_assignments
        if item["user_point_index"] is not None
    }
    assignments = []
    for local_index, (point, user_index) in enumerate(zip(points, user_indices)):
        if local_index in sample_by_index:
            sample = sample_by_index[local_index]
            assignments.append({
                "userIndex": user_index, "localPointIndex": local_index,
                "source": sample["source"], "nodeIndex": None,
                "segmentIndex": int(sample["segment_index"]),
            })
            continue
        nodes = [index for index, node in enumerate(plan.nodes) if abs(node - point).contains(0)]
        if not nodes:
            raise ValueError(f"user point {user_index} has no FlintNDE value assignment")
        assignments.append({
            "userIndex": user_index, "localPointIndex": local_index,
            "source": "node_snapshot", "nodeIndex": nodes[-1], "segmentIndex": None,
        })
    return assignments


def _point_values(result: dict[str, Any], assignments: list[dict[str, Any]], digits: int, dimension: int) -> list[dict[str, Any]]:
    """按映射合并节点快照与批量 dense 输出。"""

    dense = {
        int(item["user_point_index"]): item for item in result["sample_results"]
        if item["user_point_index"] is not None
    }
    singular = {
        int(item["user_point_index"]): item
        for item in result.get("singular_target_results", [])
        if item["user_point_index"] is not None
    }
    records = []
    for assignment in assignments:
        local_index = assignment["localPointIndex"]
        node_index = assignment["nodeIndex"]
        singular_item = singular.get(local_index)
        if singular_item is not None:
            records.append({
                "userIndex": assignment["userIndex"],
                "localPointIndex": local_index,
                "source": assignment["source"],
                "evaluationAlgorithm": None,
                "values": [
                    _value_record(value, digits)
                    for value in singular_item["values"]
                ],
            })
            continue
        item = dense.get(local_index) if node_index is None else None
        vector = item["value"] if item is not None else result["primary_snapshots"][node_index]
        records.append({
            "userIndex": assignment["userIndex"], "localPointIndex": local_index,
            "source": assignment["source"],
            "evaluationAlgorithm": None if item is None else item.get("evaluation_algorithm"),
            "values": [_acb_record(vector[row, 0], digits) for row in range(dimension)],
        })
    return records


def _value_record(value: Any, digits: int) -> dict[str, Any]:
    """把有限 Acb 或奇点状态转换为统一的用户结果记录。"""

    if isinstance(value, str):
        return {"text": value}
    return _acb_record(value, digits)


def _singular_result_difference(primary: Any, reference: Any) -> tuple[Any, str | None]:
    """比较末端奇点主/参考结果；发散状态只允许逐分量完全一致。"""

    from flint import arb

    if primary.component_classifications != reference.component_classifications:
        return arb(0), "primary/reference singular-target classifications do not agree"
    difference = arb(0)
    for primary_value, reference_value in zip(primary.values, reference.values):
        if isinstance(primary_value, str) or isinstance(reference_value, str):
            if primary_value != reference_value:
                return difference, "primary/reference singular-target divergence states do not agree"
            continue
        denominator = max(
            arb(1), abs(reference_value), key=lambda item: float(item.mid())
        )
        difference = max(difference, abs(primary_value - reference_value) / denominator)
    return difference, None


def _letters_to_rational_system(
    records: list[dict[str, Any]], dimension: int, tensor_system: Any, name: str
) -> Any:
    """从结构化张量 letters 重建局部幂对数分析所需的有理矩阵系统。"""

    from flintnde import RationalMatrixSystem, gaussian_rational, rational_function

    rows = [[rational_function(0) for _ in range(dimension)] for _ in range(dimension)]
    for record in records:
        alpha = gaussian_rational(record["alpha"])
        beta = gaussian_rational(record["beta"])
        if beta.is_zero:
            continue
        denominator = (alpha, beta)
        entries = tensor_system.residue_entries(
            _tensor_contributions_of(record["contributions"])
        )
        for (row, column), residue in entries.items():
            # record 表示 R dlog(alpha+beta*s)，局部系统必须保留导数因子 beta。
            rows[row][column] += rational_function(beta * residue, denominator)
    return RationalMatrixSystem(tuple(tuple(row) for row in rows), variable_name="s", name=name)


def _run(data: dict[str, Any]) -> dict[str, Any]:
    """在单一 Python 进程内完成全部边界和复仿射段；满足合同时整请求优先交引擎。"""

    backend_path = _resolve_backend_path(data)
    engine_result = _engine_lane(data, backend_path)
    if engine_result is not None:
        return engine_result
    from flint import acb, acb_mat, arb
    from flintnde import (
        SingularPathError, configure_working_precision, direct_user_point_path,
        evaluate_singular_target, gaussian_rational, plan_singular_target_match,
        plan_transport_path, transport_planned_path_refined,
    )

    digits = data["workingPrecisionDigits"]
    configure_working_precision(digits)
    tensor_system = _tensor_system(data)
    primary_vector, reference_vector, boundary_report = _initial_vectors(
        data, digits, tensor_system
    )
    if reference_vector is None:
        reference_vector = acb_mat(primary_vector)
    segments_out = []
    all_met = True
    primary_seconds = 0.0
    reference_seconds = 0.0
    planning_seconds = 0.0
    singularity_classifications = []
    for segment_index, segment in enumerate(data["segments"]):
        system = _letters_to_system(
            segment["letters"], data["dimension"], tensor_system,
            f"MadStree-{data['masterDigest'][:12]}-{segment_index}",
        )
        start = gaussian_rational(segment["start"]).to_acb()
        points = [gaussian_rational(point).to_acb() for point in segment["points"]]
        user_indices = list(segment["userIndices"])
        singular_positions = [
            index for index, point in enumerate(points)
            if any((point - pole).contains(0) for pole in system.singularities)
        ]
        terminal_singular_position = (
            len(points) - 1
            if singular_positions and singular_positions[-1] == len(points) - 1
            else None
        )
        intermediate_singular_positions = [
            position for position in singular_positions
            if position != terminal_singular_position
        ]
        if intermediate_singular_positions and (
            not data["pathPlanning"] or data["singularityMode"] != "singularity_jump"
        ):
            return {
                "status": "intermediateSingularPointRequiresPlannedJump", "schema": data["schema"],
                "segmentIndex": segment_index,
                "message": (
                    "intermediate singular user points require FlintNDE path planning "
                    "with singularity-jump mode"
                ),
                "singularUserIndices": [
                    user_indices[index] for index in intermediate_singular_positions
                ],
            }
        singular_position = terminal_singular_position
        regular_points = points if singular_position is None else points[:singular_position]
        regular_indices = user_indices if singular_position is None else user_indices[:singular_position]
        clock = time.perf_counter()
        point_records = []
        plan = None
        result = None
        terminal_plan = None
        terminal_result = None
        assignments = []
        met = True
        if regular_points:
            try:
                plan = (
                    plan_transport_path(
                        system, start, regular_points, singularity_mode=data["singularityMode"],
                        message_language=data["messageLanguage"],
                    )
                    if data["pathPlanning"]
                    else direct_user_point_path(
                        system, start, regular_points, message_language=data["messageLanguage"]
                    )
                )
            except SingularPathError as error:
                return {
                    "status": "singularPathRefused", "schema": data["schema"],
                    "segmentIndex": segment_index, "message": str(error),
                    "singularPathPairs": [[_point_record(a, digits), _point_record(b, digits)] for a, b in error.singular_path_pairs],
                }
            planning_seconds += time.perf_counter() - clock
            assignments = _point_assignments(plan, regular_points, regular_indices)
            result = transport_planned_path_refined(
                system, primary_vector, plan,
                reference_initial_vector=reference_vector,
                primary_order=data["primaryOrder"], reference_order=data["referenceOrder"],
                radius_fraction=0.60, target_relative_error=data["targetRelativeError"],
                certification_mode=data["certificationMode"],
            )
            primary_vector = result["primary_snapshots"][-1]
            reference_vector = result["reference_snapshots"][-1]
            primary_seconds += float(result["primary_seconds"])
            reference_seconds += float(result["reference_seconds"])
            met = result["target_relative_error_met"]
            all_met = all_met and (met is None or bool(met))
            point_records = _point_values(result, assignments, digits, data["dimension"])
            for singular_item in result.get("singular_target_results", []):
                local_index = int(singular_item["user_point_index"])
                singularity_classifications.append({
                    "segmentIndex": segment_index,
                    "userIndex": regular_indices[local_index],
                    "coordinate": _point_record(regular_points[local_index], digits),
                    **singular_item["report"],
                })
        if singular_position is not None:
            if segment_index != len(data["segments"]) - 1:
                return {
                    "status": "singularTurnContinuationUnsupported",
                    "schema": data["schema"],
                    "segmentIndex": segment_index,
                    "message": (
                        "a singular point at an affine-segment turn does not fix a unique "
                        "continuation branch"
                    ),
                    "singularUserIndex": user_indices[singular_position],
                }
            target = points[singular_position]
            target_exact = gaussian_rational(segment["points"][singular_position])
            rational_system = _letters_to_rational_system(
                segment["letters"], data["dimension"], tensor_system,
                f"MadStree-{data['masterDigest'][:12]}-{segment_index}-target",
            )
            start_for_target = (
                gaussian_rational(segment["points"][singular_position - 1])
                if regular_points else gaussian_rational(segment["start"])
            )
            match_clock = time.perf_counter()
            match_plan = plan_singular_target_match(
                rational_system, start_for_target, target_exact,
            )
            if match_plan.inserted:
                if not data["pathPlanning"]:
                    return {
                        "status": "singularTargetMatchPlanningRequired",
                        "schema": data["schema"],
                        "segmentIndex": segment_index,
                        "message": (
                            "the terminal singular target needs a hidden ordinary match point "
                            "inside its convergence disk; enable FlintNDE path planning"
                        ),
                        "singularUserIndex": user_indices[singular_position],
                        "terminalSingularMatch": match_plan.report,
                    }
                terminal_plan = plan_transport_path(
                    system,
                    start_for_target.to_acb(),
                    [match_plan.match_point],
                    singularity_mode=data["singularityMode"],
                    message_language=data["messageLanguage"],
                )
                planning_seconds += time.perf_counter() - match_clock
                terminal_result = transport_planned_path_refined(
                    system,
                    primary_vector,
                    terminal_plan,
                    reference_initial_vector=reference_vector,
                    primary_order=data["primaryOrder"],
                    reference_order=data["referenceOrder"],
                    radius_fraction=0.60,
                    target_relative_error=data["targetRelativeError"],
                    certification_mode=data["certificationMode"],
                )
                primary_vector = terminal_result["primary_snapshots"][-1]
                reference_vector = terminal_result["reference_snapshots"][-1]
                primary_seconds += float(terminal_result["primary_seconds"])
                reference_seconds += float(terminal_result["reference_seconds"])
                terminal_transport_met = terminal_result["target_relative_error_met"]
                all_met = all_met and (
                    terminal_transport_met is None or bool(terminal_transport_met)
                )
            else:
                planning_seconds += time.perf_counter() - match_clock
            primary_singular = evaluate_singular_target(
                rational_system, match_plan.match_point, primary_vector, target_exact,
                order=data["primaryOrder"], target_relative_error=data["targetRelativeError"],
            )
            reference_singular = evaluate_singular_target(
                rational_system, match_plan.match_point, reference_vector, target_exact,
                order=data["referenceOrder"], target_relative_error=data["targetRelativeError"],
            )
            singular_difference, singular_error = _singular_result_difference(
                primary_singular, reference_singular
            )
            if singular_error is not None:
                return {
                    "status": "singularTargetReferenceMismatch",
                    "schema": data["schema"],
                    "segmentIndex": segment_index,
                    "message": singular_error,
                    "singularUserIndex": user_indices[singular_position],
                    "primaryClassification": primary_singular.classification,
                    "referenceClassification": reference_singular.classification,
                    "terminalSingularMatch": match_plan.report,
                }
            singular_met = bool(singular_difference < arb(data["targetRelativeError"]))
            all_met = all_met and singular_met
            met = singular_met if met is None else bool(met) and singular_met
            point_records.append({
                "userIndex": user_indices[singular_position],
                "localPointIndex": singular_position,
                "source": "singular_target",
                "evaluationAlgorithm": None,
                "values": [_value_record(value, digits) for value in primary_singular.values],
            })
            singularity_classifications.append({
                "segmentIndex": segment_index,
                "userIndex": user_indices[singular_position],
                "coordinate": _point_record(target, digits),
                **primary_singular.report,
                "referenceClassification": reference_singular.classification,
                "primaryReferenceRelativeDifference": singular_difference.str(digits),
                "targetRelativeErrorMet": singular_met,
                "terminalSingularMatch": match_plan.report,
            })
        actual_nodes = [] if plan is None else list(plan.nodes)
        if singular_position is not None and terminal_plan is not None:
            terminal_nodes = list(terminal_plan.nodes)
            if actual_nodes and terminal_nodes and abs(actual_nodes[-1] - terminal_nodes[0]).contains(0):
                actual_nodes.extend(terminal_nodes[1:])
            else:
                actual_nodes.extend(terminal_nodes)
        plan_report = {} if plan is None else dict(plan.report)
        if singular_position is not None:
            plan_report["terminal_singular_match"] = match_plan.report
        difference_candidates = []
        if result is not None:
            difference_candidates.append(result["relative_difference_inf"])
        if singular_position is not None:
            difference_candidates.append(singular_difference)
        if terminal_result is not None:
            difference_candidates.append(terminal_result["relative_difference_inf"])
        segment_difference = max(difference_candidates) if difference_candidates else None
        point_records.sort(key=lambda item: item["userIndex"])
        segments_out.append({
            "segmentIndex": segment_index,
            "pathPlanning": data["pathPlanning"],
            "actualNodes": [_point_record(node, digits) for node in actual_nodes],
            "nodeCount": len(actual_nodes),
            "coveredSampleCount": 0 if plan is None else len(plan.sample_assignments),
            "planReport": plan_report,
            "pointAssignments": assignments,
            "pointValues": point_records,
            "relativeDifferenceInf": None if segment_difference is None else segment_difference.str(digits),
            "relativeDifferenceMidpoint": None if segment_difference is None else float(segment_difference.mid()),
            "targetRelativeErrorMet": True if met is None else bool(met),
            "certificationMode": None if result is None else result["certification_mode"],
            "continuationReady": singular_position is None,
            "endpointValues": (
                [_value_record(value, digits) for value in primary_singular.values]
                if singular_position is not None
                else [
                    _acb_record(primary_vector[row, 0], digits)
                    for row in range(data["dimension"])
                ]
            ),
        })
    message = (
        "FlintNDE planned and evaluated every MadStree affine segment in one process."
        if data["pathPlanning"] else
        "FlintNDE evaluated every MadStree affine segment using the supplied user points as nodes."
    ) if data["messageLanguage"] == "EN" else (
        "FlintNDE 已在同一进程内规划并计算全部 MadStree 复仿射段。"
        if data["pathPlanning"] else
        "FlintNDE 已把用户给定点严格作为节点，在同一进程内计算全部 MadStree 复仿射段。"
    )
    return {
        "status": "success", "schema": data["schema"],
        "workerPid": os.getpid(),
        "executionAction": "plan_and_execute" if data["pathPlanning"] else "execute_user_nodes",
        "masterDigest": data["masterDigest"], "dimension": data["dimension"],
        "backendPackagePath": str(backend_path), "messageLanguage": data["messageLanguage"],
        "message": message, "pathPlanning": data["pathPlanning"],
        "segmentCount": len(segments_out), "segments": segments_out,
        "singularityClassifications": singularity_classifications,
        "boundary": boundary_report, "planningSeconds": planning_seconds,
        "primarySeconds": primary_seconds, "referenceSeconds": reference_seconds,
        "targetRelativeErrorMet": all_met,
        "finalValues": [
            _acb_record(primary_vector[row, 0], digits)
            for row in range(data["dimension"])
        ],
    }


def _run_indexed_task(item: tuple[int, dict[str, Any]]) -> tuple[int, dict[str, Any]]:
    """在独立进程执行一个固定 ep 请求，并保留输入序号。"""

    index, data = item
    return index, _run(data)


def _validate_batch_request(value: Any) -> dict[str, Any]:
    """验证 MadStree 固定 ep 批量请求和缺省 12 并行合同。"""

    data = _exact_keys(
        value, {"schema", "parallelTaskCount", "tasks", "messageLanguage"},
        "ep batch request",
    )
    if data["schema"] != EVALUATE_BATCH_SCHEMA:
        raise ValueError("unsupported MadStree-FlintNDE ep batch schema")
    parallel_count = _integer(data["parallelTaskCount"], "parallelTaskCount")
    if data["messageLanguage"] not in {"EN", "CN"}:
        raise ValueError('messageLanguage must be "EN" or "CN"')
    if not isinstance(data["tasks"], list) or not data["tasks"]:
        raise ValueError("ep batch tasks must be a nonempty list")
    tasks = []
    for index, task in enumerate(data["tasks"]):
        item = _exact_keys(task, {"ep", "request"}, f"tasks[{index}]")
        tasks.append({"ep": item["ep"], "request": _validate_request(item["request"])})
    return {
        "schema": data["schema"],
        "parallelTaskCount": parallel_count,
        "tasks": tasks,
        "messageLanguage": data["messageLanguage"],
    }


def _run_batch(data: dict[str, Any]) -> dict[str, Any]:
    """用有界进程池运行不同 ep；每项隔离 FLINT 上下文并按输入顺序返回。"""

    tasks = data["tasks"]
    requested = data["parallelTaskCount"]
    effective = min(requested, len(tasks))
    payloads = [(index, task["request"]) for index, task in enumerate(tasks)]
    with get_context("spawn").Pool(
        processes=effective, maxtasksperchild=1
    ) as pool:
        completed = pool.map(_run_indexed_task, payloads, chunksize=1)
    results = [
        {"ep": tasks[index]["ep"], "result": result}
        for index, result in completed
    ]
    success = all(item["result"].get("status") == "success" for item in results)
    language = data["messageLanguage"]
    message = (
        f"MadStree/FlintNDE ep task pool used {effective} of requested {requested} "
        "workers (default 12); queued tasks started automatically as workers finished."
        if language == "EN" else
        f"MadStree/FlintNDE 不同 ep 任务池：请求并行数 {requested}，实际并行数 {effective}"
        "（缺省 12）；任务完成后已自动续交队列。"
    )
    return {
        "status": "success" if success else "failed",
        "schema": data["schema"],
        "message": message,
        "messageLanguage": language,
        "parallelTaskCountRequested": requested,
        "parallelTaskCountEffective": effective,
        "taskCount": len(results),
        "results": results,
    }


def _validate_series_control_request(value: Any) -> dict[str, Any]:
    """验证自适应正规化控制请求；每个阶段只接受自己的严格字段集合。"""

    if not isinstance(value, dict) or value.get("schema") != SERIES_CONTROL_SCHEMA:
        raise ValueError("unsupported MadStree-FlintNDE ep series control schema")
    action = value.get("action")
    common = {"schema", "action", "backendPackagePath", "maximumPower", "goalDigits"}
    production_base = common | {
            "leadingPower", "sampleSpacing", "validationSampleCount", "validationScale",
            "maximumSamples", "extraWorkingPrecision", "productionRound",
            "fitExtraOrder", "fitOrderIncrement", "fitMaximumRounds",
        }
    expected = {
        "production_plan": production_base,
        "fit": common | {
            "leadingPower", "workingPrecisionDigits", "points", "values",
            "validationPoints", "validationValues", "validationTolerance",
        },
    }.get(action)
    if expected is None:
        raise ValueError(f"unsupported ep series control action: {action}")
    if action == "production_plan":
        optional = {
            "samplePoints", "sampleAngleRange", "validationPoints",
            "initialInternalMaximumPower"
        }
        actual = set(value)
        if not production_base <= actual or not actual <= production_base | optional:
            raise ValueError(
                "ep series production_plan request has unexpected or missing fields"
            )
        data = dict(value)
    else:
        data = _exact_keys(value, expected, f"ep series {action} request")
    _string(data["backendPackagePath"], "backendPackagePath")
    _integer(data["maximumPower"], "maximumPower", -10**9)
    _integer(data["goalDigits"], "goalDigits")
    return data


def _fraction_text(value: Fraction) -> str:
    """把自动网格点写成 Wolfram 和 FLINT 都能精确读取的有理数字符串。"""

    return f"{value.numerator}/{value.denominator}"


def _wolfram_sample_text(value: Any, digits: int) -> str:
    """把实点或 Acb 复点写成 MadStree 可精确拉回的 Gaussian rational。"""

    from flint import acb, fmpq

    if isinstance(value, fmpq):
        return str(value)
    point = acb(value)
    real = _fraction_text(Fraction(point.real.mid().str(digits, radius=False)))
    imag = _fraction_text(Fraction(point.imag.mid().str(digits, radius=False)))
    return f"({real})+({imag})*I"


def _sample_text_acb(value: Any) -> Any:
    """读取 adapter 自己生成的 Wolfram 复数文本或原有实标量文本。"""

    from flint import acb

    if isinstance(value, str):
        match = re.fullmatch(r"\(([^()]*)\)\+\(([^()]*)\)\*I", value)
        if match is not None:
            return acb(
                match.group(1).replace("*^", "e"),
                match.group(2).replace("*^", "e"),
            )
    return acb(value)


def _acb_from_ball_record(record: Any, label: str) -> Any:
    """从 MadStree 求值结果恢复完整 Acb 球，不使用仅供显示的中点。"""

    from flint import acb
    from flintnde.core import arb_ball_from_json

    item = _exact_keys(
        record, {"real", "imag", "real_ball", "imag_ball"}, label
    )
    return acb(
        arb_ball_from_json(item["real_ball"], f"{label}.real_ball"),
        arb_ball_from_json(item["imag_ball"], f"{label}.imag_ball"),
    )


def _series_vectors(records: Any, label: str) -> tuple[Any, ...]:
    """恢复一组同维 Acb 列矢量。"""

    from flintnde.core import column_vector

    if not isinstance(records, list) or not records:
        raise ValueError(f"{label} must be a nonempty list")
    vectors = []
    for sample_index, vector in enumerate(records):
        if not isinstance(vector, list) or not vector:
            raise ValueError(f"{label}[{sample_index}] must be a nonempty vector")
        vectors.append(column_vector([
            _acb_from_ball_record(item, f"{label}[{sample_index}][{row}]")
            for row, item in enumerate(vector)
        ]))
    dimension = vectors[0].nrows()
    if any(vector.nrows() != dimension for vector in vectors):
        raise ValueError(f"{label} vectors have inconsistent dimensions")
    return tuple(vectors)


def _run_series_control(data: dict[str, Any]) -> dict[str, Any]:
    """按 MadStree 符号证书规划 ep 网格，或重构并验证 Laurent 系数。"""

    _resolve_backend_path(data)
    from flint import acb
    from flintnde import configure_working_precision, fit_sampled_series
    from flintnde.regularization import (
        AUTOMATIC,
        _resolve_plan,
    )

    action = data["action"]
    maximum_power = data["maximumPower"]
    goal_digits = data["goalDigits"]
    if action == "production_plan":
        leading = _integer(data["leadingPower"], "leadingPower", -10**9)
        production_round = _integer(data["productionRound"], "productionRound")
        fit_extra_order = _integer(data["fitExtraOrder"], "fitExtraOrder", 0)
        fit_order_increment = _integer(
            data["fitOrderIncrement"], "fitOrderIncrement", 1
        )
        fit_maximum_rounds = _integer(
            data["fitMaximumRounds"], "fitMaximumRounds", 1
        )
        if production_round > fit_maximum_rounds:
            raise ValueError("productionRound exceeds fitMaximumRounds")
        pole_depth = max(0, -leading)
        reconstruction_depth = maximum_power - leading
        automatic_initial_count = max(
            math.ceil(2.5 * reconstruction_depth + pole_depth),
            reconstruction_depth + 1 + fit_extra_order,
        )
        initial_internal_maximum = data.get("initialInternalMaximumPower")
        if initial_internal_maximum is None:
            initial_count = automatic_initial_count
        else:
            initial_internal_maximum = _integer(
                initial_internal_maximum, "initialInternalMaximumPower", maximum_power
            )
            initial_count = initial_internal_maximum - leading + 1
        target_count = initial_count + fit_order_increment * (production_round - 1)
        explicit_points = data.get("samplePoints")
        explicit_validation_points = data.get("validationPoints")
        if explicit_points is not None:
            if not isinstance(explicit_points, list) or not explicit_points:
                raise ValueError("samplePoints must be a nonempty list")
            capacity_count = len(explicit_points)
            if initial_internal_maximum is not None and initial_count > capacity_count:
                raise ValueError(
                    f"samplePoints has {capacity_count} candidates but first-round internal "
                    f"power {initial_internal_maximum} requires {initial_count}"
                )
            initial_count = min(initial_count, capacity_count)
            if initial_count < reconstruction_depth + 1:
                raise ValueError("samplePoints do not cover the requested regulator powers")
            target_count = min(
                capacity_count,
                initial_count + fit_order_increment * (production_round - 1),
            )
        else:
            capacity_count = initial_count + fit_order_increment * (fit_maximum_rounds - 1)
        maximum_samples = _integer(data["maximumSamples"], "maximumSamples")
        if target_count > maximum_samples:
            raise ValueError(
                f"incremental production needs {target_count} samples, "
                f"above maximumSamples={maximum_samples}"
            )
        capacity_count = min(capacity_count, maximum_samples)
        alpha = Decimal(pole_depth) / Decimal(4) + Decimal(goal_digits) / Decimal(
            reconstruction_depth + 1
        )
        plan = _resolve_plan(
            leading_power=leading, maximum_power=maximum_power,
            goal_digits=goal_digits,
            sample_points=(explicit_points if explicit_points is not None else AUTOMATIC),
            sample_angle_range=data.get("sampleAngleRange", AUTOMATIC),
            sample_count=capacity_count,
            base_sample=AUTOMATIC, sample_spacing=data["sampleSpacing"],
            working_precision_digits=AUTOMATIC,
            extra_working_precision=float(data["extraWorkingPrecision"]),
            transport_order=AUTOMATIC, transport_extra_order=AUTOMATIC,
            transport_sample_count=AUTOMATIC, transport_extra_sample_count=AUTOMATIC,
            validation_sample_count=_integer(
                data["validationSampleCount"], "validationSampleCount"
            ),
            validation_points=(
                explicit_validation_points
                if explicit_validation_points is not None
                else AUTOMATIC
            ),
            validation_scale=data["validationScale"],
            maximum_samples=maximum_samples,
            rationalize_sample_points=True,
        )
        return {
            "status": "success", "schema": data["schema"], "action": action,
            "points": [
                _wolfram_sample_text(item, plan.working_precision_digits)
                for item in plan.sample_arguments[:target_count]
            ],
            "validationPoints": [
                _wolfram_sample_text(item, plan.working_precision_digits)
                for item in plan.validation_arguments
            ],
            "sampleCount": target_count,
            "initialSampleCount": initial_count,
            "capacitySampleCount": capacity_count,
            "unusedCandidateCount": capacity_count - target_count,
            "internalMaximumPower": leading + target_count - 1,
            "workingPrecisionDigits": plan.working_precision_digits,
            "primaryOrder": plan.transport_order,
            "referenceOrder": plan.transport_order + plan.transport_extra_order,
            "targetRelativeError": f"1e-{goal_digits}",
            "baseSample": plan.base_sample,
            "alphaEpsilon": plan.alpha_epsilon,
            "productionRound": production_round,
            "sampleSource": plan.source,
        }
    leading = _integer(data["leadingPower"], "leadingPower", -10**9)
    precision = _integer(data["workingPrecisionDigits"], "workingPrecisionDigits")
    configure_working_precision(precision)
    values = _series_vectors(data["values"], "values")
    validation_values = _series_vectors(data["validationValues"], "validationValues")
    internal_maximum = leading + len(values) - 1
    result = fit_sampled_series(
        sample_points=[_sample_text_acb(value) for value in data["points"]],
        sample_values=values,
        maximum_power=internal_maximum, leading_power=leading,
        validation_points=[
            _sample_text_acb(value) for value in data["validationPoints"]
        ],
        validation_values=validation_values,
        validation_tolerance=None, series_parameter="ep",
    )
    tolerance = acb(data["validationTolerance"]).real
    maximum_residual = acb(
        result.diagnostics["maximum_validation_relative_residual"]
    ).real
    coefficients = {
        str(power): [_acb_record(vector[row, 0], precision) for row in range(vector.nrows())]
        for power, vector in zip(result.powers, result.coefficients)
        if power <= maximum_power
    }
    if not maximum_residual < tolerance:
        return {
            "status": "success", "schema": data["schema"], "action": action,
            "fitStatus": "retry",
            "reason": (
                f"series validation residual {maximum_residual.str(20)} is not "
                f"below {tolerance.str(20)}"
            ),
            "internalMaximumPower": internal_maximum,
            "leadingPower": result.leading_power,
            "maximumPower": result.maximum_power,
            "coefficients": coefficients,
            "diagnostics": result.diagnostics,
            "effectiveParameters": result.effective_parameters,
        }
    return {
        "status": "success", "schema": data["schema"], "action": action,
        "fitStatus": "accepted",
        "leadingPower": result.leading_power,
        "maximumPower": result.maximum_power,
        "internalMaximumPower": internal_maximum,
        "coefficients": coefficients,
        "diagnostics": result.diagnostics,
        "effectiveParameters": result.effective_parameters,
    }


def main() -> int:
    """命令行入口；失败也以 UTF-8 JSON 返回完整原因。"""

    if len(sys.argv) == 1:
        print(_citation_text())
        print("usage: flintnde_transport.py INPUT.json OUTPUT.json")
        return 0
    if len(sys.argv) != 3:
        raise SystemExit("usage: flintnde_transport.py INPUT.json OUTPUT.json")
    # 文件协议把 MadStree 引用写入响应 metadata，且不让依赖包的 import 提示进入日志。
    os.environ["FLINTNDE_SUPPRESS_CITATION_NOTICE"] = "1"
    input_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])
    try:
        raw = json.loads(input_path.read_text(encoding="utf-8"))
        if raw.get("schema") == EVALUATE_BATCH_SCHEMA:
            result = _run_batch(_validate_batch_request(raw))
        elif raw.get("schema") == SERIES_CONTROL_SCHEMA:
            result = _run_series_control(_validate_series_control_request(raw))
        else:
            result = _run(_validate_request(raw))
        exit_code = 0
    except Exception as error:  # noqa: BLE001 - CLI 边界必须序列化所有错误
        result = {
            "status": "failed", "errorType": type(error).__name__,
            "error": str(error), "traceback": traceback.format_exc(),
        }
        exit_code = 1
    result = _with_citation_metadata(result)
    output_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
