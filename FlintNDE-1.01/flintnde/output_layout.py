"""为调用 FlintNDE 的脚本建立稳定、可审计的结果目录。

调用者显式传入 ``__file__``；本模块始终在该脚本旁的 ``results/<run_name>/`` 下
组织输出，不依赖进程当前工作目录，也不向 package 安装目录写文件。固定分类目录
按需创建，路径检查阻止绝对文件名、跨目录文件名和未约定分类。
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


OUTPUT_CATEGORIES = (
    "configuration",
    "singularities",
    "frobenius",
    "transport",
    "regularization",
    "summary",
)


def validate_single_path_name(value: str, field_name: str) -> str:
    """验证单层目录名或文件名，禁止绝对路径和父目录跳转。"""

    candidate = Path(value)
    if not value or value in {".", ".."} or candidate.is_absolute() or candidate.name != value:
        raise ValueError(f"{field_name} must be one nonempty path name")
    return value


@dataclass(frozen=True)
class OutputLayout:
    """保存调用脚本、运行名和标准结果根目录，并提供受限的分类路径接口。"""

    caller_script: Path
    run_name: str
    results_root: Path
    run_root: Path

    def directory(self, category: str) -> Path:
        """返回并按需建立一个固定用途分类目录。"""

        if category not in OUTPUT_CATEGORIES:
            allowed = ", ".join(OUTPUT_CATEGORIES)
            raise ValueError(f"unknown output category {category!r}; choose one of: {allowed}")
        target = self.run_root / category
        target.mkdir(parents=True, exist_ok=True)
        return target

    def file(self, category: str, filename: str) -> Path:
        """返回分类目录中的单层文件路径；不创建空文件。"""

        safe_filename = validate_single_path_name(filename, "filename")
        return self.directory(category) / safe_filename

    def write_json(self, category: str, filename: str, payload: Any) -> Path:
        """以 UTF-8、稳定缩进把 JSON 可序列化对象写入指定分类并返回路径。"""

        target = self.file(category, filename)
        target.write_text(
            json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        return target


def initialize_output_layout(
    caller_script: str | Path,
    *,
    run_name: str | None = None,
) -> OutputLayout:
    """在调用脚本旁初始化 ``results/<run_name>/`` 并写入布局说明。

    ``caller_script`` 应直接传调用脚本的 ``__file__``。``run_name`` 缺省取脚本 stem；
    用户可为同一脚本的不同物理配置提供直观单层名称，但不能传任意路径。
    """

    script_path = Path(caller_script).absolute()
    if not script_path.is_file():
        raise ValueError(f"caller_script must identify an existing file: {script_path}")
    resolved_run_name = validate_single_path_name(run_name or script_path.stem, "run_name")
    results_root = script_path.parent / "results"
    run_root = results_root / resolved_run_name
    layout = OutputLayout(script_path, resolved_run_name, results_root, run_root)
    layout.write_json(
        "configuration",
        "output_layout.json",
        {
            "schema": "flintnde_output_layout_v1",
            "caller_script": script_path.name,
            "run_name": resolved_run_name,
            "results_directory": str(Path("results") / resolved_run_name),
            "categories": list(OUTPUT_CATEGORIES),
        },
    )
    return layout
