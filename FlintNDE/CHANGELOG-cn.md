# FlintNDE 更新日志

## 1.01

新增独立 C++/FLINT 引擎 `engine/`（MSYS2 UCRT64 中 `make` 构建出 `flintnde_engine`），在同一文件协议合同下执行 prefix 认证双链输运：双链同精度、只差阶数，参考链级数的前 primaryOrder+1 个系数即主链级数。引擎输出与 Python 侧逐位一致；Python 包 `flintnde` 的公开接口无任何增删。

新增留数内核选项 `--residues auto|sparse|dense`，缺省 `auto` 在留数矩阵非零密度不超过 0.25 时自动选稀疏跳乘内核；三种内核同精度同算术序列，输出逐位一致。

批量任务改为单进程 `std::thread` 线程池并行，各任务流式解析（`--workers` 缺省回落到请求的 `parallelTaskCount`），大请求 JSON 的解析墙钟大幅下降（基准 100 任务中位 56.8 s 降至 11.8 s）。

新增 Wolfram 直连入口 `FlintNDEEngineEvaluateEpBatch`：固定 regulator 任务直接交给引擎执行，不经 Python bridge；成功任务的 `finalValues` 与逐点值解码为 Complex。

依赖变化：Python 依赖不变，仍为 `python-flint>=0.6` 与 `sympy>=1.12`；引擎为可选组件，构建新增 MSYS2 UCRT64 工具链与 FLINT 3.x 开发库依赖，Windows 运行时需 `libflint` 运行库位于 PATH 或由 `"EngineLibraryPath"` 显式指定。

自本版起 FlintNDE 独立演进，不再与 MadStree Vendor 逐文件同步；MadStree 继续使用其 Vendor 内冻结的 FlintNDE v1.0。
