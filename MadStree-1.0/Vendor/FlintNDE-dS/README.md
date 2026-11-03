# FlintNDE Python package

发布名使用 `FlintNDE`，Python 导入名使用 `flintnde`。

`engine/` 是 C++ 认证双链引擎（本目录随 v0.21 一次性引入）：基线取独立程序包
package-FlintNDE `v1.02` 引擎（schema `flintnde_engine_*`、十进制 letter 支持），此后按
MadStree 请求合同独立扩展（`tensorSectors` 结构化 letters、embedded 认证模式等），不与独立
程序包保持同步。构建产物 `engine/build/flintnde_engine.exe` 不随源码分发，由适配器在运行时
探测使用；缺失时适配器回退 Python 认证双链。构建与契约详见「C++ 引擎」一节。

每个 Python 进程首次 `import flintnde` 时在 stderr 显示一次简洁引用提醒：
`FlintNDE package paper, arXiv identifier pending.`。该提醒不列出其它程序包或论文。

## 版本更新说明

每次修改 `pyproject.toml` 中的 FlintNDE 版本号，都必须在独立 FlintNDE 项目包根的
`更新日志.md` 顶部追加该版本增量。更新说明至少包含基线版本、新增与
修复内容、公开接口或数值 convention 变化、迁移要求、已执行验证和已知限制；条目中的版本
必须与 `pyproject.toml` 一致。Vendor 交付不复制该日志。当前 `0.1.0.dev0` 及此前已有资产
豁免，不追溯补建。建议在独立版本 branch 中完成开发和验证，但是否创建或合并 branch
完全由用户决定。

## 当前路径与多点求值合同

### 奇点目标

`evaluate_singular_target` 复用认证的局部幂对数基，在用户明确命中的有限奇点处分析
解的局部行为。可去奇点返回 Acb 数值；真实 pole 和零阶 log 发散分量返回文本
`Infinity`，并在 `SingularTargetResult.report` 中给出逐分量分类、局部阶数和样本诊断。
该接口不调用 NIntegrate，也不把绕行点值冒充奇点处值。目标位置可以是任意 exact
$\mathbb Q(i)$ 点；分类只读取完整有理矩阵的局部解，不读取变量名、零坐标或上游图信息。
专项回归把目标平移到 $2+i$，同时保留另一个有限奇点 $5+i$，并在同一个二维系统中得到
可去与真实 pole 两种分量，防止实现退化为零点或单一 letter 特判。

`plan_singular_target_match` 先以目标奇点到最近其它奇点的距离确定局部收敛圆；当前普通点
在保守圆内时直接复用，否则沿来向生成域内隐藏普通匹配点。调用方把低阶主链和高阶参考链
分别输运到同一点，再用 `evaluate_singular_target` 分别求值：主链是用户结果，参考链只检查
逐分量分类、发散状态和有限值误差。

当前实现仍面向通用单变量矩阵方程 `Y'(x)=A(x)Y(x)`。对 exact
`RationalMatrixSystem`，有限奇点、极点阶数和无穷远分类全部由包内部从矩阵元分母发现；
用户不需要先给出 dlog letters 或奇点位置。程序随后自动验证

```text
A(x) = P(x) + Sum_j R_j/(x-p_j),
```

其中 `P(x)` 可以是任意有限次数矩阵多项式。只有恒等重构、有限极点均为一阶等门禁全部
通过时，才调用多项式加简单极点的快速递推；高阶极点、一般有理矩阵和解析回调继续使用
各自的通用路线或按现有能力 fail closed。因此快速路线是内部结构特化，不是 dlog-only
公开输入合同。纯多项式系统对应空的有限极点集合，同样是合法的已认证情形。

术语上，任何经过中途节点的多点输运都可称为“折跃”；只有显式穿过奇点并用局部基连接
两侧匹配点时才称为“奇点折跃”。Python 模式只接受 `"avoid"` / `"singularity_jump"`；
Wolfram 模式只接受 `"Avoid"` / `"SingularityJump"`。

路径生成和执行是两个既有阶段，不另设 `AutoPlan` 一类开关：

- `plan_transport_path(...)` 或 `FlintNDEPlanPath[...]` 接收同一单变量复参数平面中的
  原始点序列并产生完整计划；当前节点收敛圆盘内的多个复点共享该节点的解系数，不要求
  这些点位于同一条实线段上。不同单变量拉回应分别规划，不能跨拉回共享局部系数；
- `transport_planned_path_refined(...)` 或 `FlintNDEExecutePath[...]` 只恢复并执行该计划，
  不再次调用规划器；
- 规划节点覆盖的多个用户点按节点/线段组成 evaluation bucket。点数不少于 8 的桶由
  `acb_poly.evaluate(..., algorithm="fast")` 使用子积树/余数树快速多点求值；小桶用
  iterative 算法，输出记录实际算法与点数；
- `direct_user_point_path(...)` 是关闭规划的公开入口：严格使用 `start` 后的用户点序列作为
  节点，不插点、不删点、不调用规划器；穿过奇点时 fail closed；
- 缺省为避开奇点模式。直线段命中内部奇点时返回奇点及对应线段诊断；只有显式选择
  `singularity_mode="singularity_jump"` 或 `"SingularityMode" -> "SingularityJump"`
  才允许奇点折跃。奇点折跃所选多值分支等价于某条绕行路径，必须由用户自行确认；
- `message_language="EN"|"CN"` 和 Mathematica 的
  `MessageLanguage -> "EN"|"CN"` 控制简洁提示，缺省英文。

工作精度由 `WorkingPrecisionDigits`（Python 为相应 `decimal_digits`）决定。两种接口缺省
都是 200 位：`import flintnde` 立即建立该全局 Acb 精度，`configure_working_precision()`
可重置为同一缺省；显式传入其它正整数则按用户值覆盖：

```text
working bits = ceil(WorkingPrecisionDigits log2(10)) + 32.
```

缺省 200 位请求使用 697 bit；例如显式 70、100 位请求分别使用 265、365 bit，更高十进制
精度会按同一公式继续提高 bit 数。
`OutputDigits` 只控制输出格式，不会提高内部精度。线段投影、匹配比例、旋转因子和
winding/monodromy 均在当前 Arb 精度构造，不经过 binary64 几何。两类序列化路径都保存
`planning_precision_digits`，节点和奇点折跃几何保存 Arb 中点--半径--指数；若执行请求
高于规划精度，程序拒绝执行并要求按目标精度重新规划，不能声称从低精度 JSON 补回信息。

Wolfram Language 入口已按普通程序包组织。把当前版本根加入 `$Path` 后可直接使用：

```wl
Needs["FlintNDE`"];

system = FlintNDERationalSystem[{{1/(x - 1)}}, x];
plan = FlintNDEPlanPath[
  system, 0, {1/2},
  "WorkingPrecisionDigits" -> 100,
  MessageLanguage -> "EN"
];
result = FlintNDEExecutePath[
  system, {1}, plan,
  "WorkingPrecisionDigits" -> 100
];
```

`"WorkDirectory"` 表示临时运行根本身。`Automatic` 使用当前目录下唯一的 `results_temp/`，
bridge 文件位于其 `bridge/` 子目录，文件名为短 token。Windows 下先检查所有完整路径；超过
传统 Win32/Wolfram 安全上限时返回 `RuntimePathTooLong`，不会继续启动 Python。
输入写入失败、缺少 `python-flint`、Python bridge 启动失败和零退出但无输出分别返回
`RuntimeInputWriteFailed`、`PythonFlintUnavailable`、`BridgeLaunchFailed` 和
`BridgeOutputMissing`。
bridge 使用参数列表 `RunProcess` 和显式 `ProcessDirectory`，stdout/stderr 由 Wolfram 捕获后写入
短日志；不经过 shell 命令拼接、重定向或旧 launcher fallback。

显式已知的 `P(x)+Sum R_j/(x-p_j)` 也可用 `FlintNDEPartialFractionSystem` 构造；这不是
要求用户预提取奇点，`FlintNDERationalSystem` 仍会自动发现并选择路线。

首版公开接口分为以下部分：

- `initialize_output_layout` / `OutputLayout`：在调用脚本旁建立分层结果目录；
- `GaussianRational` / `gaussian_rational`：不增加依赖的 Q(i) 精确复数公开输入；
- `RationalMatrixSystem` / `analyze_singularities`：exact Gaussian-rational 矩阵与完整奇点清单；
- `NamedPoint` / `prepare_local_expansion` / `build_adaptive_path`：统一局部调度与自描述可执行路径；
- `plan_transport_path` / `direct_user_point_path`：自动规划或严格用户节点的两种路径构造；
- `plan_singular_target_match` / `evaluate_singular_target`：末端奇点的域内隐藏匹配计划、
  有限/pole/log 分类和可序列化结果；
- `AnalyticMatrixSystem`：接收任意返回 FLINT `acb_mat` 的解析矩阵函数；
- `transport_path` / `transport_path_refined`：普通点输运与统一奇点局部基 bridge；
- `transport_frobenius_boundaries_refined`：共享局部基和普通段 Taylor 矩阵的多列 Frobenius
  初值输运，并逐列报告 refinement；
- `RegularSingularSystem` / `build_frobenius_manifest` / `build_power_log_basis`：
  exact-first 的正则奇点局部解；
- `frobenius_boundary` / `FrobeniusBoundary`：验证奇点起点的 `{a,b,C}` 领头分支；
- `exponential_boundary` / `ExponentialBoundary`：验证指数型奇点的 `{phi,a,b,C}` 可复用边界；
- `attempt_fuchsian_reduction` / `build_local_solution_basis`：高阶 pole 的 exact Lee--Moser balance
  与认证指数 sector；
- `five_term_tail_diagnostic` / `FiveTermTailDiagnostic`：对 caller 给出的标量渐近项执行可读取的
  前后五项严格阈值门禁；
- `NumericalRegularSingularSystem` / `NumericalFrobeniusOptions`：没有 exact residue 时的
  精度感知结构判别；
- `reconstruct_series_solution` / `SeriesReconstructionResult`：复用完整 NDE 输运，重构
  正规化参数的幂级数解。
- `fit_sampled_series`：对调用方已经计算的标量或列矢量样本作 Laurent/Taylor 方阵重构；
  源采样标签可与实际代入的正规化数值点分开指定。
- `flintnde.regular_point_de`：固定 regulator 后的 pole/residue 系统、Watson 边界、
  普通点输运、认证的 k=0 Frobenius 局部解与 finite part。
- `flintnde.epsilon_jet_de`：解析 epsilon jet 的普通点输运、联合 epsilon/k 局部解与
  finite-part 系数，不使用浮点 epsilon 采样拟合。

普通点默认通过 Cauchy 圆周采样和离散 Fourier 投影恢复矩阵 Taylor 系数。该接口
只要求采样圆内解析，不要求矩阵预先拆成常数项与单极点 residue。

调用脚本用 `initialize_output_layout(__file__, run_name=None)` 初始化输出。缺省目录为
`<脚本目录>/results/<脚本stem>/`；`run_name` 只能是单层直观名称。固定分类为
`configuration`、`singularities`、`frobenius`、`transport`、`regularization`、`summary`，其中只有
`configuration/output_layout.json` 在初始化时创建，其余目录按实际输出懒创建。

## 安装与调用位置

用户环境要求 Python 3.10 或更高版本、`python-flint>=0.6` 和 `sympy>=1.12`。三种安装方式为：

```powershell
# 安装发布的 wheel；用户无需自行构建
python -m pip install .\flintnde-1.0-py3-none-any.whl

# 从下载的源码普通安装
python -m pip install "path\to\package-FlintNDE\versions\FlintNDE-0.5.0"

# 开发者可编辑安装
python -m pip install -e "path\to\package-FlintNDE\versions\FlintNDE-0.5.0"
```

`pip` 会自动安装满足版本要求的 `python-flint` 与 `sympy`。安装一次后，可在任意目录 `A` 的脚本中
直接 `import flintnde`。不需要复制源码；输出位置也不取决于 package 安装目录或启动命令的
当前目录，而只取决于传给 `initialize_output_layout` 的调用脚本路径。配置文件本身不会自动
触发计算，`A` 中仍需一个 Python 调用脚本读取配置并调用公开接口。

## 正规化参数的幂级数解重构

最小调用给出完整 NDE 问题、用户需要返回的最高 regulator 幂，以及上游从符号边界和
DE 得到的最低阶证书：

```python
from flintnde import reconstruct_series_solution

series_result = reconstruct_series_solution(
    DEmatrix=DEmatrix,
    boundary=boundary,
    path=path,
    maximum_power=2,
    leading_power=-1,
    leading_power_certificate={
        "status": "certified",
        "leading_power": -1,
        "method": "symbolic-boundary-and-DE",
    },
)
```

若矩阵、边界或奇点位置依赖 `ep`，用工厂保持每个自动生成的有理样本在构造矩阵时仍是
exact rational：

```python
def DEmatrix(ep):
    # 自动 rationalization 时 ep 是 FLINT fmpq。
    return build_system_at(ep)

def boundary(ep):
    return build_ordinary_initial_vector_at(ep)

def path(ep, system_at_ep):
    return build_path_at(ep, system_at_ep)
```

固定的 `AnalyticMatrixSystem`/`RationalMatrixSystem`、边界和路径也可直接传入。
普通点边界是 `acb_mat` 或标量列表；正则奇点边界使用下文的 `frobenius_boundary`。
`boundary(ep)` 可以按每个 regulator 样本返回其中任一种格式。

完整可选项为：

```python
series_result = reconstruct_series_solution(
    DEmatrix=DEmatrix,
    boundary=boundary,
    path=path,
    maximum_power=2,
    leading_power=-1,
    leading_power_certificate={
        "status": "certified",
        "leading_power": -1,
        "method": "symbolic-boundary-and-DE",
    },
    series_parameter="ep",
    goal_digits=30,
    sample_points="automatic",
    sample_angle_range="automatic",
    sample_count="automatic",
    base_sample="automatic",
    sample_spacing=0.01,
    working_precision_digits="automatic",
    extra_working_precision=0.0,
    transport_order="automatic",
    transport_extra_order="automatic",
    transport_sample_count="automatic",
    transport_extra_sample_count="automatic",
    radius_fraction=0.60,
    guard_bits=32,
    validation_sample_count=2,
    validation_points="automatic",
    validation_scale=0.5,
    validation_tolerance="automatic",
    fit_extra_order=2,
    initial_internal_maximum_power="automatic",
    fit_order_increment=2,
    fit_max_rounds=3,
    maximum_samples=100,
    rationalize_sample_points=True,
    parallel_task_count=12,
    output_layout=None,
    result_name="series_reconstruction",
)
```

`parallel_task_count` 控制同时运行的独立固定 `ep` NDE 任务，缺省为 `12`。程序实际启动
`min(parallel_task_count, ep任务数)` 个 worker；任务更多时，任一 worker 完成后自动续交
下一个任务，不需要用户手动分批。因为 python-flint 的精度与线程 context 是进程级全局状态，
这里使用独立进程而不是 Python 线程。并行模式下 `DEmatrix(ep)`、`boundary(ep)` 和
`path(ep,system)` 应定义为模块顶层函数，并在 worker 内构造 FLINT 对象。固定
`RationalMatrixSystem`、普通边界向量和 `AdaptivePath` 会由程序结构化序列化并在 worker
恢复；任意 `AnalyticMatrixSystem` 实例或局部闭包不能可靠跨 Windows worker 传输，需改为
模块顶层 factory，或在调试时显式设置 `parallel_task_count=1`。该选项与 `ctx.threads` 不同：前者并行不同 `ep` 任务，后者只给
单个进程内部分 FLINT 算术内核设置最大线程数。

也可直接对任意模块顶层固定 `ep` 任务使用：

```python
from flintnde import run_ep_tasks

batch = run_ep_tasks(ep_values, solve_ep, parallel_task_count=12)
```

完整可运行示例见仓库级 `examples/ep_parallel.py`。

`maximum_power=m` 与整数 `leading_power=n_min` 都是必填参数。FlintNDE 只消费矩阵 NDE，
任意 Python callable 不携带足以证明 Laurent 支撑的符号结构，因此不再用数值样本猜测
最低阶。调用方必须同时提供严格证书
`{"status":"certified","leading_power":n_min,"method":"..."}`；字段不全、状态错误、
`automatic` 或证书整数不一致都在任何 NDE 样本前 fail closed。MadStree 会从自己的实际
符号边界条件和 dlog DE 自动生成证书，MadStree 用户无需手填。随后定义

```text
P = max(0, -n_min)
K = m - n_min
N_fit = max(ceil(5 K / 2 + P), K + 1)
alpha_ep = P / 4 + goal_digits / (K + 1)
ep0 = 10^(-alpha_ep)
p0 = max(ceil((N_fit + P) alpha_ep), 30)
working_precision_digits = max(200, ceil(2 (1 + extra_working_precision) p0))
```

其中 200 是缺省自动规划下限；若用户显式给出 `working_precision_digits`，则直接采用用户值。
自动路线首轮使用全部 `N_fit` 个生产点做 Acb 方阵插值，内部拟合到
`n_min + N_fit - 1`，只把 `n_min..maximum_power` 返回用户；不使用最小二乘或伪逆。
`transport_order=4*p0`，`transport_extra_order=ceil(max(50,p0/5))`，后者作为参考链增加的
阶数。`transport_sample_count` 和 `transport_extra_sample_count` 分别覆盖两条 Cauchy--DFT
链的圆周点数；`"automatic"` 继续使用基础输运的 `max(32,2*order)`。

显式 `sample_points` 是按给定顺序消费的生产候选池，而不是要求首轮全部使用。若显式给出
`initial_internal_maximum_power=q`，首轮恰取前 `q-n_min+1` 个候选点；缺省仍按上述自动公式，
但候选池较短且足以覆盖用户要求时使用整个候选池。验证失败后每轮按
`fit_order_increment` 从剩余候选中增加点并复用所有旧值；候选池耗尽时绝不生成池外取值，
但会返回当前最佳系数并设置 `precision_target_met=False`、
`precision_failure_reason="candidate_pool_exhausted"`。`sample_count` 与显式候选池同时给出时
必须等于候选总数。

`sample_angle_range=(theta_min,theta_max)` 以弧度指定自动复采样的开角区间，不能与显式
`sample_points` 同时使用。程序在区间内部均匀选择最多三条射线并循环分配生产点，不取边界、
也不刻意贴近边界；模长仍为原自动 `ep0*(1+j*sample_spacing)`，`ep0` 继续由 pole 深度、
目标位数和拟合阶数决定，没有用户模长上限参数。自动验证点保持对应角度，只把模长乘
`validation_scale`。缺省 `"automatic"` 严格保留原正实轴网格。普通 Python 浮点可用于角度
边界控制；需要更高边界精度时传十进制字符串。

独立验证点缺省位于生产网格的 `validation_scale` 倍尺度；也可用字符串、`fmpq` 或 `acb`
序列显式传入 `validation_points`。验证点完全不参与插值；NDE 主/参考链以及
返回截断级数都以 `validation_tolerance` 检查，缺省为 `10^(-goal_digits)`。已有方阵拟合结果
但精度仍未达标时，程序发出 `RuntimeWarning` 并返回当前最佳系数；调用方必须检查
`effective_parameters["precision_target_met"]` 或同名 diagnostics 字段。候选池不足、达到
`maximum_samples` 和达到轮数上限分别有结构化失败原因。显式 `sample_points` 应使用字符串、`fmpq` 或 `acb`，不接受会掩盖
输入精度的 Python `float/complex`；显式 `validation_points` 必须与整个生产候选池分离。若要严格限制 regulator
取值范围，应同时显式给出 `sample_points` 和 `validation_points`。`rationalize_sample_points=False` 会让自动生产和
验证点以 Acb 而非 `fmpq` 传给工厂。

结果的 `powers` 与 `coefficients` 一一对应；`coefficient(power)` 取单阶系数矢量，
`evaluate(ep)` 计算返回的截断级数。`sample_values`、`validation_values`、
`effective_parameters` 和 `diagnostics` 保存生产值、验证值、实际生效配置、每次 NDE refinement
与残差。提供 `output_layout` 时，JSON 写入调用脚本旁
`results/<run>/regularization/<result_name>_series_reconstruction.json`。

若目标量不是单条 NDE 的终点矢量，而是调用方已组合出的 Wronskian 或其它标量，可直接
复用相同的 Acb 方阵拟合层。下面以实轴奇点上方的 lateral 样本为例：`sample_points`
保存正的 `delta` 标签，`regularization_points` 才是级数变量 `+delta`。下方 lateral
使用同一组正标签，并把正规化点改为 `-delta`。

```python
from flintnde import fit_sampled_series

fit_plus = fit_sampled_series(
    sample_points=("0.006", "0.0041", "0.0027", "0.0013"),
    sample_values=(w_plus_1, w_plus_2, w_plus_3, w_plus_4),
    regularization_points=("0.006", "0.0041", "0.0027", "0.0013"),
    maximum_power=3,
    leading_power=0,
    validation_points=("0.0008", "0.00045"),
    validation_values=(w_plus_validation_1, w_plus_validation_2),
    validation_regularization_points=("0.0008", "0.00045"),
    validation_tolerance="1e-20",
    series_parameter="signed_iw_delta",
)
w_plus_limit = fit_plus.coefficient(0)[0, 0]
```

`regularization_points` 与 `validation_regularization_points` 都可缺省；缺省时分别直接
使用 `sample_points` 与 `validation_points`，输入格式、系数和摘要字段保持原接口行为。
显式提供时，两组点必须分别一一对应、非零且互异，验证正规化点还必须与生产正规化点
分离。该接口不重算 NDE，也不改变 FLINT 全局精度；全部生产正规化点一次进入方阵
Vandermonde 求解，不使用最小二乘或伪逆，验证点只用于残差和 fail-closed 门禁。

## 包内 regular-point 与 epsilon-jet DE

固定 regulator 的矩阵系统从 `flintnde.regular_point_de` 导入。它接收已经完成
convention 数值代入的 pole/residue payload，提供 Watson 边界、普通点 Taylor 输运、
由 exact manifest 认证的 `k=0` 局部基和有限部读取。解析 regulator jet 使用同一模块的
中性 `regulator` 接口。

```python
import flintnde.regular_point_de as regular_de
import flintnde.epsilon_jet_de as epsilon_de
```

`flintnde.epsilon_jet_de` 在此基础上处理截断的 epsilon 系数矢量。epsilon 只作为解析
jet 逐阶进入矩阵、边界和局部递推；普通点输运及联合 epsilon/k Frobenius 解都不需要
多个浮点 epsilon 样本。两模块均属于 FlintNDE 包，调用方只把包根加入导入路径并核对
实际模块文件位于该包根之下。
## exact Q(i) 输入、奇点清单与统一路径

自动奇点发现不从数值黑箱回归分母。矩阵元应写成 exact rational coefficient records；系数
顺序恒为常数项、一次项、二次项等：

```python
from flintnde import *

layout = initialize_output_layout(__file__, run_name="physical_point_01")
system = RationalMatrixSystem(
    (
        (rational_function([0, "1+I"], ["-2-3*I", 1]), 0),
        (0, rational_function({"real": "1/3", "imag": "-2/5"}, [1, 0, 1])),
    ),
    variable_name="s",
    name="example-system",
)
inventory = analyze_singularities(system, output_layout=layout)
path = build_adaptive_path(
    system,
    NamedPoint("left_boundary", 2),
    NamedPoint("right_match", 3),
    path_name="left_to_right",
    max_step_over_radius=0.45,
    output_layout=layout,
)
local = prepare_local_expansion(
    system,
    NamedPoint("left_boundary", 0),
    order=32,
    radius_fraction=0.60,
    sample_count=96,
    output_layout=layout,
)
```

`numerator` 和 `denominator` 中每个 exact 系数可使用以下任一格式：整数、FLINT `fmpq`、
分数/十进制字符串、`"a+b*I"` 字符串、`(real, imag)` 二元组，或
`{"real": ..., "imag": ...}` 字典。`gaussian_rational(value)` 可单独验证并标准化这些输入，
返回 `GaussianRational`。Python `float`、`complex`、`arb` 和 Acb ball 会被 exact 路线拒绝；
它们只能进入显式 numerical 接口，不能靠提高 FLINT 工作精度恢复已经丢失的输入位数。
由于 `rational_function((c0, c1))` 的 tuple 本身表示按幂次排列的多项式系数，单个复常数若用
二元组，应写成 `rational_function([("a", "b")])` 或先调用
`rational_function(gaussian_rational(("a", "b")))`，避免把实部、虚部误读成常数项和一次项。
若一个矩阵元由多个分式项组成，用
`rational_function(...) + rational_function(...)` 先形成完整矩阵元，程序会 exact 合并约分，
再收集约分后分母。一般 `AnalyticMatrixSystem` 的回调无法从有限采样可靠恢复全部 pole，仍须
显式给 `singularities`。

一个实际 16 维 complex-exact kECep 系统的维护交叉检查脚本位于仓库根相对路径
`package-FlintNDE/check_01_nde_npackage_16x16/check_01_nde_npackage_16x16.py`。它通过通用
`RationalMatrixSystem over Q(i)(k)` 构造共同系统并比较完整端点向量；运行时结果不作为
发布产物保留。当前正式独立验证位于仓库根相对路径 `package-FlintNDE/independent-validation/`。

Q(i)[s] 分母能严格分裂出的奇点会保存 `location_exact`，并可直接进入 exact 局部调度。
若一个 square-free 因子同时含零根和不在 Q(i) 中的代数根，奇点发现会先精确除出 `s`；
因此边界起点 `s=0` 仍保存 `location_exact="0"` 与因子 `("0","1")`，不会因同组其它根
需要 Acb 隔离而丢失 exact 身份。
一般代数根仍可通过 Acb ball 隔离并进入奇点清单，但 `regular_singular_system_at` 当前只在
奇点属于 Q(i) 时构造 exact Frobenius 数据。indicial polynomial 也必须完全分裂于 Q(i)；
否则 exact gate fail closed，且不会自动转用阈值分类。只有用户原始输入本来就是浮点数据时，
才应显式构造 `NumericalRegularSingularSystem`。

奇点 pole 阶由完整矩阵的所有矩阵元共同决定，不能只检查对角元。若任一约分后矩阵元在
某点有二阶或更高 pole，清单记录为 `non_fuchsian_input_basis`。这只说明输入基底不是
Fuchsian；局部调度随后按三层路线处理：

1. 逐轮由最高两阶 Laurent 矩阵构造 exact Lee--Moser projector，应用
   `I=(1-Q+Q/z)J`；每步要求 `Q^2=Q` 且完整矩阵 pole 阶严格下降，并保存 ordered
   projectors、逐步 pole 阶和正逆系统 exact round-trip；
2. 若仍有高阶 pole，尝试在 Q(i) 中 exact 对角化 leading matrix，并只接受所有高阶 Laurent
   矩阵同为对角、不同指数 signature 之间完整有理矩阵元严格为零的 sector；
3. 抽掉每个 sector 的标量指数后，残余系统必须至多为 simple pole，否则 fail closed；

在第二条路线中 `Phi'(z)=sum_{m>=2} c_m z**(-m)`，所以
`Phi(z)=-sum_{m>=2} c_m/(m-1) z**(-(m-1))`。特别地，`A[-2]=k` 给出
`exp(-k/z)`。`k=0` 只表示没有该标量指数，不能把 defective/nilpotent 高阶块变成正则奇点。
4. 若严格解耦失败，但 pole 阶恰为二、`A_-2` 在 Q(i) 中有单重互异谱，则构造
`exp(-k/z) z**rho Sum[v_n z**n]` 的 exact 形式递推。每次固定求到用户指定的 `N` 阶，并额外
生成 `N+1,...,N+5` 五阶用于块比和相对 refinement；仍同时报告最小项位置作为诊断，但不据此
提前截断。该方法仅可作为奇点起点；内部/终点需要 Stokes
connection，因而 fail closed。

当前 Moser 路线在 exact Q(i)(x)、有限局部变量和无 ramification 的 meromorphic gauge
范围内工作。需要重复/defective formal block、代数数域扩张、ramification 或 Stokes
矩阵的系统继续拒绝，也可由用户显式给普通绕行路径。

### 计算前能力预检

`build_adaptive_path_plan(system,start,target,...)` 对路径端点和全部内部奇点执行与正式建路
相同的局部分类，不做数值输运。返回的 `PathPlan.continuation_ready` 是能力门禁；`False`
时 `messages` 和各 segment method 说明是 defective/repeated formal block、ramification、
代数扩域、缺失 Stokes connection 或其它未认证类型。普通 simple pole 也必须通过该门禁：
非 Q(i) 奇点位置以及不受支持的 exact indicial/Jordan/resonance 谱均返回 `False`。此时调用者只能改用普通绕行点，或
自行提供本包外的局部连接数据。正式 `build_adaptive_path` 对同一未认证内部点/终点抛出
`LocalReductionError`，不会退回普通点 Cauchy--DFT，也不会把提高工作精度当作结构证明。

start-only 的单重二阶-pole 形式渐近是唯一例外：它只负责从给定奇点边界初始化到同一
Stokes sector 内首个普通匹配点；相同局部基一旦位于路径内部或终点，`continuation_ready`
即为 `False`。

`exp(-k/z)` 本身不决定剩余级数是否收敛。程序抽掉完整的已认证标量高阶部分后，会重新检查
residual 系统；residual 至多为 simple pole 时，后续是收敛的 Frobenius power-log 级数，
可在最近其他奇点给出的半径内按 `max_step_over_radius` 选匹配点。若 residual 仍有高阶项，
一般只能形成 sectorial/Gevrey 渐近级数。上述单重二阶-pole 路线不提供 Stokes 连接；最近
其他奇点距离只是首匹配点参考尺度，不是收敛半径。自动首匹配点对分支 `k_i` 使用最近指数
根差 `Delta_i=min_{j!=i}|k_i-k_j|` 和晚期项估算 `n_min≈Delta_i/|z|`。缺省
`formal_minimum_order_factor=3`，所以建议距离为
`min_i Delta_i/(3*N)`，令请求阶数 `N` 约为预计最小项阶数的三分之一；仅当
`max_step_over_radius*R` 或路径区间给出更小距离时再缩短步长。其它渐近结构拒绝，不会把离
奇点较远的点当作可靠求值点。

形式级数令 `T_n=v_n z**n`。完成 `N+5` 阶后报告
`r5=norm_inf(sum(T[N+1:N+6]))/norm_inf(sum(T[N-4:N+1]))`，即先做五个矢量项之和再取
无穷范数；同时报告 `norm_inf(S[N+5]-S[N])/norm_inf(S[N+5])`。希望 `r5<1`；若不满足，
程序仍返回结果，但发出醒目 warning 并在 `formal_accuracy_issues` 中记录。即使 `r5<1`，实际
数值也始终保留供用户判断下降是否足够。

当下游已有标量级数项而需要硬门禁时，调用
`five_term_tail_diagnostic(terms,N,threshold="0.1")`。它使用不同且明确的标量定义
`abs(sum(T[N+1:N+6]))/sum(abs(T[N-4:N+1]))`；结果对象公开
`next_five_sum_absolute`、`previous_five_absolute_sum`、`next_over_previous`、`threshold`
和 `passed`。比较为严格 `<`，等于阈值不通过；前五项绝对值之和包含零时直接抛错。

奇点清单写入 `results/<run>/singularities/singularity_inventory.json`；若路径含 `inf`，还会
保存反演变量的 `singularity_inventory_sinv.json`。命名路径写入
`transport/<path_name>_path.json`，逐段记录方法、步长、控制收敛半径和两者比值。`inf` 只接受
精确字符串 `"inf"`，反演后矩阵为 `-sinv**(-2) A(1/sinv)`。

`path` 的运行时类型就是 `list[acb]`，可直接使用：

```python
snapshots, segment_reports, elapsed = transport_path(
    system, initial_vector, path,
    order=32, sample_count=96, radius_fraction=0.60,
    target_relative_error="1e-15",
)
```

`transport_path_refined(..., target_relative_error="1e-15")` 还会把主链/参考链末点相对差与
该阈值比较，返回 `target_relative_error_met`。未达到时只警告并保留结果；外层幂级数重构仍按
自身的样本精度门禁拒绝不合格样本。

`primary_order` 是生产阶数：用户可见末点、dense samples、保存点及正规化拟合输入均来自
primary chain。`reference_order` 必须更高，但只用于与主链比较和认证精度；明确命名的
`reference_snapshots`/`referenceFinalVector` 仅作为诊断。跨多个调用段续算时，两条链分别
传播自己的末点，不在段界把参考链状态重新作为主链初值。

`transport_path` / `transport_path_refined` 已支持把正则奇点保存成可复用的
`{"resultType":"frobenius_boundary","result":{"terms":[{"a":...,"b":...,"C":[...]}]}}`，
也支持把已认证指数型奇点保存为 `resultType: "exponential_boundary"` 及 `{phi,a,b,C}` terms。
以下限制只属于多列批量加速入口；它不表示 FlintNDE 不能保存奇点边界常数。

仓库级 `examples/` 另提供 `ep_parallel.py` 与 `ep_parallel_mathematica.wl`：分别演示 Python
和 Wolfram 入口在缺省上限 12 下有界并行不同固定 `ep`，并与闭式 `2^ep` 比较。MMA 使用
`FlintNDEEvaluateEpBatch[jobs, ParallelTaskCount -> 12]`；实际并行数及输入同序结果写入返回值。
原有 `qnm_2x2.py`、`regular_singular_save.py` 和 `exponential_boundary_save.py` 分别覆盖 QNM、
`{a,b,C}` 与 `{phi,a,b,C}` 保存合同。

多个 Frobenius 初值属于同一 exact `RationalMatrixSystem` 时，可用

```python
batch = transport_frobenius_boundaries_refined(
    system,
    boundaries,
    path,
    primary_order=24,
    reference_order=32,
    target_relative_error="1e-12",
)
```

`boundaries` 是非空 `FrobeniusBoundary` 列表；返回的每个 snapshot 是以这些初值为列的
Acb 矩阵。主链与参考链各共享一次局部基和每段 Taylor 矩阵，但
`relative_differences_inf`、`target_relative_error_met` 与 boundary report 仍按输入列分别给出。
当前批量入口只接受正则奇点启动后全部 transition 均为 `ordinary_taylor` 的路径，并且自身
不写 `save` 文件；需要保存普通点、保存奇点边界或使用其它局部 bridge 时，调用方应
逐个调用 `transport_path_refined`。MadStree adapter 会按该逐分支路线执行。

### 路径保存点

保存点没有名称字段。用户只在坐标后加 `"save"`：

```python
path = [(z0, "save"), z1, (z2, "save")]
snapshots, reports, elapsed = transport_path(
    system,
    initial_vector,
    path,
    order=32,
    save_output_directory=output_directory,  # 可省略；缺省为 Path.cwd()
)
```

`build_adaptive_path(system, (z0, "save"), (z1, "save"))` 以及带保存标签的
`detour_points` 使用同一格式。内部 `NamedPoint` 只服务路径分类与奇点诊断，其名称不进入保存合同。

程序到达一个保存点就立即写 `flintnde_save_001.json`、`flintnde_save_002.json` 等文件；
整条路径完成后再写 `flintnde_save_points.json`。后续步骤失败不删除已完成文件，也不写完整
汇总。`transport_path_refined` 只保存 primary chain；reference chain 不写用户保存点。
普通点的 `resultType` 是
`ordinary_vector`，结果含坐标和 Acb 列向量；正则奇点的 `resultType` 是
`frobenius_boundary`，结果含可再次交给 `frobenius_boundary` 的 `{a,b,C}` terms。后者是
局部边界数据，不是奇点处的 `Y(z_*)`。严格解耦指数型奇点的 `resultType` 是
`exponential_boundary`，结果可直接交给 `exponential_boundary`；显式带 `"save"` 的中间
奇点在 bridge 入射侧反解局部常数后立即写出。Lee--Moser 后尚无原基逆局部 jet 的情形和
需要 Stokes connection 的中间/终点仍拒绝执行；start-only formal 分支只允许保存起点，并
在输出中标记 `continuationReady: false`。

## 指数型奇点的 `{phi,a,b,C}` 边界

对 FlintNDE 已认证的指数局部解，令局部变量为有限点的 `z=s-s0` 或无穷远的
`z=sinv=1/s`。可移植边界格式表示
`exp(phi(z)) * z**a * log(z)**b * C`：

```python
from flintnde import exponential_boundary

boundary = exponential_boundary([
    {
        "phi": [
            {"power": -2, "coefficient": "3/2"},
            {"power": -1, "coefficient": "-1"},
        ],
        "a": "2/3",
        "b": 0,
        "C": [1, 0],
    }
])
```

`phi` 只允许负整数幂，系数、`a` 和 `C` 必须是 exact Q(i)；同幂项会合并，零项会删除。
`phi=[]` 表示同一高阶-pole 系统中该 sector 的 `exp(phi)=1`，用于完整保存可能同时含零/非零
指数签名的解。
程序从 DE 的高阶 Laurent 矩阵独立推导指数 sector，并逐项验证输入 `phi`，因此它不是用户
可任意选择的新参数。指数型奇点只接受 `{phi,a,b,C}`；不会从 `{a,b,C}` 的 `C` 静默猜测
指数签名。保存输出始终显式写出 `phi`，便于跨进程复用和审计。若要保存路径中的本性奇点，可写
`detour_points=((s0,"save"),)`；该坐标不是绕行普通点，而是明确要求保存的 singular checkpoint。
只有 `continuation_ready=True` 的指数局部基可这样作为中间点。

## 奇点起点的 `{a,b,C}` 边界

若路径起点是奇点，奇点本身通常没有有限数值。用户给出局部变量 `z` 下
`z**a * log(z)**b * C` 的最高 log 领头项：有限奇点使用 `z=s-s0`，无穷远使用
`z=sinv=1/s`。推荐输入格式为：

```python
from flintnde import frobenius_boundary

boundary = frobenius_boundary([
    {"a": "0", "b": 1, "C": [1, 0]},
    {"a": "0", "b": 0, "C": [2, 0]},
])

snapshots, reports, elapsed = transport_path(
    system,
    boundary,
    path,                 # build_adaptive_path 返回的奇点起点路径
    order=32,
    sample_count=96,
    radius_fraction=0.60,
)
```

每项也可写成 `(a, b, C)` 三元组；多个记录表示线性叠加。`a` 和 `C` 的每个分量必须是
exact Q(i) 输入，`b` 必须是非负 Python 整数，`C` 的长度必须等于微分方程维数且不能为零。
输入 `float`、`complex` 或 Acb ball 会在输运前拒绝。

程序按 exact manifest 验证：

- 半单 residue：`a` 必须是 indicial root，`b=0`，`C` 必须属于该 root 的 exact 本征子空间；
  整数差 resonance 在更高阶产生的 log 由 recurrence 自动加入。
- 单根 Jordan residue：`C` 是该分支最高 `log**b` 的系数，必须满足 `N C=0` 且属于
  `N**b` 的像；`b` 不得超过 exact 最大 log 次数。程序用 exact RREF 取自由变量为零的
  规范广义本征向量，随后由 `exp(N log(z))` 自动补齐较低 log 次数。需要额外低 log
  齐次分支时，再增加相同 `a`、较小 `b` 的记录。
- Moser 或指数型高阶 pole：Moser 路线生成降阶 power-log 基的 exact jet，乘回完整
  Laurent gauge 后直接在用户原积分基验证 `a,b,C`；它不会只对 leading `C` 应用逆变换，
  因而不会漏掉相邻阶分量。指数 sector 由 exact 变换后的 `C` 自动识别，不新增用户参数；
  单项 `C` 若跨多个指数 sectors 会被拒绝，必须拆成多项。

验证通过后，程序在路径的第一个普通匹配点求值统一局部基，再进入原有普通点输运。
`reports[0]` 的方法名记录实际局部路线，保存原始 terms、规范基常数、
工作变量、分支约定和匹配点。奇点起点传普通有限列向量、普通点起点传 `{a,b,C}`、错误的
`a/b/C` 或 numerical Frobenius manifest 都会 fail closed。

正规化参数重构直接复用同一协议：

```python
def boundary(ep):
    return frobenius_boundary([
        {"a": 1, "b": 0, "C": [1 + ep]},
    ])

series_result = reconstruct_series_solution(
    DEmatrix=DEmatrix,
    boundary=boundary,
    path=path,
    maximum_power=2,
    leading_power=0,
    leading_power_certificate={
        "status": "certified",
        "leading_power": 0,
        "method": "symbolic-frobenius-boundary-and-DE",
    },
)
```

每个生产和验证样本都会重新验证对应系统的 indicial/log 结构；参考链的边界初始化
审计记录保存在 `diagnostics` 的各样本 `boundary_initialization` 字段中。

这里 `convergence_radius` 是由奇点位置决定的收敛半径；`max_step_over_radius=0.45` 只表示
每一步不超过该半径的 `0.45` 倍。后续 `transport_path` 的 `radius_fraction` 应大于
`max_step_over_radius`，使实际步长严格落在 Cauchy 采样圆内。

若直线或某条折线路段遇到内部奇点，`build_adaptive_path` 在缺省
`singularity_mode="avoid"` 下返回结构化拒绝。用户可显式给出普通绕行点：

```python
path = build_adaptive_path(
    system,
    NamedPoint("left", -1),
    NamedPoint("right", 1),
    detour_points=(NamedPoint("upper_detour", 0.5j),),
    max_step_over_radius=0.45,
)
```

库函数不调用交互式 `input()`，也不会在缺省模式下静默穿越奇点。只有显式设置
`singularity_mode="singularity_jump"` 时，内部奇点才由
左/右普通匹配点之间的统一局部基 bridge 表示；程序先在入射侧反解局部基常数，再在
出射侧求值，并提醒用户确认与某一绕行同伦类等价的多值分支。起点或终点本身为奇点时，
返回列表从/到奇点收敛圆内的普通匹配点；奇点局部边界值由
`prepare_local_expansion` 计算。

返回对象 `AdaptivePath` 是 `list[acb]` 的子类，无需格式转换；额外只读诊断包括
`step_reports`、`internal_singularities`、`start_classification`、`target_classification`、
`working_variable`、`infinity_transformation` 和
`max_step_over_convergence_radius`。`step_reports` 按实际路径顺序列出每一步的长度、控制
半径、`radius_owner` 与 `step_over_convergence_radius`；进入奇点时 `radius_owner` 是
`singularity`，不会误用前一普通点半径。相同内容写入路径 JSON 的 `steps`、
`step_over_convergence_radius_list` 和 `actual_max_step_over_convergence_radius`。含 `inf`
时工作点自动位于 `sinv`，把原 `RationalMatrixSystem` 直接交给主输运即可。

## C++ 引擎

`engine/` 是 C++ 认证双链引擎（认证双链输运、`auto`/`sparse`/`dense` 三种留数内核模式、
批量流式解析、单进程 worker 线程池），基线取独立程序包 package-FlintNDE `v1.02` 引擎的
一次性冻结复制，此后按 MadStree 请求合同独立扩展，不与独立程序包同步。请求 schema 为
`flintnde_engine_evaluate_v1` / `flintnde_engine_ep_batch_v1` /
`flintnde_engine_ep_series_control_v1`；MadStree 适配器当前只使用第一个。扩展内容：

- 结构化 letters 重建：请求可携带 `tensorSectors` 与每条 letter 的
  `{alpha, beta, contributions}` 结构化贡献，引擎在 C++ 内按精确 Q(i) 高斯有理重建张量积
  留数（槽轴链、局部基、U 与逆、共轭块缓存与稀疏累加，移植自 `flintnde/tensor_dlog.py`），
  接既有精确极点合并路径；与逐极点稠密 residue 输入逐位一致，请求不再需要传输稠密留数矩阵。
- embedded 认证模式：接受 `certificationMode: "embedded"`；段生效模式＝请求 embedded 且
  段点数为一时 embedded（单参考链＋前缀中点快照），否则强制升 certified 双链，镜像 Python
  车道的升级规则。
- 生效模式上报：每段响应报告实际生效的 `certificationMode`，并以
  `planReport.engine = "flintnde_engine"` 标记所走车道（Python 车道无该标记）。

适配器（MadStree 侧 `Backend/flintnde_transport.py`）在请求满足合同时把整请求交引擎：
`pathPlanning=false`、`singularityMode="avoid"`、有限边界、段起点与全部用户点不与任何
letter 极点（`-alpha/beta`，`beta≠0`）重合，且探测到可执行文件。任一条不满足即回退
Python 认证双链，结果语义不变；引擎运行期错误 fail-closed 上抛，不静默回退（如
`alpha=0` 的退化 letter——letter 在段锚点为零——与缺 `tensorSectors` 的结构化 letter 均被
引擎拒绝并原样上报）。ep 批量请求保持 Python 进程池，池内每任务独立按同一合同尝试引擎
车道。三顶点 zd 一类正则奇点边界与开启自动规划的请求按合同始终走 Python 车道。

可执行文件探测与环境变量：

- 缺省探测 `engine/build/flintnde_engine.exe`；构建产物不随源码分发，缺失即回退 Python。
- `MADSTREE_FLINTNDE_ENGINE`：`python` 或空串强制关闭引擎车道；可执行文件路径覆盖探测；
  指向不存在的文件时报错。
- `FLINTNDE_ENGINE_LIBPATH`：libflint 运行库目录；未设置时依次探测
  `F:/codinglib/msys64/ucrt64/bin`、`C:/msys64/ucrt64/bin`，命中则前置 PATH。

构建：MSYS2 UCRT64 shell 内 `cd engine && make`（g++ `-std=c++17 -O3`，FLINT 3.x）。Git
Bash 内构建需先
`PATH=/f/codinglib/msys64/ucrt64/bin:/c/WINDOWS/System32/downlevel:$PATH`（cc1plus 的传递
依赖如 `libgmp-10.dll` 在 ucrt64/bin），且 `-flto` 在 MSYS2 shell 外不可用
（`make LDFLAGS=` 覆盖）。

车道对表：`test/test_engine_lane.py` 在真实张量请求上覆盖结构化 letter 与稠密形式逐位
一致、引擎与 Python 车道容差级一致、合同拒绝、embedded/certified 生效模式上报、
pathPlanning 与环境变量的回退行为；引擎未构建时引擎组跳过，资格与回退组仍运行。

## 依赖与奇点门禁

运行依赖只有 `python-flint>=0.6`。普通点 Cauchy--DFT 重建不调用符号求导，也不对解做
多项式回归。正则奇点缺省读取 exact Q(i) residue，以两张同维 FLINT `fmpq_mat` 精确判断
indicial roots、Jordan 结构和整数差共振；当前 exact gate 对非 Q(i) 谱 fail closed，且
失败后绝不自动降级为浮点判别。

若只有浮点 residue，使用 `NumericalRegularSingularSystem`。缺省相对零斩杀线为

```text
dimension * 10^(-precision_digits)
```

再乘 residue/regular coefficient 的最大矩阵元尺度。Python `float`/`complex` 自动按
15 位可靠十进制输入处理，即使计算工作精度更高也不会提高该值；高精度十进制或 Acb
记录应通过 `input_precision_digits` 或 `precision_digits` 声明可靠输入位数。无法判断
输入精度且用户未声明时 fail closed。`relative_zero_tolerance` 可由用户覆盖。每次浮点判别都会
发出 warning，并在 manifest 中记录 `matrix_scale`、相对/绝对 cutoff、被判零数量，
以及最大被判零量和最小保留量各自的绝对值与相对 `matrix_scale` 的比值。即使用户
覆盖阈值，这些量也会完整出现；零矩阵尺度的相对比值明确记为 `null`。
无法稳定分类的 mixed-root defective 系统会阻断，要求提高精度或提供 exact 输入。

一旦 exact 或数值 gate 判定存在 log，局部解必定使用 power-log recurrence。程序包
没有在奇点用纯多项式回归拟合解的入口。`reconstruct_series_solution` 只在多个完整 NDE
终点值之间重构 regulator 的 Laurent 系数，与局部解是否含 log 无关。

运行测试：

```powershell
python -m unittest discover -s tests -v
```

2026-11-03 fresh 结果为 `unittest discover 195/195`（v0.21 引擎接入＋精度自适应
清扫移植后）；Wolfram `Needs["FlintNDE`"]` 端到端检查为 `26/26`（2026-08-25 口径，
此后 vendor Wolfram 侧零改动）。验证覆盖一般有理矩阵、任意次数
多项式加简单极点的内部特化、缺省避开奇点、显式奇点折跃、严格消息语言、计划序列化
精度、执行期不重规划、浅层目录和 Python 前置路径/写入门禁，并覆盖 fast multipoint、
严格用户节点、符号最低阶证书门禁、增量扩阶与缓存复用。本目录 Python 侧以独立程序包
发布版 1.0 的 Python 侧为母本：core 与 routing/regularization/singularity_jump/
singular_targets/asymptotics/savepoints 等十一个模块与独立发展版 `v1.02` 逐字节一致，
其余模块保留 MadStree 特化（known-roots 证书、block 留数装配、张量重建等，
`block_operators.py`/`tensor_dlog.py` 为本目录专有）；与独立程序包不保持逐文件同步，
版本更新说明只在独立包根目录集中保存，不属于 Vendor 交付文件。

完整安装方式、普通点输运、正则奇点、显式数值路线、正规化重建及全部公开参数见
`../../Documentation/FlintNDE.pdf` 的 “Installation and public interface” 一节。
