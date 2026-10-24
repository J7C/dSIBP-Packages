# MadStree

English version: [README-en.md](README-en.md)。

更新日志：[CHANGELOG-cn.md](CHANGELOG-cn.md)；英文版 [CHANGELOG-en.md](CHANGELOG-en.md)。

`MadStree` 是一个直接使用 time-integral tensor formula 的 Wolfram Language 程序包。输入可以是有序 dS 树拓扑，也可以是不积分 loop momentum 的纯 time-only incidence graph；输出包括 contact-reachable sectors、同序主积分、逐步或完整迭代约化、block-triangular dlog connection 和自动边界证书。主算法不生成一般 IBP 方程组，也不调用 Kira。

当前缺省 convention 为

```text
time power: (-tau)^A
nu: 原样作为 Hankel order，不取绝对值
h(nu,0;z): z^(sNu nu) H_nu(z)
formulaNu: sNu nu，缺省 sNu=+1
NuConvention: "Positive"
```

只使用函数名 `h`，且 `h(nu,0;z)=z^(sNu nu) H_nu(z)`：输入的 `nu` 原样作为 Hankel order，不取绝对值。基准公式写在 `sNu=+1` 的 `+nu` 形式上，公式中凡出现阶的位置统一读取 `formulaNu := sNu nu`。`NuConvention -> "Negative"` 选择 2401.00129 的负 prefactor convention，此时 `formulaNu=-nu`；两套公式的差别只是 `formulaNu` 换号，Hankel order 始终不变。

`NuConvention` 只在 `MSInitTree` 初始化时选择。context 建立后，公式矩阵、递推与 dlog 都固定读取该值，不提供逐调用修改 convention 的入口。

## 奇点用户点输出

`MSEvaluatePath` 仍把最大连续复仿射段交给 FlintNDE。候选奇点只由该段完整解析 dlog DE
拉回后的 exact 有理矩阵/letters 决定；变量名、坐标是否为 0、图名和顶点数都不参与判定。
当用户点确实命中有限奇点时，局部幂对数解再逐分量区分可去奇点、真实 pole 和零阶 log
发散：可去奇点返回 Acb 数值，发散分量写为 `{"text":"Infinity"}`。每次命中都在顶层结果的
`singularityClassifications` 表中列出坐标、用户点序号、局部基阶数、分量分类和判定样本；
该表与 `pointResults` 中的处理方式一一对应。分类与继续输运使用同一通用 FlintNDE 路线，
不包含针对 `k0`、`k2=0` 或某个 topology 的生产分支。

## 加载

```wl
packageRoot = "/absolute/path/to/MadStree";
AppendTo[$Path, packageRoot];
Needs["MadStree`"];
```

首次成功加载会显示一次引用提醒。在 Notebook 中 arXiv 号为可点击链接；headless kernel
保留完整 URL。MadStree 只列出 [arXiv:2401.00129](https://arxiv.org/abs/2401.00129)、
[arXiv:2411.03088](https://arxiv.org/abs/2411.03088) 和 MadStree package paper
（arXiv identifier pending）。

## 最小流程

```wl
spec = <|
  "vertices" -> {
    <|"id" -> 1, "externalLegEnergy" -> k1, "timePower" -> a1, "vertexType" -> "+"|>,
    <|"id" -> 2, "externalLegEnergy" -> k2, "timePower" -> a2, "vertexType" -> "+"|>
  },
  "lines" -> {
    <|"type" -> "massless", "endpoints" -> {1, 2},
      "momentum" -> q, "nu" -> 1/2|>
  }
|>;

context = MSInitTree[spec];
masters = MSMasterIntegrals[context];
masterDefinitions = KeyTake[#, {"integral", "normalization", "bareIntegral", "definition"}] & /@ masters;
topKey = First[context["sectorOrder"]];
matrices = MSFormulaMatrices[context, topKey];
de = MSDLogDE[context];

formulaData = MSFormulaData[context];
written = MSWriteFormulaArtifacts[
  context,
  TimePowerRules -> Automatic
];
written["outputDirectory"]
written["files", "dlogDE"]

targetRules = {k1 -> 9 I, k2 -> 3 I, q -> 1, a1 -> 1, a2 -> 1};
chart = MSBoundaryChartCertificate[context, targetRules];
```

`MSIntegral[sectorKey,timeShifts,stateBits]` 表示 normalized master
$J_s(\mathbf n;\mathbf a)$；同一组三项也唯一标识裸积分 $I_s(\mathbf n;\mathbf a)$。
`MSMasterIntegrals` 直接给出每个 master 的精确 `normalization` 和
`J_s=normalization I_s`；单个 shifted integral 可用
`MSIntegralDefinition[MSIntegral[sectorKey,timeShifts,stateBits],context]` 查询。
massive Full 传播子收缩时，该 normalization 含唯一的 contour 因子：省略共同动量幂后，
`++` subsector 为 `-4 I/Pi`，`--` subsector 为 `+4 I/Pi`；该符号不重复进入 dlog contact block。
二态因子要么是 massive 传播子在其中一个顶点上的二维态，要么是整条 massless Full 传播子
共享的 quotient 二态；其逐位顺序由 `MSSlotRegistry` 给出。`++/--` 的 h 组合只在手册统一定义，
不在每条输出中展开。`MSSectors[context]` 返回 contact-reachable sectors 的确定顺序（键、
被收缩传播子与 master 数），与 `context["sectorOrder"]` 一致；`MSContextQ[obj]` 判定对象
是否为本包合法 context，是全部公开入口共用的输入门禁。

`"vertexType" -> "+"|"-"` 是顶点的 Schwinger--Keldysh 轮廓支，必须在每个顶点上
显式给出。传播子不再输入 `skType`、`sigma` 或顶点符号；程序从 `"endpoints"` 指向的
顶点自动推导 Full/Cross/External、`++/--/+-/-+`、整体符号和缺省 Hankel branches。
公开传播子 `"type"` 只接受 `"massive"` 或 `"massless"`。传播子不输入 ID；程序始终按
`"lines"` 输入顺序在内部编号。只有显式 `thetaBundles` 需要引用传播子时，才使用这一位置编号。

`"externalLegEnergy"` 不是“顶点自身具有的能量”，而是附着在该顶点的外腿无 theta
指数参数。`vertexType="+"` 唯一给出 `Exp[-I externalLegEnergy tau]`，`"-"` 唯一给出
`Exp[+I externalLegEnergy tau]`。SK contour sign 也由同一 `vertexType` 派生，但在内部与
外腿指数符号分字段保存，传播子 contact map 不读取第二套用户符号。对 `+` 支的正阻尼坐标
`K>0`，应取 `externalLegEnergy=+I K`；对 `-` 支则取 `-I K`。

上述 massless 图可直接生成完整多 sector Frobenius 边界，并通过唯一数值入口输运到一个或多个用户点：

```wl
result = MSEvaluatePath[
  context,
  {{k1, k2}, {9 I, 3 I}},
  ParameterRules -> {q -> 1, a1 -> 1, a2 -> 1},
  FlintNDEPathPlanning -> True,
  BoundaryScale -> 4,
  RankOrder -> {1, 2},
  PythonExecutable -> "...",
  WorkingPrecision -> 200,
  TransportOrder -> 80,
  ReferenceTransportOrder -> 104,
  TargetRelativeError -> "1e-20",
  MessageLanguage -> "EN"
];
```

第二个参数始终叫 `pointSequence`。首行是有序可跑动坐标符号，后续是等宽坐标值行；单点只有
一行值，多点继续追加值行。固定、不可偏导参数只在 `ParameterRules` 中一次给出。表头与
`ParameterRules` 必须互斥，二者合起来必须覆盖解析 DE 和边界需要的全部符号；缺失或替换后
仍非数值时，程序在创建 runtime 和启动 Python/FlintNDE 前返回明确 `Failure`，其中列出问题
符号并说明 NDE 未启动。MadStree 按输入顺序划分最大连续复仿射单变量段 `x(s)=x0+s v`，每段
只拉回一次并把全部 exact 复参数点交给 FlintNDE；MadStree 不规划节点或绕行。普通值行缺省
返回，`{{values...},"tmp"}` 是不进入最终保存结果的途经点；其它字符串标签一律拒绝。点结果
携带 `coordinate`、`value`、`status` 和 `userIndex`。

外能量可以在用户扫描中保持固定并由 `ParameterRules` 给出。由于自动无穷远边界会改变全部
外能量，程序在内部把这些固定外能量补入路径坐标：从真实的边界匹配点输运到首个用户点，
随后在各用户点保持给定值。返回的 `coordinateSymbols` 和逐点 `coordinate` 因而可能比输入
表头多出这些外能量；返回的 `parameterRules` 不再含已补入坐标的外能量。用户也可以直接把
全部外能量列进表头，并在各行重复固定值。内动量、质量和时间幂次参数的处理不变。

`TransportOrder -> 80` 是生产阶数，返回的节点值、dense 用户点、段末值及后续正规化拟合均取
低阶 primary chain。`ReferenceTransportOrder -> 104` 只运行高阶 reference chain，并用两链
差值检查 `TargetRelativeError`。跨多个复仿射段时两条链分别续传自己的末点，高阶值不会替换
用户结果，也不会在下一段成为主链初值。

共同解析 regulator `ep` 的 Laurent 极限使用：

```wl
epSeries = MSReconstructEpSeries[
  context,
  ep,
  {{k1, k2, k3}, {9 I, 3 I, 5 I}},
  ParameterRules -> {q12 -> 1, q23 -> 2,
    a1 -> 1 + ep, a2 -> 1 + ep, a3 -> 1 + ep},
  MaximumEpPower -> 0,
  EpGoalDigits -> 20,
  ParallelTaskCount -> 12
];
```

`MaximumEpPower` 是用户需要返回的最高 `ep` 幂；缺省 `0` 表示需要 pole（若存在）及有限项。
缺省时用户不提供任何 `ep` 取值。任何数值 NDE 启动前，程序先从实际符号边界条件与 dlog DE
认证最低整数幂：检查 `ep=0` 处 DE 无负 Laurent 阶、定义积分边界解析且物理分支合并后
存在可证明非零的最低系数；正规化参数进入未认证路径坐标或结构无法证明时 fail closed。
输入覆盖和 boundary chart 的数值预检在固定非零 exact 点 `ep=1/1009` 进行，只用于确认除
regulator 外的参数已经给全并且 generic chart 可构造。这个 probe 不进入生产/验证点池、不参与
拟合，也不决定最低幂；DE、完整物理边界与共振记录的正式 Laurent valuation 仍严格取 `ep=0`。
随后按最低幂、最高幂和 `EpGoalDigits` 自动决定生产点、`ep` 尺度、工作精度、输运阶数及
不参与拟合的独立验证点。内部拟合缺省比用户最高阶多两阶；验证失败时每轮再增加两阶，
只计算新增生产点，既有生产点和验证点分别缓存复用。最多三轮仍不满足原误差门槛时返回
当前最佳系数，并用 `status -> "computed_with_warning"`、`precisionTargetMet -> False` 和
结构化原因明确撤销精度认证。
自动工作精度取 200 位与自适应估计的较大者。
`ParallelTaskCount` 缺省为 12，控制每批独立 `ep` 后端进程；超出部分自动续交。它不同于
python-flint 单进程内的 `ctx.threads`。

若需把所有取值限制在用户认可范围，可显式给出 `EpSamplePoints -> {...}` 与
`EpValidationPoints -> {...}`。前者是有序生产候选池，可以多于首轮所需点数；设置
`EpInitialInternalMaximumPower -> q` 后，首轮只取拟合到 `ep^q` 所需的前缀，验证失败才从
剩余候选中增量取点并复用旧值。验证点从不参与拟合且必须与整个生产候选池分离；候选耗尽
不会自动越出用户范围，结果原因是 `candidate_pool_exhausted`。

另一种自动入口是 `EpSampleAngleRange -> {thetaMin,thetaMax}`，单位为弧度并按开区间解释，
不能与显式 `EpSamplePoints` 同时使用。程序在内部均匀选择最多三条射线，不取边界；模长仍由
pole 深度、`EpGoalDigits` 和拟合阶数自动决定，没有模长上限选项。自动验证点保持角度并缩小
模长。为满足 MadStree exact dlog 拉回合同，adapter 把高精度生成的复点转为 exact
Gaussian-rational 后再求值。上述选项缺省均为 `Automatic`，因此缺省调用及后端输入 schema
不变。

`FlintNDEPathPlanning -> True` 让 FlintNDE 在每段内部规划节点。落在同一节点收敛圆盘内的用户点组成一个 evaluation bucket，并用该节点保存的向量级数做快速多点求值；点数不少于 8 的桶使用子积树/余数树，小桶使用 iterative 算法。`False` 则严格把用户点依次作为节点，不插点、不删点、不调用规划器。不同复仿射段不共享局部系数，因此没有多变量高维 Taylor 球。若一段全部 dlog letters 为常量，拉回连接为零并正常输运。

`SingularityMode -> "Automatic"` 是缺省：用户显式给出的奇点不被静默绕开。FlintNDE 把精确奇点从普通节点链移除，用同一局部基覆盖收敛域内双侧用户点，真实发散返回文本 `Infinity`，并从出射普通点继续；末端奇点按目标奇点到最近其它奇点的距离自动选择域内隐藏普通匹配点。`"Avoid"` 显式拒绝穿越；`"SingularityJump"` 显式选择局部穿越分支，其多值性等价于某条绕行路径，用户必须确认。奇点处非共线转向和关闭规划后的中间奇点 fail closed。`MessageLanguage -> "EN"|"CN"` 严格区分大小写，缺省英文。`MSBoundaryData` 与 `MSEvaluatePath` 的 `WorkingPrecision` 缺省为 200 位；显式给值时直接覆盖。内部 bit 数为 `ceil(WorkingPrecision*log2(10))+32`。

`de["masters"]` 是 DE 行列和边界向量顺序的唯一 authority。`MSBoundaryData` 对已闭合的 tree/time-only context 统一生成 nested curve、拉回完整 dlog connection、求 residue，并按 ancestor sectors 解每个 master 的 indicial leading vector。`MSEvaluatePath` 把 exact `{a,b,C}` 分支、全部复仿射段和用户选项放入同一 UTF-8 请求；FlintNDE 在同一进程内完成边界初始化和后续有限段输运。生产代码不调用定义积分或 `NIntegrate`；这些只允许作为独立验证 oracle。

通用入口不按图 id、顶点数、master 数或论文数值点分派。只有 dlog 未闭合、rank/chart 未通过 normal-crossing 证书、late-time 指数不衰减、拉回系统不是 exact regular singular，或目标/anchor 落在 DE letter 上时才结构化 fail closed。

一条传播子的两个顶点同属一个轮廓支时，公开 `"massless"` 传播子在内部派生为 `masslessFull`，
并且唯一使用整条传播子共享的 quotient 二态 slot。两个顶点各自的二维态构成的四态乘积只用于
手册推导 quotient，不是程序输入或第二套生产表示；当前接口没有表示选择选项。

每个 sector 的 `sectorKey` 是按初始化时 root `lines` 顺序排列的定长字符串：`0` 表示对应传播子已收缩，`1` 表示未收缩；不能 contact shrink 的传播子位始终为 `1`。top 是全 `1` 字符串。前导零是身份的一部分，禁止把 key 当整数。`context["sectorKeySchema"]` 给出 `rootLineOrder`、位宽和位语义；例如 `{e1,e2,e3,e4}` 中收缩 `{e1,e2,e4}` 得到 `"0010"`。master 的完整身份是 `MSIntegral[sectorKey,timeShifts,stateBits]`，所以不同 subsector 即使局部 shifts/bits 相同也不是同一积分。`MSInitTree` 的 `sectorIdentityCertificate` 会一次性检查 contraction set、sector key、完整 master、global index 和按完整 master 顺序计算的 SHA-256 digest；碰撞时拒绝初始化。massless top-to-sub dlog 的 contact 原子按 event 中实际被收缩的 massless 传播子数乘 `(-1)^N`，simultaneous contact 不使用图专用符号表。

`MSBoundaryChartCertificate[context,targetRules]` 在边界计算前构造
`1/K[sigma[j]]=Product[x[r],{r,j,V}]` 的 nested blow-up，逐 sector 固化所有 theta，并检查完整 dlog letters、normalization、坐标 Jacobian 与 shifted-contact 新分母都是 coordinate monomial 乘边界非零 unit。`RankOrder -> All` 检查全部 root-time strict charts；证书不通过时 `MSBoundaryData` 返回 `BoundaryChartNotCertified`。证书与边界记录中的局部量使用两个公开符号：`MSBlowupCoordinate[i]` 是嵌套 blow-up 的第 `i` 个局部坐标，`MSDampingEnergy[v]` 是顶点 `v` 的正阻尼能量 `K_v`。

## 共同 theta 与 time-only 圈图

同一当前 component pair 上的 Full 传播子自动组成共同-theta bundle。程序只生成非空奇数子集 event，系数为 `2^(1-Length[selected])`；一个 event 同时删除所选传播子并只合并一次顶点。sector 由 event BFS 生成，不是 Full 传播子的幂集。合并后未选的 massless Full 传播子只保留 even state，massive Full 传播子使用 `10 -> 01` equal-time canonical；raw tensor matrix 通过保存的 projection/embedding 投到真实 master space。

纯 time-only 圈图使用专用入口：

```wl
cycleContext = MSInitTimeGraph[<|
  "vertices" -> {...},
  "lines" -> {...}
|>];
```

它只允许固定传播子动量的时间积分，不读取或生成 loop momenta、ISP、Landau/threshold 数据。active 的自环传播子不再触发 contact；多条收缩路径到同一 vertex partition/contracted set 被 canonical 成一个 sector。

## 单顶点函数族

单顶点 family 不需要伪造 topology。紧凑输入沿用参考代码的 `ki/nui` 习惯：

```wl
vertexContext = MSInitVertexFamily[<|
  "ki" -> {k0, k1, k2},
  "nui" -> {a0, nu1, nu2},
  "hankelBranches" -> {1, 2}
|>];

vertexDE = MSDLogDE[vertexContext];
singlePointSequence = {{k0}, {3 I}};
vertexValue = MSEvaluatePath[
  vertexContext,
  singlePointSequence,
  ParameterRules -> {k1 -> 1, k2 -> 2, a0 -> 2, nu1 -> 1/2, nu2 -> 3/2}
];
```

这里 `First[ki]` 是附着外腿的无 theta 指数参数，`First[nui]` 是 `(-tau)` 的基准幂；其余位置一一给出 h block 的动量和原样取用的 Hankel 阶 `nu`。显式模型改用
`externalLegEnergy/timePower/hBlocks/exponentialBlocks/vertexType/normalization`，其中每个 h block
给 `id/momentum/nu/hankelBranch`，每个额外指数 block 给 `id/momentum/exponentType`。这套初始化
schema 与树图的 `vertices/lines` 不同，因为它定义的是单个顶点上的函数积而不是 incidence graph；
初始化完成后，主积分、递推、DE、边界及上述 `MSEvaluatePath` 数值接口完全相同。纯指数 block
是 `1x1`，不增加状态位；每个 h block 是 `2x2`。初始化时选定的 `NuConvention` 此后固定，不能在递推、DE 或数值计算时另行覆盖。

## 直接约化

```wl
reduction = MSReduce[
  2 MSIntegral[First[vertexContext["sectorOrder"]], {1}, {0, 0}] -
  3 MSIntegral[First[vertexContext["sectorOrder"]], {-1}, {1, 0}],
  vertexContext,
  MasterBasis -> Automatic
];
```

`MSReduce` 只处理固定初始化 context 内的合法 `MSIntegral[sectorKey,timeShifts,stateBits]`：`sectorKey` 必须属于 contact-reachable sector DAG，`timeShifts` 是该 sector 每个 component 的整数平移，`stateBits` 按 2401.00129 的二进制顺序遍历该 sector 全部二维 slots。它可线性处理有限多个这样的对象，也可直接接受 `List`、规则 `Table`、空列表和 ragged 嵌套列表；返回的 `result`、逐元素 `coefficientVector/masterRules`、residual 和 remaining shifts 保持输入形状，并共享同一次调用的递归 memo。任一元素失败时会报告其精确位置。它不支持新增传播子幂、新离散指标、context 外 sector 或另一个函数族。约化沿 shift/contact DAG 递归到零 shift masters，并缓存重复子问题。标量返回的 `masterBasis`、`coefficientVector` 和 `masterRules` 严格同序；`result` 是显式主积分线性组合；`nonMasterResidual` 与 `remainingShiftedIntegrals` 必须同时为空才是完全约化。`MasterBasis` 可给 context 全部 masters 的任意排列，但缺项、重复或额外对象都会拒绝。`singularLayers` 逐步保存方向、分量和分母：正 shift 向零的 raising 只报告 `M0` energy letters，负 shift 向零的 lowering 只报告 `M1` 本征值。约化使用初始化时固定的顶点时间幂：要代入数值或正规化幂（如 `a1=-1+ep`），以这些值作为 `timePower` 重新初始化 context 后约化（Example 06 演示），或对符号结果事后代入；`MSReduce` 本身没有替换选项。

单步入口为 `MSRecurrenceStep[integral, component, context]`：对选定 component 把一个非零 shift
移向零一步，省略 component 时自动取第一个非零 shift，返回 `reducedOneStep`、`nearerShift`
和该步实际被反演的分母；`MSReduce` 即其 memo 化的完整入口。contact 注入结构可用
`MSContactMaps[context, sectorKey]` 单独读取：每个 contact 事件一条记录，含按分量的注入矩阵
（`matricesByComponent`）、目标 sector 与零/升幂两档目标 time shifts（`R0TargetShift`、
`R1TargetShiftsByComponent`）、normalization 比值和被收缩传播子列表；省略 sectorKey 或传
`All` 时按 `context["sectorOrder"]` 返回全部 sector。

初始化会用已知 h 两支小时间领头幂次，为每个 reachable component 保存时间-IBP 安全 shift
下界。零 shift master 已命中 Gamma 非正整数时立即拒绝；正整数起点仍可初始化，但实际约化
跨过安全下界时才拒绝，并列出原始顶点、sector、分支、Gamma 参数和所请求 shift。加入通用
符号解析正规化参数后，程序保留该参数并继续符号计算。

若 top-to-sub 的 `R^(1)` 落在 subsector 非零 shift，`MSDLogDE` 会逐列调用同一公式递推，生成 shifted subsector 到全局同序 masters 的 reduction matrix，再右复合 contact block。每个 contact block 保存 `shiftReductionRecords`、residual、remaining shifts 和 singular layers；任何一列未闭合时 dlog 状态为 `contactShiftReductionFailed`，不会输出伪闭合 connection。

公式层不构造一般大矩阵逆。massive 顶点态的固定 `sigma2` 使用论文的 `2x2` 变换及其
显式逆，massless Full 传播子共享 quotient slot 的固定 `sigma1` 使用自逆 Hadamard；全 sector 只作这些
局部矩阵的 Kronecker 积。`M1` 在 state-bit basis 逐对角元取倒数，`M0` 在对应局部张量基逐
energy letter 取倒数。若未来同一个 slot 出现多个不能由同一局部变换对角化的组分，这条快速
路径必须 fail closed。

FlintNDE 数值后端作为版本内嵌副本放在版本目录内，不跟随独立包的后续版本；缺省位置只在
一个相对路径变量中定义：

```wl
MSFlintNDEConfiguration[]
(* relativePath -> Vendor/FlintNDE *)

MSSetFlintNDERelativePath[FileNameJoin[{"new-name", "code", "package"}]]
```

路径始终相对当前 MadStree 版本目录；Vendor 此后改名或移动时只改这一项。`MSNumericalSystem` 保留为用户自行提供边界向量时的低层入口。

MMA 自动调用的运行根缺省是调用脚本旁唯一的 `results_temp/`：

```text
results_temp/
  nde/    # 当前请求输入与日志
  cache/  # 按请求和后端源码身份命名的成功 JSON cache
```

`MSRuntimeDirectory -> path` 的显式路径就是运行根本身；绝对路径直接使用，相对路径相对当前调用脚本目录解析，不依赖进程工作目录，也不会再次追加 `results_temp`。Windows 下在目录创建和 Python 启动前检查完整 cache/input/log 路径；超过安全上限时独立返回 `RuntimePathTooLong` 及实际长度。`RuntimeInputWriteFailed`、`PythonFlintUnavailable`、`FlintNDELaunchFailed`、`FlintNDEOutputMissing` 和 `FlintNDEOutputInvalid` 保持互斥故障边界。成功结果缓存键同时包含请求内容以及 `Backend/flintnde_transport.py`、Vendor `pyproject.toml` 和排序后的 `flintnde/*.py` 的 SHA-256，因此当前版本原位修复后不会误复用与源码身份不匹配的计划。Python 按缺省规则在相应 package 旁建立 `__pycache__/`，该目录及 `*.pyc` 已由 Git 忽略。除这些调用侧缓存外，程序包源码目录不接收运行产物。
MadStree adapter 与 Vendor Wolfram bridge 均通过参数列表 `RunProcess` 启动 Python，不经过
shell `Run`、命令引号拼接或重定向，也不提供重试 fallback；连续调用的 Windows DLL 初始化
回归由当前数值门禁和 Vendor `25/25` 接口测试共同覆盖。

解析 `de = MSDLogDE[context]` 本身是正式用户结果。纯符号工作流可显式执行
`MSWriteFormulaArtifacts[context]`，它的两个专用选项是 `TimePowerRules`（只特化写出的递推
元数据，不改解析 DE）与 `MSOutputDirectory`（重定向产物目录）；`MSEvaluatePath` 与 `MSReconstructEpSeries` 则在启动任何
Python/FlintNDE 数值任务前自动写出或复用同一 context/调用根的解析资产。写出失败时数值 NDE
不启动；成功返回值的 `"formulaArtifacts"` 给出正式目录与全部文件路径。缺省写入：

```text
results/
  madstree_formula/run-<UUID>/
    masters.wl
    recurrence_metadata.wl
    dlog_de.wl
    manifest.wl
```

`dlog_de.wl` 始终保存完整解析 DE；`ParameterRules` 只生成求值函数内的局部副本，不覆盖或
回写正式公式。`manifest.wl` 记录版本、case、convention、sector/master 摘要和全部实际路径。
同一 Wolfram 会话中重复计算同一 context 时复用仍完整存在的这组文件；文件被删除则下次调用
重新生成。四个文件未全部写成时返回 `FormulaArtifactWriteFailed`，不会声称成功或启动后端。

普通用户点结果保存在 `MSEvaluatePath` 返回值的 `"saved"` 记录中；`"tmp"` 点只参与输运。可用

```wl
MSExportEvaluationData[
  result,
  MSOutputDirectory -> outputDirectory,
  ExportFormats -> {"CSV", "JSON"}
]
```

把已保存普通点导出为 CSV/JSON。缺省目录是
`results/madstree_evaluation/run-<UUID>/`，文件固定为 `evaluation_data.csv` 和
`evaluation_data.json`。任一请求格式写出失败时返回 `EvaluationExportFailed`，不会因另一格式
已写出而声称整体 `"written"`；路径过长仍单独返回 `RuntimePathTooLong`。临时点不进入普通点
导出。`ExportFormats` 可取 `Automatic`、单个 `"CSV"`/`"JSON"` 或两者组成的列表；JSON 中
本来就是字符串的精度/相对差字段原样 round-trip，不会再包一层字面引号。`results/` 只放用户正式结果，`results_temp/` 只放请求、日志和 cache，
`test/results_test/` 只放开发测试产物；程序包源码目录不接收运行产物。
完整公式、massless `4 -> 2` quotient、contact shift 与 top-to-sub dlog 推导见 [Documentation/用户手册.pdf](Documentation/用户手册.pdf)；全部公开符号的一句话总表见该手册「公开接口总表」一节。

## Examples

- [01_massless_full_edge.wl](Examples/01_massless_full_edge.wl)：massless quotient、主积分、递推、dlog 和自动边界/数值入口（含用户自定义有限边界与批量多点求值、CSV/JSON 导出）。
- [02_single_vertex_family.wl](Examples/02_single_vertex_family.wl)：单顶点函数族的紧凑/显式定义、主积分、约化，以及统一 `pointSequence` 下的单点与多点数值互检。
- [03_time_only_cycle_chart.wl](Examples/03_time_only_cycle_chart.wl)：time-only 圈图、共同 theta、contact sector 与全部 strict-rank chart。
- [04_three_vertex_tree.wl](Examples/04_three_vertex_tree.wl)：三顶点 massless 树图（+++ 顶点结构），后接批量多点求值与 CSV/JSON 导出。
- [05_massive_three_vertex_tree.wl](Examples/05_massive_three_vertex_tree.wl)：三顶点 massive 树图
  （+++ 顶点结构、非半整数 nu），含 normalized/bare 定义自检：逐 master 比对
  `MSIntegralDefinition` 与 `MSMasterIntegrals` 记录的 integral/bareIntegral/normalization 字段，
  并校验 shifted 积分读取所在 sector 的 normalization。同一 context 下先计算并导出原非零点列，再以第二次独立
  `MSEvaluatePath` 调用计算中间顶点外腿能量 `k2=0` 的四点列；两组不混合规划或结果。
- [06_massless_three_vertex_ep_regularization.wl](Examples/06_massless_three_vertex_ep_regularization.wl)：
  三顶点 massless 树的时间幂正规化 `a1=-1+ep`、`a2=a3=1+ep`，`ep=0` 处时间积分真正发散；
  数值 NDE 前先做两项符号自检：`MSReduce` 一次作用于两个相对主积分幂偏移 `-1` 的积分之和并核对
  status/residual/remaining/singularLayers；以 `a1=-1+ep, a2=a3=1+ep` 代入 `timePower` 重新初始化
  context 后重复同一约化（系数显式带 `ep`）。
  之后由符号边界与 DE 认证最低幂为 `-1`，再选择生产/验证点并一并提取 `ep^-1` 极点系数与有限项。
- [07_zero_external_leg_energy.wl](Examples/07_zero_external_leg_energy.wl)：无物理外腿的顶点省略
  `externalLegEnergy` 或显式输入 0；私有辅助能量保留在解析 dlog 中并自动输运到物理零点。

七个 examples 已在删除既有 `results/` 与 `results_temp/` 后从当前版本路径全部 fresh 运行并退出 `0`；
Example 06 自适应检查为 `17/17`，使用 4 个生产点、2 个独立验证点，缺省请求 12 并按实际
点数自动使用 6 个并行任务。

## 开发与验证说明

本程序包的记号、功能与对外接口由作者设计，[Documentation/用户手册.pdf](Documentation/用户手册.pdf)
记录的推导亦由作者给出，实现由作者指导的多个大语言模型 agent 程序完成，每条公式都对照推导审阅过。
可靠性由作者设计的验证任务确立；下面只列出交叉验证计算结果的任务及其参照数据来源。

- **单顶点函数族与两顶点 massive 关联函数**：精确微分方程、边界系数与数值输运结果与
  [arXiv:2411.03088](https://arxiv.org/abs/2411.03088) 的已知结果比较，后者由独立的多变量
  超几何构造给出。
- **三顶点 massive 树图**：参照数据由一个 agent 程序仅依据此前手写的两顶点代码，沿独立路线
  推广到三顶点得到；其输出先与 [arXiv:2309.10849](https://arxiv.org/abs/2309.10849) 发表的
  八分支结果核对通过，再用于检验本程序包，因此参照数据与本程序包不共享推导。
- **IBP 关系与微分方程**：与 `dSIBP` 程序包（in preparation，与本包一同发布于
  [github.com/J7C/dSIBP-Packages](https://github.com/J7C/dSIBP-Packages)）比较。`dSIBP` 采样
  naive IBP 种子并求解所得线性约化系统，本程序包直接套用推导所得的通用公式、不再次求解 IBP，
  因此一致性检验的是闭式公式本身。

## 当前边界

- 已支持共同-theta odd-subset simultaneous contact、coincident massless/massive quotient、event-reachable sector BFS 与纯 time-only 圈图；这不等于支持 loop-momentum integration。
- 生产边界不按图名或 master 数分派。单顶点 massiveExternal 也走同一通用路线（2411.03088 Sec.3.3 显式级数仅保留为测试对照基准）；已闭合的 tree/time-only context 统一由 sector DAG、component/slot metadata、normalization 和 strict time rank 生成 nested curve、完整 dlog pullback residue 与 ancestor-sector leading system。公式/dlog/chart 未闭合、late-time 指数不衰减或拉回系统不是 exact regular singular 时结构化 fail closed；有限点定义积分只用于独立验证。
- FlintNDE 输运要求拉回 connection 属于 exact `Q(i)(s)` 或 `Q(i)(t)`。本轮 fresh 开发检查
  覆盖单顶点、pure massless、massive Full 和 mixed 三顶点，并统一使用当前 generic boundary producer。
- FlintNDE 保留 exact Lee--Moser、高阶 pole 与奇点折跃能力，并增加按节点覆盖桶的 fast multipoint evaluation 和公开严格用户节点入口。未认证 high-pole、需要 ramification、代数扩域或一般 Stokes connection 的内部点与终点继续 fail closed。

当前 fresh 回归包括 core `54/54`、公式 artifact `23/23`、点序列/规划 `29/29`、
normalized/bare integral 定义 `20/20`、二维
900 点复平面分组 `8/8`、
单顶点 NDE `13/13`、runtime/export `10/10`、辅助外腿 `6/6` 及 Examples 01--07 `7/7`；
23 个 Wolfram 开发脚本、Python adapter `17/17` 和 Vendor FlintNDE `193/193` 均通过。
