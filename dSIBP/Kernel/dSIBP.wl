(* ::Package:: *)

(* ::Chapter:: *)
(*公开上下文与接口声明*)

BeginPackage["dSIBP`"];

$dSIBPVersion::usage = "$dSIBPVersion 给出当前加载的 dSIBP package 发布号。";
AuditLevel::usage = "AuditLevel 控制同源数据是否重跑完整开发证书；\"standard\" 复用 sealed producer 状态，\"full\" 重算内容 hash、coverage、representation 或 residual 审计。";
KiraRequireCompleteSystem::usage = "KiraRequireCompleteSystem 指定 Kira 导出是否要求 completeSystemQ；正式约化缺省为 True。";

J::usage = "J 是统一积分 Head。full 模式使用 J[aList,linePacks,ispList]；timeOnly 公开接口只使用 J[sectorKey,timeShifts,stateBits]，其中 sectorKey 是 root propagator 顺序的定长 0/1 字符串。";
sp::usage = "sp[p,q] 是用户输入使用的对称标量积。";
qq::usage = "qq[i,j] 是内部圈动量标量积坐标。";
qk::usage = "qk[i,j] 是内部圈动量与外动量标量积坐标。";
kk::usage = "kk[i,j] 是内部外动量标量积坐标。";
a::usage = "a[v] 是 makeBaseIntegral 为顶点 v 建立的时间幂次整数指标。";
b::usage = "b[e] 是未缩并传播子 e 的分母幂次整数指标。";
bS::usage = "bS[e] 是缩并传播子 e 的分母幂次整数指标。";
n::usage = "n[e,...] 是 full-line 的二元离散态指标。";
ispN::usage = "ispN[i] 是第 i 个 ISP 坐标的整数幂次指标，定义零点固定为 0；正值为 numerator 幂，用户显式给出的负值表示该坐标的额外 denominator，package 不阻断。";
a0::usage = "a0[v] 是顶点时间幂次零点；物理幂次为 a+a0。";
b0::usage = "b0[e] 是未缩并传播子分母幂次零点；物理分母幂次为 b+b0。";
bS0::usage = "bS0[e] 是缩并传播子分母幂次零点；物理分母幂次为 bS+bS0。";
rho::usage = "rho[i] 是第 i 个 ISP 的内部坐标符号。";
dim::usage = "dim 是圈动量积分的空间维数参数。";
ke::usage = "ke[i] 是与动量向量无关的独立顶点相位能量参数；实际无圈动量模长由 independentExternalMomenta 依次绑定为 sE1,sE2,...。";
kEpower::usage = "kEpower[be1,be2,...] 是 sector metadata 中按 independentExternalMomenta 稳定顺序保存的无圈模长幂向量；当前参数表达式由同一 metadata 的 kEParameterExpressions 单独给出。";
tau::usage = "tau[v] 是 rep2Integrand 输出中的顶点共形时间。";
xi::usage = "xi[e] 是 rep2Integrand 输出中的第 e 条线的动量模。";
Hh::usage = "Hh[block] 是 rep2Integrand 使用的惰性传播子 building-block 包装。";
MassiveBlock::usage = "MassiveBlock[...] 是 massive line 的惰性 integrand block。";
MasslessBlock::usage = "MasslessBlock[...] 是同分支 massless line 的惰性 integrand block。";
MasslessCrossBlock::usage = "MasslessCrossBlock[...] 是异分支 massless line 的惰性 integrand block。";
Tuserweight::usage = "Tuserweight[id] 是 Kira user-defined system 结果中的积分编号 token。";
userMI::usage = "userMI[i] 是用户定义主积分坐标的稳定公开 token；其物理内容始终由同源 J 线性组合给出。";

dtau::usage = "dtau[vertex,expr] 生成指定顶点的时间 IBP；三参数形式接受 parsed topology 或 DSInit context。";
dqq::usage = "dqq[dLoop,vectorLoop,expr] 生成圈动量沿圈动量方向的 IBP；四参数形式接受 parsed topology 或 DSInit context。";
dqk::usage = "dqk[dLoop,vectorExternal,expr] 生成圈动量沿外动量方向的 IBP；四参数形式接受 parsed topology 或 DSInit context。";
ds::usage = "ds[expr,var] 对 exact 初始化 context 的外部变量 var 求总导数；三参数形式接受 parsed topology 或 DSInit context。缺省 var 为 loop Gram 根号 ssij 或显式独立无圈模长 sE1,sE2,...，并同时作用于积分指标、sector prefactor 和显式动力学系数。";
rep2innerform::usage = "rep2innerform[expr] 把用户 sp/ssij/sEe 表示转换为当前 topology 的内部坐标；双参数形式接受 parsed topology 或 DSInit context。一般混合或过完备坐标没有唯一反向映射时返回 $Failed。";
rep2outform::usage = "rep2outform[expr] 把内部标量积坐标按当前规则转换为用户 sp/ssij/sEe 表示；双参数形式接受 parsed topology 或 DSInit context。";
rep2Integrand::usage = "rep2Integrand[expr] 把统一 J 表示展开为用于核对的形式 integrand；双参数形式接受 parsed topology 或 DSInit context。";
symmetry::usage = "symmetry[expr,topo] 一次应用 topology 的内建、用户和 tadpole 对称性规则。";
repSymmetry0::usage = "repSymmetry0[topo] 返回 topology 输入的原始用户对称性规则。";
repIterative0::usage = "repIterative0 保存最近一次 tree 单步迭代生成的原始替换规则。";
repIterative::usage = "repIterative[expr,end] 把 tree 积分迭代约化到各顶点的目标时间幂次；sector-tagged treeLinearData 输入会保持 sector 身份并返回同结构约化结果，end 缺省为全零。";
DSTreeSeeds::usage = "DSTreeSeeds 直接生成带 sector/contact 审计的 pure-time/tree 种子；loop time-IBP 投影只保留为交叉验证。";
DSTreeNaiveIBP::usage = "DSTreeNaiveIBP[context,masters] 把 loop time-IBP 投影成 sector-tagged tree 线性系统，并在指定有序 master basis 下直接求解全部一步升幂对象；masters 缺省取 DSTreeDLogDE 的同序归一化 masters。";
DSTreeNaiveDE::usage = "DSTreeNaiveDE[context,variables,masters] 通过 loop 顶点相位导数投影、h 的传播子动量模长导数和 DSTreeNaiveIBP 约化构造 tree 微分方程；结果保持指定 master 顺序和 normalization。";
DSTreeDLogDE::usage = "DSTreeDLogDE[data] 返回 tree vertex-family 的 dlog 微分方程、同序 master 列表和 letters；DSInit context 输入使用 direct pure-time contact selectors 组装全部可达 sector 的 block-triangular connection、normalization 审计与同序 tagged masters。严格表示与 dlog residual 证书只在 AuditLevel->\"full\" 执行。";

DSInit::usage = "DSInit[input,opts] 验证 topology/ISP、初始化完整 contact-reachable sector，并可写出版本化 init metadata。";
DSInfo::usage = "DSInfo[] 返回当前初始化的简要信息；DSInfo[context,\"Full\"] 返回完整初始化 Association。";
DSKinematics::usage = "DSKinematics[input,rules] 返回 topology 的缺省动力学变量提案、全部必需模长覆盖、从属 binding、可复制的参数重定义格式，以及给定规则的秩、零空间、完备性和可逆性审计；rules 缺省读取 input 或使用自动提案。";
DSParameterNotation::usage = "DSParameterNotation[context] 返回圈外 Gram 根号、独立无圈模长、全部必需模长覆盖、当前用户变量规则及 DSRedefineParameters 的可复制示例；无参数形式读取当前 context。";
DSRedefineParameters::usage = "DSRedefineParameters[context,rules] 用新的完整动力学变量规则重新初始化并返回新 context；rules 左端写 baseCoordinateOrder 中的 sp[原始动量,...]，右端写自定义参数表达式，不写 ssij->custom；DSRedefineParameters[rules] 只在成功后更新当前 context。";
DSSeeds::usage = "DSSeeds[context,opts] 生成全部 contact-reachable sector、全部离散态的符号 canonical general templates；连续指标域由 DSGenerateIBP 指定，数值系数由 DSLinear 的 CoefficientRules 指定。";
DSAllSeeds::usage = "DSAllSeeds[seedData] 取出 DSSeeds 返回的一维 allSeeds 模板列表；DSAllSeeds[] 取最近一次成功生成的模板。模板已遍历 n_i=0,1 并执行 EOM/canonical，连续指标仍为 general。";
DSSeedGroups::usage = "DSSeedGroups[seedData] 返回 DSSeeds 按 {sectorKey,ibpClass,generator} 排列的嵌套 seedGroups；每个二级组包含该 IBP 算符的全部离散态模板。无参数形式读取最近结果。";
DSSeedGroupMetadata::usage = "DSSeedGroupMetadata[seedData] 返回与 DSSeedGroups 同序的 sector、IBP 类型、生成元、模板 ordinal 和计数。无参数形式读取最近结果。";
DSMetaSeedRange::usage = "DSMetaSeedRange[seeds,{indices...}] 按用户给定的 seeds 外层结构初始化逐组 seed shift metadata；flat seeds 作为一组，nested seeds 的每个顶层元素作为一组。积分参数先统一 Flatten，再由 Variables 审计实际指标。声明多余或遗漏只 warning，metadata 仍按实际完整集合建立；再次初始化覆盖旧状态。DSMetaSeedRange[] 返回当前状态。";
DSGenerateIBP::usage = "DSGenerateIBP[seeds,{min,max}] 把统一目标积分指标包络按已初始化的逐输入组 shift metadata 反推成 seed 点域；DSGenerateIBP[seeds,{index,min,max},...] 要求 root 指标 exact cover。用户可显式给负 ISP 下界；自动反推不会把 ISP seed 下界降到用户 target 下界以下。先筛 parity，再代入数值点生成 IBP；n_i 不得再次撒点。";
DSLinear::usage = "DSLinear[seedData,context,opts] 把 canonical seeds 转换为 backend-neutral linearData。";
DSReorderIntegrals::usage = "DSReorderIntegrals[linearData,order] 按用户给定的 J 或现有积分 ID 显式重排 linearData；其结果的 integralList 是后续 plan/export/import 的唯一积分顺序。";
DSUserMI::usage = "DSUserMI[linearData,expressions,spec] 构造并附加有序 userMI basis；expressions 可为单个 J 或齐次 J 线性组合。DSUserMI[data] 或 DSUserMI[data,key] 查询同源可逆映射、秩和 backend ID。";
DSKiraPlan::usage = "DSKiraPlan[linearData,spec] 按 linearData 的既定积分顺序生成 preReduction/formal 两阶段 Kira 计划；formal 计划先解析构造 active basis 一阶导数及最小 target closure。";
DSKiraExport::usage = "DSKiraExport[linearData,opts] 或 DSKiraExport[kiraPlan] 序列化 Kira 基础输入和同源 manifest；不会启动 Kira。缺省禁止数值化微分变量，只有已冻结解析导数闭包的 formal plan 可显式选择 postDerivative 数值阶段。";
DSKiraImport::usage = "DSKiraImport[path,context,opts] 导入并验证完整 Kira reduction、master 顺序和积分双向映射；完成日志只作诊断，实际 artifact identity 与结构闭合是硬边界。";
DSDE::usage = "DSDE[reductionData,variables,opts] 用 ds 和 reduction rules 构造保持 master 顺序的微分方程矩阵。";
DSScaleCheck::usage = "DSScaleCheck[deData,spec,opts] 检查约化后的 Euler/标度关系；证书范围由 manifest 实际固定的 scaling variables 决定，返回 symbolic、fixedSection 或 exactPoint。";
DSMessagesOn::usage = "DSMessagesOn[] 开启 info、progress 和 warning 提醒。";
DSMessagesOff::usage = "DSMessagesOff[] 关闭可选提醒；fatal error 始终保留。";
DSMessagesQ::usage = "DSMessagesQ[] 返回可选提醒当前是否开启。";
DSPublicAPI::usage = "DSPublicAPI[] 返回当前版本需要用户掌握的公开函数、逻辑分组及各高层函数缺省选项；供手册附录和 example 覆盖门禁共同读取。";

WriteInitializationFiles::usage = "WriteInitializationFiles 是 DSInit 的选项；缺省 False。";
InitializationDirectory::usage = "InitializationDirectory 指定 init metadata 目录；缺省 Automatic，解析到调用脚本同目录的 init/。";
GenerateDerivativeMetadata::usage = "GenerateDerivativeMetadata 控制 DSInit 是否预生成独立变量微分算符 metadata；缺省 False。";
OverwriteInitialization::usage = "OverwriteInitialization 控制是否允许覆盖已有但输入哈希不一致的 init metadata；缺省 False。";
RegisterAsCurrent::usage = "RegisterAsCurrent 控制 DSInit 是否注册为无参公开接口的当前 context；缺省 True。";
ProgressReporting::usage = "ProgressReporting 控制单次高层调用的阶段进度；Automatic 跟随全局消息开关。";
KinematicRules::usage = "KinematicRules 是 DSInit 的可选动力学变量替换规则；缺省 Automatic 使用 input 中的 kinematicRules，若仍未给出则采用 DSKinematics 的缺省提案。";

PrecomputeShrinkSectorMetadata::usage = "PrecomputeShrinkSectorMetadata 控制底层是否预枚举 contact-reachable sector metadata。";
GenerateShrinkSectors::usage = "GenerateShrinkSectors 控制 canonical seed 是否生成 contact-reachable shrink sectors；高层完整工作流缺省 True。";
KiraOrdering::usage = "KiraOrdering 指定 backend-neutral linearData 的 Kira 排序约定。";
CoefficientRules::usage = "CoefficientRules 指定 linearData 层的小规模系数替换规则；DE 工作流不得直接或间接消去微分变量。";
LinearSystemMode::usage = "LinearSystemMode 选择 DSLinear 的 \"symbolic\" 或 \"numeric\" 模式。";
ExportKira::usage = "ExportKira 是底层组合工作流的导出开关；DSKiraExport 本身不运行 Kira。";
OutputDirectory::usage = "OutputDirectory 指定 serializer 输出目录；None 表示只返回内存数据。";
KiraCoefficientRules::usage = "KiraCoefficientRules 指定 Kira 导出前的系数规则；配置 KiraActiveBasis 时规则不得触及其 derivativeVariables。";
KiraTargetIntegrals::usage = "KiraTargetIntegrals 指定 Kira job 的目标积分。";
KiraActiveBasis::usage = "KiraActiveBasis 为 DSKiraExport 指定有序 active basis 线性组合、名称和导数变量；缺省 None。";
KiraNumericStage::usage = "KiraNumericStage 选择 \"symbolic\"（缺省，DE 变量不得数值化）或 \"postDerivative\"（只在 active basis 的解析一阶导数及 target closure 已构造后允许定点数值 reduction）。";
KiraJobOptions::usage = "KiraJobOptions 指定仅用于生成 Kira job 文件的选项 Association。";
KiraReductionFile::usage = "KiraReductionFile 指定 DSKiraImport 读取的 reduction 规则文件；缺省为 results/Tuserweight/kira_list.m。";
KiraMasterFile::usage = "KiraMasterFile 指定 DSKiraImport 读取的有序 master 文件；缺省为 results/Tuserweight/masters。";
KiraCompletionFile::usage = "KiraCompletionFile 指定 DSKiraImport 检查的完成日志；缺省为 kira.log。";
KiraCompletionPatterns::usage = "KiraCompletionPatterns 指定完成日志必须匹配的字符串或 RegularExpression 列表。";
ScalingRelation::usage = "ScalingRelation 指定 DSScaleCheck 使用的 \"Custom\"、\"LoopTopology\" 或 \"PureMassiveBubble\" 标度关系。";
ScalingVariables::usage = "ScalingVariables 指定 Euler 算符中的变量顺序。";
ScalingWeights::usage = "ScalingWeights 指定 Euler 算符中各变量的系数；ssij 与独立无圈模长 sEi 都是动量一次量，缺省物理权重为 1。";
ScalingDegrees::usage = "ScalingDegrees 指定各 master 的预期齐次次数；PureMassiveBubble 可设 Automatic。";

(* ::Chapter:: *)
(*私有实现加载*)

Begin["`Private`"];

$dSIBPPackageRoot = DirectoryName[DirectoryName[$InputFileName]];
$dSIBPVersion = "1.0";

Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "Core", "TopologyKinematics.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "Core", "LoopCore.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "Core", "Context.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "Core", "Metadata.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "IBP", "Loop.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "IBP", "Tree.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "IBP", "PureTime.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "IBP", "GenerateIBP.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "Backends", "KiraExport.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "Backends", "UserMI.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "Backends", "KiraPlan.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "Backends", "KiraImport.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "DE", "BuildDE.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "DE", "Scaling.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "Tree", "VertexFamily.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "Core", "KinematicCoordinates.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "Core", "SectorModel.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "Core", "TimeOnlyRepresentation.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "Core", "RelationsParity.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "Core", "ParameterInterface.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "Tree", "PublicBoundary.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{$dSIBPPackageRoot, "Kernel", "Core", "BilingualMessages.wl"}], CharacterEncoding -> "UTF-8"];


(* ::Chapter:: *)
(*加载完整性与引用提醒*)

(* 只有关键公开入口都获得实现后才宣布 package 可用；模块中途失败时直接停止，
   避免成功提示或引用提醒掩盖不完整加载。 *)
dsPackageLoadCompleteQ[] := And @@ (
    Length[DownValues[#]] > 0 & /@ {
      DSInit, DSSeeds, DSGenerateIBP, DSLinear, DSReorderIntegrals, DSUserMI,
      DSKiraPlan, DSKiraExport, DSKiraImport, DSDE, DSScaleCheck,
      DSTreeSeeds, DSTreeNaiveIBP, DSTreeNaiveDE, DSTreeDLogDE
      }
    );

If[! TrueQ[dsPackageLoadCompleteQ[]],
 Print["dSIBP package loading failed: one or more required modules are incomplete."];
 Abort[]
 ];

$dSIBPCitationEntries = {
   <|
    "authors" -> "Jiaqi Chen and Bo Feng",
    "title" -> "Towards Systematic Evaluation of de Sitter Correlators via Generalized Integration-By-Parts Relations",
    "arXiv" -> "2401.00129",
    "url" -> "https://arxiv.org/abs/2401.00129"
    |>,
   <|
    "authors" -> "Jiaqi Chen, Bo Feng and Yi-Xiao Tao",
    "title" -> "Multivariate hypergeometric solutions of cosmological (dS) correlators by d log-form differential equations",
    "arXiv" -> "2411.03088",
    "url" -> "https://arxiv.org/abs/2411.03088"
    |>,
   <|
    "authors" -> "Jiaqi Chen, Bo Feng, Zhehan Qin and Yi-Xiao Tao",
    "title" -> "Loop integrals in de Sitter spacetime: The parity-split IBP system and d log-form differential equations",
    "arXiv" -> "2604.14549",
    "url" -> "https://arxiv.org/abs/2604.14549"
    |>
   };

(* Notebook 使用真正的 Hyperlink；终端保留完整 URL，便于日志直接访问。 *)
dsCitationNotebookLine[index_Integer, entry_Association] := Row[{
    ToString[index] <> ". ", entry["authors"], ", \"", entry["title"], "\", ",
    Hyperlink["arXiv:" <> entry["arXiv"], entry["url"]], "."
    }];

dsCitationTextLine[index_Integer, entry_Association] := StringJoin[
   ToString[index], ". ", entry["authors"], ", \"", entry["title"], "\", arXiv:",
   entry["arXiv"], ", ", entry["url"], "."
   ];

If[! TrueQ[$dSIBPCitationReminderShown],
 $dSIBPCitationReminderShown = True;
 If[TrueQ[$Notebooks],
  Print@Column@Join[
     {"To support dSIBP, please acknowledge our work when you use dSIBP in your research. You can do so by citing the following papers:"},
     MapIndexed[dsCitationNotebookLine[First[#2], #1] &, $dSIBPCitationEntries],
     {"4. dSIBP package paper, arXiv identifier pending."}
     ],
  Print["To support dSIBP, please acknowledge our work when you use dSIBP in your research. You can do so by citing the following papers:"];
  MapIndexed[Print[dsCitationTextLine[First[#2], #1]] &, $dSIBPCitationEntries];
  Print["4. dSIBP package paper, arXiv identifier pending."]
  ]
 ];

End[];
EndPackage[];
