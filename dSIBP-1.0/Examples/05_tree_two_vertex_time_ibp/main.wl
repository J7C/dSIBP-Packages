(* ::Package:: *)
(* pure-time 示例：公开积分使用 J[sectorKey,timeShifts,stateBits]。两顶点 massive-only family
   演示 naive tree IBP/DE 路线；atomic massless family 演示 quotient canonical 与 naive 路线的
   PendingRederivation 门禁。 *)

(* ::Chapter:: *)
(*加载标准 package*)

exampleDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{exampleDir, "..", "load_current_package.wl"}], CharacterEncoding -> "UTF-8"];


(* ::Chapter:: *)
(*详细物理输入*)

(* 图论圈数为零，因此缺省就是 timeOnly。这里仍显式写出模式，作为用户输入模板。
   p12 是 lineInput 中实际出现的独立无圈动量；程序只为其模长建立 sE1。 *)
treeCaseInput = <|
   "name" -> "treeTwoVertexPlusPlus",
   "vertices" -> {
     <|"id" -> v1, "vertexType" -> "+", "externalLegEnergy" -> E1|>,
     <|"id" -> v2, "vertexType" -> "+", "externalLegEnergy" -> E2|>
     },
   "lines" -> {
     <|"id" -> 1, "endpoints" -> {v1, v2}, "momentum" -> p12,
       "nu" -> nu12, "massType" -> "massive"|>
     },
   "loopMomenta" -> {},
   "loopExternalMomenta" -> {},
   "independentExternalMomenta" -> {p12},
   "ibpMode" -> "timeOnly",
   "ispData" -> {},
   "zeroPointRules" -> {a0[v1] -> alpha1, a0[v2] -> alpha2, b0[1] -> beta12},
   "symmetryRules" -> {}
   |>;

(* masslessT 在 Private relation 层约到二维物理基；公开 stateBits 中整边共享一位。
   DSSeeds 与 naive tree 结果都使用同一个 time-only J，不需要私有状态选项。 *)
atomicMasslessInput = <|
   "name" -> "atomicMasslessTimeOnly",
   "vertices" -> {
     <|"id" -> u1, "vertexType" -> "+", "externalLegEnergy" -> EM1|>,
     <|"id" -> u2, "vertexType" -> "+", "externalLegEnergy" -> EM2|>
     },
   "lines" -> {
     <|"id" -> 1, "endpoints" -> {u1, u2}, "momentum" -> pM,
       "massType" -> "massless"|>
     },
   "loopMomenta" -> {},
   "loopExternalMomenta" -> {},
   "independentExternalMomenta" -> {pM},
   "ibpMode" -> "timeOnly",
   "ispData" -> {},
   "zeroPointRules" -> {a0[u1] -> gamma1, a0[u2] -> gamma2, b0[1] -> deltaM},
   "symmetryRules" -> {}
   |>;


(* ::Chapter:: *)
(*缺省选项*)

(* 本例只消费内存 context，不在 example 目录保存初始化工作树。 *)
treeInitOptions = {
   WriteInitializationFiles -> False,
   GenerateDerivativeMetadata -> False,
   ProgressReporting -> True
   };


(* ::Chapter:: *)
(*原生 pure-time seeds 与 linearData*)

DSMessagesOn[];
treeContext = DSInit[treeCaseInput, Sequence @@ treeInitOptions];
treeSeedBatch = DSSeeds[treeContext, ProgressReporting -> True];
treeAllSeeds = DSAllSeeds[treeSeedBatch];
treeGeneratedIBP = DSGenerateIBP[treeAllSeeds, {-1, 1}];
treeLinearData = DSLinear[treeGeneratedIBP, treeContext, ProgressReporting -> True];

atomicContext = DSInit[
   atomicMasslessInput,
   RegisterAsCurrent -> True,
   WriteInitializationFiles -> False,
   GenerateDerivativeMetadata -> False,
   ProgressReporting -> True
   ];
atomicSeedBatch = DSSeeds[
   atomicContext,
   ProgressReporting -> True
   ];
atomicAllSeeds = DSAllSeeds[atomicSeedBatch];
atomicGeneratedIBP = DSGenerateIBP[atomicAllSeeds, {-1, 1}];
atomicLinearData = DSLinear[atomicGeneratedIBP, atomicContext, ProgressReporting -> True];
atomicStates = DeleteDuplicates@Flatten[
    Lookup[atomicAllSeeds, "discreteRules", {}]
    ];
atomicSectors = Sort@Lookup[Lookup[atomicSeedBatch, "sectorMetadataList", {}], "sectorKey", {}];


(* ::Chapter:: *)
(*naive tree IBP 与 DE*)

(* masters 缺省 Automatic：由内部 dlog master 判定给出同序归一化 masters，
   再对主积分求全时间导数并 LinearSolve 得到 naive IBP 约化关系。
   masters 是可选位置参数，显式写出 Automatic 以免选项规则被绑进该参数。 *)
treeNaiveIBP = DSTreeNaiveIBP[treeContext, Automatic, ProgressReporting -> True];
treeDEVariables = {E1, E2, sE1};

(* 把公开 naive IBP 结果直接传回 DSTreeNaiveDE；只在表示审计通过时跨越公开边界。 *)
treeNaiveDE = DSTreeNaiveDE[treeNaiveIBP, treeDEVariables, ProgressReporting -> True];

(* massless quotient 尚未在三槽基上重新推导，naive 路线返回 PendingRederivation 门禁。 *)
atomicNaiveStatus = Lookup[
   DSTreeNaiveIBP[atomicContext, Automatic, ProgressReporting -> True],
   {"status", "reason"},
   Missing["NotAvailable"]
   ];

(* 结构一致性检查：DE 与 IBP 共享同一批 master、矩阵为方阵且 DE 闭合（无残留对象）。
   它捕获的是两个消费路径的实现分叉，不是整体符号约定。
   符号约定的绝对锚在 check_tree_normalization_closure（缺省 Positive 记法下对角块
   -(1-2 nu1)/sE1）与 independent-benchmark 对 arXiv:2411.03088 两顶点 DE 的数值复核
   （含 off-diagonal top-to-child col5）里，不在本 example。 *)
treeNaiveMasterIntegrals = Lookup[treeNaiveIBP["masters"], "integral", {}];
treeNaiveMasterCount = Lookup[treeNaiveIBP, "masterCount", 0];
treeNaiveDEConsistent = TrueQ[
   Lookup[treeNaiveDE, "masterCount", -1] === treeNaiveMasterCount &&
    And @@ (Dimensions[#] === {treeNaiveMasterCount, treeNaiveMasterCount} & /@
       Values[Lookup[treeNaiveDE, "matrices", <||>]]) &&
    Flatten[Values[Lookup[treeNaiveDE, "residualObjects", <||>]]] === {}
   ];

summary = <|
   "version" -> dSIBP`$dSIBPVersion,
   "initStatus" -> Lookup[treeContext, "status", "missing"],
   "seedRepresentation" -> Lookup[treeSeedBatch, "representation", Missing["representation"]],
   "linearRepresentation" -> Lookup[treeLinearData, "representation", Missing["representation"]],
   "naiveIBPStatus" -> Lookup[treeNaiveIBP, "status", "missing"],
   "naiveMasterCount" -> treeNaiveMasterCount,
   "naiveMasterSource" -> Lookup[treeNaiveIBP, "masterSource", None],
   "naiveMasterRepresentation" -> MatchQ[treeNaiveMasterIntegrals, {J[_String, _, _] ..}],
   "masters" -> treeNaiveMasterIntegrals,
   "naiveDEStatus" -> Lookup[treeNaiveDE, "status", "missing"],
   "naiveDEConsistent" -> treeNaiveDEConsistent,
   "atomicInitStatus" -> Lookup[atomicContext, "status", "missing"],
   "atomicSeedStatus" -> Lookup[atomicSeedBatch, "status", "missing"],
   "atomicRepresentation" -> Lookup[atomicSeedBatch, "representation", Missing["representation"]],
   "atomicStates" -> atomicStates,
   "atomicSectors" -> atomicSectors,
   "atomicLinearStatus" -> Lookup[atomicLinearData, "status", "missing"],
   "atomicNaiveStatus" -> atomicNaiveStatus
   |>;

Print[InputForm[summary]];
If[! And[
    ToString[summary["version"]] === runtimeVersionToken,
    summary["initStatus"] === "initialized",
    summary["seedRepresentation"] === "J[sectorKey,timeShifts,stateBits]",
    summary["naiveIBPStatus"] === "solved",
    summary["naiveMasterCount"] > 0,
    summary["naiveMasterSource"] === "automaticDLog",
    summary["naiveMasterRepresentation"],
    summary["naiveDEStatus"] === "generated",
    summary["naiveDEConsistent"],
    summary["atomicInitStatus"] === "initialized",
    summary["atomicSeedStatus"] === "generated",
    summary["atomicRepresentation"] === "J[sectorKey,timeShifts,stateBits]",
    Sort[summary["atomicStates"]] === Sort[{
       n[1, 1] -> 0, n[1, 1] -> 1, n[1, 2] -> 0
       }],
    summary["atomicSectors"] === {"0", "1"},
    summary["atomicLinearStatus"] === "generated",
    atomicNaiveStatus === {"PendingRederivation", "masslessQuotientFormulaNotCertified"}
    ], Exit[1]];