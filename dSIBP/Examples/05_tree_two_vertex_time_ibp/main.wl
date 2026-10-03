(* ::Package:: *)
(* pure-time 示例：公开积分使用 J[sectorKey,timeShifts,stateBits]。两顶点 massive-only family 演示
   Private 论文公式适配；atomic massless family 演示 quotient canonical 与公式路线的明确边界。 *)

(* ::Chapter:: *)
(*加载标准 package*)

exampleDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{exampleDir, "..", "load_current_package.wl"}], CharacterEncoding -> "UTF-8"];


(* ::Chapter:: *)
(*详细物理输入*)

(* 图论圈数为零，因此缺省就是 timeOnly。这里仍显式写出模式，作为用户输入模板。
   p12 是 lineInput 中实际出现的独立无圈动量；程序只为其模长建立 sE1。 *)
treeCaseInput = <|
   "name" -> "020TreeTwoVertexPlusPlus",
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

(* masslessFull 在 Private relation 层约到二维物理基；公开 stateBits 中整边共享一位。
   DSSeeds 与 tree 结果都使用同一个 time-only J，不需要私有状态选项。 *)
atomicMasslessInput = <|
   "name" -> "020AtomicMasslessTimeOnly",
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

(* 缺省：repIterative 按当前 sector 的 active 顶点把全部 a 约到 0。 *)
treeReductionEndpoint = Automatic;


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

selectedIntegral = J["1", {1, 0}, {1, 0}];
selectedSeed = DSTreeSeeds[v1, selectedIntegral, treeContext];


(* ::Chapter:: *)
(*迭代约化、naive IBP/DE 与公式型 dlog DE*)

treeDLog = DSTreeDLogDE[treeContext];
treeTarget = First[treeDLog["masters"]]["integral"];
treeReduction = repIterative[treeTarget, treeReductionEndpoint, treeContext];

treeDEVariables = {E1, E2, sE1};
treeNaiveIBP = DSTreeNaiveIBP[treeContext, treeDLog["masters"], ProgressReporting -> True];
treeNaiveDE = DSTreeNaiveDE[treeNaiveIBP, treeDEVariables, ProgressReporting -> True];
treeDEResiduals = Association@Table[
    variable -> (Together /@ Flatten[treeNaiveDE["matrices", variable] - D[treeDLog["omega"], variable]]),
    {variable, treeDEVariables}
    ];
treeDERoutesAgree = And @@ Flatten[Map[TrueQ[# === 0] &, Values[treeDEResiduals], {2}]];
atomicFormulaStatus = Lookup[
   DSTreeDLogDE[atomicContext],
   {"status", "reason"},
   Missing["NotAvailable"]
   ];

summary = <|
   "version" -> dSIBP`$dSIBPVersion,
   "initStatus" -> Lookup[treeContext, "status", "missing"],
   "seedRepresentation" -> Lookup[treeSeedBatch, "representation", Missing["representation"]],
   "linearRepresentation" -> Lookup[treeLinearData, "representation", Missing["representation"]],
   "selectedSeedRoute" -> Lookup[selectedSeed, "generationRoute", Missing["route"]],
   "iterativeReductionFreeOfFailure" -> FreeQ[treeReduction, $Failed],
   "dlogStatus" -> Lookup[treeDLog, "status", "missing"],
   "naiveIBPStatus" -> Lookup[treeNaiveIBP, "status", "missing"],
   "naiveDEStatus" -> Lookup[treeNaiveDE, "status", "missing"],
   "deRoutesAgree" -> treeDERoutesAgree,
   "masters" -> Lookup[treeDLog, "masters", {}],
   "atomicInitStatus" -> Lookup[atomicContext, "status", "missing"],
   "atomicSeedStatus" -> Lookup[atomicSeedBatch, "status", "missing"],
   "atomicRepresentation" -> Lookup[atomicSeedBatch, "representation", Missing["representation"]],
   "atomicStates" -> atomicStates,
   "atomicSectors" -> atomicSectors,
   "atomicLinearStatus" -> Lookup[atomicLinearData, "status", "missing"],
   "atomicFormulaStatus" -> atomicFormulaStatus
   |>;

Print[InputForm[summary]];
If[! And[
    ToString[summary["version"]] === currentVersion,
    summary["initStatus"] === "initialized",
    summary["seedRepresentation"] === "J[sectorKey,timeShifts,stateBits]",
    summary["selectedSeedRoute"] === "directPureTime",
    summary["iterativeReductionFreeOfFailure"],
    summary["dlogStatus"] === "generated",
    summary["naiveIBPStatus"] === "solved",
    summary["naiveDEStatus"] === "generated",
    summary["deRoutesAgree"],
    summary["atomicInitStatus"] === "initialized",
    summary["atomicSeedStatus"] === "generated",
    summary["atomicRepresentation"] === "J[sectorKey,timeShifts,stateBits]",
    Sort[summary["atomicStates"]] === Sort[{
       n[1, 1] -> 0, n[1, 1] -> 1, n[1, 2] -> 0
       }],
    summary["atomicSectors"] === {"0", "1"},
    summary["atomicLinearStatus"] === "generated",
    atomicFormulaStatus === {"PendingRederivation", "masslessQuotientFormulaNotCertified"}
    ], Exit[1]];
