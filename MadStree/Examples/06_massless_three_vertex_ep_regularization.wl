(* ::Package:: *)

(***
文件：06_massless_three_vertex_ep_regularization.wl
用途：演示真实无质量三顶点树图的时间幂解析正规化 a1=-1+ep、a2=a3=1+ep；
      顶点 1 的幂取 -1+ep，使 eps=0 处时间积分真正发散，Laurent 展开含单一
      1/ep 极点：请求返回到 ep^0 时，极点系数 c_{-1} 与有限部 c_0 一并返回。
      数值阶段之前先做两项符号检查：任取两个相对主积分 tau 幂偏移 -1 的积分，
      其和由 MSReduce 一次（而非分别）约化并检查闭合；以 a1=-1+ep、a2=a3=1+ep
      作为 timePower 重新初始化后作代入约化，结果显式携带 ep 依赖。
      程序在数值 NDE 前从符号边界与 dlog DE 认证最低幂 n_min=-1，
      再选择 exact ep 点并重构 Laurent 系数。用户也可提供冗余生产候选池、独立验证点
      和首轮内部最高幂；程序只按每轮所需点数消费候选前缀并复用既有结果。本例缺省
      选择开角域 {-Pi/3,Pi/3}，模长仍由目标精度自动决定。
并行：ParallelTaskCount 缺省为 12；生产与独立验证任务都由有界进程池自动调度，
      任务多于并行上限时自动续交，不需要用户提供或手动分批 ep 取值。批量基准可在
      wolframscript 命令末尾传入一个正整数覆盖本例请求值，普通用户直接修改参数行即可。
输出：运行产物由缺省目录合同写入本 example 目录下 results_temp，不写入程序包源码模块目录。
***)

(* ::Chapter:: *)
(*程序包加载与用户参数*)

exampleDirectory = DirectoryName[$InputFileName];
packageRoot = DirectoryName[exampleDirectory];
AppendTo[$Path, FileNameJoin[{packageRoot, "Kernel"}]];
Needs["MadStree`"];

(* 用户可修改：缺省并行上限是 12；命令行末尾正整数只用于重复基准。 *)
epParallelOverrideText = If[
  ListQ[$ScriptCommandLine] && Length[$ScriptCommandLine] >= 2,
  Last[$ScriptCommandLine],
  ""
];
epParallelOverrideQ = StringQ[epParallelOverrideText] &&
  StringMatchQ[epParallelOverrideText, DigitCharacter ..] &&
  ToExpression[epParallelOverrideText] >= 1;
epParallelTaskCount = If[
  epParallelOverrideQ,
  ToExpression[epParallelOverrideText],
  12
];
maximumEpPower = 0;
epGoalDigits = 12;
useCustomEpGrid = False;
useComplexEpAngleRange = True;
epSampleAngleRange = {-Pi/3, Pi/3};
epProductionCandidates = {
  10/10^13, 9/10^13, 8/10^13, 7/10^13, 6/10^13,
  5/10^13, 4/10^13, 3/10^13, 2/10^13, 1/10^13
};
epValidationCandidates = {9/10^14, 8/10^14};
epInitialInternalMaximumPower = 2;


(* ::Chapter:: *)
(*无质量三顶点 family 与单顶点降幂 ep 正规化*)

treeSpec = <|
  "vertices" -> {
    <|"id" -> 1, "externalLegEnergy" -> k1, "timePower" -> a1, "vertexType" -> "+"|>,
    <|"id" -> 2, "externalLegEnergy" -> k2, "timePower" -> a2, "vertexType" -> "+"|>,
    <|"id" -> 3, "externalLegEnergy" -> k3, "timePower" -> a3, "vertexType" -> "+"|>
  },
  "lines" -> {
    <|"type" -> "massless", "endpoints" -> {1, 2},
      "momentum" -> q12, "nu" -> 1/2|>,
    <|"type" -> "massless", "endpoints" -> {2, 3},
      "momentum" -> q23, "nu" -> 1/2|>
  }
|>;

context = MSInitTree[treeSpec];
formulaArtifacts = MSWriteFormulaArtifacts[context];

(* ::Chapter:: *)
(*数值前的符号检查：组合迭代约化与代入约化*)

topKey = First[context["sectorOrder"]];

(* 任取两个相对主积分 tau 幂偏移 -1 的积分；MSReduce 一次作用于它们的和，
   而不是分别约化。 *)
shiftedCombination =
  MSIntegral[topKey, {-1, 0, 0}, {0, 0}] +
  MSIntegral[topKey, {0, -1, 0}, {1, 0}];
combinationReduction = MSReduce[shiftedCombination, context];
combinationReduction["status"]
combinationReduction["coefficientVector"]
combinationReduction["singularLayers"]

If[! (AssociationQ[combinationReduction] &&
      combinationReduction["status"] === "reduced" &&
      TrueQ[combinationReduction["nonMasterResidual"] === 0] &&
      combinationReduction["remainingShiftedIntegrals"] === {}),
  Print["Example failed at MSReduce: the combined reduction did not close."];
  Exit[1]
];

(* 代入正规化幂后的约化：约化使用 context 初始化时固定的时间幂，把
   a1=-1+ep、a2=a3=1+ep 直接作为 timePower 重新初始化同一拓扑，
   同一组合的约化结果即显式携带 ep 依赖与奇面。 *)
regulatedSpec = treeSpec /. {a1 -> -1 + ep, a2 -> 1 + ep, a3 -> 1 + ep};
regulatedContext = MSInitTree[regulatedSpec];
regulatedReduction = MSReduce[shiftedCombination, regulatedContext];
regulatedReduction["status"]
regulatedReduction["coefficientVector"]
regulatedReduction["singularLayers"]

If[! (AssociationQ[regulatedReduction] &&
      regulatedReduction["status"] === "reduced" &&
      TrueQ[regulatedReduction["nonMasterResidual"] === 0] &&
      regulatedReduction["remainingShiftedIntegrals"] === {}),
  Print["Example failed at MSReduce: the regulated reduction did not close."];
  Exit[1]
];

pointSequence = {{k1, k2, k3}, {9 I, 3 I, 5 I}};
parameterRules = {
  q12 -> 1, q23 -> 2,
  a1 -> -1 + ep, a2 -> 1 + ep, a3 -> 1 + ep
};


(* ::Chapter:: *)
(*自适应正规化级数重构*)

{epSeriesWallSeconds, epSeries} = AbsoluteTiming[MSReconstructEpSeries[
  context,
  ep,
  pointSequence,
  ParameterRules -> parameterRules,
  MaximumEpPower -> maximumEpPower,
  EpGoalDigits -> epGoalDigits,
  EpSamplePoints -> If[useCustomEpGrid, epProductionCandidates, Automatic],
  EpSampleAngleRange -> If[
    useComplexEpAngleRange && ! useCustomEpGrid,
    epSampleAngleRange,
    Automatic
  ],
  EpValidationPoints -> If[useCustomEpGrid, epValidationCandidates, Automatic],
  EpInitialInternalMaximumPower -> If[
    useCustomEpGrid,
    epInitialInternalMaximumPower,
    Automatic
  ],
  ParallelTaskCount -> epParallelTaskCount,
  FlintNDEPathPlanning -> True,
  BoundaryScale -> 4,
  BoundarySeriesOrder -> 24,
  RankOrder -> {1, 2, 3},
  MessageLanguage -> "CN",
  MSRuntimeDirectory -> Automatic
]];

If[Head[epSeries] === Failure,
  Print["Example failed at MSReconstructEpSeries: ", InputForm[epSeries]];
  Exit[1]
];


(* ::Chapter:: *)
(*结果与验收*)

checks = <|
  "threeVertices" -> Length[treeSpec["vertices"]] === 3,
  "twoMasslessPropagators" -> Length[treeSpec["lines"]] === 2 &&
    AllTrue[treeSpec["lines"], #["type"] === "massless" &] &&
    AllTrue[context["lines"], #["type"] === "masslessFull" &],
  "oneVertexLoweredEpRegularization" ->
    ({a1, a2, a3} /. parameterRules) === {-1 + ep, 1 + ep, 1 + ep},
  "analyticFormulaWritten" -> Lookup[formulaArtifacts, "status", None] === "generated" &&
    FileExistsQ[formulaArtifacts["files", "dlogDE"]] &&
    epSeries["formulaArtifacts", "outputDirectory"] === formulaArtifacts["outputDirectory"] &&
    Get[epSeries["formulaArtifacts", "files", "dlogDE"]] === MSDLogDE[context],
  "defaultParallelCountDocumented" ->
    MemberQ[Options[MSReconstructEpSeries], ParallelTaskCount -> 12] &&
    MemberQ[Options[MSReconstructEpSeries], MaximumEpPower -> 0] &&
    MemberQ[Options[MSReconstructEpSeries], EpSamplePoints -> Automatic] &&
    MemberQ[Options[MSReconstructEpSeries], EpSampleAngleRange -> Automatic] &&
    MemberQ[Options[MSReconstructEpSeries], EpValidationPoints -> Automatic] &&
    MemberQ[Options[MSReconstructEpSeries],
      EpInitialInternalMaximumPower -> Automatic],
  "requestedFinitePart" -> maximumEpPower === 0 &&
    epSeries["maximumPower"] === 0 && KeyExistsQ[epSeries["coefficients"], 0],
  "symbolicSupportCertifiedBeforeNDE" ->
    epSeries["laurentSupportCertificate", "status"] === "certified" &&
    epSeries["laurentSupportCertificate", "method"] ===
      "symbolicBoundaryAndRegularDLogValuation",
  "certifiedSinglePoleLeadingPower" -> epSeries["leadingPower"] === -1 &&
    epSeries["laurentSupportCertificate", "leadingPower"] === -1,
  "poleCoefficientReturned" -> KeyExistsQ[epSeries["coefficients"], -1] &&
    ! AllTrue[Lookup[epSeries["coefficients"], -1], TrueQ[# === 0] &],
  "regularDEAndMeromorphicBoundary" ->
    TrueQ[epSeries["laurentSupportCertificate",
      "differentialEquationRegularAtEpZeroQ"]] &&
    TrueQ[epSeries["laurentSupportCertificate",
      "boundaryFiniteLaurentAtEpZeroQ"]] &&
    epSeries["laurentSupportCertificate", "boundaryDefinitionAnalyticAtEpZeroQ"] ===
      False,
  "automaticProductionGrid" -> Length[epSeries["productionEpValues"]] >=
    epSeries["maximumPower"] - epSeries["leadingPower"] + 1,
  "openComplexAngleRange" -> If[
    useComplexEpAngleRange && ! useCustomEpGrid,
    epSeries["sampleAngleRange"] === epSampleAngleRange &&
      AllTrue[
        Arg /@ N[epSeries["productionEpValues"], 40],
        TrueQ[N[First[epSampleAngleRange], 40] < # <
          N[Last[epSampleAngleRange], 40]] &
      ] &&
      Length[DeleteDuplicates[Round[
        Arg /@ N[epSeries["productionEpValues"], 40], 10^-20]]] <= 3,
    epSeries["sampleAngleRange"] === Automatic
  ],
  "customCandidatePoolRespected" -> If[
    useCustomEpGrid,
    SubsetQ[epSeries["productionEpValues"], epProductionCandidates] &&
      epSeries["productionEpCandidateValues"] === epProductionCandidates,
    epSeries["productionEpCandidateValues"] === Automatic
  ],
  "independentValidationGrid" -> Length[epSeries["validationEpValues"]] >= 1 &&
    Intersection[epSeries["productionEpValues"], epSeries["validationEpValues"]] === {},
  "validationPassed" -> TrueQ[epSeries["precisionTargetMet"]] &&
    AllTrue[epSeries["validation"], TrueQ[#["passed"]] &],
  "effectiveParallelCount" ->
    epSeries["parallelTaskCountEffective"] === Min[
      epParallelTaskCount,
      Length[epSeries["productionEpValues"]] + Length[epSeries["validationEpValues"]]
    ],
  "allTargetsMet" -> AllTrue[
    Lookup[epSeries["pointEvaluations"], "evaluation"],
    TrueQ[#["flintNDE", "targetRelativeErrorMet"]] &
  ]
|>;

Print["symbolically certified ep powers: ", epSeries["leadingPower"], " through ",
  epSeries["maximumPower"]];
Print["Laurent support certificate: ",
  InputForm[KeyDrop[epSeries["laurentSupportCertificate"], "boundaryRecords"]]];
Print["automatic production ep values: ", epSeries["productionEpValues"]];
Print["independent validation ep values: ", epSeries["validationEpValues"]];
Print["parallel requested/effective: ",
  epSeries["parallelTaskCountRequested"], "/", epSeries["parallelTaskCountEffective"]];
Print["adaptive reconstruction wall seconds: ", epSeriesWallSeconds];
(* 展示前把绝对值低于 10^-20 的分量归零：本次重构相对残差约 10^-25、物理分量量级
   10^-5 以上，低于 10^-20 的分量是数值噪声（数学上为零），不是物理信号；返回记录
   epSeries 本身保持未截断的高精度原值。 *)
displayChopThreshold = 10^-20;
poleVector = Lookup[epSeries["poleCoefficients"], -1];
Print["pole coefficients: ",
  InputForm[Chop[N[poleVector], displayChopThreshold]]];
Print["finite part: ",
  InputForm[Chop[N[epSeries["finitePart"]], displayChopThreshold]]];
Print["maximum validation relative residual: ",
  epSeries["maximumValidationRelativeResidual"]];
Print["returned coefficient precision: ",
  Min[Precision /@ Flatten[Values[epSeries["coefficients"]]]]];
Print["checks: ", Count[Values[checks], True], "/", Length[checks]];

If[! And @@ Values[checks],
  Print["failed checks: ", Keys@Select[checks, Not@TrueQ[#] &]];
  Exit[1]
];
Print["Example PASSED: single-pole ep Laurent reconstruction completed."];
