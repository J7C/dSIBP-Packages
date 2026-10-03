(* ::Package:: *)

(***
文件：06_massless_three_vertex_ep_regularization.wl
用途：演示真实无质量三顶点树图的共同时间幂解析正规化 a1=a2=a3=1+ep；
      用户只指定需要返回到 ep^0，程序在数值 NDE 前从符号边界与 dlog DE 认证最低幂，
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
(*无质量三顶点 family 与共同 ep 正规化*)

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
pointSequence = {{k1, k2, k3}, {9 I, 3 I, 5 I}};
parameterRules = {
  q12 -> 1, q23 -> 2,
  a1 -> 1 + ep, a2 -> 1 + ep, a3 -> 1 + ep
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
  "sharedEpRegularization" ->
    ({a1, a2, a3} /. parameterRules) === ConstantArray[1 + ep, 3],
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
  "automaticLeadingPowerZero" -> epSeries["leadingPower"] === 0 &&
    epSeries["laurentSupportCertificate", "leadingPower"] === 0,
  "regularDEAndAnalyticBoundary" ->
    TrueQ[epSeries["laurentSupportCertificate",
      "differentialEquationRegularAtEpZeroQ"]] &&
    TrueQ[epSeries["laurentSupportCertificate",
      "boundaryFiniteLaurentAtEpZeroQ"]] &&
    TrueQ[epSeries["laurentSupportCertificate",
      "boundaryDefinitionAnalyticAtEpZeroQ"]],
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
Print["pole coefficients: ", InputForm[epSeries["poleCoefficients"]]];
Print["finite part: ", InputForm[epSeries["finitePart"]]];
Print["maximum validation relative residual: ",
  epSeries["maximumValidationRelativeResidual"]];
Print["checks: ", Count[Values[checks], True], "/", Length[checks]];

If[! And @@ Values[checks],
  Print["failed checks: ", Keys@Select[checks, Not@TrueQ[#] &]];
  Exit[1]
];
Print["Example PASSED: adaptive shared-ep Laurent reconstruction completed."];
