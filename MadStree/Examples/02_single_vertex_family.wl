(* ::Package:: *)

(***
文件：02_single_vertex_family.wl
用途：演示单顶点函数族的紧凑/显式定义、解析公式保存、主积分与约化，以及统一 pointSequence 表格下的单点和多点数值计算。
核心逻辑：两个初始化 schema 构造同一 context；ParameterRules 只给一次固定参数，单点和多点只相差 pointSequence 的值行数。
运行：可在 Mathematica 前端逐节执行，或用 wolframscript -file 运行整个文件。
***)

(* ::Chapter:: *)
(*加载 MadStree*)

packageRoot = DirectoryName[DirectoryName[$InputFileName]];
AppendTo[$Path, FileNameJoin[{packageRoot, "Kernel"}]];
Needs["MadStree`"];


(* ::Chapter:: *)
(*定义单顶点函数族*)

(* 紧凑模型沿用论文/参考代码的 ki、nui 排列：首项给外腿指数参数和时间幂，
   后续同位置元素给各 h block 的动量与 nu。 *)
compactContext = MSInitVertexFamily[<|
  "ki" -> {k0, q},
  "nui" -> {a, nu},
  "hankelBranches" -> {2},
  "vertexType" -> "+"
|>];

(* 显式模型表达同一个函数族，适合需要逐 block 命名的输入。 *)
explicitContext = MSInitVertexFamily[<|
  "externalLegEnergy" -> k0,
  "timePower" -> a,
  "hBlocks" -> {
    <|"id" -> h1, "momentum" -> q, "nu" -> nu, "hankelBranch" -> 2|>
  },
  "exponentialBlocks" -> {},
  "vertexType" -> "+",
  "normalization" -> 1
|>];

topKey = First[compactContext["sectorOrder"]];
masters = MSMasterIntegrals[compactContext];
dlogDE = MSDLogDE[compactContext];
formulaArtifacts = MSWriteFormulaArtifacts[compactContext];

Lookup[masters, "integral"]
dlogDE["omegaPotential"] // MatrixForm
formulaArtifacts["files", "dlogDE"]


(* ::Chapter:: *)
(*约化一个标量线性组合*)

shiftedIntegral = MSIntegral[topKey, {1}, {0}];
shiftedCombination =
  2 shiftedIntegral - 3 MSIntegral[topKey, {-1}, {1}];
reduction = MSReduce[shiftedCombination, compactContext];

reduction["status"]
reduction["result"]


(* ::Chapter:: *)
(*统一的单点与多点数值接口*)

(* ParameterRules 只替换一次；pointSequence 首行固定坐标列，后续只写数值。 *)
parameterRules = {q -> 1, nu -> 1/2, a -> 2};
singlePointSequence = {
  {k0},
  {3 I}
};

(* 多点不切换接口，只在同一个表头之后继续追加等宽数值行。 *)
multipointSequence = {{k0}, {3 I}, {7 I/2}, {4 I}};

commonNumericalOptions = {
  ParameterRules -> parameterRules,
  FlintNDEPathPlanning -> True,
  BoundaryScale -> 4,
  WorkingPrecision -> 32,
  TransportOrder -> 64,
  ReferenceTransportOrder -> 88,
  TargetRelativeError -> "1e-16",
  MessageLanguage -> "CN"
};

singleEvaluation = MSEvaluatePath[
  compactContext,
  singlePointSequence,
  Sequence @@ commonNumericalOptions
];

multipointEvaluation = MSEvaluatePath[
  compactContext,
  multipointSequence,
  Sequence @@ commonNumericalOptions
];

singlePointRecord = If[
  AssociationQ[singleEvaluation],
  First[Lookup[singleEvaluation, "pointResults", {}], Missing["PointResultAbsent"]],
  Missing["EvaluationFailed"]
];
multipointFirstRecord = If[
  AssociationQ[multipointEvaluation],
  First[Lookup[multipointEvaluation, "pointResults", {}], Missing["PointResultAbsent"]],
  Missing["EvaluationFailed"]
];
commonPointRelativeDifference = If[
  AssociationQ[singlePointRecord] && AssociationQ[multipointFirstRecord],
  Max[Abs[singlePointRecord["value"] - multipointFirstRecord["value"]]]/
    Max[Max[Abs[singlePointRecord["value"]]], 10^-30],
  Infinity
];

singlePointRecord
Lookup[multipointEvaluation["pointResults"], {"coordinate", "status", "userIndex"}]
commonPointRelativeDifference


(* ::Chapter:: *)
(*验收门禁*)

exampleChecks = <|
  "compactContext" -> MSContextQ[compactContext],
  "explicitContext" -> MSContextQ[explicitContext],
  "analyticFormulaWritten" -> Lookup[formulaArtifacts, "status", None] === "generated" &&
    FileExistsQ[formulaArtifacts["files", "dlogDE"]],
  "schemaEquivalentMasters" ->
    Lookup[MSMasterIntegrals[compactContext], "integral"] ===
      Lookup[MSMasterIntegrals[explicitContext], "integral"],
  "schemaEquivalentDE" ->
    MSDLogDE[compactContext]["connection"] === MSDLogDE[explicitContext]["connection"],
  "scalarCombinationReduced" -> Lookup[reduction, "status", None] === "reduced" &&
    reduction["nonMasterResidual"] === 0 &&
    reduction["remainingShiftedIntegrals"] === {},
  "singleComputed" -> Lookup[singleEvaluation, "status", None] === "computed",
  "multipointComputed" -> Lookup[multipointEvaluation, "status", None] === "computed",
  "singlePointCount" -> Length[Lookup[singleEvaluation, "pointResults", {}]] === 1,
  "multipointCount" -> Length[Lookup[multipointEvaluation, "pointResults", {}]] === 3,
  "commonPointCoordinate" ->
    Lookup[singlePointRecord, "coordinate", None] ===
      Lookup[multipointFirstRecord, "coordinate", None],
  "commonPointAgreement" -> TrueQ[commonPointRelativeDifference < 10^-20]
|>;

If[! And @@ Values[exampleChecks],
  Print["Example FAILED: ", InputForm[Select[exampleChecks, Not]]];
  Exit[1]
];
Print["Example PASSED: ", Count[Values[exampleChecks], True], "/", Length[exampleChecks],
  " checks; single and multipoint calls use the same pointSequence schema."];
