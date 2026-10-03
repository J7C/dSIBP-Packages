(* ::Package:: *)

(***
文件：05_massive_three_vertex_tree.wl
用途：演示三顶点、两条 massive 传播子的树图从拓扑初始化到主积分、递推、dlog DE 和多点数值求值。
数值分组：原非零点列与中间顶点外腿能量 k2=0 的点列共用同一 context，但分别规划、求值和导出。
结构：1 --(q12)-- 2 --(q23)-- 3；传播子按输入顺序在内部编号为 1、2。
nu 约定：取 nu12=3/4、nu23=1/3，避免半整数引起的表示退化。
运行：可在 Mathematica 前端逐节执行，或用 wolframscript -file 运行整个文件。
***)

(* ::Chapter:: *)
(* Load MadStree *)

packageRoot = DirectoryName[DirectoryName[$InputFileName]];
AppendTo[$Path, FileNameJoin[{packageRoot, "Kernel"}]];
Needs["MadStree`"];


(* ::Chapter:: *)
(* Define ordered tree topology (+++ contour vertices; internal type is derived as massiveFull) *)

treeSpec = <|
  "vertices" -> {
    <|"id" -> 1, "externalLegEnergy" -> k1, "timePower" -> a1, "vertexType" -> "+"|>,
    <|"id" -> 2, "externalLegEnergy" -> k2, "timePower" -> a2, "vertexType" -> "+"|>,
    <|"id" -> 3, "externalLegEnergy" -> k3, "timePower" -> a3, "vertexType" -> "+"|>
  },
  "lines" -> {
    <|"type" -> "massive", "endpoints" -> {1, 2},
      "momentum" -> q12, "nu" -> nu12|>,
    <|"type" -> "massive", "endpoints" -> {2, 3},
      "momentum" -> q23, "nu" -> nu23|>
  }
|>;

context = MSInitTree[treeSpec];
topKey = First[context["sectorOrder"]];
MSSectors[context]


(* ::Chapter:: *)
(* Direct formula results *)

masters = MSMasterIntegrals[context];
topMatrices = MSFormulaMatrices[context, topKey];
contactMaps = MSContactMaps[context, topKey];
dlogDE = MSDLogDE[context];
formulaArtifacts = MSWriteFormulaArtifacts[context];

(* 每条记录直接显示 normalized master、裸积分标签和精确 normalization。 *)
KeyTake[#, {"integral", "normalization", "bareIntegral", "definition"}] & /@ masters

dlogDE["omegaPotential"] // MatrixForm
formulaArtifacts["files", "dlogDE"]


(* ::Chapter:: *)
(* Iterative reduction and automatic numerical boundary *)

shiftedIntegral = MSIntegral[topKey, {1, 0, 0}, {0, 0, 0, 0}];
shiftedDefinition = MSIntegralDefinition[shiftedIntegral, context];
reduction = MSReduce[shiftedIntegral, context];
topStates = Take[First[MSSectors[context]]["stateOrder"], 2];
shiftedTable = Table[
  MSIntegral[topKey, {shift, 0, 0}, stateBits],
  {shift, {1, 2}},
  {stateBits, topStates}
];
tableReduction = MSReduce[shiftedTable, context];
shiftedDefinition

reduction["result"]
tableReduction["coefficientVector"]

numericalTemplate = MSNumericalSystem[dlogDE];
numericalTemplate["status"]

targetRules = {
  k1 -> 9 I, k2 -> 3 I, k3 -> 5 I,
  q12 -> 1, q23 -> 2,
  nu12 -> 3/4, nu23 -> 1/3,
  a1 -> 1, a2 -> 1, a3 -> 1
};
parameterRules = {
  k3 -> 5 I, q12 -> 1, q23 -> 2,
  nu12 -> 3/4, nu23 -> 1/3,
  a1 -> 1, a2 -> 1, a3 -> 1
};
singlePointSequence = {{k1, k2}, {9 I, 3 I}};
boundary = MSBoundaryData[
  context,
  targetRules,
  BoundaryScale -> 4,
  WorkingPrecision -> 32
];

targetValue = MSEvaluatePath[
  context,
  singlePointSequence,
  ParameterRules -> parameterRules,
  BoundaryScale -> 4,
  WorkingPrecision -> 32,
  TransportOrder -> 72,
  ReferenceTransportOrder -> 96,
  TargetRelativeError -> "1e-14"
];

(* Failure gate: stop immediately if the transport failed, so the batch
   section (which reuses targetValue["values"]) does not cascade on a
   Failure. The package has already printed the backend diagnostic message. *)
If[Head[targetValue] === Failure,
  Print["Example failed at MSEvaluatePath: ", targetValue];
  Exit[1]
];

targetValue["values"]
Lookup[targetValue["flintNDE", "segments"], "relativeDifferenceInf"]


(* ::Chapter:: *)
(* Multipoint evaluation and export *)

batchParameterRules = parameterRules;
pointSequence = {
  {k1, k2},
  {8 I, 2 I},
  {7 I, 4 I},
  {6 I, 3 I},
  {9 I, 5 I}
};

batchEvaluation = MSEvaluatePath[
  context,
  pointSequence,
  ParameterRules -> batchParameterRules,
  FlintNDEPathPlanning -> True,
  BoundaryScale -> 4,
  WorkingPrecision -> 32,
  TransportOrder -> 72,
  ReferenceTransportOrder -> 96,
  TargetRelativeError -> "1e-14"
];
batchPointResults = If[
  AssociationQ[batchEvaluation],
  Lookup[batchEvaluation, "pointResults", {}],
  {}
];

Lookup[batchPointResults, "value"]
Lookup[batchPointResults, "status"]

batchExport = If[
  Head[batchEvaluation] === Failure,
  batchEvaluation,
  MSExportEvaluationData[
    batchEvaluation,
    MSOutputDirectory -> "results/madstree_evaluation_05",
    SignificantDigits -> 16
  ]
];
batchExport


(* ::Chapter:: *)
(*中间顶点外腿能量为零的独立点列*)

(* 这组点不与原 pointSequence 混合：坐标 k2 虽为零，是否为 DE 奇点仍由完整
   top/contact-sector letters 决定，不能按坐标名称或零值预判。 *)
zeroMiddleEnergyPointSequence = {
  {k1, k2},
  {8 I, 0},
  {7 I, 0},
  {6 I, 0},
  {9 I, 0}
};
zeroMiddleEnergyCoordinateRules =
  Thread[First[zeroMiddleEnergyPointSequence] -> #] & /@
    Rest[zeroMiddleEnergyPointSequence];
zeroMiddleEnergyLetterValues =
  Simplify[dlogDE["letters"] /. Join[#, parameterRules]] & /@
    zeroMiddleEnergyCoordinateRules;
zeroMiddleEnergyZeroLetters =
  Select[#, TrueQ[Simplify[#] === 0] &] & /@ zeroMiddleEnergyLetterValues;

zeroMiddleEnergyEvaluation = MSEvaluatePath[
  context,
  zeroMiddleEnergyPointSequence,
  ParameterRules -> parameterRules,
  FlintNDEPathPlanning -> True,
  BoundaryScale -> 4,
  WorkingPrecision -> 32,
  TransportOrder -> 72,
  ReferenceTransportOrder -> 96,
  TargetRelativeError -> "1e-14"
];
zeroMiddleEnergyPointResults = If[
  AssociationQ[zeroMiddleEnergyEvaluation],
  Lookup[zeroMiddleEnergyEvaluation, "pointResults", {}],
  {}
];
zeroMiddleEnergyClassifications = If[
  AssociationQ[zeroMiddleEnergyEvaluation],
  Lookup[zeroMiddleEnergyEvaluation, "singularityClassifications", {}],
  Missing["EvaluationFailed"]
];

Lookup[zeroMiddleEnergyPointResults, {"coordinate", "status", "value"}]
zeroMiddleEnergyZeroLetters
zeroMiddleEnergyClassifications

zeroMiddleEnergyExport = If[
  Head[zeroMiddleEnergyEvaluation] === Failure,
  zeroMiddleEnergyEvaluation,
  MSExportEvaluationData[
    zeroMiddleEnergyEvaluation,
    MSOutputDirectory -> "results/madstree_evaluation_05_zero_middle_energy",
    SignificantDigits -> 16
  ]
];
zeroMiddleEnergyExport


(* ::Chapter:: *)
(* Failure gate *)

(* Exit non-zero when any stage produced a Failure, so a fresh run cannot
   report success while hiding errors. *)
exampleGateFailures = Select[
  {
    formulaArtifacts,
    tableReduction,
    targetValue,
    batchEvaluation,
    batchExport,
    zeroMiddleEnergyEvaluation,
    zeroMiddleEnergyExport
  },
  Head[#] === Failure &
];
If[exampleGateFailures =!= {},
  Print["Example FAILED: ", Length[exampleGateFailures], " result(s) are Failure objects."];
  Scan[Print["  ", #] &, exampleGateFailures];
  Exit[1]
];
If[
  Lookup[tableReduction, "status", None] =!= "reduced" ||
    tableReduction["nonMasterResidual"] =!= ConstantArray[0, {2, 2}] ||
    tableReduction["remainingShiftedIntegrals"] =!= ConstantArray[{}, {2, 2}],
  Print["Example FAILED: Table reduction did not close element by element."];
  Exit[1]
];
zeroMiddleEnergyChecks = <|
  "fourSeparatePhysicalPoints" -> Length[zeroMiddleEnergyPointResults] === 4,
  "middleExternalLegEnergyIsZero" ->
    And @@ (Last[#] === 0 & /@ Rest[zeroMiddleEnergyPointSequence]),
  "noDLogLetterVanishes" ->
    And @@ (# === {} & /@ zeroMiddleEnergyZeroLetters),
  "allPointsSaved" ->
    Lookup[zeroMiddleEnergyPointResults, "status"] === ConstantArray["saved", 4],
  "allPointsMeetReferencePrecision" ->
    And @@ Lookup[
      Lookup[zeroMiddleEnergyPointResults, "flintNDE", <||>],
      "targetRelativeErrorMet",
      False
    ],
  "noUserPointMisclassifiedAsSingular" -> zeroMiddleEnergyClassifications === {}
|>;
If[! And @@ Values[zeroMiddleEnergyChecks],
  Print["Example FAILED: ", InputForm[Select[zeroMiddleEnergyChecks, Not]]];
  Exit[1]
];
Print[
  "Example PASSED: all checked results are non-Failure; zero-middle-energy checks ",
  Count[Values[zeroMiddleEnergyChecks], True], "/", Length[zeroMiddleEnergyChecks], "."
]
