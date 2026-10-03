(* ::Package:: *)

(***
File: 01_massless_full_edge.wl
Purpose: Demonstrates a single theta-carrying massless full edge from topology initialization to master integrals, recurrence and the dlog DE, followed by the unique two-phase single-point and multipoint numerical workflow.
Run: execute section by section in the Mathematica front end, or run the whole file with wolframscript -file.
***)

(* ::Chapter:: *)
(* Load MadStree *)

packageRoot = DirectoryName[DirectoryName[$InputFileName]];
AppendTo[$Path, FileNameJoin[{packageRoot, "Kernel"}]];
Needs["MadStree`"];


(* ::Chapter:: *)
(* Define ordered tree topology *)

treeSpec = <|
  "vertices" -> {
    <|"id" -> 1, "externalLegEnergy" -> k1, "timePower" -> a1, "vertexType" -> "+"|>,
    <|"id" -> 2, "externalLegEnergy" -> k2, "timePower" -> a2, "vertexType" -> "+"|>
  },
  "lines" -> {
    <|"type" -> "massless", "endpoints" -> {1, 2},
      "momentum" -> q, "nu" -> 1/2|>
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

Lookup[masters, "integral"]

dlogDE["omegaPotential"] // MatrixForm
formulaArtifacts["files", "dlogDE"]


(* ::Chapter:: *)
(* Iterative reduction and automatic numerical boundary *)

shiftedIntegral = MSIntegral[topKey, {1, 0}, {0}];
reduction = MSReduce[shiftedIntegral, context];
reduction["result"]

numericalTemplate = MSNumericalSystem[dlogDE];
numericalTemplate["status"]

targetRules = {k1 -> 9 I, k2 -> 3 I, q -> 1, a1 -> 1, a2 -> 1};
parameterRules = {q -> 1, a1 -> 1, a2 -> 1};
singlePointSequence = {{k1, k2}, {9 I, 3 I}};
boundary = MSBoundaryData[
  context,
  targetRules,
  BoundaryScale -> 4,
  WorkingPrecision -> 40
];

targetValue = MSEvaluatePath[
  context,
  singlePointSequence,
  ParameterRules -> parameterRules,
  BoundaryScale -> 4,
  WorkingPrecision -> 40,
  TransportOrder -> 80,
  ReferenceTransportOrder -> 104,
  TargetRelativeError -> "1e-20"
];

(* Failure gate: stop immediately if the transport failed, so later sections
   do not cascade on a Failure. The package has already printed the backend
   diagnostic message. *)
If[Head[targetValue] === Failure,
  Print["Example failed at MSEvaluatePath: ", targetValue];
  Exit[1]
];

targetValue["values"]
Lookup[targetValue["flintNDE", "segments"], "relativeDifferenceInf"]


(* ::Chapter:: *)
(* Multipoint evaluation and export *)

(* Fixed parameters shared by all points; each row gives {k1,k2}. Bare points
   are saved by default, so the execution result directly contains all values. *)
batchParameterRules = {q -> 1, a1 -> 1, a2 -> 1};
pointSequence = {
  {k1, k2},
  {8 I, 2 I},
  {7 I, 4 I},
  {6 I, 5 I}
};

batchEvaluation = MSEvaluatePath[
  context,
  pointSequence,
  ParameterRules -> batchParameterRules,
  FlintNDEPathPlanning -> True,
  BoundaryScale -> 4,
  WorkingPrecision -> 40,
  TransportOrder -> 80,
  ReferenceTransportOrder -> 104,
  TargetRelativeError -> "1e-20"
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
    MSOutputDirectory -> "results/madstree_evaluation",
    SignificantDigits -> 16
  ]
];
batchExport

(* ::Chapter:: *)
(* Failure gate *)

(* Exit non-zero when any stage produced a Failure, so a fresh run cannot
   report success while hiding errors. *)
exampleGateFailures = Select[
  {formulaArtifacts, targetValue, batchEvaluation, batchExport},
  Head[#] === Failure &
];
If[exampleGateFailures =!= {},
  Print["Example FAILED: ", Length[exampleGateFailures], " result(s) are Failure objects."];
  Scan[Print["  ", #] &, exampleGateFailures];
  Exit[1]
];
Print["Example PASSED: all checked results are non-Failure."]
