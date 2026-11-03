(* ::Package:: *)

(***
File: 03_time_only_cycle_chart.wl
Purpose: Demonstrates time-only cycle initialization, the common-theta contact sector, the dlog DE and all strict time-rank chart certificates.
Structure: representative public calls only; validation assertions and private helper checks are not part of the example.
Run: execute section by section in the Mathematica front end, or run the whole file with wolframscript -file.
***)

(* ::Chapter:: *)
(* Load MadStree *)

packageRoot = DirectoryName[DirectoryName[$InputFileName]];
AppendTo[$Path, FileNameJoin[{packageRoot, "Kernel"}]];
Needs["MadStree`"];


(* ::Chapter:: *)
(* Initialize the triangle time-only graph *)

spec = <|
  "vertices" -> {
    <|"id" -> 1, "externalLegEnergy" -> k1, "timePower" -> a1, "vertexType" -> "+"|>,
    <|"id" -> 2, "externalLegEnergy" -> k2, "timePower" -> a2, "vertexType" -> "+"|>,
    <|"id" -> 3, "externalLegEnergy" -> k3, "timePower" -> a3, "vertexType" -> "+"|>
  },
  "lines" -> {
    <|"type" -> "massless", "endpoints" -> {1, 2}, "momentum" -> q12, "nu" -> 1/2|>,
    <|"type" -> "massless", "endpoints" -> {2, 3}, "momentum" -> q23, "nu" -> 1/2|>,
    <|"type" -> "massless", "endpoints" -> {3, 1}, "momentum" -> q31, "nu" -> 1/2|>
  }
|>;

context = MSInitTimeGraph[spec];
sectors = MSSectors[context];
de = MSDLogDE[context];
formulaArtifacts = MSWriteFormulaArtifacts[context];

Lookup[sectors, {"sectorKey", "contractedLineIds", "masterCount"}]
de["dlogStatus"]
formulaArtifacts["files", "dlogDE"]


(* ::Chapter:: *)
(* Build all strict-rank chart certificates *)

targetRules = {
  k1 -> 11 I, k2 -> 7 I, k3 -> 5 I,
  a1 -> 0, a2 -> 0, a3 -> 0,
  q12 -> 2, q23 -> 3, q31 -> 4
};

certificate = MSBoundaryChartCertificate[
  context,
  targetRules,
  RankOrder -> All
];

Lookup[certificate["charts"], {"rankOrder", "normalCrossingQ"}]


(* ::Chapter:: *)
(* Generate the Frobenius boundary of one ordered chart *)

boundary = MSBoundaryData[
  context,
  targetRules,
  BoundaryScale -> 3,
  BoundarySeriesOrder -> 12,
  RankOrder -> {1, 2, 3},
  WorkingPrecision -> 30
];

boundary["method"]
Length[boundary["leadingBranches"]]


(* ::Chapter:: *)
(* Failure gate *)

(* Exit non-zero when any stage produced a Failure, so a fresh run cannot
   report success while hiding errors. *)
exampleGateFailures = Select[
  {context, de, formulaArtifacts, certificate, boundary},
  Head[#] === Failure &
];
If[exampleGateFailures =!= {},
  Print["Example FAILED: ", Length[exampleGateFailures], " result(s) are Failure objects."];
  Scan[Print["  ", #] &, exampleGateFailures];
  Exit[1]
];
Print["Example PASSED: all checked results are non-Failure."]
