(* ::Package:: *)
(***
File: kira_result_summary.wl
Purpose: Readable record of the external Kira reduction, DE and scaling outcome for the mixed
bubble+tree example.
Boundary: No reduction table, log, database or cache is kept here.
***)


(* ::Chapter:: *)
(*Formal result record*)

<|
  "schema" -> "dsibp_example_kira_result_summary_v1",
  "status" -> "passed",
  "case" -> "mixBubbleTree",
  "packageVersion" -> "1.0",
  "branch" -> "+++",
  "kira" -> <|
    "version" -> "2.3 (Git: 2.3-7-geb541f9)",
    "runtime" -> "WSL",
    "parallelConfig" -> "w10*1",
    "wallTimeSeconds" -> 160.11,
    "exitStatus" -> 0
  |>,
  "masterIDs" -> Range[81],
  "masterCount" -> 81,
  "targetCount" -> 676,
  "selectedEquationCount" -> 63410,
  "unreducedCount" -> 0,
  "deVariables" -> {loopScale, legScale1, legScale2, E1, E2, E3},
  "deDimensions" -> ConstantArray[{81, 81}, 6],
  "sourceIdentities" -> {81, 81},
  "scalingCertificateScope" -> "exactPoint",
  "scalingSymbolicQ" -> False,
  "matrixResidualNonzeroCount" -> 0,
  "sourceResidualNonzeroCount" -> 0,
  "ordinaryExportedIntegralsCheckedForBridgePack" -> 310688,
  "fixedTreeBridgeExcludedFromParity" -> True,
  "intermediateArtifactsRetainedInRepository" -> False
|>
