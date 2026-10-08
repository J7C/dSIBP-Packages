(* ::Package:: *)
(***
File: kira_result_summary.wl
Purpose: Readable record of the external Kira reduction, DE and scaling outcome for the pure
massive bubble example.
Boundary: No reduction table, log, database or cache is kept here.
***)


(* ::Chapter:: *)
(*Formal result record*)

<|
  "schema" -> "dsibp_example_kira_result_summary_v1",
  "status" -> "passed",
  "case" -> "pureMassiveBubble",
  "packageVersion" -> "1.0",
  "kira" -> <|
    "version" -> "2.3 (Git: 2.3-7-geb541f9)",
    "runtime" -> "WSL",
    "parallelConfig" -> "w10*1",
    "wallTimeSeconds" -> 109.09,
    "exitStatus" -> 0
  |>,
  "masterIDs" -> Range[19],
  "masterCount" -> 19,
  "targetCount" -> 300,
  "selectedEquationCount" -> 2179,
  "unreducedCount" -> 0,
  "deVariables" -> {ss11, P0},
  "deDimensions" -> {{19, 19}, {19, 19}},
  "scalingCertificateScope" -> "symbolic",
  "sourceIntegralIdentities" -> {53, 53},
  "sourceActiveBasisIdentities" -> {19, 19},
  "referenceMatrixEqualCounts" -> <|
    "P0" -> {361, 361}, "ip0" -> {361, 361}, "ks" -> {361, 361}
  |>,
  "intermediateArtifactsRetainedInRepository" -> False
|>
