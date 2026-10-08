(* ::Package:: *)
(***
File: kira_input_summary.wl
Purpose: Readable record of the formal Kira input built by the pure massive bubble example.
Boundary: The full input is generated in a workspace outside this repository and handed to an
external Kira. This file holds no equation tables and no run paths, and it is not package input.
***)


(* ::Chapter:: *)
(*Formal input record*)

<|
  "schema" -> "dsibp_example_kira_input_summary_v1",
  "case" -> "pureMassiveBubble",
  "packageVersion" -> "1.0",
  "executionBoundary" -> <|
    "wolframInputDirectory" -> "external workspace selected by DSIBP_KIRA_WORKSPACE",
    "kiraRuntime" -> "WSL",
    "kiraRunInsideExampleDirectory" -> False,
    "intermediateArtifactsRetainedInRepository" -> False
  |>,
  "topology" -> <|
    "branch" -> "--",
    "massiveCycleLineCount" -> 2,
    "loopCount" -> 1,
    "parity" -> "both cycle propagators are even"
  |>,
  "parameterRules" -> {
    dim -> 37/11, nu -> 7/13, etaNu -> 23/17,
    analyticRegulator -> 0
  },
  "numericStage" -> "postDerivative",
  "seedTemplateCount" -> 88,
  "canonicalEquationCount" -> 14986,
  "formalEquationCount" -> 15004,
  "integralCount" -> 5728,
  "activeMasterCount" -> 19,
  "formalTargetCount" -> 300,
  "retainedInputFiles" -> {
    "main.wl", "family_conventions.wl", "reference_user_mi_basis.wl",
    "dlog_basis.wl"
  }
|>
