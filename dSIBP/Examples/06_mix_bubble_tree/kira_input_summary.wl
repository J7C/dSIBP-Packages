(* ::Package:: *)
(***
File: kira_input_summary.wl
Purpose: Readable record of the formal Kira input built by the mixed bubble+tree example.
Boundary: The full input is generated in a workspace outside this repository and handed to an
external Kira. This file holds no equation tables and no run paths, and it is not package input.
***)


(* ::Chapter:: *)
(*Formal input record*)

<|
  "schema" -> "dsibp_example_kira_input_summary_v1",
  "case" -> "mixBubbleTree",
  "packageVersion" -> "1.0",
  "executionBoundary" -> <|
    "wolframInputDirectory" -> "external validation workspace",
    "kiraRuntime" -> "WSL",
    "kiraRunInsideExampleDirectory" -> False,
    "intermediateArtifactsRetainedInRepository" -> False
  |>,
  "branch" -> "+++",
  "parityConstraints" -> {
    b[1] + n[1, 1] + n[1, 2] -> 0,
    b[2] + n[2, 1] + n[2, 2] -> 0
  },
  "fixedTreeBridgeLine" -> 3,
  "fixedTreeBridgeExcludedFromParity" -> True,
  "formalEnvelope" -> {
    {a[v1], -1, 7}, {a[v2], -1, 7}, {a[v3], -1, 6},
    {b[1], -2, 7}, {b[2], -2, 7}
  },
  "exactPointRules" -> {
    dim -> 37/11, nu1 -> 7/13,
    alpha1 -> 17/19, alpha2 -> 23/29, alpha3 -> 31/37,
    loopScale -> 43/17, legScale1 -> 47/19, legScale2 -> 53/23,
    E1 -> 29/13, E2 -> 31/17, E3 -> 37/19
  },
  "derivativeVariables" -> {loopScale, legScale1, legScale2, E1, E2, E3},
  "seedTemplateCount" -> 178,
  "canonicalEquationCount" -> 818217,
  "formalEquationCount" -> 818297,
  "integralCount" -> 310769,
  "activeMasterCount" -> 81,
  "derivativeTargetCount" -> 595,
  "formalTargetCount" -> 676
|>
