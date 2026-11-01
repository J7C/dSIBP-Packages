(* ::Package:: *)
(***
文件：kira_input_summary.wl
用途：保存 mixed bubble+tree 在 dSIBP 029 上 fresh Kira 输入的轻量摘要。
边界：完整输入在仓库外生成并交给 WSL Kira；本文件不含方程表、运行路径或 reduction table。
***)


(* ::Chapter:: *)
(*正式输入摘要*)

<|
  "schema" -> "dsibp_example_kira_input_summary_v1",
  "case" -> "mixBubbleTree",
  "packageVersion" -> "1.0",
  "evidenceRelease" -> "029",
  "uncutEquivalentOptionsForCurrentPackage" -> {
    SeedEnvelopeCut -> None, TotalPowerLimit -> None
  },
  "executionBoundary" -> <|
    "wolframInputDirectory" -> "external validation workspace",
    "kiraRuntime" -> "WSL",
    "kiraRunInsideExampleDirectory" -> False,
    "intermediateArtifactsRetainedInRepository" -> False
  |>,
  "branch" -> "+++",
  "nuConvention" -> "Positive",
  "parityConstraints" -> {
    b[1] + n[1, 1] + n[1, 2] -> 0,
    b[2] + n[2, 1] + n[2, 2] -> 0
  },
  "fixedTreeBridgeLine" -> 3,
  "fixedTreeBridgeExcludedFromParity" -> True,
  "formalEnvelope" -> {
    {a[v1], 0, 5}, {a[v2], 0, 5}, {a[v3], 0, 0},
    {b[1], -4, 7}, {b[2], -4, 7}
  },
  "targetRule" -> "903 fixed integrals: independent a[v3]=0, nonnegative continuous slots, total absolute power <= 5",
  "exactPointRules" -> {
    dim -> 37/11, nu1 -> 7/13,
    alpha1 -> 17/19, alpha2 -> 23/29, alpha3 -> 31/37,
    loopScale -> 43/17, legScale1 -> 47/19, legScale2 -> 53/23,
    E1 -> 29/13, E2 -> 31/17, E3 -> 37/19
  },
  "symbolicCoefficientVariables" -> {beta1, beta2},
  "seedTemplateCount" -> 167,
  "canonicalEquationCount" -> 68418,
  "integralCount" -> 35412,
  "activeMasterCount" -> 84,
  "formalTargetCount" -> 903,
  "a3TargetExtension" -> <|
    "seedEnvelopeA3" -> {-1, 1},
    "targetA3Values" -> {0, 1},
    "targetCount" -> 1332,
    "masterCount" -> 84
  |>
|>
