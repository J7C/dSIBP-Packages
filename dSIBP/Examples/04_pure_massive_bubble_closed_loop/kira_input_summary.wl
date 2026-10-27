(* ::Package:: *)
(***
文件：kira_input_summary.wl
用途：保存 pure massive bubble 正式 Kira 输入的轻量、可读摘要。
边界：完整输入在仓库外工作区生成并交给 WSL Kira；本文件不含方程表或运行路径，也不作为程序包输入。
***)


(* ::Chapter:: *)
(*正式输入摘要*)

<|
  "schema" -> "dsibp_example_kira_input_summary_v1",
  "case" -> "pureMassiveBubble",
  "packageVersion" -> "1.0",
  "evidenceRelease" -> "025.0",
  "generationOptions" -> {
    SeedEnvelopeCut -> Automatic, TotalPowerLimit -> None
  },
  "targetEnvelope" -> {
    {a[v1], 0, 3}, {a[v2], 0, 3}, {b[1], -3, 5}, {b[2], -3, 5}
  },
  "packageSHA256" ->
    "6B043F0B1C0BEB6687D03A5D993896EA894481A617D8251D014F325F6733EA82",
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
  "canonicalEquationCount" -> 2442,
  "formalEquationCount" -> 2460,
  "integralCount" -> 1317,
  "activeMasterCount" -> 19,
  "formalTargetCount" -> 300,
  "retainedInputFiles" -> {
    "main.wl", "family_conventions.wl", "reference_user_mi_basis.wl",
    "dlog_basis.wl"
  }
|>
