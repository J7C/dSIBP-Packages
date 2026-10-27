(* ::Package:: *)
(***
文件：kira_result_summary.wl
用途：保存 mixed bubble+tree 的 fresh WSL Kira、DE 与 scaling 最终摘要。
边界：本文件不保留 reduction table、日志、database 或 cache；正式证据见独立检验报告。
***)


(* ::Chapter:: *)
(*正式结果摘要*)

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
  "formalValidationManifestSHA256" ->
    "B3C9ED791EB2C0247A114FA2AC29B526A860DB77FD9B1C6921E8AA412A7A7106",
  "formalValidationReductionSHA256" ->
    "7D66B032E78EB724B7EA643092CEAA1DCACA81AA3A6C2968D362E8127FD9C5BD",
  "intermediateArtifactsRetainedInRepository" -> False
|>
