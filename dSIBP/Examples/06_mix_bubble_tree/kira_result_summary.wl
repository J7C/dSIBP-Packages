(* ::Package:: *)
(***
文件：kira_result_summary.wl
用途：保存 mixed bubble+tree 在 dSIBP 029 上 fresh WSL Kira 的最终轻量摘要。
边界：不保留 reduction table、日志、database 或 cache；完整台账见版本外 example-test-ledger/029_example6.md。
***)


(* ::Chapter:: *)
(*正式结果摘要*)

<|
  "schema" -> "dsibp_example_kira_result_summary_v1",
  "status" -> "passed",
  "case" -> "mixBubbleTree",
  "packageVersion" -> "1.0",
  "evidenceRelease" -> "029",
  "branch" -> "+++",
  "kira" -> <|
    "version" -> "2.3 (Git: 2.3)",
    "runtime" -> "WSL",
    "parallelConfig" -> "10 jobs",
    "stageWallTimeSeconds" -> 253.923,
    "exitStatus" -> 0
  |>,
  "masterCount" -> 84,
  "targetCount" -> 903,
  "canonicalEquationCount" -> 68418,
  "integralCount" -> 35412,
  "selectedEquationCount" -> 15874,
  "unreducedCount" -> 0,
  "mastersBySector" -> <|
    "top" -> 44, "e3" -> 22, "e2" -> 8,
    "e1_e2" -> 4, "e2_e3" -> 4, "e1_e2_e3" -> 2
  |>,
  "mastersByTotalAbsolutePower" -> <|"0" -> 9, "1" -> 33, "2" -> 42|>,
  "masterListFile" -> "master_integrals.wl",
  "a3SeedExpansionKeepsSameMasters" -> True,
  "a3TargetExtension" -> <|
    "targetCount" -> 1332,
    "reducedTargetCount" -> 1248,
    "masterCount" -> 84,
    "sameMasterIntegralsAsA3ZeroTarget" -> True
  |>,
  "deStatus" -> "notRun",
  "scalingStatus" -> "notRun",
  "reductionSHA256" -> "2D9DB31CF05FA483781181373BDA9DBFDC7121B0E2C2672D2917596EED823871",
  "mastersSHA256" -> "87B39955410DDE0A28A0A805645722A9E00E37D4B701B9378C357FCAAD60224B",
  "intermediateArtifactsRetainedInRepository" -> False
|>
