(* ::Package:: *)
(***
文件：kira_result_summary.wl
用途：保存 pure massive bubble 的 fresh WSL Kira、DE 与 scaling 最终摘要。
边界：本文件不保留 reduction table、日志、database 或 cache；正式证据见独立检验报告。
***)


(* ::Chapter:: *)
(*正式结果摘要*)

<|
  "schema" -> "dsibp_example_kira_result_summary_v1",
  "status" -> "passed",
  "case" -> "pureMassiveBubble",
  "packageVersion" -> "1.0",
  "evidenceRelease" -> "030",
  "generationOptions" -> {
    SeedEnvelopeCut -> Automatic, TotalPowerLimit -> None
  },
  "targetEnvelope" -> {
    {a[v1], 0, 3}, {a[v2], 0, 3}, {b[1], -3, 5}, {b[2], -3, 5}
  },
  "kira" -> <|
    "version" -> "2.3 (Git: 2.3)",
    "runtime" -> "WSL",
    "parallelConfig" -> "w10*1",
    "wallTimeSeconds" -> 6,
    "exitStatus" -> 0
  |>,
  "masterIDs" -> Range[19],
  "masterCount" -> 19,
  "backendMasterIDs" -> Join[Range[19], {151, 152, 154, 155, 156, 166, 167, 168, 169, 173, 178,
    182, 183, 189, 190, 234, 235, 236, 239, 240, 244, 247, 248, 256, 307, 308, 317, 318,
    319, 323, 324, 330, 331, 340, 342, 345, 346, 347, 348, 349, 352, 353, 354, 355, 356,
    358, 360, 366, 367, 372, 377, 378, 383, 391, 392, 393, 394, 700, 701, 702, 703, 704,
    705, 706, 707, 714, 716, 736, 738, 920, 921, 922, 923, 928, 929, 931, 932, 942, 943,
    944, 945}],
  "globalMasterCount" -> 100,
  "targetCount" -> 300,
  "selectedEquationCount" -> 795,
  "unreducedCount" -> 0,
  "deVariables" -> {ss11, P0},
  "deDimensions" -> {{19, 19}, {19, 19}},
  "postKiraChecks" -> {23, 23},
  "scalingCertificateScope" -> "symbolic",
  "sourceIntegralIdentities" -> {53, 53},
  "sourceActiveBasisIdentities" -> {19, 19},
  "referenceMatrixEqualCounts" -> <|
    "P0" -> {361, 361}, "ip0" -> {361, 361}, "ks" -> {361, 361}
  |>,
  "formalValidationManifestSHA256" ->
    "405127690D0FCEC621E951543697CF33745B76A4C1440A674265119F2697284D",
  "formalValidationReductionSHA256" ->
    "E837D5AE39BAE249400BF80515D018161166525AA89AE2C8F45414FDEBAFD3E0",
  "intermediateArtifactsRetainedInRepository" -> False
|>
