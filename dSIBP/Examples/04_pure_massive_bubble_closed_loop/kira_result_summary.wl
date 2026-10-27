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
  "evidenceRelease" -> "025.0",
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
    "wallTimeSeconds" -> 5,
    "exitStatus" -> 0
  |>,
  "masterIDs" -> Range[19],
  "masterCount" -> 19,
  "backendMasterIDs" -> Join[Range[19], {888, 925, 927, 1145, 1147}],
  "globalMasterCount" -> 24,
  "targetCount" -> 300,
  "selectedEquationCount" -> 795,
  "unreducedCount" -> 0,
  "deVariables" -> {ss11, P0},
  "deDimensions" -> {{19, 19}, {19, 19}},
  "postKiraChecks" -> {23, 23},
  "scalingCertificateScope" -> "symbolic",
  "sourceIntegralIdentities" -> {53, 53},
  "sourceActiveBasisIdentities" -> {19, 19},
  "sourceIdentityEvidenceRelease" -> "024.0",
  "referenceMatrixEqualCounts" -> <|
    "P0" -> {361, 361}, "ip0" -> {361, 361}, "ks" -> {361, 361}
  |>,
  "formalValidationManifestSHA256" ->
    "8BEF56D1361F01086468A1EAF7B26B733BDFA73E3E5C4C92B8281D87021542F8",
  "formalValidationReductionSHA256" ->
    "57943724084180761C66FC756C0E38714624784B7877D97A2B3051188DF491FF",
  "intermediateArtifactsRetainedInRepository" -> False
|>
