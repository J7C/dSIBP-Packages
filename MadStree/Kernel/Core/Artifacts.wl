(* ::Package:: *)

(***
File: Artifacts.wl
Purpose: Aggregates and writes all-sector master integrals, recurrence metadata and the dlog DE of a MadStree context.
Scope: This module only serializes the results of the existing formula producers; it does not generate IBP equations or derive differential operators.
***)

(* ::Chapter:: *)
(* All-sector formula data *)

(* timePowerRules are substituted uniformly only in the final formula data; the sector/master identities frozen in the context remain unchanged. *)
msFormulaTimePowerSymbols[context_?MSContextQ] := DeleteDuplicates@Cases[
  Lookup[context["vertices"], "timePower", {}],
  symbol_Symbol /; Context[symbol] =!= "System`",
  Infinity
];


msFormulaTimePowerRules[Automatic] := {};
msFormulaTimePowerRules[rules_List] /; And @@ (MatchQ[#, _Rule | _RuleDelayed] & /@ rules) := rules;
msFormulaTimePowerRules[value_] := msFailure[
  "TimePowerRulesRequired",
  <|"expected" -> "Automatic or a list of rules", "actual" -> HoldForm[value]|>
];


Options[MSFormulaData] = {TimePowerRules -> Automatic};


MSFormulaData[context_?MSContextQ, OptionsPattern[]] := Module[
  {rules, recurrenceMetadata, de, rawData, substitutedData},
  rules = msFormulaTimePowerRules[OptionValue[TimePowerRules]];
  If[Head[rules] === Failure, Return[rules]];
  recurrenceMetadata = Association@Table[
    sectorKey -> MSFormulaMatrices[context, sectorKey],
    {sectorKey, context["sectorOrder"]}
  ];
  If[! FreeQ[Values[recurrenceMetadata], _Failure],
    Return[msFailure["FormulaMatrixGenerationFailed", <|"sectors" -> Keys[recurrenceMetadata]|>]]
  ];
  de = MSDLogDE[context];
  If[Lookup[de, "dlogStatus", None] =!= "certifiedByFormulaChecks",
    Return[msFailure[
      "CertifiedDLogRequired",
      <|"dlogStatus" -> Lookup[de, "dlogStatus", Missing["Absent"]]|>
    ]]
  ];
  rawData = <|
    "status" -> "generated",
    "version" -> $MadStreeVersion,
    "caseName" -> Lookup[context, "caseName", "MadStreeContext"],
    "nuConvention" -> context["convention", "nuConvention"],
    "sectorOrder" -> context["sectorOrder"],
    "sectorKeySchema" -> context["sectorKeySchema"],
    "masters" -> context["masters"],
    "masterDigest" -> context["masterDigest"],
    "timePowerSymbols" -> msFormulaTimePowerSymbols[context],
    "timePowerRules" -> rules,
    "recurrenceMetadata" -> recurrenceMetadata,
    "dlogDE" -> de
  |>;
  substitutedData = rawData /. rules;
  Join[
    substitutedData,
    <|
      (* dlog_de.wl 始终保存未代数值的解析 DE；数值参数只允许在求值器的局部副本中使用。 *)
      "dlogDE" -> de,
      "timePowerSymbols" -> rawData["timePowerSymbols"],
      "timePowerRules" -> rules,
      "artifactDigest" -> IntegerString[
        Hash[Join[substitutedData, <|"dlogDE" -> de|>], "SHA256"], 16, 64
      ]
    |>
  ]
];


MSFormulaData[___] := msFailure[
  "InitializedContextRequired",
  <|"function" -> "MSFormulaData"|>
];


(* ::Chapter:: *)
(* Calling directory and file writing *)

$MSFormulaArtifactRegistry = <||>;


msFormulaRuntimeDirectory[] := msRuntimeDirectory[];


msFormulaAbsolutePathQ[path_String] := msAbsolutePathQ[path];


msResolveFormulaOutputDirectory[Automatic] := FileNameJoin[{
  msFormulaRuntimeDirectory[],
  "results",
  "madstree_formula",
  "run-" <> CreateUUID[]
}];


msResolveFormulaOutputDirectory[path_String] := ExpandFileName[
  If[msFormulaAbsolutePathQ[path], path, FileNameJoin[{msFormulaRuntimeDirectory[], path}]]
];


msResolveFormulaOutputDirectory[value_] := msFailure[
  "OutputDirectoryRequired",
  <|"expected" -> "Automatic or a path string", "actual" -> HoldForm[value]|>
];


(* 同一 context 在同一调用根只保留一份有效解析资产引用；用户删除文件后会在下一次数值调用前重建。 *)
msFormulaArtifactRegistryKey[context_?MSContextQ] := IntegerString[
  Hash[{ExpandFileName[msFormulaRuntimeDirectory[]], context}, "SHA256"],
  16,
  64
];


msFormulaArtifactFilesExistQ[data_] := AssociationQ[data] &&
  AssociationQ[Lookup[data, "files", None]] &&
  AllTrue[Values[data["files"]], FileExistsQ];


msFormulaArtifactReference[data_Association] := KeyTake[
  data,
  {"outputDirectory", "files", "artifactDigest"}
];


Options[MSWriteFormulaArtifacts] = {
  TimePowerRules -> Automatic,
  MSOutputDirectory -> Automatic
};


MSWriteFormulaArtifacts[context_?MSContextQ, OptionsPattern[]] := Module[
  {data, outputDirectory, paths, manifestPath, manifest, result, registryKey,
   pathFailure, writeResult, failedFiles},
  data = MSFormulaData[
    context,
    TimePowerRules -> OptionValue[TimePowerRules]
  ];
  If[Head[data] === Failure, Return[data]];
  outputDirectory = msResolveFormulaOutputDirectory[OptionValue[MSOutputDirectory]];
  If[Head[outputDirectory] === Failure, Return[outputDirectory]];
  paths = <|
    "masters" -> FileNameJoin[{outputDirectory, "masters.wl"}],
    "recurrenceMetadata" -> FileNameJoin[{outputDirectory, "recurrence_metadata.wl"}],
    "dlogDE" -> FileNameJoin[{outputDirectory, "dlog_de.wl"}]
  |>;
  manifestPath = FileNameJoin[{outputDirectory, "manifest.wl"}];
  pathFailure = msRuntimePathLengthFailure[
    Join[{outputDirectory, manifestPath}, Values[paths]]
  ];
  If[Head[pathFailure] === Failure, Return[pathFailure]];
  If[msEnsureDirectory[outputDirectory] === $Failed,
    Return[msFailure["OutputDirectoryCreationFailed", <|"path" -> outputDirectory|>]]
  ];
  manifest = KeyTake[
    data,
    {
      "status", "version", "caseName", "nuConvention", "sectorOrder",
      "sectorKeySchema", "masterDigest", "timePowerSymbols", "timePowerRules", "artifactDigest"
    }
  ];
  manifest = Join[manifest, <|
    "outputDirectory" -> outputDirectory,
    "files" -> Join[paths, <|"manifest" -> manifestPath|>],
    "allSectorQ" -> True,
    "dlogStatus" -> data["dlogDE", "dlogStatus"]
  |>];
  writeResult = Quiet@Check[
    Put[data["masters"], paths["masters"]];
    Put[data["recurrenceMetadata"], paths["recurrenceMetadata"]];
    Put[data["dlogDE"], paths["dlogDE"]];
    Put[manifest, manifestPath];
    "written",
    $Failed
  ];
  failedFiles = Select[Join[Values[paths], {manifestPath}], ! FileExistsQ[#] &];
  If[writeResult === $Failed || failedFiles =!= {},
    Return[msFailure["FormulaArtifactWriteFailed", <|
      "outputDirectory" -> outputDirectory,
      "writtenFiles" -> Select[Join[Values[paths], {manifestPath}], FileExistsQ],
      "failedFiles" -> failedFiles
    |>]]
  ];
  result = Join[
    data,
    <|
      "outputDirectory" -> outputDirectory,
      "files" -> Join[paths, <|"manifest" -> manifestPath|>],
      "manifest" -> manifest
    |>
  ];
  registryKey = msFormulaArtifactRegistryKey[context];
  AssociateTo[$MSFormulaArtifactRegistry, registryKey -> result];
  result
];


MSWriteFormulaArtifacts[___] := msFailure[
  "InitializedContextRequired",
  <|"function" -> "MSWriteFormulaArtifacts"|>
];


(* 数值入口使用本 helper 强制先取得完整解析公式；写出失败直接向上传递 Failure。 *)
msEnsureFormulaArtifacts[context_?MSContextQ] := Module[
  {registryKey, cached, written},
  registryKey = msFormulaArtifactRegistryKey[context];
  cached = Lookup[$MSFormulaArtifactRegistry, registryKey, Missing["NotWritten"]];
  If[msFormulaArtifactFilesExistQ[cached], Return[cached]];
  written = MSWriteFormulaArtifacts[context];
  If[Head[written] === Failure, Return[written]];
  written
];
