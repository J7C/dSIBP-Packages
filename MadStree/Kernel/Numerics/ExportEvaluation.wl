(* ::Package:: *)

(***
File: ExportEvaluation.wl
Purpose: Exports multi-point numerical evaluation data (point coordinates and master values) to CSV and JSON
         for downstream plotting tools.
Scope: Consumes only point/result records from MSEvaluatePath; it does not
       rerun transports or guess physical conventions.
***)

(* ::Chapter:: *)
(* Public export entry *)

Options[MSExportEvaluationData] = {
  MSOutputDirectory -> Automatic,
  ExportFormats -> Automatic,
  SignificantDigits -> 16
};


(* ::Chapter:: *)
(* Path and record helpers *)

msExportResolveDirectory[Automatic] := FileNameJoin[{
  msRuntimeDirectory[],
  "results",
  "madstree_evaluation",
  "run-" <> CreateUUID[]
}];

msExportResolveDirectory[path_String] := ExpandFileName[
  If[msAbsolutePathQ[path],
    path,
    FileNameJoin[{msRuntimeDirectory[], path}]
  ]
];

msExportResolveDirectory[value_] := msFailure[
  "OutputDirectoryRequired",
  <|"expected" -> "Automatic or a path string", "actual" -> HoldForm[value]|>
];


msExportComplexRecord[value_?NumericQ, digits_Integer] := <|
  "real" -> N[Re[value], digits],
  "imag" -> N[Im[value], digits]
|>;

msExportComplexRecord[value_String, digits_Integer] := <|
  "real" -> value,
  "imag" -> ""
 |>;

msExportComplexRecord[value_, digits_Integer] := <|"real" -> "NonNumeric", "imag" -> "NonNumeric"|>;

msExportRelativeDifference[result_Association] := Replace[
  Lookup[
    Lookup[result, "flintNDE", <||>],
    "relativeDifferenceInf",
    Missing["Unavailable"]
  ],
  {
    text_String :> text,
    value_ :> ToString[value, InputForm]
  }
];


(* ::Chapter:: *)
(* Public export *)

MSExportEvaluationData[
  evaluation_Association,
  opts : OptionsPattern[]
] := Module[
  {pointRecords, points, results, outputDirectory, formats, formatOption, digits, variableSymbols,
   masterCount, header, rows, jsonPayload, csvFile, jsonFile, statuses,
    relativeDifferences, metQ, files, writeResults, failedFiles, writtenFiles,
    pathFailure},
  If[Lookup[evaluation, "status", None] =!= "computed",
    Return[msFailure[
      "ComputedEvaluationRequired",
      <|"status" -> Lookup[evaluation, "status", Missing["Absent"]]|>
    ]]
  ];
  (* 只导出已保存普通点；临时点不冒充用户数值结果。 *)
  pointRecords = Select[
    Lookup[evaluation, "pointResults", {}],
    AssociationQ[#] && Lookup[#, "status", None] === "saved" &&
      ListQ[Lookup[#, "value", None]] &
  ];
  If[pointRecords === {},
    Return[msFailure["SavedPointResultsRequired", <||>]]
  ];
  points = Lookup[pointRecords, "coordinate"];
  results = Map[
    <|
      "values" -> #["value"],
      "status" -> #["status"],
      "flintNDE" -> Lookup[#, "flintNDE", <||>]
    |> &,
    pointRecords
  ];

  outputDirectory = msExportResolveDirectory[OptionValue[MSOutputDirectory]];
  If[Head[outputDirectory] === Failure, Return[outputDirectory]];

  digits = OptionValue[SignificantDigits];
  If[! IntegerQ[digits] || digits < 2,
    Return[msFailure["InvalidSignificantDigits", <|"value" -> digits|>]]
  ];
  (* Replace 的 RuleDelayed 显式绑定列表值，避免 Switch 的命名 pattern 留成未求值符号。 *)
  formatOption = OptionValue[ExportFormats];
  formats = Replace[formatOption, {
    Automatic -> {"CSV", "JSON"},
    "CSV" -> {"CSV"},
    "JSON" -> {"JSON"},
    value_List :> DeleteDuplicates@Select[value, MemberQ[{"CSV", "JSON"}, #] &],
    _ -> {"CSV"}
  }];
  If[formats === {},
    Return[msFailure["InvalidExportFormats", <|"value" -> OptionValue[ExportFormats]|>]]
  ];

  variableSymbols = DeleteDuplicates[Flatten[
    Cases[#, Rule[left_Symbol, _] :> left, Infinity] & /@ points
  ]];
  masterCount = Length[Lookup[First[results], "values"]];
  statuses = Lookup[results, "status"];
  relativeDifferences = msExportRelativeDifference /@ results;
  metQ = TrueQ /@ Lookup[Lookup[results, "flintNDE", <||>], "targetRelativeErrorMet"];

  rows = MapThread[
    Function[{point, result},
      Join[
        Flatten[
          Map[
            Function[symbol,
              With[{record = msExportComplexRecord[symbol /. point, digits]},
                {record["real"], record["imag"]}
              ]
            ],
            variableSymbols
          ]
        ],
        Flatten[
          Map[
            Function[value,
              With[{record = msExportComplexRecord[value, digits]},
                {record["real"], record["imag"]}
              ]
            ],
            result["values"]
          ]
        ],
        {ToString[Lookup[result, "status", Missing["Unavailable"]]],
         msExportRelativeDifference[result],
         ToString[TrueQ[Lookup[Lookup[result, "flintNDE", <||>], "targetRelativeErrorMet"]]]}
      ]
    ],
    {points, results}
  ];

  header = Join[
    Flatten[Map[Function[symbol, {SymbolName[symbol] <> "_re", SymbolName[symbol] <> "_im"}], variableSymbols]],
    Flatten[Table[{("M" <> ToString[index] <> "_re"), ("M" <> ToString[index] <> "_im")}, {index, masterCount}]],
    {"status", "relativeDifferenceInf", "targetRelativeErrorMet"}
  ];

  csvFile = FileNameJoin[{outputDirectory, "evaluation_data.csv"}];
  jsonFile = FileNameJoin[{outputDirectory, "evaluation_data.json"}];
  files = Join[
    If[MemberQ[formats, "CSV"], {csvFile}, {}],
    If[MemberQ[formats, "JSON"], {jsonFile}, {}]
  ];
  pathFailure = msRuntimePathLengthFailure[Prepend[files, outputDirectory]];
  If[Head[pathFailure] === Failure, Return[pathFailure]];
  If[msEnsureDirectory[outputDirectory] === $Failed,
    Return[msFailure["OutputDirectoryCreationFailed", <|"path" -> outputDirectory|>]]
  ];

  jsonPayload = <|
    "schema" -> "madstree_evaluation_data_v1",
    "pointCount" -> Length[points],
    "masterCount" -> masterCount,
    "digits" -> digits,
    "variables" -> (SymbolName /@ variableSymbols),
    "points" -> Map[
      Function[point,
        AssociationThread[
          SymbolName /@ variableSymbols,
          Map[msExportComplexRecord[# /. point, digits] &, variableSymbols]
        ]
      ],
      points
    ],
    "values" -> Map[
      Function[result, Map[msExportComplexRecord[#, digits] &, result["values"]]],
      results
    ],
    "statuses" -> statuses,
    "relativeDifferenceInf" -> relativeDifferences,
    "targetRelativeErrorMet" -> metQ
  |>;
  writeResults = Association@DeleteCases[{
    If[MemberQ[formats, "CSV"],
      csvFile -> Quiet@Check[Export[csvFile, Prepend[rows, header], "CSV"], $Failed],
      Nothing
    ],
    If[MemberQ[formats, "JSON"],
      jsonFile -> Quiet@Check[Export[jsonFile, jsonPayload, "JSON"], $Failed],
      Nothing
    ]
  }, Nothing];
  failedFiles = Select[
    Keys[writeResults],
    Lookup[writeResults, #] === $Failed || ! FileExistsQ[#] &
  ];
  writtenFiles = Select[Keys[writeResults], FileExistsQ];
  If[failedFiles =!= {},
    Return[msFailure["EvaluationExportFailed", <|
      "outputDirectory" -> outputDirectory,
      "writtenFiles" -> writtenFiles,
      "failedFiles" -> failedFiles,
      "pathLengths" -> AssociationMap[StringLength, Keys[writeResults]]
    |>]]
  ];

  <|
    "status" -> "written",
    "outputDirectory" -> outputDirectory,
    "files" -> writtenFiles,
    "formats" -> formats,
    "pointCount" -> Length[points],
    "masterCount" -> masterCount
  |>
];

MSExportEvaluationData[___] := msFailure[
  "ComputedEvaluationRequired",
  <|"function" -> "MSExportEvaluationData"|>
];
