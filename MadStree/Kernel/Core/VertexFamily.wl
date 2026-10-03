(* ::Package:: *)

(***
File: VertexFamily.wl
Purpose: Initializes single-vertex function families from reference-code ki/nui data or explicit building-block associations.
Scope: This entry does not require a topology; internally it only reuses the certified sector/matrix producers and does not generate spurious internal lines.
***)

(* ::Chapter:: *)
(* Input normalization *)

Options[MSInitVertexFamily] = {NuConvention -> "Positive"};


(* 紧凑输入与显式 block 输入是两种现行数据模型；每种模型只接受自己的字段集合。 *)
$msVertexFamilyCompactKeys = {
  "ki", "nui", "hankelBranches", "exponentialBlocks",
  "vertexType", "normalization"
};
$msVertexFamilyExplicitKeys = {
  "externalLegEnergy", "timePower", "hBlocks", "exponentialBlocks",
  "vertexType", "normalization"
};
$msVertexFamilyHBlockKeys = {"id", "momentum", "nu", "hankelBranch"};
$msVertexFamilyExponentialBlockKeys = {"id", "momentum", "exponentType"};


msVertexFamilyInputSchemaIssues[spec_Association] := Module[
  {compactQ, allowedKeys, fields},
  compactQ = KeyExistsQ[spec, "ki"] || KeyExistsQ[spec, "nui"];
  allowedKeys = If[compactQ, $msVertexFamilyCompactKeys, $msVertexFamilyExplicitKeys];
  fields = Complement[Keys[spec], allowedKeys];
  If[fields === {}, {}, {<|
    "code" -> "unknownVertexFamilyFields",
    "inputMode" -> If[compactQ, "compact", "explicit"],
    "fields" -> fields
  |>}]
];


msVertexFamilyBlockSchemaIssues[blocks_List, allowedKeys_List, kind_String] :=
  MapIndexed[
    Function[{block, index},
      With[{fields = Complement[Keys[block], allowedKeys]},
        If[fields === {}, Nothing, <|
          "code" -> "unknownVertexFamilyBlockFields",
          "blockKind" -> kind,
          "position" -> First[index],
          "fields" -> fields
        |>]
      ]
    ],
    blocks
  ];

msVertexFamilyCompactInput[spec_Association] := Module[
  {ki = Lookup[spec, "ki", Missing["ki"]],
   nui = Lookup[spec, "nui", Missing["nui"]], branches, blockCount},
  If[! ListQ[ki] || ! ListQ[nui] || Length[ki] < 1 || Length[ki] =!= Length[nui],
    Return[msFailure["VertexFamilyKiNuiLengths", <|"ki" -> ki, "nui" -> nui|>]]
  ];
  blockCount = Length[ki] - 1;
  branches = Lookup[spec, "hankelBranches", ConstantArray[1, blockCount]];
  If[! ListQ[branches] || Length[branches] =!= blockCount ||
     ! And @@ (MemberQ[{1, 2}, #] & /@ branches),
    Return[msFailure[
      "VertexFamilyHankelBranches",
      <|"expectedLength" -> blockCount, "actual" -> branches|>
    ]]
  ];
  <|
    "externalLegEnergy" -> First[ki],
    "timePower" -> First[nui],
    "hBlocks" -> MapThread[
      <|"momentum" -> #1, "nu" -> #2, "hankelBranch" -> #3|> &,
      {Rest[ki], Rest[nui], branches}
    ],
    "exponentialBlocks" -> Lookup[spec, "exponentialBlocks", {}],
    "vertexType" -> Lookup[spec, "vertexType", "+"],
    "normalization" -> Lookup[spec, "normalization", 1]
  |>
];

msVertexFamilyAssociationInput[spec_Association] := <|
  "externalLegEnergy" -> Lookup[spec, "externalLegEnergy", Missing["externalLegEnergy"]],
  "timePower" -> Lookup[spec, "timePower", Missing["timePower"]],
  "hBlocks" -> Lookup[spec, "hBlocks", {}],
  "exponentialBlocks" -> Lookup[spec, "exponentialBlocks", {}],
  "vertexType" -> Lookup[spec, "vertexType", "+"],
  "normalization" -> Lookup[spec, "normalization", 1]
|>;

msNormalizeVertexFamilyInput[spec_Association] := If[
  KeyExistsQ[spec, "ki"] || KeyExistsQ[spec, "nui"],
  msVertexFamilyCompactInput[spec],
  msVertexFamilyAssociationInput[spec]
];

msNormalizeVertexFamilyHBlock[block_Association, position_Integer] := <|
  "id" -> Lookup[block, "id", Symbol["h" <> ToString[position]]],
  "momentum" -> Lookup[block, "momentum", Missing["momentum", position]],
  "nu" -> Lookup[block, "nu", Missing["nu", position]],
  "hankelBranch" -> Lookup[block, "hankelBranch", 1]
|>;


(* exponentType 是额外指数块进入 effective external-leg energy 的字面代数号，不是顶点轮廓支。 *)
msExponentTypeSign["+"] := 1;
msExponentTypeSign["-"] := -1;
msExponentTypeSign[other_] := Missing["InvalidExponentType", other];


msNormalizeVertexFamilyExponentialBlock[block_Association, position_Integer] := <|
  "id" -> Lookup[block, "id", Symbol["exp" <> ToString[position]]],
  "momentum" -> Lookup[block, "momentum", Missing["momentum", position]],
  "exponentType" -> Lookup[block, "exponentType", "+"],
  "exponentSign" -> msExponentTypeSign[Lookup[block, "exponentType", "+"]]
|>;


(* ::Chapter:: *)
(* Dedicated public initialization *)

MSInitVertexFamily[spec_Association, OptionsPattern[]] := Module[
  {normalized, rawHBlocks, rawExponentialBlocks, hBlocks, exponentialBlocks, issues = {},
   effectiveExternalLegEnergy, treeSpec, context, vertexId = vertexFamilyRoot, inputIssues,
   nuConvention = OptionValue[NuConvention]},
  inputIssues = msVertexFamilyInputSchemaIssues[spec];
  If[inputIssues =!= {},
    Message[MSInitVertexFamily::badinput, inputIssues];
    Return[msFailure["InvalidVertexFamilyInputFields", <|"issues" -> inputIssues|>]]
  ];
  normalized = msNormalizeVertexFamilyInput[spec];
  If[Head[normalized] === Failure, Return[normalized]];
  rawHBlocks = normalized["hBlocks"];
  rawExponentialBlocks = normalized["exponentialBlocks"];
  If[! ListQ[rawHBlocks] || ! And @@ (AssociationQ /@ rawHBlocks),
    AppendTo[issues, "hBlocks must be a list of associations"]
  ];
  If[! ListQ[rawExponentialBlocks] || ! And @@ (AssociationQ /@ rawExponentialBlocks),
    AppendTo[issues, "exponentialBlocks must be a list of associations"]
  ];
  If[issues === {},
    issues = Join[
      msVertexFamilyBlockSchemaIssues[
        rawHBlocks, $msVertexFamilyHBlockKeys, "h"
      ],
      msVertexFamilyBlockSchemaIssues[
        rawExponentialBlocks, $msVertexFamilyExponentialBlockKeys, "exponential"
      ]
    ]
  ];
  If[issues =!= {},
    Message[MSInitVertexFamily::badinput, issues];
    Return[msFailure["InvalidVertexFamilyInput", <|"issues" -> issues|>]]
  ];
  hBlocks = MapIndexed[msNormalizeVertexFamilyHBlock[#1, First[#2]] &, rawHBlocks];
  exponentialBlocks = MapIndexed[
    msNormalizeVertexFamilyExponentialBlock[#1, First[#2]] &,
    rawExponentialBlocks
  ];
  If[Head[normalized["externalLegEnergy"]] === Missing,
    AppendTo[issues, "missing externalLegEnergy or ki[[1]]"]
  ];
  If[Head[normalized["timePower"]] === Missing, AppendTo[issues, "missing timePower or nui[[1]]"]];
  If[AnyTrue[
      hBlocks,
      Head[#["momentum"]] === Missing || Head[#["nu"]] === Missing ||
        ! MemberQ[{1, 2}, #["hankelBranch"]] &
    ],
    AppendTo[issues, "every h block requires momentum, nu, and Hankel branch 1 or 2"]
  ];
  If[AnyTrue[
      exponentialBlocks,
      Head[#["momentum"]] === Missing || ! MemberQ[{"+", "-"}, #["exponentType"]] &
    ],
    AppendTo[issues, "every exponential block requires momentum and exponentType \"+\" or \"-\""]
  ];
  If[! MemberQ[{"+", "-"}, normalized["vertexType"]],
    AppendTo[issues, "vertexType must be \"+\" or \"-\""]
  ];
  If[issues =!= {},
    Message[MSInitVertexFamily::badinput, issues];
    Return[msFailure["InvalidVertexFamilyInput", <|"issues" -> issues|>]]
  ];
  effectiveExternalLegEnergy = Simplify[
    normalized["externalLegEnergy"] + Total[(#["exponentSign"] #["momentum"]) & /@ exponentialBlocks]
  ];
  treeSpec = <|
    "vertices" -> {
      <|
        "id" -> vertexId,
        "externalLegEnergy" -> effectiveExternalLegEnergy,
        "timePower" -> normalized["timePower"],
        "vertexType" -> normalized["vertexType"]
      |>
    },
    "lines" -> Map[
      <|
        "type" -> "massive",
        "endpoints" -> {vertexId},
        "momentum" -> #["momentum"],
        "nu" -> #["nu"],
        "hankelBranches" -> {#["hankelBranch"]}
      |> &,
      hBlocks
    ],
    "normalization" -> normalized["normalization"]
  |>;
  context = MSInitTree[treeSpec, NuConvention -> nuConvention];
  If[Head[context] === Failure, Return[context]];
  Join[
    context,
    <|
      "contextKind" -> "vertexFamily",
      "vertexFamily" -> <|
        "inputConvention" -> "ki={k0,k1,...}; nui={timePower,nu1,...}",
        "baseExternalLegEnergy" -> normalized["externalLegEnergy"],
        "effectiveExternalLegEnergy" -> effectiveExternalLegEnergy,
        "timePower" -> normalized["timePower"],
        "phaseBlock" -> <|
          "kind" -> "phaseExponent", "dimension" -> 1,
          "momentum" -> normalized["externalLegEnergy"], "exponentType" -> "+", "exponentSign" -> 1
        |>,
        "exponentialBlocks" -> Map[
          Join[#, <|"kind" -> "phaseExponent", "dimension" -> 1|>] &,
          exponentialBlocks
        ],
        "hBlocks" -> Map[
          Join[#, <|"kind" -> "massiveEndpoint", "dimension" -> 2|>] &,
          hBlocks
        ]
      |>
    |>
  ]
];

MSInitVertexFamily[other_, ___] := (
  Message[MSInitVertexFamily::badinput, HoldForm[other]];
  msFailure["MalformedVertexFamilySpec", <|"input" -> HoldForm[other]|>]
);
