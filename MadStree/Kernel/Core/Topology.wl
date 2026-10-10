(* ::Package:: *)

(***
文件：Topology.wl
用途：规范化并验证 MadStree 的有序树图或纯时间关联图输入。
范围：用户只在顶点上指定轮廓支；传播子的轮廓类型与符号由端点派生。massless Full 线固定使用
共享二态 quotient，不读取其它表示选择；本模块不处理圈动量、ISP 或标量积闭合。
***)

(* ::Chapter:: *)
(* Input normalization *)

$msSupportedInputLineTypes = {"massive", "massless"};


(* Association 键顺序和额外字段不影响物理输入；这里只检查后续公式必需的字段是否齐全。 *)
$msRequiredTreeSpecKeys = {"vertices", "lines"};
$msRequiredVertexInputKeys = {"id", "timePower", "vertexType"};
$msRequiredLineInputKeys = {"type", "endpoints", "momentum"};


msInputSchemaIssues[spec_Association, vertices_List, lines_List] := Module[
  {missingRootFields, missingVertexFields, missingLineFields, issues = {}},
  missingRootFields = msAssociationKeysMissing[spec, $msRequiredTreeSpecKeys];
  If[missingRootFields =!= {},
    AppendTo[issues, <|"code" -> "missingTreeFields", "fields" -> missingRootFields|>]
  ];
  Do[
    missingVertexFields = msAssociationKeysMissing[vertex, $msRequiredVertexInputKeys];
    If[missingVertexFields =!= {},
      AppendTo[issues, <|
        "code" -> "missingVertexFields", "position" -> position,
        "fields" -> missingVertexFields
      |>]
    ],
    {position, Length[vertices]}, {vertex, {vertices[[position]]}}
  ];
  Do[
    missingLineFields = msAssociationKeysMissing[line, $msRequiredLineInputKeys];
    If[missingLineFields =!= {},
      AppendTo[issues, <|
        "code" -> "missingLineFields", "position" -> position,
        "fields" -> missingLineFields
      |>]
    ],
    {position, Length[lines]}, {line, {lines[[position]]}}
  ];
  issues
];

msAssociationKeysMissing[association_Association, keys_List] := Select[keys, ! KeyExistsQ[association, #] &];

(* 顶点 + 对应 Exp[-I E tau]，顶点 - 对应 Exp[+I E tau]。该符号直接进入
   component energy，因此只允许从 vertexType 在此处派生。 *)
msVertexTypeSign["+"] := -1;
msVertexTypeSign["-"] := 1;
msVertexTypeSign[other_] := Missing["InvalidVertexType", other];

(* SK 轮廓权重与外腿指数相位是不同合同：+/- 轮廓本身仍对应 +1/-1，
   传播子 contact map 只能读取这一层。 *)
msVertexTypeContourSign["+"] := 1;
msVertexTypeContourSign["-"] := -1;
msVertexTypeContourSign[other_] := Missing["InvalidVertexType", other];

msFullContourSign["++"] := 1;
msFullContourSign["--"] := -1;
msFullContourSign[_] := Missing["NotFullContourType"];

msDerivedLineClass[endpointTypes_List] := Switch[
  Length[endpointTypes],
  1, "External",
  2, If[SameQ @@ endpointTypes, "Full", "Cross"],
  _, Missing["InvalidEndpointCount", Length[endpointTypes]]
];

(* Branch numbering follows HankelH[1|2,nu,z]. The branch of a full propagator is determined dynamically by the time ordering. *)
msDefaultHankelBranches["massiveCross", "+-"] := {2, 1};
msDefaultHankelBranches["massiveCross", "-+"] := {1, 2};
msDefaultHankelBranches["massiveExternal", "+"] := {2};
msDefaultHankelBranches["massiveExternal", "-"] := {1};
msDefaultHankelBranches["massiveFull", "++"] := "timeOrdered++";
msDefaultHankelBranches["massiveFull", "--"] := "timeOrdered--";
msDefaultHankelBranches["masslessFull", "++"] := "timeOrdered++";
msDefaultHankelBranches["masslessFull", "--"] := "timeOrdered--";
msDefaultHankelBranches[_, _] := Missing["NotHankelLine"];

msAuxiliaryExternalLegEnergySymbol[position_Integer] := Symbol[
  "MadStree`Private`msAuxiliaryExternalLegEnergy" <> IntegerString[position]
];


msNormalizeVertex[vertex_Association, position_Integer] := Module[
  {inputEnergy, auxiliaryQ, effectiveEnergy},
  inputEnergy = Lookup[vertex, "externalLegEnergy", 0];
  auxiliaryQ = TrueQ[inputEnergy === 0];
  effectiveEnergy = If[
    auxiliaryQ,
    msAuxiliaryExternalLegEnergySymbol[position],
    inputEnergy
  ];
  <|
    "id" -> Lookup[vertex, "id", Missing["VertexId", position]],
    "position" -> position,
    "externalLegEnergy" -> effectiveEnergy,
    "externalLegEnergyTarget" -> inputEnergy,
    "auxiliaryExternalLegEnergyQ" -> auxiliaryQ,
    "timePower" -> Lookup[vertex, "timePower", 0],
    "vertexType" -> Lookup[vertex, "vertexType", Missing["VertexType", position]],
    "contourSign" -> msVertexTypeContourSign[Lookup[vertex, "vertexType", Missing["VertexType", position]]],
    "vertexSign" -> msVertexTypeSign[Lookup[vertex, "vertexType", Missing["VertexType", position]]]
  |>
];


msAuxiliaryExternalLegEnergyRules[context_?MSContextQ] := Cases[
  context["vertices"],
  vertex_Association /; TrueQ[vertex["auxiliaryExternalLegEnergyQ"]] :>
    Rule[vertex["externalLegEnergy"], vertex["externalLegEnergyTarget"]]
];


msAuxiliaryExternalLegEnergySymbols[context_?MSContextQ] :=
  First /@ msAuxiliaryExternalLegEnergyRules[context];

msNormalizeLine[
  line_Association,
  position_Integer,
  nuConvention_String,
  vertices_List
] := Module[
  {inputType, type, endpoints, endpointVertices, endpointTypes, endpointSigns, lineClass,
   momentum, nuMagnitude, hankelOrder, formulaNu, contourType, fullContourSign,
   contactRawPower, pinchNormalization, functionSystem, hankelBranches, vertexById},
  inputType = Lookup[line, "type", Missing["LineType", position]];
  endpoints = Lookup[line, "endpoints", {}];
  vertexById = AssociationThread[Lookup[vertices, "id"] -> vertices];
  endpointVertices = Lookup[vertexById, endpoints, Missing["UnknownEndpoint"]];
  endpointTypes = Lookup[endpointVertices, "vertexType", Missing["UnknownEndpoint"]];
  endpointSigns = msVertexTypeContourSign /@ endpointTypes;
  lineClass = msDerivedLineClass[endpointTypes];
  type = If[
    MemberQ[$msSupportedInputLineTypes, inputType] && StringQ[lineClass],
    inputType <> lineClass,
    Missing["DerivedLineType", {inputType, lineClass}]
  ];
  momentum = Lookup[line, "momentum", Missing["Momentum", position]];
  nuMagnitude = Lookup[line, "nu", If[inputType === "massless", 1/2, Missing["NuMagnitude", position]]];
  hankelOrder = nuMagnitude;
  formulaNu = If[nuConvention === "Positive", nuMagnitude, -nuMagnitude];
  contourType = If[And @@ StringQ /@ endpointTypes, StringJoin[endpointTypes], Missing["ContourType"]];
  fullContourSign = msFullContourSign[contourType];
  contactRawPower = Lookup[
    line,
    "contactRawPower",
    Switch[type, "massiveFull", 2 formulaNu - 1, "masslessFull", 0, _, Missing["NoContact"]]
  ];
  pinchNormalization = Lookup[
    line,
    "pinchNormalization",
    Switch[
      type,
      (* MadStree 两个端点使用共同 Hankel 阶，因此不重复加入 paper endpoint basis 的
         Boltzmann 因子；但 ++/-- 收缩定义的物理 child master 相差 fullContourSign。
         pinch 幂次取 +nu 基准形状 k^(2 formulaNu-1)（等价于论文的 k^(-2 nu_ref-1)）。 *)
      "massiveFull", -(4 I/Pi) fullContourSign momentum^(2 formulaNu - 1),
      "masslessFull", 1,
      _, 1
    ]
  ];
  functionSystem = Which[
    type === "masslessFull", "masslessQuotient",
    StringStartsQ[ToString[type], "massive"], "h",
    True, "exponential"
  ];
  hankelBranches = Lookup[line, "hankelBranches", msDefaultHankelBranches[type, contourType]];
  <|
    "id" -> position,
    "position" -> position,
    "inputType" -> inputType,
    "type" -> type,
    "endpoints" -> endpoints,
    "momentum" -> momentum,
    "nu" -> nuMagnitude,
    "hankelOrder" -> hankelOrder,
    (* hPrefactorPower 是 h=z^(s_nu nu) 的带号 prefactor 幂，外部记录与 Sectors 的小 t 幂表依赖它；
       新语义下 formulaNu=s_nu nu 与它数值恒等，此处仍保留独立字段不变。 *)
    "hPrefactorPower" -> If[nuConvention === "Positive", nuMagnitude, -nuMagnitude],
    "formulaNu" -> formulaNu,
    "nuConvention" -> nuConvention,
    "contourType" -> contourType,
    "fullContourSign" -> fullContourSign,
    "endpointSigns" -> endpointSigns,
    "hankelBranches" -> hankelBranches,
    "functionSystem" -> functionSystem,
    "contactRawPower" -> contactRawPower,
    "pinchNormalization" -> pinchNormalization
  |>
];


(* ::Section:: *)
(* Correctness boundaries *)

msInternalLineQ[line_Association] := Length[line["endpoints"]] === 2;
msFullLineQ[line_Association] := MemberQ[{"massiveFull", "masslessFull"}, line["type"]];

msNormalizeContactBundles[spec_Association, lines_List] := Module[
  {fullLines, fullIds, explicit, explicitGroups, usedIds, remaining, automaticGroups, groups,
   lineById, invalidGroups},
  fullLines = Select[lines, msFullLineQ];
  fullIds = Lookup[fullLines, "id", {}];
  explicit = Lookup[spec, "thetaBundles", Automatic];
  explicitGroups = Replace[
    explicit,
    {
      Automatic -> {},
      items_List :> Replace[
        items,
        {
          item_Association :> Lookup[item, "lineIds", Missing["LineIds"]],
          item_List :> item,
          item_ :> Missing["MalformedBundle", item]
        },
        {1}
      ],
      other_ :> {Missing["MalformedBundles", other]}
    }
  ];
  If[MemberQ[explicitGroups, _Missing, Infinity],
    Return[msFailure["MalformedContactBundles", <|"value" -> explicit|>]]
  ];
  usedIds = Flatten[explicitGroups];
  If[! DuplicateFreeQ[usedIds] || ! SubsetQ[fullIds, usedIds],
    Return[msFailure[
      "InvalidContactBundleLines",
      <|"fullLineIds" -> fullIds, "explicitLineIds" -> usedIds|>
    ]]
  ];
  lineById = AssociationThread[Lookup[lines, "id"] -> lines];
  invalidGroups = Select[
    explicitGroups,
    # === {} || Length[DeleteDuplicates[Sort /@ Lookup[lineById /@ #, "endpoints"]]] =!= 1 &
  ];
  If[invalidGroups =!= {},
    Return[msFailure[
      "ContactBundleMustShareThetaArgument",
      <|"groups" -> invalidGroups|>
    ]]
  ];
  remaining = Select[fullLines, ! MemberQ[usedIds, # ["id"]] &];
  automaticGroups = Lookup[#, "id"] & /@ GatherBy[remaining, Sort[# ["endpoints"]] &];
  groups = Join[explicitGroups, automaticGroups];
  MapIndexed[
    <|
      "bundleId" -> "thetaBundle:" <> ToString[First[#2]],
      "lineIds" -> #1,
      "rootEndpoints" -> If[#1 === {}, {}, Sort[lineById[First[#1]]["endpoints"]]],
      "source" -> If[MemberQ[explicitGroups, #1], "explicit", "inferredFromCommonEndpoints"]
    |> &,
    groups
  ]
];


msValidateTreeInput[spec_Association, vertices_List, lines_List] := Module[
  {issues = {}, vertexIds, internalLines, graph, expectedEndpointCount, graphMode,
   connectedQ, edgeCount, bundleData, expectedBranches},
  vertexIds = Lookup[vertices, "id"];
  If[MemberQ[vertexIds, _Missing] || ! DuplicateFreeQ[vertexIds],
    AppendTo[issues, <|"code" -> "vertexIdsMustBePresentAndUnique", "ids" -> vertexIds|>]
  ];
  Do[
    If[! MemberQ[{"+", "-"}, vertex["vertexType"]],
      AppendTo[issues, <|"code" -> "vertexTypeMustBePlusOrMinus", "vertex" -> vertex["id"],
        "value" -> vertex["vertexType"]|>]
    ],
    {vertex, vertices}
  ];
  Do[
    If[! MemberQ[$msSupportedInputLineTypes, line["inputType"]],
      AppendTo[issues, <|"code" -> "unsupportedLineType", "line" -> line["id"],
        "type" -> line["inputType"], "allowed" -> $msSupportedInputLineTypes|>]
    ];
    expectedEndpointCount = Length[line["endpoints"]];
    If[! MemberQ[{1, 2}, expectedEndpointCount],
      AppendTo[issues, <|"code" -> "lineMustHaveOneOrTwoEndpoints", "line" -> line["id"],
        "actual" -> expectedEndpointCount|>]
    ];
    If[! SubsetQ[vertexIds, line["endpoints"]],
      AppendTo[issues, <|"code" -> "unknownEndpoint", "line" -> line["id"], "endpoints" -> line["endpoints"]|>]
    ];
    If[MemberQ[{"massiveFull", "massiveCross", "massiveExternal"}, line["type"]] && Head[line["nu"]] === Missing,
      AppendTo[issues, <|"code" -> "missingNuMagnitude", "line" -> line["id"]|>]
    ];
    If[Head[line["momentum"]] === Missing,
      AppendTo[issues, <|"code" -> "missingMomentum", "line" -> line["id"]|>]
    ];
    If[StringStartsQ[line["type"], "massive"],
      expectedBranches = If[line["type"] === "massiveExternal", 1, 2];
      If[! StringQ[line["hankelBranches"]] &&
         (! ListQ[line["hankelBranches"]] || Length[line["hankelBranches"]] =!= expectedBranches ||
          ! And @@ (MemberQ[{1, 2}, #] & /@ line["hankelBranches"])),
        AppendTo[issues, <|"code" -> "invalidHankelBranches", "line" -> line["id"],
          "value" -> line["hankelBranches"]|>]
      ]
    ],
    {line, lines}
  ];
  internalLines = Select[lines, msInternalLineQ];
  graphMode = Lookup[spec, "graphMode", "Tree"];
  If[! MemberQ[{"Tree", "TimeOnly"}, graphMode],
    AppendTo[issues, <|"code" -> "unknownGraphMode", "value" -> graphMode|>]
  ];
  If[vertexIds =!= {} && FreeQ[vertexIds, _Missing],
    graph = Graph[vertexIds, UndirectedEdge @@@ Lookup[internalLines, "endpoints"]];
    connectedQ = TrueQ[ConnectedGraphQ[graph]];
    edgeCount = Length[internalLines];
    If[Length[vertexIds] > 1 && ! connectedQ,
      AppendTo[issues, <|"code" -> "internalGraphMustBeConnected", "edgeCount" -> edgeCount|>]
    ];
    If[graphMode === "Tree" && Length[vertexIds] > 1 && edgeCount =!= Length[vertexIds] - 1,
      AppendTo[issues, <|"code" -> "internalGraphMustBeATree", "edgeCount" -> edgeCount,
        "vertexCount" -> Length[vertexIds]|>]
    ]
  ];
  bundleData = msNormalizeContactBundles[spec, lines];
  If[Head[bundleData] === Failure,
    AppendTo[issues, <|"code" -> First[bundleData], "data" -> Last[bundleData]|>]
  ];
  issues
];
