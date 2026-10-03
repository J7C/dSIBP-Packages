(* ::Package:: *)

(* ::Chapter:: *)
(*初始化与 metadata 序列化*)

Options[DSInit] = {
   WriteInitializationFiles -> False,
   InitializationDirectory -> Automatic,
   GenerateDerivativeMetadata -> False,
   OverwriteInitialization -> False,
   RegisterAsCurrent -> True,
   ProgressReporting -> Automatic,
   KinematicRules -> Automatic
   };

DSInit::badinput = "DSInit 输入不是有效的 topology Association，或 ISP/动量坐标不闭合。";
DSInit::sectorincomplete = "无法完整初始化 contact-reachable sectors：`1`。";
DSInit::initconflict = "初始化目录 `1` 已含不同输入哈希或未知文件；如确认覆盖，请显式设置 OverwriteInitialization -> True。";
DSInit::writefailed = "初始化 metadata 写入失败：`1`。";
DSInfo::noinit = "当前没有已注册的 DSInit context。";
DSInfo::badcontext = "给定对象不是有效的 DSInit context。";

dsInputHash[input_Association] := IntegerString[Hash[HoldComplete[input], "SHA256"], 16, 64];

dsCallerDirectory[] := Which[
   StringQ[$InputFileName] && $InputFileName =!= "", DirectoryName[$InputFileName],
   TrueQ[$Notebooks], With[{directory = Quiet[NotebookDirectory[]]}, If[StringQ[directory], directory, Directory[]]],
   True, Directory[]
   ];

dsResolveInitializationDirectory[Automatic] := FileNameJoin[{dsCallerDirectory[], "init"}];
dsAbsolutePathQ[path_String] := StringStartsQ[path, "/"] || StringStartsQ[path, "\\\\"] ||
   (StringLength[path] >= 3 && StringMatchQ[StringTake[path, 1], LetterCharacter] &&
     StringTake[path, {2, 2}] === ":" && MemberQ[{"\\", "/"}, StringTake[path, {3, 3}]]);
dsResolveInitializationDirectory[path_String] := ExpandFileName[
   If[dsAbsolutePathQ[path], path, FileNameJoin[{dsCallerDirectory[], path}]]
   ];
dsResolveInitializationDirectory[_] := $Failed;

dsDerivativeMetadata[topo_Association, progressSetting_] := Module[{generators, operators},
   generators = makeIndependentVariableDerivativeGenerators[topo];
   operators = dsProgressMap[
     "正在生成微分算符 / Building differential operators",
     generators,
     Function[generator,
      <|
       "variable" -> generator["variable"],
       "userVariable" -> generator["userVariable"],
       "kind" -> generator["kind"],
        "decomposition" -> Switch[generator["kind"],
          "externalInvariant",
          makeExternalInvariantDerivativeDecomposition[topo, generator["variable"]],
          "kinematicCoordinate",
          <|
           "status" -> "chainRuleAdapter",
           "atomicCoordinates" -> Lookup[kinematicAtomicDerivativeData[topo], "inputExpression", {}],
           "atomicJacobian" -> Lookup[generator, "atomicJacobian", {}]
           |>,
          _,
          Missing["DirectExternalLegEnergyDerivative"]
          ]
       |>
      ],
     progressSetting
     ];
   <|"status" -> If[FreeQ[operators, $Failed], "generated", "failed"], "variableCount" -> Length[generators], "operators" -> operators|>
   ];

(* 内部 issue code 保持稳定供程序判断；公开文本只解释原因和下一步，不暴露驼峰状态名。 *)
dsInitializationIssueText[issue_Association] := Module[{code, details, sentence, detailText},
   code = Lookup[issue, "code", "unknownInitializationIssue"];
   details = KeyDrop[issue, {"severity", "code"}];
   sentence = Switch[code,
     Alternatives[
      "missingRequiredCaseKeys", "malformedVertices", "verticesMissingRequiredKeys",
      "malformedLines", "linesMissingRequiredKeys", "malformedLineEndpoints",
      "malformedLoopMomenta", "malformedLoopExternalMomenta",
      "malformedIndependentExternalMomenta", "malformedIBPMode", "malformedSymmetryRules",
      "malformedISPData", "ispDataMissingRequiredKeys", "topologyParseFailed"
      ],
     "拓扑输入缺少必需字段，或字段的类型和结构不正确。请按当前 topology 格式补齐并修正输入后重新初始化。 " <>
      "The topology input is missing required fields or contains fields with an invalid type or structure. Complete and correct the current topology schema, then initialize again.",

     Alternatives[
      "malformedEndpoints", "unknownEndpointVertices", "lineEndpointNotInVertexData",
      "activeVertexIdsNotInVertexData", "fixedAVertexValuesNotInVertexData",
      "malformedExtLegs", "extLegVertexNotInVertexData", "malformedExternalLegEnergies",
      "sectorExternalLegEnergyNotInVertices"
      ],
     "线或外腿引用了无效顶点，或端点结构不完整。请让每个端点和外腿都引用已声明的顶点。 " <>
      "A line or external leg refers to an invalid vertex, or its endpoint data are incomplete. Make every endpoint and external leg refer to a declared vertex.",

     Alternatives[
      "missingFunctionSystemKeys", "functionSystemTNot2x2", "functionSystemTSingular",
      "functionSystemWInconsistent", "functionSystemWTInconsistent",
      "functionSystemTraceInconsistent", "functionSystemDerivativeNotFiniteLaurent",
      "functionSystemWTNotFiniteCompatibleLaurent", "invalidLineFunctionSystem"
      ],
     "传播子的函数系统不完整、不可逆或不满足声明的微分关系。请修正完整 functionSystem 后重新初始化。 " <>
      "A propagator function system is incomplete, singular, or inconsistent with its declared differential relations. Correct the complete functionSystem and initialize again.",

     Alternatives[
      "invalidIBPMode", "loopMomentumCountMismatch", "loopRoutingRankMismatch",
      "unsupportedLoopRoutingCoefficients", "loopRoutingOutsideCycleSpace",
      "bridgeCarriesLoopMomentum", "nonlinearLoopMomentumRouting"
      ],
     "圈动量声明与图的圈空间或路由不一致。请给出独立圈动量，并让 bridge 不携带圈动量。 " <>
      "The loop-momentum declaration is inconsistent with the graph cycle space or routing. Supply independent loop momenta and keep bridge lines free of loop momentum.",

     Alternatives[
      "undercompleteLoopExternalMomenta", "undercompleteIndependentExternalMomenta",
      "incompleteKinematicCoordinates"
      ],
     "动量或动力学坐标欠完备，无法唯一表示当前函数族。请补齐缺失的独立方向或坐标后重新初始化。 " <>
      "The momentum or kinematic coordinates are undercomplete and cannot represent the current family uniquely. Add the missing independent directions or coordinates, then initialize again.",

     Alternatives[
      "overcompleteLoopExternalMomenta", "overcompleteIndependentExternalMomenta",
      "overcompleteKinematicCoordinates"
      ],
     "动量或动力学坐标过完备。符号 IBP 可以继续，但唯一微分算符和后端导出已禁用；请改用一组独立坐标以启用这些功能。 " <>
      "The momentum or kinematic coordinates are overcomplete. Symbolic IBP remains available, but unique differential operators and backend export are disabled; choose an independent coordinate set to enable them.",

     Alternatives[
      "duplicateLoopMomenta", "duplicateExternalMomenta", "duplicateExternalLegMomenta",
      "loopExternalMomentumOverlap", "externalLegMomentumOverlap"
      ],
     "动量列表含重复项或不同角色之间有重叠。请删除重复项，并分别声明圈动量、进入圈传播子的外动量和独立外腿动量。 " <>
      "The momentum lists contain duplicates or overlap between roles. Remove duplicates and declare loop momenta, external momenta entering loop propagators, and independent external-leg momenta separately.",

     Alternatives[
      "duplicateLineIds", "duplicateVertexIds", "duplicateISPNames", "unknownVertexSigns",
      "unknownPackTypes", "unknownLineMassTypes", "unknownLineSKTypes", "unknownLineStates"
      ],
     "拓扑标识或离散类型无效。请使用唯一标识，并把顶点、质量和线状态改为当前接口允许的取值。 " <>
      "A topology identifier or discrete type is invalid. Use unique identifiers and values allowed by the current interface for vertex, mass, and line states.",

     Alternatives[
      "invalidExternalLegEnergyMomentumDependence", "undeclaredMomentumVariables",
      "nonLinearLineMomenta", "nonLinearScalarProductArguments"
      ],
     "动量或外腿能量表达式使用了未声明变量或非线性路由。请只使用已声明动量的线性组合，并把标量积写在这些组合上。 " <>
      "A momentum or external-leg-energy expression uses undeclared variables or nonlinear routing. Use only linear combinations of declared momenta and form scalar products from those combinations.",

     Alternatives[
      "malformedISPRangeSpecs", "invalidZeroPointRules", "unsupportedISPExpressions",
      "insufficientISPData", "scalarProductCoordinateCountMismatch",
      "scalarProductCoordinateSolveFailed"
      ],
     "ISP、零点规则或标量积坐标不能闭合当前函数族。请修正 ISP 表达式和范围，并提供足够的独立坐标。 " <>
      "The ISP data, zero-point rules, or scalar-product coordinates do not close the current family. Correct the ISP expressions and ranges and provide enough independent coordinates.",

     "unsupportedSeedFeatures",
     "当前输入含尚未支持的 seed 功能，因此相关 producer 已禁用。请移除这些功能或改用当前支持的函数族设置。 " <>
      "The input requests seed features that are not yet supported, so the affected producer is disabled. Remove those features or use a currently supported family setup.",

     _,
     "初始化检查发现无法归类的输入问题。请查看 validationReport[\"issues\"] 中保存的结构化详情，修正输入后重新初始化。 " <>
      "Initialization found an input problem that could not be classified for display. Inspect the structured details saved in validationReport[\"issues\"], correct the input, and initialize again."
     ];
   detailText = If[
     details === <||>,
     "",
     " 相关输入：" <> ToString[details, InputForm] <> ". Relevant input: " <> ToString[details, InputForm] <> "."
     ];
   sentence <> detailText
   ];

dsReadExistingManifest[path_String] := If[FileExistsQ[path], Quiet[Check[Get[path], $Failed]], Missing["NoManifest"]];


(* 初始化 metadata 固定为 UTF-8/LF 的可读 InputForm，避免 Windows Put 产生 CRLF 与行尾空格。 *)
dsWriteInitializationExpression[expr_, path_String] := Module[{text},
   text = ToString[expr, InputForm, PageWidth -> 120] <> "\n";
   text = StringReplace[text, RegularExpression["[ \\t]+(?=\\n|$)"] -> ""];
   writeKiraUTF8LFText[path, text]
   ];


dsInitializationConflictQ[directory_String, inputHash_String, overwriteQ_] := Module[{manifestPath, manifest, knownFiles},
   If[TrueQ[overwriteQ] || ! DirectoryQ[directory], Return[False]];
   manifestPath = FileNameJoin[{directory, "manifest.wl"}];
   manifest = dsReadExistingManifest[manifestPath];
   knownFiles = FileExistsQ /@ (FileNameJoin[{directory, #}] & /@ {"topology.wl", "sectors.wl", "conventions.wl", "derivatives.wl"});
   Which[
    AssociationQ[manifest], Lookup[manifest, "inputHash", Missing["inputHash"]] =!= inputHash,
    Head[manifest] === Missing && ! Or @@ knownFiles, False,
    True, True
    ]
   ];

dsWriteInitializationFiles[context_Association, directory_String, overwriteQ_] := Module[
   {manifestPath, fileData, paths, manifest, writeResult},
   If[dsInitializationConflictQ[directory, context["inputHash"], overwriteQ], Return[<|"status" -> "conflict", "directory" -> directory|>]];
   Quiet[CreateDirectory[directory, CreateIntermediateDirectories -> True]];
   fileData = <|
     "topology.wl" -> context["topologyData"],
     "sectors.wl" -> context["sectors"],
     "conventions.wl" -> context["conventions"]
     |>;
   If[AssociationQ[context["derivatives"]], AssociateTo[fileData, "derivatives.wl" -> context["derivatives"]]];
   paths = AssociationMap[FileNameJoin[{directory, #}] &, Keys[fileData]];
   writeResult = Quiet[Check[KeyValueMap[(dsWriteInitializationExpression[#2, paths[#1]]; #1) &, fileData], $Failed]];
   If[writeResult === $Failed, Return[<|"status" -> "failed", "directory" -> directory|>]];
   manifest = <|
     "status" -> "initialized",
     "packageVersion" -> context["packageVersion"],
     "inputHash" -> context["inputHash"],
     "caseName" -> context["caseName"],
     "generatedAt" -> DateString[{"ISODate", "T", "Time", "TimeZone"}],
      "files" -> Map[FileNameTake, paths],
     "sectorCount" -> Length[context["sectors"]],
     "derivativeMetadataQ" -> AssociationQ[context["derivatives"]]
     |>;
   manifestPath = FileNameJoin[{directory, "manifest.wl"}];
   If[Quiet[Check[dsWriteInitializationExpression[manifest, manifestPath]; True, False]] =!= True, Return[<|"status" -> "failed", "directory" -> directory|>]];
   <|"status" -> "written", "directory" -> directory, "manifest" -> manifestPath, "files" -> Append[paths, "manifest.wl" -> manifestPath]|>
   ];

(* ::Section:: *)
(*时间 IBP 的小 t 整数端点门禁*)

(* 只把可以 exact 证明的整数归入门禁；含通用符号解析正规化参数的表达式保持未判定。 *)
dsProvablyIntegerQ[value_] := Module[{simplified = Quiet[FullSimplify[value]]},
   IntegerQ[simplified] || TrueQ[Quiet[FullSimplify[Element[simplified, Integers]]]]
   ];


(* h/H 的两支领头幂次来自固定局部解；state=1 表示对函数自变量求一次导数。 *)
dsKnownSmallTLeadingPowerTable[line_Association] := Module[
   {compiled, variable, p, q, nuValue, mode},
   compiled = Lookup[line, "compiledFunctionSystem", <||>];
   variable = Lookup[compiled, "variable", Missing["Variable"]];
   p = Lookup[compiled, "P", Missing["P"]];
   q = Lookup[compiled, "Q", Missing["Q"]];
   nuValue = Lookup[line, "nu", nu];
   mode = Which[
     Head[variable] =!= Missing &&
       functionSystemZeroQ[p - (2 nuValue + 1)/variable] &&
       functionSystemZeroQ[q - 1],
       "h",
     Head[variable] =!= Missing &&
       functionSystemZeroQ[p - 1/variable] &&
       functionSystemZeroQ[q - (1 - nuValue^2/variable^2)],
       "H",
     True,
       Missing["UnsupportedSmallTFunctionSystem"]
     ];
   Switch[mode,
    "h", <|0 -> Simplify /@ {-2 nuValue, 0}, 1 -> Simplify /@ {-2 nuValue - 1, 1}|>,
    "H", <|0 -> Simplify /@ {-nuValue, nuValue}, 1 -> Simplify /@ {-nuValue - 1, nuValue - 1}|>,
    _, mode
    ]
   ];


dsEndpointSmallTChoices[line_Association, endpointSlot_Integer] := Module[{table},
   table = dsKnownSmallTLeadingPowerTable[line];
   If[Head[table] === Missing, Return[table]];
   Flatten[
    Table[
     <|
      "lineId" -> line["id"],
      "endpointSlot" -> endpointSlot,
      "state" -> state,
      "branch" -> branch,
      "power" -> table[state][[branch]]
      |>,
     {state, {0, 1}},
     {branch, {1, 2}}
     ],
    1
    ]
   ];


dsSectorVertexOriginalIds[topo_Association, vertexId_] := Module[{representatives},
   representatives = Lookup[
     topo,
     "sectorVertexRepresentativeMap",
     AssociationThread[topo["vertexIds"] -> topo["vertexIds"]]
     ];
   Select[Keys[representatives], representatives[#] === vertexId &]
   ];


dsTimeIBPVertexCertificate[topo_Association, sectorKey_, vertexId_] := Module[
   {endpointData, unsupported, choiceTables, selections, branchRecords, integerRecords},
   endpointData = Flatten[
    Table[
     If[
      Lookup[line, "massType", "massive"] === "massive" &&
       Lookup[line, "packType", ""] =!= "shrunk",
      Map[
       Function[endpointSlot,
        <|
         "lineId" -> line["id"],
         "endpointSlot" -> endpointSlot,
         "choices" -> dsEndpointSmallTChoices[line, endpointSlot]
         |>
        ],
       Flatten@Position[line["endpoints"], vertexId]
       ],
      {}
      ],
     {line, topo["lines"]}
     ],
    1
    ];
   unsupported = Select[endpointData, Head[Lookup[#, "choices", None]] === Missing &];
   If[unsupported =!= {},
    Return[<|
      "status" -> "leadingPowersUnsupported",
      "sectorKey" -> sectorKey,
      "vertexId" -> vertexId,
      "vertexIds" -> dsSectorVertexOriginalIds[topo, vertexId],
      "unsupportedEndpoints" -> KeyDrop[unsupported, "choices"]
      |>]
    ];
   choiceTables = Lookup[endpointData, "choices", {}];
   selections = If[choiceTables === {}, {{}}, Tuples[choiceTables]];
   branchRecords = MapIndexed[
     Function[{selection, ordinal},
      <|
       "branchOrdinal" -> First[ordinal],
       "localChoices" -> selection,
       "alpha" -> Quiet[FullSimplify[
         vertexZeroPoint[topo, vertexId] + 1 + Total[Lookup[selection, "power", {}]]
         ]]
       |>
      ],
     selections
     ];
   integerRecords = Select[branchRecords, dsProvablyIntegerQ[Lookup[#, "alpha"]] &];
   <|
    "status" -> If[integerRecords === {}, "accepted", "requiresRegularization"],
    "sectorKey" -> sectorKey,
    "vertexId" -> vertexId,
    "vertexIds" -> dsSectorVertexOriginalIds[topo, vertexId],
    "endpointCount" -> Length[endpointData],
    "branchCount" -> Length[branchRecords],
    "branchRecords" -> branchRecords,
    "integerBranchRecords" -> integerRecords
    |>
   ];


dsTimeIBPDomainCertificate[topologyData_Association] := Module[
   {metadataList, sectorTopologies, records, integerIssues, unsupported},
   metadataList = Lookup[topologyData, "sectorMetadataList", {}];
   sectorTopologies = sectorTopologyForMetadata018[topologyData, #] & /@ metadataList;
   records = Flatten[
    MapThread[
     Function[{sectorTopology, metadata},
      dsTimeIBPVertexCertificate[
        sectorTopology,
        Lookup[metadata, "sectorKey", Missing["SectorKey"]],
        #
        ] & /@ Lookup[sectorTopology, "activeVertexIds", sectorTopology["vertexIds"]]
      ],
     {sectorTopologies, metadataList}
     ],
    1
    ];
   integerIssues = Select[records, Lookup[#, "status", None] === "requiresRegularization" &];
   unsupported = Select[records, Lookup[#, "status", None] === "leadingPowersUnsupported" &];
   <|
    "status" -> Which[
      integerIssues =!= {}, "requiresRegularization",
      unsupported =!= {}, "partiallyClassified",
      True, "certified"
      ],
    "vertexRecords" -> records,
    "integerIssues" -> integerIssues,
    "unsupportedRecords" -> unsupported
    |>
   ];


dsTimeIBPIntegerIssueSummary[issues_List] := Map[
   <|
     "sectorKey" -> Lookup[#, "sectorKey"],
     "vertexIds" -> Lookup[#, "vertexIds", {}],
     "triggeredBranches" -> Map[
       KeyTake[#, {"branchOrdinal", "localChoices", "alpha"}] &,
       Lookup[#, "integerBranchRecords", {}]
       ]
     |> &,
   issues
   ];


DSInit[input_Association, OptionsPattern[]] := Module[
   {progress = OptionValue[ProgressReporting], topologyData, validation, subsetSummary, derivatives, context,
     inputHash, initDirectory, writeResult = <|"status" -> "notRequested"|>, warnings, errors,
     effectiveInput, kinematicAudit, declarationAudit, guide, guideText, parityDataList,
     parityRequestedQ, parityUsableQ, parityFailures, timeIBPDomainCertificate},
   effectiveInput = If[
     OptionValue[KinematicRules] === Automatic,
     input,
     Join[input, <|"kinematicRules" -> OptionValue[KinematicRules]|>]
     ];
   inputHash = dsInputHash[effectiveInput];
   topologyData = dsStageRun[
     "初始化 topology、ISP 与完整 sector metadata / Initializing topology, ISP, and complete sector metadata",
     makeTopologyData[
       effectiveInput,
      PrecomputeShrinkSectorMetadata -> True
      ],
     progress
     ];
   validation = Lookup[topologyData, "validationReport", <|"errorCount" -> 1, "issues" -> {}|>];
   kinematicAudit = Lookup[topologyData, "kinematicCoordinateAudit", <||>];
   declarationAudit = Lookup[topologyData, "momentumDeclarationAudit", <||>];
   guide = kinematicParameterRedefinitionGuide[kinematicAudit];
   guideText = If[
     StringQ[Lookup[guide, "commandExample", None]],
     Lookup[guide, "commandExample", ""],
     Lookup[guide, "defaultBehavior", ""]
     ];
   dsInfoPrint[
     "动量角色：loopMomenta " <> ToString[Lookup[topologyData, "loopMomenta", {}], InputForm] <>
      "；loopExternalMomenta " <> ToString[Lookup[topologyData, "loopExternalMomenta", {}], InputForm] <>
      "；effectiveLoopExternalMomenta " <> ToString[Lookup[topologyData, "effectiveLoopExternalMomenta", {}], InputForm] <>
      "；independentExternalMomenta " <> ToString[Lookup[topologyData, "independentExternalMomenta", {}], InputForm] <>
      "；实际需要的 loop 方向 " <> ToString[Lookup[declarationAudit, "requiredLoopExternalDirections", {}], InputForm] <>
      "；实际需要的无圈模长 " <> ToString[Lookup[declarationAudit, "requiredIndependentMomentumMagnitudes", {}], InputForm] <>
      ". Momentum roles: loopMomenta " <> ToString[Lookup[topologyData, "loopMomenta", {}], InputForm] <>
      "; loopExternalMomenta " <> ToString[Lookup[topologyData, "loopExternalMomenta", {}], InputForm] <>
      "; effectiveLoopExternalMomenta " <> ToString[Lookup[topologyData, "effectiveLoopExternalMomenta", {}], InputForm] <>
      "; independentExternalMomenta " <> ToString[Lookup[topologyData, "independentExternalMomenta", {}], InputForm] <>
      "; required loop directions " <> ToString[Lookup[declarationAudit, "requiredLoopExternalDirections", {}], InputForm] <>
      "; required loop-free magnitudes " <> ToString[Lookup[declarationAudit, "requiredIndependentMomentumMagnitudes", {}], InputForm],
     progress
     ];
   dsInfoPrint[
     "动力学变量选择：" <> ToString[Lookup[kinematicAudit, "status", "unknown"]] <>
      "；缺省规则 " <> ToString[Lookup[kinematicAudit, "defaultRules", {}], InputForm] <>
      "；当前规则 " <> ToString[Lookup[kinematicAudit, "selectedRules", {}], InputForm] <>
      "；从属模长绑定 " <> ToString[Lookup[kinematicAudit, "dependentMagnitudeBindings", {}], InputForm] <>
      ". Kinematic-variable selection: " <> ToString[Lookup[kinematicAudit, "status", "unknown"]] <>
      "; default rules " <> ToString[Lookup[kinematicAudit, "defaultRules", {}], InputForm] <>
      "; selected rules " <> ToString[Lookup[kinematicAudit, "selectedRules", {}], InputForm] <>
      "; dependent magnitude bindings " <> ToString[Lookup[kinematicAudit, "dependentMagnitudeBindings", {}], InputForm],
     progress
     ];
   dsInfoPrint[
     "必需模长的参数覆盖 " <> ToString[kinematicRequiredMagnitudeCoverage[topologyData], InputForm] <>
      ". Parameter coverage of required magnitudes " <> ToString[kinematicRequiredMagnitudeCoverage[topologyData], InputForm],
     progress
     ];
   dsInfoPrint[
     "参数可保持缺省，也可复制以下格式重定义：" <> guideText <>
      "。规则左端写原始 sp[...]，不要写 ssij/sEi -> custom；右端写自定义参数表达式。 " <>
      "Keep the default parameters or copy this form to redefine them: " <> guideText <>
      ". Put the original sp[...] on the left, not ssij/sEi -> custom, and the custom parameter expression on the right.",
     progress
     ];
   If[Lookup[topologyData, "status", None] === "invalidInput" || topologyValidationErrorQ[validation],
    errors = Select[Lookup[validation, "issues", {}], Lookup[#, "severity", ""] === "error" &];
    Scan[dsErrorPrint[dsInitializationIssueText[#]] &, errors];
    Message[DSInit::badinput]; dsErrorPrint["topology/ISP 初始化失败；上述详情同时保存在 validationReport[\"issues\"]。 Topology/ISP initialization failed; the details above are also stored in validationReport[\"issues\"]."];
    Return[dsFailedInitializationData[
      "invalidInputOrTopology",
      <|"inputHash" -> inputHash, "topologyData" -> topologyData, "validationReport" -> validation|>
      ]]
    ];
   subsetSummary = Lookup[topologyData, "precomputedShrinkSectorSummary", <||>];
   If[Lookup[subsetSummary, "status", "missing"] =!= "generated" || ! TrueQ[Lookup[subsetSummary, "completeCoverageQ", False]],
    Message[DSInit::sectorincomplete, subsetSummary]; dsErrorPrint["contact-reachable sector 未完整初始化。 Contact-reachable sectors were not initialized completely."];
    Return[dsFailedInitializationData[
      "incompleteSectorMetadata",
      <|"inputHash" -> inputHash, "topologyData" -> topologyData|>
      ]]
    ];
   timeIBPDomainCertificate = dsTimeIBPDomainCertificate[topologyData];
   If[Lookup[timeIBPDomainCertificate, "status", None] === "requiresRegularization",
    dsErrorPrint[
     "这个函数族的小 t 领头幂会在整数指标移位下触及发散的时间积分边界，因此初始化已停止。" <>
      " 请在这些顶点加入通用符号解析正规化参数后重试：" <>
      ToString[dsTimeIBPIntegerIssueSummary[timeIBPDomainCertificate["integerIssues"]], InputForm] <>
      ". The small-t leading powers of this family reach divergent time-integral boundaries under integer index shifts, so initialization was stopped." <>
      " Add a symbolic analytic regulator at the listed vertices and try again: " <>
      ToString[dsTimeIBPIntegerIssueSummary[timeIBPDomainCertificate["integerIssues"]], InputForm]
     ];
    Return[dsFailedInitializationData[
      "timeIBPRequiresRegularization",
      <|
       "inputHash" -> inputHash,
       "topologyData" -> topologyData,
       "timeIBPDomainCertificate" -> timeIBPDomainCertificate
       |>
      ]]
    ];
   If[Lookup[timeIBPDomainCertificate, "status", None] === "partiallyClassified",
    dsWarningPrint[
     "至少一个自定义函数系统没有可直接查表的小 t 领头幂，因此本次初始化只能对已识别的 h/H 端点执行整数边界检查。" <>
      " 如果自定义系统可能在 t=0 发散，请先加入解析正规化。" <>
      " At least one custom function system has no tabulated small-t leading powers, so this initialization checks integer boundaries only for recognized h/H endpoints." <>
      " Add an analytic regulator first if the custom system may diverge at t=0.",
     progress
     ]
    ];
   parityDataList = Lookup[Lookup[topologyData, "sectorMetadataList", {}], "parityData", {}];
   parityRequestedQ = Lookup[topologyData, "parityConstraints", {}] =!= {};
   parityUsableQ = parityDataList =!= {} && And @@ Lookup[parityDataList, "parityUsableQ", {False}];
   parityFailures = Select[
     parityDataList,
     Lookup[#, "status", "disabled"] === "disabled" &&
       Lookup[#, "reason", "noParityConstraints"] =!= "noParityConstraints" &
     ];
   topologyData = Join[topologyData, <|
      "capabilities" -> Join[
        Lookup[topologyData, "capabilities", <||>],
        <|"parityUsableQ" -> parityUsableQ|>
        ]
      |>];
   If[parityRequestedQ && parityFailures =!= {},
    dsErrorPrint[
     "显式 parity 约束无法在全部 sector 上定义，初始化已拒绝：" <>
      ToString[parityFailures, InputForm] <>
      ". Explicit parity constraints are not defined on every sector; initialization was rejected: " <>
      ToString[parityFailures, InputForm]
     ];
    Return[dsFailedInitializationData[
      "invalidParityConstraints",
      <|"inputHash" -> inputHash, "topologyData" -> topologyData,
        "parityFailures" -> parityFailures|>
      ]]
    ];
    If[! parityRequestedQ && ! parityUsableQ,
     dsWarningPrint[
      "当前函数系统没有已证明可运输的 parity 闭合，parity 筛选已禁用；普通 IBP 仍可继续。 " <>
       "The current function system has no proved transportable parity closure; parity filtering is disabled, while ordinary IBP remains available.",
      progress
      ]
     ];
    warnings = Select[Lookup[validation, "issues", {}], Lookup[#, "severity", ""] === "warning" &];
   Scan[dsWarningPrint[dsInitializationIssueText[#], progress] &, warnings];
   derivatives = If[
     TrueQ[OptionValue[GenerateDerivativeMetadata]] && dsTopologyCapabilityQ[topologyData, "derivativeUsableQ"],
     dsDerivativeMetadata[topologyData, progress],
     Missing["NotGenerated"]
     ];
   context = <|
     "status" -> "initialized",
     "packageVersion" -> $dSIBPVersion,
     "inputHash" -> inputHash,
     "caseName" -> Lookup[topologyData, "name", "unnamed"],
      "input" -> effectiveInput,
     "topology" -> topologyData,
     "topologyData" -> topologyData,
     "capabilities" -> Lookup[topologyData, "capabilities", <||>],
     "sectors" -> Lookup[topologyData, "sectorMetadataList", {}],
     "conventions" -> dsConventionMetadata[topologyData],
     "derivatives" -> derivatives,
     "validationReport" -> validation,
     "timeIBPDomainCertificate" -> timeIBPDomainCertificate,
     "loopTreeProjectionConvention" -> <|
       "targetAZeroPointBecomesTreeNu0" -> True,
       "removedLineZeroPointsBecomeExplicitEnergyPowers" -> True,
       "relativePhysicalPowerNormalization" -> True,
       "unsafePowerExpand" -> False
       |>,
     "initializationWrite" -> writeResult
     |>;
   If[TrueQ[OptionValue[WriteInitializationFiles]],
    initDirectory = dsResolveInitializationDirectory[OptionValue[InitializationDirectory]];
    If[initDirectory === $Failed,
     Message[DSInit::writefailed, OptionValue[InitializationDirectory]];
     dsWarningPrint["InitializationDirectory 无效；内存数学 context 仍然有效。 InitializationDirectory is invalid; the in-memory mathematical context remains valid."];
     writeResult = <|"status" -> "failed", "reason" -> "invalidInitializationDirectory",
       "requestedDirectory" -> OptionValue[InitializationDirectory]|>,
     writeResult = dsWriteInitializationFiles[context, initDirectory, OptionValue[OverwriteInitialization]];
     If[writeResult["status"] === "conflict",
      Message[DSInit::initconflict, initDirectory];
      dsWarningPrint["已有初始化信息与当前输入不一致，未覆盖；内存数学 context 仍然有效。 Existing initialization data do not match the current input and were not overwritten; the in-memory mathematical context remains valid."]
      ];
     If[! MemberQ[{"written", "conflict"}, writeResult["status"]],
      Message[DSInit::writefailed, initDirectory];
      dsWarningPrint["初始化文件未完整写入；内存数学 context 仍然有效。 Initialization files were not written completely; the in-memory mathematical context remains valid."]
      ]
     ];
    context = Join[context, <|"initializationWrite" -> writeResult|>]
    ];
   If[TrueQ[OptionValue[RegisterAsCurrent]],
    $dSIBPCurrentContext = context;
    setIBPTopologyContext[context["topology"]]
    ];
   dsInfoPrint[
    "初始化完成：" <> context["caseName"] <> "，sector " <> ToString[Length[context["sectors"]]] <> "/" <> ToString[Length[context["sectors"]]] <>
     ". Initialization completed: " <> context["caseName"] <> ", sectors " <> ToString[Length[context["sectors"]]] <> "/" <> ToString[Length[context["sectors"]]],
    progress
    ];
   context
   ];

DSInit[input_, OptionsPattern[]] := (
   Message[DSInit::badinput];
   dsErrorPrint["DSInit 需要 Association 输入。 DSInit requires an Association input."];
   dsFailedInitializationData["inputNotAssociation", <|"input" -> HoldForm[input]|>]
   );

DSInfo[] := Module[{context = dsResolveContext[Automatic]},
   If[Head[context] === Missing, Message[DSInfo::noinit]; Return[<|"status" -> "notInitialized"|>]];
   DSInfo[context]
   ];

DSInfo[context_Association] := Module[{resolved = dsResolveContext[context]},
   If[Head[resolved] === Missing, Message[DSInfo::badcontext]; Return[<|"status" -> "invalidContext"|>]];
   Join[<|"status" -> "initialized"|>, dsContextSummary[resolved], <|
     "sectorCount" -> Length[resolved["sectors"]],
     "independentVariables" -> Lookup[resolved["conventions"], "independentVariables", {}],
     "initializationWrite" -> Lookup[resolved, "initializationWrite", <||>]
     |>]
   ];

DSInfo[context_Association, "Full"] := Module[{resolved = dsResolveContext[context]},
   If[Head[resolved] === Missing, Message[DSInfo::badcontext]; <|"status" -> "invalidContext"|>, resolved]
   ];
