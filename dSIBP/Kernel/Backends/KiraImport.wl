(* ::Package:: *)

(* ::Chapter:: *)
(*Kira 结果取回*)

Options[DSKiraImport] = {
   KiraReductionFile -> Automatic,
   KiraMasterFile -> Automatic,
   KiraCompletionFile -> Automatic,
   KiraCompletionPatterns -> Automatic,
   ProgressReporting -> Automatic
   };

DSKiraImport::badpath = "Kira workspace 路径不存在或不是目录：`1`。";
DSKiraImport::missing = "Kira 结果缺少必需文件：`1`。";
DSKiraImport::incomplete = "Kira 完成日志没有成功标记：`1`。";
DSKiraImport::mismatch = "Kira 结果与当前 export/context 不一致：`1`。";
DSKiraImport::invalid = "Kira reduction 数据未通过完整性检查：`1`。";

dsKiraResolveFile[root_String, Automatic, relative : {__String}] := FileNameJoin[Prepend[relative, root]];
dsKiraResolveFile[root_String, path_String, _List] := ExpandFileName[If[dsAbsolutePathQ[path], path, FileNameJoin[{root, path}]]];
dsKiraResolveFile[_String, _, _List] := $Failed;

dsKiraCompletionPatterns[Automatic] := {
   RegularExpression["(?i)kira[^\\n]*(?:finished|completed)[^\\n]*(?:success|successfully)"],
   RegularExpression["(?i)all jobs[^\\n]*(?:finished|completed)[^\\n]*(?:success|successfully)"],
   RegularExpression["(?i)unreduced integrals:\\s*0\\.?"],
   "Kira finished successfully"
   };
dsKiraCompletionPatterns[patterns_List] := patterns;
dsKiraCompletionPatterns[pattern_] := {pattern};

dsKiraCompletionQ[text_String, patterns_] := AnyTrue[
   dsKiraCompletionPatterns[patterns],
   Function[pattern, StringContainsQ[text, pattern] || StringMatchQ[text, ___ ~~ pattern ~~ ___]]
   ];

dsKiraReadExpression[path_String] := Quiet[Check[Get[path], $Failed]];

dsKiraMasterIDs[text_String] := DeleteDuplicates @ ToExpression @ StringCases[
   text,
   StartOfLine ~~ WhitespaceCharacter ... ~~ id : NumberString ~~ WhitespaceCharacter ... ~~ "#" :> id
   ];

dsRuleListQ[rules_] := ListQ[rules] && And @@ (MatchQ[Unevaluated[#], _Rule | _RuleDelayed] & /@ rules);

dsKiraIntegralTokenQ[expr_] := MatchQ[expr, _J | dsTreeToken[_String, J[_List]]];

dsJToIDPairs[rules_List] := Select[
   Cases[rules, HoldPattern[Rule[integral_, id_Integer]] :> {integral, id}],
   dsKiraIntegralTokenQ[First[#]] &
   ];
dsIDToJPairs[rules_List] := Select[
   Cases[rules, HoldPattern[Rule[Tuserweight[id_Integer], integral_]] :> {id, integral}],
   dsKiraIntegralTokenQ[Last[#]] &
   ];

(* 双向文本相等还不够：同步重复的积分或 ID 也会相等，但并不是可逆编号。 *)
dsKiraInverseMapQ[jToKira_List, kiraToJ_List] := Module[{forward, backward},
   forward = SortBy[dsJToIDPairs[jToKira], Last];
   backward = SortBy[Reverse /@ dsIDToJPairs[kiraToJ], Last];
   forward === backward &&
    Length[forward] === Length[jToKira] && Length[backward] === Length[kiraToJ] &&
    DuplicateFreeQ[First /@ forward] && DuplicateFreeQ[Last /@ forward]
   ];

dsKiraRuleIDs[rules_List] := <|
   "lhs" -> DeleteDuplicates @ Cases[rules, HoldPattern[Tuserweight[id_Integer] -> _] :> id],
   "rhs" -> DeleteDuplicates @ Cases[Last /@ rules, Tuserweight[id_Integer] :> id, Infinity],
   "all" -> DeleteDuplicates @ Cases[rules, Tuserweight[id_Integer] :> id, Infinity]
   |>;

dsKiraCoefficientVariables[rules_List] := Module[{rhs},
   rhs = Last /@ rules /. Tuserweight[_Integer] -> 1;
   DeleteDuplicates[Variables[Together /@ rhs]]
   ];

dsKiraBackendVariableNames[rules_List] := DeleteDuplicates[
   SymbolName /@ Select[dsKiraCoefficientVariables[rules], Head[#] === Symbol &]
   ];

dsKiraBackendRestoreRules[variableMap_List, imaginaryUnit_] := Join[
   Map[
    Function[item,
     With[{name = Lookup[item, "backend"], original = Lookup[item, "original"]},
      HoldPattern[s_Symbol /; SymbolName[s] === name] :> original
      ]
     ],
    variableMap
    ],
   If[StringQ[imaginaryUnit],
    {With[{name = imaginaryUnit}, HoldPattern[s_Symbol /; SymbolName[s] === name] :> I]},
    {}
    ]
   ];

dsKiraRestoreBackendCoefficients[rules_List, variableMap_List, imaginaryUnit_] :=
   rules /. dsKiraBackendRestoreRules[variableMap, imaginaryUnit];


(* ::Section::Closed:: *)
(*Kira 内部虚轴运动学坐标恢复*)

dsKiraBackendKinematicConventionDataQ[data_Association] := Module[
   {status, records, physical, backend, names, forward, backward},
   status = Lookup[data, "status", "notRequired"];
   If[status === "notRequired", Return[True]];
   If[status =!= "configured" ||
     Lookup[data, "scope", Missing["scope"]] =!= "KiraBackendOnly" ||
     Lookup[data, "physicalToBackendConvention", Missing["convention"]] =!=
       "physicalKinematic == -I backendKinematic" ||
     Lookup[data, "derivativeConvention", Missing["derivativeConvention"]] =!=
       "D[physicalKinematic] == I D[backendKinematic]" ||
     Lookup[data, "eulerConvention", Missing["eulerConvention"]] =!=
       "physicalKinematic D[physicalKinematic] == backendKinematic D[backendKinematic]",
    Return[False]
    ];
   records = Lookup[data, "records", Missing["records"]];
   If[! ListQ[records] || records === {} || ! And @@ (AssociationQ /@ records), Return[False]];
   physical = Lookup[records, "physical", {}];
   backend = Lookup[records, "backend", {}];
   names = Lookup[records, "backendName", {}];
   forward = Lookup[records, "physicalToBackendRule", {}];
   backward = Lookup[records, "backendToPhysicalRule", {}];
   And[
    And @@ (Head[#] === Symbol & /@ physical),
    And @@ (Head[#] === Symbol & /@ backend),
    And @@ (StringQ /@ names),
    DuplicateFreeQ[physical], DuplicateFreeQ[backend], DuplicateFreeQ[names],
    backend === (kiraBackendSymbol /@ names),
    forward === MapThread[Rule[#1, -I #2] &, {physical, backend}],
    backward === MapThread[Rule[#1, I #2] &, {backend, physical}],
    Lookup[data, "physicalToBackendRules", Missing["forward"]] === forward,
    Lookup[data, "backendToPhysicalRules", Missing["backward"]] === backward,
    And @@ (# === I & /@ Lookup[records, "physicalDerivativeFromBackendFactor", {}]),
    And @@ (# === -I & /@ Lookup[records, "backendDerivativeFromPhysicalFactor", {}]),
    And @@ (TrueQ /@ Lookup[records, "eulerOperatorInvariantQ", {}])
    ]
   ];


dsKiraBackendKinematicRuleDataQ[data_Association, topo_Association, convention_Association] := Module[
   {userRules, expected, keys},
   userRules = Lookup[data, "userRules", Missing["userRules"]];
   If[! dsRuleListQ[userRules], Return[False]];
   expected = kiraKinematicCoefficientRulesData[userRules, topo, convention];
   keys = {
     "status", "userRules", "normalizedPhysicalRules", "backendRules",
     "physicalPointRules", "kinematicNumericRules", "coordinateConsistencyQ"
     };
   AssociationQ[expected] && KeyTake[data, keys] === KeyTake[expected, keys]
   ];


dsKiraRestorePhysicalKinematicVariables[rules_List, data_Association] := Module[{restoreRules},
   If[Lookup[data, "status", "notRequired"] =!= "configured", Return[rules]];
   restoreRules = Map[
     Function[record,
      With[{name = Lookup[record, "backendName"], physical = Lookup[record, "physical"]},
       HoldPattern[s_Symbol /; SymbolName[Unevaluated[s]] === name] :> I physical
       ]
      ],
     Lookup[data, "records", {}]
     ];
   rules /. restoreRules
   ];


(* ::Section::Closed:: *)
(*Gaussian 相位规范恢复*)

dsKiraGaussianPhaseGaugeDataQ[data_Association, integralCount_Integer] := Module[
   {status, convention, phaseRules, ids, phases},
   status = Lookup[data, "status", "notApplicable"];
   If[status === "notApplicable", Return[True]];
   If[! MemberQ[{"applied", "notRequired"}, status], Return[False]];
   convention = Lookup[data, "physicalToBackendConvention", Missing["convention"]];
   phaseRules = Lookup[data, "integralPhaseRules", Missing["integralPhaseRules"]];
   If[convention =!= "J[id] == I^phase[id] Kira[id]" || ! dsRuleListQ[phaseRules], Return[False]];
   ids = First /@ phaseRules;
   phases = Last /@ phaseRules;
   Lookup[data, "integralCount", Missing["integralCount"]] === integralCount &&
    ids === Range[integralCount] && And @@ (MemberQ[{0, 1}, #] & /@ phases) &&
    DuplicateFreeQ[ids] && TrueQ[Lookup[data, "conflictCount", 0] === 0]
   ];


dsKiraRestorePhysicalPhaseGauge[rules_List, data_Association] := Module[
   {status = Lookup[data, "status", "notApplicable"], phaseByID},
   If[status === "notApplicable", Return[rules]];
   phaseByID = Association[Lookup[data, "integralPhaseRules", {}]];
   rules /. HoldPattern[Rule[Tuserweight[lhsID_Integer], rhs_]] :>
     Rule[
      Tuserweight[lhsID],
      Expand[I^Lookup[phaseByID, lhsID, 0] (rhs /.
          Tuserweight[rhsID_Integer] :> I^(-Lookup[phaseByID, rhsID, 0]) Tuserweight[rhsID])]
      ]
   ];


dsKiraContextMatchQ[manifest_Association, context_Association] := And[
   Lookup[Lookup[manifest, "context", <||>], "inputHash", Missing["inputHash"]] === Lookup[context, "inputHash", Missing["contextHash"]],
   Lookup[manifest, "zeroPointRules", Missing["zeroPointRules"]] === Lookup[context["topology"], "zeroPointRules", Missing["contextZeroPointRules"]],
   Lookup[manifest, "symmetryRules", Missing["symmetryRules"]] === Lookup[context["topology"], "symmetryRules", Missing["contextSymmetryRules"]],
   Lookup[manifest, "tadpoleSymmetryData", Missing["tadpoleSymmetryData"]] === dsStableTadpoleSymmetryData[Lookup[context["topology"], "tadpoleSymmetryData", tadpoleSymmetryData[context["topology"]]]],
   Lookup[manifest, "loopTreeProjectionConvention", Missing["projectionConvention"]] === Lookup[context, "loopTreeProjectionConvention", Missing["contextProjectionConvention"]]
   ];


dsKiraExportFileDigestsQ[identity_Association, workspace_String] := Module[
   {root = ExpandFileName[workspace], records, paths},
   records = Lookup[identity, "exportFiles", Missing["exportFiles"]];
   If[! ListQ[records] || records === {}, Return[False]];
   paths = FileNameJoin[{root, Lookup[#, "path", ""]}] & /@ records;
   And @@ MapThread[
     Function[{record, path},
      StringStartsQ[ExpandFileName[path], root] && FileExistsQ[path] &&
       IntegerString[FileHash[path, "SHA256"], 16, 64] === Lookup[record, "sha256", Missing["sha256"]]
      ],
     {records, paths}
     ]
   ];


dsKiraArtifactIdentityQ[
   manifest_Association,
   workspace_String,
   repJ2Kira_List
   ] := Module[{identity, payload, contract, serializedKiraInputDigestQ},
   identity = Lookup[manifest, "artifactIdentity", Missing["artifactIdentity"]];
   contract = Lookup[manifest, "linearArtifactContract", <||>];
   If[! AssociationQ[identity] || ! AssociationQ[contract], Return[False]];
   payload = KeyDrop[identity, "exportContentDigest"];
   serializedKiraInputDigestQ = TrueQ[
     Lookup[identity, "identityVersion", Missing["identityVersion"]] === 2 &&
      Lookup[identity, "serializedKiraInputDigest", Missing["serializedKiraInputDigest"]] ===
       dsKiraExportFileDigest[Lookup[identity, "exportFiles", {}], "userSystem/ibp.kira"]
     ];
   TrueQ[
    Lookup[identity, "exportContentDigest", Missing["exportContentDigest"]] === dsKiraExpressionDigest[payload] &&
     serializedKiraInputDigestQ &&
     Lookup[identity, "linearSourceDigest", Missing["linearSourceDigest"]] === Lookup[contract, "sourceDigest", Missing["sourceDigest"]] &&
     Lookup[identity, "integralMapDigest", Missing["integralMapDigest"]] === dsKiraExpressionDigest[Lookup[manifest, "integralRules", {}]] &&
     Lookup[identity, "integralMapDigest", Missing["integralMapDigest"]] === dsKiraExpressionDigest[repJ2Kira] &&
     Lookup[identity, "targetDigest", Missing["targetDigest"]] === dsKiraExpressionDigest[Lookup[manifest, "targetIntegralIDs", {}]] &&
     Lookup[identity, "coefficientRulesDigest", Missing["coefficientRulesDigest"]] === dsKiraExpressionDigest[Lookup[manifest, "coefficientRulesApplied", {}]] &&
     Lookup[identity, "activeBasisDigest", Missing["activeBasisDigest"]] === dsKiraExpressionDigest[Lookup[manifest, "activeBasis", <||>]] &&
     dsKiraExportFileDigestsQ[identity, workspace]
    ]
   ];


(* ::Section::Closed:: *)
(*Active basis manifest 门禁*)

dsKiraActiveBasisData[manifest_Association] := Lookup[manifest, "activeBasis", <|"status" -> "disabled", "count" -> 0|>];

dsKiraUserMIDataQ[userData_Association, activeData_Association] := Module[
   {count, expressions, activeIndices, tokens, activeTokens, backendIDs, backendTokens, payload},
   count = Lookup[userData, "count", -1];
   expressions = Lookup[userData, "expressions", {}];
   activeIndices = Lookup[userData, "activeIndices", {}];
   tokens = Lookup[userData, "tokens", {}];
   activeTokens = Lookup[userData, "activeTokens", {}];
   backendIDs = Lookup[userData, "backendIDs", {}];
   backendTokens = Lookup[userData, "backendTokens", {}];
   payload = KeyDrop[userData, {
      "mappingDigest", "backendIDs", "backendTokens", "activeBackendIDs",
      "activeBackendTokens", "userMIToBackendRules", "backendToUserMIRules"
      }];
   Lookup[userData, "status", "failed"] === "configured" &&
    count === Lookup[activeData, "count", -2] &&
    expressions === Lookup[activeData, "expressions", Missing["expressions"]] &&
    activeIndices === Lookup[activeData, "activeIndices", Missing["activeIndices"]] &&
    tokens === (userMI /@ Range[count]) && activeTokens === tokens[[activeIndices]] &&
    backendIDs === Lookup[activeData, "ids", Missing["ids"]] &&
    backendTokens === (Tuserweight /@ backendIDs) &&
    Lookup[userData, "activeBackendIDs", {}] === Lookup[activeData, "activeIDs", Missing["activeIDs"]] &&
    Lookup[userData, "activeBackendTokens", {}] === (Tuserweight /@ Lookup[activeData, "activeIDs", {}]) &&
    Lookup[userData, "userMIToBackendRules", {}] === Thread[tokens -> backendTokens] &&
    Lookup[userData, "backendToUserMIRules", {}] === Thread[backendTokens -> tokens] &&
    TrueQ[Lookup[userData, "reversibleQ", False]] &&
    Lookup[userData, "rank", -1] === count &&
    Lookup[userData, "mappingDigest", Missing["mappingDigest"]] === dsKiraExpressionDigest[payload]
   ];


dsKiraUserMIDataQ[None, _Association] := True;
dsKiraUserMIDataQ[_, _Association] := False;


dsKiraActiveBasisDataQ[data_Association] := Module[
   {status, count, activeCount, names, expressions, ids, tokens, activeIndices, activeNames,
    activeExpressions, activeIDs, activeTokens, auxiliaryIDs, variables, targetIDs, userData},
   status = Lookup[data, "status", "disabled"];
   If[status === "disabled", Return[True]];
   count = Lookup[data, "count", -1];
   activeCount = Lookup[data, "activeCount", -1];
   names = Lookup[data, "names", {}];
   expressions = Lookup[data, "expressions", {}];
   ids = Lookup[data, "ids", {}];
   tokens = Lookup[data, "tokens", {}];
   activeIndices = Lookup[data, "activeIndices", {}];
   activeNames = Lookup[data, "activeNames", {}];
   activeExpressions = Lookup[data, "activeExpressions", {}];
   activeIDs = Lookup[data, "activeIDs", {}];
   activeTokens = Lookup[data, "activeTokens", {}];
   auxiliaryIDs = Lookup[data, "auxiliaryIDs", {}];
   variables = Lookup[data, "derivativeVariables", {}];
   targetIDs = Lookup[data, "targetIntegralIDs", {}];
   userData = Lookup[data, "userMI", None];
   status === "configured" && IntegerQ[count] && count > 0 && IntegerQ[activeCount] && activeCount > 0 &&
    Length[names] === count && Length[expressions] === count && Length[ids] === count && Length[tokens] === count &&
    And @@ (StringQ[#] && # =!= "" & /@ names) && DuplicateFreeQ[names] &&
    ids === Range[count] && tokens === (Tuserweight /@ ids) &&
    Length[activeIndices] === activeCount && DuplicateFreeQ[activeIndices] &&
    And @@ (IntegerQ[#] && 1 <= # <= count & /@ activeIndices) &&
    activeIDs === ids[[activeIndices]] && activeNames === names[[activeIndices]] &&
    activeExpressions === expressions[[activeIndices]] && activeTokens === (Tuserweight /@ activeIDs) &&
    auxiliaryIDs === Complement[ids, activeIDs] && variables =!= {} &&
    DuplicateFreeQ[targetIDs] && Complement[activeIDs, targetIDs] === {} &&
    dsKiraUserMIDataQ[userData, data]
   ];

dsKiraBackendMasterObject[id_Integer, activeIDs_List, idToJ_Association] := If[
   MemberQ[activeIDs, id],
   Tuserweight[id],
   Lookup[idToJ, id, Missing["unrecognizedBackendMasterID", id]]
   ];

DSKiraImport[root_String, context_: Automatic, OptionsPattern[]] := Module[
   {resolved, workspace, files, missingFiles, manifest, repJ2Kira, repKira2J, reductionRulesBackend, reductionRulesRaw, masterText,
    completionText, completionQ, requiredFiles, mapQ, manifestMapQ, contextQ, artifactIdentityQ, masterIDs, mapPairs, idToJ, mapIDs,
    dummyID, targetIDs, ruleIDs, completeTargetsQ, rhsMastersQ, coefficientVariables, allowedCoefficientVariables,
    coefficientQ, coefficientVariableMap, backendImaginaryUnit, backendVariableNames, allowedBackendVariableNames,
     backendKinematicConvention, backendKinematicConventionQ, backendKinematicRuleData,
     backendKinematicRuleDataQ, physicalCoefficientRulesQ,
    gaussianPhaseGauge, gaussianPhaseGaugeQ,
    backendCoefficientQ, activeData, activeDataQ, activeQ, relationIDs, activeIDs, auxiliaryIDs, activeExpressions, activeTokens, activeUserMITokens,
    activeMasterOrderQ, auxiliaryNotMastersQ, recognizedIDs, backendMasters, boundaryMasterIDs, boundaryMasters,
     activeBackendToUserMIRules, checks, diagnostics, issues, reductionRules, masters, masterTokens, returnedMasterIDs, progress = OptionValue[ProgressReporting]},
   resolved = dsResolveContext[context];
   If[Head[resolved] === Missing,
    Message[DSKiraImport::mismatch, "Kira import 缺少同源 DSInit context，无法验证结果归属；请先成功调用 DSInit，或把该 context 显式传给 DSKiraImport。 Kira import has no matching DSInit context, so result provenance cannot be verified; run DSInit successfully first or pass the context explicitly to DSKiraImport."]; dsErrorPrint["Kira import 需要同源 DSInit context。 Kira import requires the matching DSInit context."]; Return[<|"status" -> "failed", "reason" -> "missingContext"|>]
    ];
   workspace = ExpandFileName[root];
   If[! DirectoryQ[workspace],
    Message[DSKiraImport::badpath, workspace]; dsErrorPrint["Kira workspace 不存在。 The Kira workspace does not exist."]; Return[<|"status" -> "failed", "reason" -> "invalidWorkspace", "workspace" -> workspace|>]
    ];
   files = <|
     "manifest" -> FileNameJoin[{workspace, "dsibp-export-manifest.wl"}],
     "repJ2Kira" -> FileNameJoin[{workspace, "result", "repJ2kira.m"}],
     "repKira2J" -> FileNameJoin[{workspace, "result", "repkira2J.m"}],
      "reduction" -> dsKiraResolveFile[workspace, OptionValue[KiraReductionFile], {"results", "Tuserweight", "kira_list.m"}],
      "masters" -> dsKiraResolveFile[workspace, OptionValue[KiraMasterFile], {"results", "Tuserweight", "masters"}],
     "completion" -> dsKiraResolveFile[workspace, OptionValue[KiraCompletionFile], {"kira.log"}]
     |>;
   requiredFiles = KeyDrop[files, "completion"];
   missingFiles = Select[Normal[requiredFiles], ! StringQ[Last[#]] || ! FileExistsQ[Last[#]] &];
   If[missingFiles =!= {},
    Message[DSKiraImport::missing, missingFiles]; dsErrorPrint["完整 Kira 结果文件不足，未导入。 Required Kira result files are missing, so no result was imported."]; Return[<|"status" -> "failed", "reason" -> "missingFiles", "workspace" -> workspace, "files" -> files, "missingFiles" -> missingFiles|>]
    ];
   {manifest, repJ2Kira, repKira2J, reductionRulesBackend} = dsStageRun[
     "读取 Kira manifest、映射与 reduction / Reading the Kira manifest, maps, and reduction",
     dsKiraReadExpression /@ Lookup[files, {"manifest", "repJ2Kira", "repKira2J", "reduction"}],
     progress
     ];
   masterText = Import[files["masters"], "Text"];
   completionText = If[StringQ[files["completion"]] && FileExistsQ[files["completion"]],
     Import[files["completion"], "Text"], Missing["CompletionLogNotFound"]];
   completionQ = StringQ[completionText] && dsKiraCompletionQ[completionText, OptionValue[KiraCompletionPatterns]];
   If[! TrueQ[completionQ],
    dsWarningPrint["Kira 日志未确认成功完成；将以 reduction、targets、masters、映射和系数域的结构闭合作为硬边界。 The Kira log does not confirm successful completion; structural reduction, target, master, map, and coefficient-domain closure remains authoritative."]
    ];
   If[! AssociationQ[manifest] || ! dsRuleListQ[repJ2Kira] || ! dsRuleListQ[repKira2J] || ! dsRuleListQ[reductionRulesBackend],
    Message[DSKiraImport::invalid, "Kira manifest、映射或 reduction 文件不是预期的 Wolfram 表达式：manifest 必须是 Association，repJ2kira 与 repkira2J 必须是规则列表，reduction 必须是规则列表；请核对 workspace 文件后重新导入。 The Kira manifest, map, or reduction file is not the expected Wolfram expression: the manifest must be an Association and repJ2kira, repkira2J, and the reduction must be rule lists; check the workspace files and import again."]; dsErrorPrint["Kira 文件不是预期的 Wolfram 表达式。 The Kira files are not the expected Wolfram expressions."]; Return[<|"status" -> "failed", "reason" -> "malformedExpressions", "workspace" -> workspace|>]
    ];
   coefficientVariableMap = Lookup[manifest, "coefficientVariableMap", {}];
   backendImaginaryUnit = Lookup[manifest, "backendImaginaryUnit", None];
   reductionRulesRaw = dsKiraRestoreBackendCoefficients[
     reductionRulesBackend,
     coefficientVariableMap,
     backendImaginaryUnit
     ];
    backendKinematicConvention = Lookup[manifest, "backendKinematicConvention", <|"status" -> "notRequired"|>];
    backendKinematicConventionQ = AssociationQ[backendKinematicConvention] &&
      dsKiraBackendKinematicConventionDataQ[backendKinematicConvention];
    backendKinematicRuleData = Lookup[manifest, "backendKinematicRuleData", Missing["backendKinematicRuleData"]];
    backendKinematicRuleDataQ = AssociationQ[backendKinematicRuleData] &&
      dsKiraBackendKinematicRuleDataQ[backendKinematicRuleData, resolved["topology"], backendKinematicConvention];
    physicalCoefficientRulesQ = backendKinematicRuleDataQ &&
      Lookup[manifest, "physicalCoefficientRulesApplied", Missing["physicalCoefficientRulesApplied"]] ===
        Lookup[backendKinematicRuleData, "physicalPointRules", Missing["physicalPointRules"]];
    If[backendKinematicConventionQ,
     reductionRulesRaw = dsKiraRestorePhysicalKinematicVariables[reductionRulesRaw, backendKinematicConvention]
    ];
   gaussianPhaseGauge = Lookup[manifest, "gaussianPhaseGauge", <|"status" -> "notApplicable"|>];
   gaussianPhaseGaugeQ = AssociationQ[gaussianPhaseGauge] &&
     dsKiraGaussianPhaseGaugeDataQ[gaussianPhaseGauge, Lookup[manifest, "integralCount", -1]];
   If[gaussianPhaseGaugeQ,
    reductionRulesRaw = dsKiraRestorePhysicalPhaseGauge[reductionRulesRaw, gaussianPhaseGauge]
    ];
   mapQ = dsKiraInverseMapQ[repJ2Kira, repKira2J];
   manifestMapQ = SortBy[dsJToIDPairs[Lookup[manifest, "integralRules", {}]], Last] === SortBy[dsJToIDPairs[repJ2Kira], Last];
   contextQ = dsKiraContextMatchQ[manifest, resolved];
   artifactIdentityQ = dsKiraArtifactIdentityQ[manifest, workspace, repJ2Kira];
   masterIDs = dsKiraMasterIDs[masterText];
   mapPairs = dsIDToJPairs[repKira2J];
   idToJ = Association[Rule @@@ mapPairs];
   mapIDs = Keys[idToJ];
   activeData = dsKiraActiveBasisData[manifest];
   activeDataQ = AssociationQ[activeData] && dsKiraActiveBasisDataQ[activeData];
   activeQ = TrueQ[activeDataQ] && Lookup[activeData, "status", "disabled"] === "configured";
   relationIDs = If[activeQ, Lookup[activeData, "ids", {}], {}];
   activeIDs = If[activeQ, Lookup[activeData, "activeIDs", {}], {}];
   auxiliaryIDs = If[activeQ, Lookup[activeData, "auxiliaryIDs", {}], {}];
    activeExpressions = If[activeQ, Lookup[activeData, "activeExpressions", {}], {}];
    activeTokens = If[activeQ, Lookup[activeData, "activeTokens", {}], {}];
    activeBackendToUserMIRules = If[
      activeQ,
      Lookup[Lookup[activeData, "userMI", <||>], "backendToUserMIRules", {}],
      {}
      ];
   activeUserMITokens = If[activeQ,
     Lookup[Lookup[activeData, "userMI", <||>], "activeTokens", activeTokens],
     {}
     ];
   activeMasterOrderQ = ! activeQ || Select[masterIDs, MemberQ[activeIDs, #] &] === activeIDs;
   auxiliaryNotMastersQ = ! activeQ || Intersection[auxiliaryIDs, masterIDs] === {};
   recognizedIDs = Join[mapIDs, relationIDs];
   dummyID = Lookup[manifest, "numericDummyIntegralId", None];
   targetIDs = DeleteCases[Lookup[manifest, "targetIntegralIDs", mapIDs], dummyID];
   ruleIDs = dsKiraRuleIDs[reductionRulesBackend];
   completeTargetsQ = Complement[targetIDs, Union[ruleIDs["lhs"], masterIDs]] === {};
   rhsMastersQ = Complement[ruleIDs["rhs"], masterIDs] === {};
   coefficientVariables = dsKiraCoefficientVariables[reductionRulesRaw];
   allowedCoefficientVariables = Lookup[manifest, "coefficientVariables", {}];
   coefficientQ = Complement[coefficientVariables, allowedCoefficientVariables] === {};
   backendVariableNames = dsKiraBackendVariableNames[reductionRulesBackend];
   allowedBackendVariableNames = Lookup[manifest, "backendCoefficientVariables", {}];
   backendCoefficientQ = If[coefficientVariableMap === {} && allowedBackendVariableNames === {},
     True,
     Complement[backendVariableNames, allowedBackendVariableNames] === {}
     ];
   diagnostics = <|"completionMarker" -> completionQ|>;
   checks = <|
     "inverseIntegralMaps" -> mapQ,
     "manifestIntegralMap" -> manifestMapQ,
     "contextConventionMatch" -> contextQ,
     "exportArtifactIdentity" -> artifactIdentityQ,
     "activeBasisManifest" -> activeDataQ,
     "nonemptyMasterOrder" -> (masterIDs =!= {}),
     "masterIDsRecognized" -> (Complement[masterIDs, recognizedIDs] === {}),
     "activeBasisIDsAreMasters" -> (! activeQ || Complement[activeIDs, masterIDs] === {}),
     "activeBasisMasterOrder" -> activeMasterOrderQ,
     "auxiliaryBasisIDsNotMasters" -> auxiliaryNotMastersQ,
     "allReductionIDsRecognized" -> (Complement[ruleIDs["all"], Append[recognizedIDs, dummyID]] === {}),
     "completeTargetCoverage" -> completeTargetsQ,
     "rhsContainsOnlyMasters" -> rhsMastersQ,
      "backendCoefficientVariablesRecognized" -> backendCoefficientQ,
      "backendKinematicConventionManifest" -> backendKinematicConventionQ,
      "backendKinematicRuleDataManifest" -> backendKinematicRuleDataQ,
      "physicalCoefficientRulesConsistent" -> physicalCoefficientRulesQ,
      "physicalKinematicVariablesRestored" -> backendKinematicConventionQ,
     "gaussianPhaseGaugeManifest" -> gaussianPhaseGaugeQ,
     "physicalIntegralPhaseRestored" -> gaussianPhaseGaugeQ,
     "coefficientVariablesRecognized" -> coefficientQ
     |>;
   issues = Keys @ Select[checks, ! TrueQ[#] &];
   If[issues =!= {},
    Message[DSKiraImport::mismatch, StringRiffle[dsKiraCheckSentence /@ issues, "；"]]; dsErrorPrint["Kira 结果未通过同源性/完整性门禁。 The Kira results failed the provenance or completeness gate."]; Return[<|"status" -> "failed", "reason" -> "validationFailed", "workspace" -> workspace, "files" -> files, "validationReport" -> <|"checks" -> checks, "issues" -> issues|>|>]
    ];
   backendMasters = dsKiraBackendMasterObject[#, relationIDs, idToJ] & /@ masterIDs;
   boundaryMasterIDs = Complement[masterIDs, relationIDs];
   boundaryMasters = Lookup[idToJ, boundaryMasterIDs, {}];
   masters = If[activeQ, activeExpressions, backendMasters];
   masterTokens = If[activeQ, activeUserMITokens, masters];
   returnedMasterIDs = If[activeQ, activeIDs, masterIDs];
    reductionRules = reductionRulesRaw /. Normal[Association[repKira2J]];
    If[activeQ,
     reductionRules = reductionRules /. activeBackendToUserMIRules
     ];
    If[Cases[reductionRules, Tuserweight[_Integer], Infinity] =!= {},
     Message[DSKiraImport::invalid, "约化右端在 userMI 恢复之后仍残留 active basis 的 backend token；请使用同一次导出的文件重新运行，或核对 active basis 与 userMI 映射。 Reduction right-hand sides still contain active-basis backend tokens after userMI restoration; rerun with files from the same export, or verify the active-basis and userMI maps."];
     dsErrorPrint["Active basis 的 Kira token 未能恢复为 userMI。 Active-basis Kira tokens were not restored to userMI."];
     Return[<|"status" -> "failed", "reason" -> "activeBasisTokenRestorationFailed", "workspace" -> workspace|>]
     ];
   <|
    "status" -> "imported",
    "workspace" -> workspace,
    "files" -> files,
    "reductionRules" -> reductionRules,
    "backendReductionRules" -> reductionRulesBackend,
    "masters" -> masters,
    "masterTokens" -> masterTokens,
    "backendMasterTokens" -> If[activeQ, activeTokens, backendMasters],
    "masterIDs" -> returnedMasterIDs,
    "backendMasters" -> backendMasters,
    "backendMasterIDs" -> masterIDs,
    "boundaryMasters" -> boundaryMasters,
    "boundaryMasterIDs" -> boundaryMasterIDs,
    "activeBasis" -> activeData,
    "integralMap" -> <|"JToKira" -> repJ2Kira, "KiraToJ" -> repKira2J|>,
    "coefficientVariables" -> coefficientVariables,
    "backendCoefficientVariables" -> backendVariableNames,
    "coefficientVariableMap" -> coefficientVariableMap,
    "backendKinematicConvention" -> backendKinematicConvention,
    "gaussianPhaseGauge" -> gaussianPhaseGauge,
    "sourceManifest" -> manifest,
    "context" -> resolved,
    "validationReport" -> <|"status" -> "passed", "checks" -> checks,
      "diagnostics" -> diagnostics,
      "warnings" -> If[completionQ, {}, {"completionMarkerMissing"}], "issues" -> {}|>
    |>
   ];

DSKiraImport[root_, context_: Automatic, OptionsPattern[]] := (Message[DSKiraImport::badpath, root]; dsErrorPrint["DSKiraImport 的第一个参数必须是目录字符串。 The first DSKiraImport argument must be a directory string."]; <|"status" -> "failed", "reason" -> "workspaceNotString"|>);
