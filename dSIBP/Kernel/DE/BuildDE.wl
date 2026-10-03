(* ::Package:: *)

(* ::Chapter:: *)
(*微分方程构造*)

(* DSDE 只消费经 KiraImport 验证的 reduction data；不会从不完整日志猜测 master 或规则。 *)

Options[DSDE] = {
   OutputDirectory -> None,
   ProgressReporting -> Automatic
   };

DSDE::badreduction = "DSDE 只接受 DSKiraImport 验证通过的 reductionData。";
DSDE::badvars = "微分变量必须是当前 family 初始化的外部独立变量：`1`。";
DSDE::writefailed = "DE 结果写入失败：`1`。";


(* ::Section::Closed:: *)
(*Kira 内部虚轴运动学坐标的微分接口*)

(* DSDE 的公开矩阵仍表示物理 D_k。这里同时给出直接 backend 坐标的
   D_ik=-I D_k 结果，并保存 D_k=I D_ik 与 Euler 不变量供 scaling 审计。 *)
dsDEBackendKinematicDerivativeView[matrices_Association, sources_Association, variables_List, convention_Association] := Module[
   {records, activeRecords, backendVariables, backendMatrices, backendSources},
   If[Lookup[convention, "status", "notRequired"] =!= "configured",
    Return[<|"status" -> "notRequired", "scope" -> "KiraBackendOnly"|>]
    ];
   records = Lookup[convention, "records", {}];
   activeRecords = Select[records, MemberQ[variables, Lookup[#, "physical"]] &];
   backendVariables = Lookup[activeRecords, "backend", {}];
   backendMatrices = AssociationThread[
     backendVariables,
     (-I Lookup[matrices, Lookup[#, "physical"]]) & /@ activeRecords
     ];
   backendSources = AssociationThread[
     backendVariables,
     (-I Lookup[sources, Lookup[#, "physical"]]) & /@ activeRecords
     ];
   <|
    "status" -> If[activeRecords === {}, "notRequired", "generated"],
    "scope" -> "KiraBackendOnly",
    "backendVariables" -> backendVariables,
    "matrices" -> backendMatrices,
    "sources" -> backendSources,
    "backendFromPhysicalDerivativeFactor" -> -I,
    "physicalFromBackendDerivativeFactor" -> I,
     "ordinaryDerivativeConvention" -> "D[physicalKinematic] == I D[backendKinematic]",
     "eulerConvention" -> "physicalKinematic D[physicalKinematic] == backendKinematic D[backendKinematic]",
    "records" -> activeRecords
    |>
   ];

(* DSKiraImport 已验证每条 reduction 的右端只含 masters，因此一次替换就是完整约化；
   若外部结果违反该合同，后续 residual gate 会拒绝 DE，而不是用任意迭代次数掩盖问题。 *)
dsReduceExpression[expr_, rules_List] := Expand[expr /. rules];

(* Kira 关系可以使用内部 kk/ISP 坐标；公开 DE 必须只含 family 声明的外部不变量。 *)
dsDEReducedExpressionToUser[expr_, context_Association] := Module[{topo = context["topology"]},
   Expand[scalarProductInternalToUser[expr /. internalISPToUserRules[topo], topo]]
   ];

dsDEMasterDecomposition[expr_, masterTokens_List] := Module[
   {coefficientTokens, tokenExpr, coefficients, source, residualIntegrals, residualBackendTokens, residualObjects},
   coefficientTokens = Array[Unique["dsMaster$"] &, Length[masterTokens]];
   tokenExpr = expr /. Thread[masterTokens -> coefficientTokens];
   coefficients = Coefficient[tokenExpr, #] & /@ coefficientTokens;
   source = Expand[tokenExpr - coefficients.coefficientTokens];
   residualIntegrals = DeleteDuplicates[Cases[source, _J, Infinity]];
   residualBackendTokens = DeleteDuplicates[Cases[source, Tuserweight[_Integer], Infinity]];
   residualObjects = Join[residualIntegrals, residualBackendTokens];
   <|
    "coefficients" -> coefficients,
    "source" -> source,
    "residualIntegrals" -> residualIntegrals,
    "residualBackendTokens" -> residualBackendTokens,
    "residualObjects" -> residualObjects,
    "closedQ" -> (residualObjects === {})
    |>
   ];

(* active basis 可同时含 top 与 lower-sector J。公开 ds 已按每个 J 的 shrunk pack 选择一次
   sector topology，并同时处理显式系数、sector prefactor 与 canonical；这里不得预先再缩一次。 *)
dsSectorAwareDerivative[expr_, variable_, context_Association] := Module[
   {rootTopology, result},
   rootTopology = Lookup[context, "topology", Missing["NoRootTopology"]];
   If[! parsedTopologyQ[rootTopology], Return[$Failed]];
   result = Check[ds[expr, variable, rootTopology], $Failed];
   If[result === $Failed, $Failed, Expand[result]]
   ];

dsDEVariableData[variable_, masterDefinitions_List, masterTokens_List, rules_List, parameterRules_List, context_Association, progress_] := Module[
   {raw, reduced, decompositions},
   raw = dsProgressMap[
     "正在构造 " <> ToString[variable, InputForm] <> " 导数 / Building " <> ToString[variable, InputForm] <> " derivatives",
     masterDefinitions,
     Function[master, dsSectorAwareDerivative[master, variable, context] /. parameterRules],
     progress
     ];
   If[MemberQ[raw, $Failed], Return[<|"status" -> "failed", "variable" -> variable, "reason" -> "dsFailed", "rawDerivatives" -> raw|>]];
   reduced = dsProgressMap[
     "正在约化 " <> ToString[variable, InputForm] <> " 导数 / Reducing " <> ToString[variable, InputForm] <> " derivatives",
     raw,
     Function[expr, dsReduceExpression[expr, rules]],
     progress
     ];
   reduced = dsDEReducedExpressionToUser[#, context] & /@ reduced;
   decompositions = dsDEMasterDecomposition[#, masterTokens] & /@ reduced;
   <|
    "status" -> If[And @@ Lookup[decompositions, "closedQ", False], "generated", "notClosed"],
    "variable" -> variable,
    "matrix" -> Lookup[decompositions, "coefficients", {}],
    "source" -> Lookup[decompositions, "source", {}],
    "rawDerivatives" -> raw,
    "reducedDerivatives" -> reduced,
    "residualIntegrals" -> DeleteDuplicates[Flatten[Lookup[decompositions, "residualIntegrals", {}]]],
    "residualBackendTokens" -> DeleteDuplicates[Flatten[Lookup[decompositions, "residualBackendTokens", {}]]],
    "residualObjects" -> DeleteDuplicates[Flatten[Lookup[decompositions, "residualObjects", {}]]]
    |>
   ];

dsWriteDEResult[data_Association, directory_String] := Module[{paths, compact},
   Quiet[CreateDirectory[directory, CreateIntermediateDirectories -> True]];
   paths = <|
     "masters" -> FileNameJoin[{directory, "masters.wl"}],
     "de" -> FileNameJoin[{directory, "de.wl"}],
     "manifest" -> FileNameJoin[{directory, "manifest.wl"}]
     |>;
   compact = KeyDrop[data, {"variableData", "context"}];
   If[Quiet[Check[
       Put[data["masters"], paths["masters"]];
       Put[KeyTake[data, {"status", "variables", "matrices", "sources", "residualIntegrals", "equationConvention"}], paths["de"]];
       Put[Join[compact, <|"files" -> AssociationMap[FileNameTake, paths]|>], paths["manifest"]];
       True,
       False
       ]],
    <|"status" -> "written", "directory" -> directory, "files" -> paths|>,
    <|"status" -> "failed", "directory" -> directory|>
    ]
   ];

DSDE[reductionData_Association, variables_: Automatic, OptionsPattern[]] := Module[
   {context, masters, masterTokens, rules, resolvedVariables, allowedVariables, badVariables,
    sourceManifest, kiraPlan, postDerivativeRules, physicalPostDerivativeRules, parameterRules, variableRecords,
    variableData, status, matrices, sources, backendKinematicConvention, backendDerivativeView,
    result, outputDirectory = OptionValue[OutputDirectory], writeResult, notClosedVariables, residualObjectCount},
   If[Lookup[reductionData, "status", "missing"] =!= "imported" ||
     Lookup[Lookup[reductionData, "validationReport", <||>], "status", "missing"] =!= "passed",
    Message[DSDE::badreduction]; dsErrorPrint["reductionData 未经 DSKiraImport 完整验证。 reductionData has not passed complete DSKiraImport validation."]; Return[<|"status" -> "failed", "reason" -> "unvalidatedReductionData"|>]
    ];
   context = Lookup[reductionData, "context", Missing["context"]];
   If[! dsContextQ[context], Message[DSDE::badreduction]; Return[<|"status" -> "failed", "reason" -> "missingContext"|>]];
   If[! dsContextCapabilityQ[context, "derivativeUsableQ"],
    Message[DSDE::badreduction];
    dsErrorPrint["当前参数声明不支持唯一微分算符。 The current parameter declaration does not define unique differential operators."]; Return[<|
      "status" -> "failed", "reason" -> "derivativeCapabilityGate",
      "capabilities" -> dsContextCapabilities[context]
      |>]
    ];
   masters = reductionData["masters"];
   masterTokens = Lookup[reductionData, "masterTokens", masters];
   If[! ListQ[masters] || ! ListQ[masterTokens] || Length[masters] =!= Length[masterTokens] || masters === {},
    Message[DSDE::badreduction]; Return[<|"status" -> "failed", "reason" -> "invalidMasterDefinitionsOrTokens"|>]
    ];
   rules = reductionData["reductionRules"];
   sourceManifest = Lookup[reductionData, "sourceManifest", <||>];
   kiraPlan = Lookup[sourceManifest, "kiraPlan", <||>];
   postDerivativeRules = If[
     Lookup[kiraPlan, "numericStage", "symbolic"] === "postDerivative",
     Lookup[kiraPlan, "coefficientRules", {}],
     {}
     ];
   physicalPostDerivativeRules = If[postDerivativeRules === {},
     {},
     Lookup[sourceManifest, "physicalCoefficientRulesApplied", postDerivativeRules]
     ];
   (* fixed-rational export 在解析导数 closure 冻结后才数值化；DSDE 必须把同一规则
      同时用于内部原子和公开坐标，否则 h EOM 或系数乘积法则会重新引入参数。 *)
   parameterRules = DeleteDuplicatesBy[
     Join[
      If[ListQ[physicalPostDerivativeRules],
       normalizeCoefficientRulesForTopology[physicalPostDerivativeRules, context["topology"]], {}],
      If[ListQ[physicalPostDerivativeRules], physicalPostDerivativeRules, {}]
      ],
     First
     ];
   resolvedVariables = dsDEResolveVariables[variables, context];
   allowedVariables = dsDEResolveVariables[Automatic, context];
   badVariables = Complement[resolvedVariables, allowedVariables];
   If[badVariables =!= {},
    Message[DSDE::badvars, badVariables]; dsErrorPrint["DSDE 变量不属于当前 family 的外部表示。 The DSDE variables are not external coordinates of the current family."]; Return[<|"status" -> "failed", "reason" -> "invalidVariables", "badVariables" -> badVariables, "allowedVariables" -> allowedVariables|>]
    ];
   variableRecords = dsProgressMap[
     "正在生成微分方程 / Building differential equations",
     resolvedVariables,
     Function[variable, dsDEVariableData[variable, masters, masterTokens, rules, parameterRules, context, OptionValue[ProgressReporting]]],
     OptionValue[ProgressReporting]
     ];
   variableData = AssociationThread[resolvedVariables, variableRecords];
   status = Which[
     AnyTrue[variableRecords, Lookup[#, "status", "failed"] === "failed" &], "failed",
     AnyTrue[variableRecords, Lookup[#, "status", "notClosed"] === "notClosed" &], "notClosed",
     True, "generated"
     ];
    matrices = AssociationThread[resolvedVariables, Lookup[variableRecords, "matrix", {}]];
    sources = AssociationThread[resolvedVariables, Lookup[variableRecords, "source", {}]];
    If[status === "notClosed",
     notClosedVariables = Keys@Select[variableData, Lookup[#, "status", "failed"] === "notClosed" &];
     residualObjectCount = Total[Length /@ Lookup[variableData, notClosedVariables, <||>][[All, "residualObjects"]]];
     Message[DSDE::notclosed, notClosedVariables, residualObjectCount];
     dsErrorPrint[
      "部分导数仍含未约化积分或后端主积分标记，因此没有生成可用的闭合微分方程。请为返回结果中 residualObjects 列出的对象补齐约化规则。 " <>
       "Some derivatives still contain unreduced integrals or backend master tokens, so no usable closed differential equation was generated. Add reduction rules for the objects listed in residualObjects."
      ]
     ];
    backendKinematicConvention = Lookup[sourceManifest, "backendKinematicConvention", <|"status" -> "notRequired"|>];
   backendDerivativeView = dsDEBackendKinematicDerivativeView[
     matrices,
     sources,
     resolvedVariables,
     backendKinematicConvention
     ];
   result = <|
     "status" -> status,
     "masters" -> masters,
     "masterTokens" -> masterTokens,
     "masterCount" -> Length[masters],
     "variables" -> resolvedVariables,
     "matrices" -> matrices,
     "sources" -> sources,
     "backendKinematicDerivativeView" -> backendDerivativeView,
     "residualIntegrals" -> AssociationThread[resolvedVariables, Lookup[variableRecords, "residualIntegrals", {}]],
     "residualBackendTokens" -> AssociationThread[resolvedVariables, Lookup[variableRecords, "residualBackendTokens", {}]],
     "residualObjects" -> AssociationThread[resolvedVariables, Lookup[variableRecords, "residualObjects", {}]],
     "variableData" -> variableData,
     "sourceManifest" -> reductionData["sourceManifest"],
     "activeBasis" -> Lookup[reductionData, "activeBasis", <|"status" -> "disabled", "count" -> 0|>],
     "parameterRulesApplied" -> parameterRules,
     "postDerivativeRulesApplied" -> postDerivativeRules,
     "physicalPostDerivativeRulesApplied" -> physicalPostDerivativeRules,
     "context" -> context,
     "equationConvention" -> "D[masters,var] == matrices[var].masters + sources[var]",
     "reductionValidationReport" -> reductionData["validationReport"]
     |>;
   writeResult = If[StringQ[outputDirectory], dsWriteDEResult[result, ExpandFileName[outputDirectory]], <|"status" -> "notRequested"|>];
   If[Lookup[writeResult, "status", "failed"] === "failed", Message[DSDE::writefailed, outputDirectory]; dsErrorPrint["DE 文件未写出。 DE files were not written."]];
   Join[result, <|"writeResult" -> writeResult|>]
   ];

DSDE[reductionData_, variables_: Automatic, OptionsPattern[]] := (Message[DSDE::badreduction]; dsErrorPrint["DSDE 输入必须是 reductionData Association。 DSDE input must be a reductionData Association."]; <|"status" -> "failed", "reason" -> "inputNotAssociation"|>);


(* ::Chapter:: *)
(*Naive tree IBP 微分方程*)

(* tree 原生导数直接作用于 Exp[I k0_v tau_v] 与 h building block。
   tree J 使用统一的 (-tau)^a convention，故顶点相位贡献为
   +I D[k0_v,x] J[...,a_v+1,...]；contact 合并不再另乘 (-1)^Delta a。
   sector master 的 normalization N 最后另按乘积法则求导。 *)

Options[DSTreeNaiveDE] = {AuditLevel -> "standard", ProgressReporting -> Automatic};

DSTreeNaiveDE::badibp = "DSTreeNaiveDE 需要 DSTreeNaiveIBP 成功返回的数据或合法 DSInit context。";
DSTreeNaiveDE::badvars = "tree 微分变量必须是当前 family 初始化的外部独立变量：`1`。";


dsTreeZeroTokenTerms[expr_] := If[TrueQ[Expand[expr] === 0], {}, dsTreeTokenTerms[expr]];


dsTreeLineMomentumMagnitudeDerivative[int_J, variable_, family_Association] := Module[
   {packs = First[int], terms = {}, vertex, leg, magnitudeDerivative, state, newPacks, shiftedIntegral},
   Do[
    vertex = family["vertices"][[vertexIndex]];
    Do[
     leg = vertex["massiveLegs"][[legIndex]];
     magnitudeDerivative = D[leg["momentumMagnitude"], variable];
     If[! TrueQ[magnitudeDerivative === 0],
      state = packs[[vertexIndex, 1 + legIndex]];
      If[! MemberQ[{0, 1}, state], Return[$Failed]];
      newPacks = ReplacePart[packs, {vertexIndex, 1} -> packs[[vertexIndex, 1]] + 1];
      newPacks = ReplacePart[newPacks, {vertexIndex, 1 + legIndex} -> 1 - state];
      shiftedIntegral = J[newPacks];
      AppendTo[terms,
        magnitudeDerivative If[
         state === 0,
         -dsTreeToken[family["sector"], shiftedIntegral],
         dsTreeToken[family["sector"], shiftedIntegral] -
          (2 leg["nu"] + 1)/leg["momentumMagnitude"] dsTreeToken[family["sector"], int]
         ]
       ]
      ],
     {legIndex, Length[vertex["massiveLegs"]]}
     ],
    {vertexIndex, Length[family["vertices"]]}
    ];
   Expand[Total[terms]]
   ];


(* loop 投影只保留为单向 oracle，不参与生产 DE。 *)
dsTreePhaseDerivativeProjectionOracle[loopIntegral_J, variable_, family_Association, rootTopology_Association] := Module[
   {internalVariable, loopDerivative, projectedData, expression},
   internalVariable = scalarProductInputToInternal[variable, family["topology"]];
   loopDerivative = directExternalLegEnergyVariableDerivativeSeed[family["topology"], loopIntegral, internalVariable];
   If[TrueQ[Expand[loopDerivative] === 0],
    Return[<|"status" -> "generated", "loopDerivative" -> 0,
      "projectedData" -> <|"status" -> "generated", "terms" -> {}, "termCount" -> 0|>,
      "internalExpression" -> 0|>]
    ];
   projectedData = dsTreeLinearData[<|"loopSeed" -> loopDerivative|>, rootTopology, loopIntegral];
   If[Lookup[projectedData, "status", "failed"] =!= "generated",
    Return[<|"status" -> "failed", "reason" -> "phaseDerivativeProjectionFailed",
      "loopDerivative" -> loopDerivative, "projectedData" -> projectedData|>]
    ];
   expression = scalarProductInternalToUser[dsTreeTokenExpression[projectedData], rootTopology];
   <|"status" -> "generated", "loopDerivative" -> loopDerivative,
    "projectedData" -> projectedData, "internalExpression" -> expression|>
   ];


dsTreePhaseDerivativeDirect[int_J, variable_, family_Association, rootTopology_Association] := Module[
   {packs = First[int], terms, publicEnergy, derivative, shiftedPacks, shiftedIntegral},
   terms = Table[
     publicEnergy = scalarProductInternalToUser[family["vertices"][[vertexIndex, "signedExternalLegEnergy"]], rootTopology];
     derivative = D[publicEnergy, variable];
     If[TrueQ[derivative === 0],
      0,
      shiftedPacks = ReplacePart[packs, {vertexIndex, 1} -> packs[[vertexIndex, 1]] + 1];
      shiftedIntegral = J[shiftedPacks];
      I derivative dsTreeToken[family["sector"], shiftedIntegral]
      ],
     {vertexIndex, Length[family["vertices"]]}
     ];
   <|
    "status" -> "generated",
    "generationRoute" -> "directTreePhase",
    "terms" -> dsTreeZeroTokenTerms[Expand[Total[terms]]],
    "internalExpression" -> Expand[Total[terms]]
    |>
   ];


dsTreeNaiveAllowedVariables[context_Association, familyContext_Association] := DeleteDuplicates@Join[
   dsDEResolveVariables[Automatic, context],
   Variables[Cases[
     familyContext["families"],
     leg_Association /; KeyExistsQ[leg, "nu"] && KeyExistsQ[leg, "momentumMagnitude"] :> leg["momentumMagnitude"],
     Infinity
     ]]
   ];


dsTreeNaiveMasterDerivative[master_Association, variable_, familyContext_Association, context_Association] := Module[
   {family, rootTopology, phaseData, lineDerivative, bareToken, bareDerivative,
    physicalRecord, physicalLogDerivative, normalizedDerivative, publicTerms},
   family = dsTreeFamilyBySector[master["sectorKey"], familyContext];
   If[Head[family] === Missing, Return[<|"status" -> "failed", "reason" -> "unknownSector"|>]];
   rootTopology = context["topology"];
   phaseData = dsTreePhaseDerivativeDirect[master["integral"], variable, family, rootTopology];
   If[Lookup[phaseData, "status", "failed"] =!= "generated", Return[phaseData]];
   lineDerivative = dsTreeLineMomentumMagnitudeDerivative[master["integral"], variable, family];
   If[lineDerivative === $Failed, Return[<|"status" -> "failed", "reason" -> "lineMomentumMagnitudeDerivativeFailed"|>]];
   physicalRecord = sectorPrefactorRecordForKey018[context, master["sectorKey"]];
   If[! AssociationQ[physicalRecord],
    Return[<|"status" -> "failed", "reason" -> "physicalSectorPrefactorMissing",
      "sectorKey" -> master["sectorKey"], "record" -> physicalRecord|>]
    ];
   physicalLogDerivative = sectorPrefactorLogDerivative018[
     <|"sectorPrefactorData" -> physicalRecord["physicalSectorPrefactorData"]|>,
     variable
     ];
   If[physicalLogDerivative === $Failed,
    Return[<|"status" -> "failed", "reason" -> "physicalSectorPrefactorDerivativeFailed",
      "sectorKey" -> master["sectorKey"], "variable" -> variable|>]
    ];
   bareToken = dsTreeToken[master["sectorKey"], master["integral"]];
   bareDerivative = Expand[
     phaseData["internalExpression"] + lineDerivative + physicalLogDerivative bareToken
     ];
   normalizedDerivative = Expand[
      D[master["coefficient"], variable] bareToken + master["coefficient"] bareDerivative
     ];
   publicTerms = dsTreeZeroTokenTerms[normalizedDerivative];
   <|
    "status" -> If[publicTerms === $Failed, "failed", "generated"],
    "master" -> master,
    "variable" -> variable,
    "phaseDerivativeRoute" -> phaseData["generationRoute"],
    "directPhaseDerivativeTerms" -> phaseData["terms"],
    "lineMomentumMagnitudeDerivativeTerms" -> dsTreeZeroTokenTerms[lineDerivative],
    "selectorCoefficientDerivative" -> D[master["coefficient"], variable],
    "physicalSectorPrefactorLogDerivative" -> physicalLogDerivative,
    "completeNormalizationLogDerivative" -> Cancel[
      D[master["coefficient"], variable]/master["coefficient"] + physicalLogDerivative
      ],
    "rawTerms" -> publicTerms,
    "internalExpression" -> normalizedDerivative
    |>
   ];


dsTreeNaiveVariableData[variable_, ibpData_Association, familyContext_Association, context_Association, progress_] := Module[
   {masters, derivativeRecords, rules, reduced, masterTokens, coefficientTokens, normalizedMasterRules,
    tokenExpressions, coefficients, residuals, residualTokens, rows},
   masters = ibpData["masters"];
   derivativeRecords = dsProgressMap[
     "正在构造 naive tree " <> ToString[variable, InputForm] <> " 导数 / Building naive-tree " <> ToString[variable, InputForm] <> " derivatives",
     masters,
     Function[master, dsTreeNaiveMasterDerivative[master, variable, familyContext, context]],
     progress
     ];
   If[AnyTrue[derivativeRecords, Lookup[#, "status", "failed"] =!= "generated" &],
    Return[<|"status" -> "failed", "variable" -> variable, "reason" -> "derivativeGenerationFailed",
      "derivativeRecords" -> derivativeRecords|>]
    ];
   rules = dsTreeInternalReductionRules[ibpData];
   reduced = Expand[Lookup[derivativeRecords, "internalExpression"] /. rules];
   masterTokens = dsTreeToken[Lookup[#, "sectorKey"], Lookup[#, "integral"]] & /@ masters;
   coefficientTokens = Array[Unique["dsTreeMaster$"] &, Length[masters]];
   normalizedMasterRules = MapThread[#1 -> #2/#3 &,
     {masterTokens, coefficientTokens, Lookup[masters, "coefficient"]}];
   tokenExpressions = Expand[reduced /. normalizedMasterRules];
   coefficients = Table[Coefficient[tokenExpressions[[row]], coefficientTokens[[column]]],
     {row, Length[masters]}, {column, Length[masters]}];
   residuals = Expand[tokenExpressions - coefficients . coefficientTokens];
   residualTokens = DeleteDuplicates[Cases[residuals, _dsTreeToken, Infinity]];
   rows = MapThread[
     <|"master" -> #1, "coefficients" -> #2, "source" -> #3|> &,
     {masters, coefficients, residuals}
     ];
   <|
    "status" -> If[residualTokens === {} && And @@ (TrueQ[# === 0] & /@ residuals), "generated", "notClosed"],
    "variable" -> variable,
    "matrix" -> coefficients,
    "source" -> residuals,
    "rows" -> rows,
    "derivativeRecords" -> (KeyDrop[#, "internalExpression"] & /@ derivativeRecords),
    "residualObjects" -> residualTokens
    |>
   ];


dsTreeNaiveDEFromIBP[ibpData_Association, variables_, OptionsPattern[DSTreeNaiveDE]] := Module[
   {context, familyContext, resolvedVariables, allowedVariables, badVariables, variableRecords, variableData, status},
   If[Lookup[ibpData, "status", "failed"] =!= "solved" || ! ListQ[Lookup[ibpData, "masters", None]],
    Message[DSTreeNaiveDE::badibp]; Return[<|"status" -> "failed", "reason" -> "invalidNaiveIBPData"|>]
    ];
   context = Lookup[ibpData, "context", Missing["context"]];
   If[! dsContextQ[context], Message[DSTreeNaiveDE::badibp]; Return[<|"status" -> "failed", "reason" -> "missingContext"|>]];
   familyContext = dsTreeFamilyContext[context];
   allowedVariables = dsTreeNaiveAllowedVariables[context, familyContext];
   resolvedVariables = If[variables === Automatic, allowedVariables, dsDEResolveVariables[variables, context]];
   badVariables = Complement[resolvedVariables, allowedVariables];
   If[badVariables =!= {},
    Message[DSTreeNaiveDE::badvars, badVariables];
    Return[<|"status" -> "failed", "reason" -> "invalidVariables", "badVariables" -> badVariables,
      "allowedVariables" -> allowedVariables|>]
    ];
   variableRecords = dsProgressMap[
     "正在生成 naive tree 微分方程 / Building naive-tree differential equations",
     resolvedVariables,
     Function[variable, dsTreeNaiveVariableData[
       variable, ibpData, familyContext, context, OptionValue[ProgressReporting]
       ]],
     OptionValue[ProgressReporting]
     ];
   variableData = AssociationThread[resolvedVariables, variableRecords];
   status = Which[
     AnyTrue[variableRecords, Lookup[#, "status", "failed"] === "failed" &], "failed",
     AnyTrue[variableRecords, Lookup[#, "status", "failed"] === "notClosed" &], "notClosed",
     True, "generated"
     ];
   <|
    "status" -> status,
    "masters" -> ibpData["masters"],
    "masterCount" -> Length[ibpData["masters"]],
    "variables" -> resolvedVariables,
    "matrices" -> AssociationThread[resolvedVariables, Lookup[variableRecords, "matrix", {}]],
    "sources" -> AssociationThread[resolvedVariables, Lookup[variableRecords, "source", {}]],
    "residualObjects" -> AssociationThread[resolvedVariables, Lookup[variableRecords, "residualObjects", {}]],
    "variableData" -> variableData,
    "naiveIBP" -> KeyDrop[ibpData, {"seedRecords", "context"}],
    "context" -> context,
    "equationConvention" -> "D[normalized tagged masters,var] == matrices[var].normalized tagged masters + sources[var]",
    "derivativeRoute" -> "direct tree phase derivative + direct h line momentum-magnitude derivative -> direct tree dtau reduction",
    "formulaDLogUsedQ" -> False
    |>
   ];


dsTreeNaiveDERaw018[context_Association, variables_: Automatic, masters_: Automatic, OptionsPattern[DSTreeNaiveDE]] /; dsContextQ[context] := Module[
   {ibpData},
   If[treeFormulaMasslessPendingQ018[context],
    Return[treeFormulaPendingRederivation018["DSTreeNaiveDE", context]]
    ];
   ibpData = dsTreeNaiveIBPRaw018[context, masters,
     AuditLevel -> OptionValue[AuditLevel], ProgressReporting -> OptionValue[ProgressReporting]];
   If[Lookup[ibpData, "status", "failed"] =!= "solved", ibpData,
    dsTreeNaiveDEFromIBP[ibpData, variables, ProgressReporting -> OptionValue[ProgressReporting]]]
   ];


DSTreeNaiveDE[_, ___] := (Message[DSTreeNaiveDE::badibp]; <|"status" -> "failed", "reason" -> "invalidContextOrIBPData"|>);
