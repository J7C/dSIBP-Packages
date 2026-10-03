(* ::Package:: *)

(* 本文件按 linearData 已冻结的积分顺序构造预约化 targets 和解析 derivative closure。
   package 只生成计划与输入，不运行 Kira。 *)


(* ::Chapter:: *)
(*Kira 两阶段 reduction 计划*)

DSKiraPlan::badlinear = "DSKiraPlan 需要 DSLinear 返回的 backend-neutral linearData。 DSKiraPlan requires backend-neutral linearData returned by DSLinear.";
DSKiraPlan::badspec = "Kira 计划配置无效：`1`。 The Kira plan specification is invalid: `1`.";
DSKiraPlan::badstage = "stage 只允许 \"preReduction\" 或 \"formal\"，收到 `1`。 stage must be \"preReduction\" or \"formal\"; received `1`.";
DSKiraPlan::badbasis = "formal 计划需要可闭合的 activeBasis：`1`。 A formal plan requires a closed activeBasis: `1`.";
DSKiraPlan::incomplete = "formal 计划只接受 completeSystemQ=True 的 linearData。 A formal plan requires linearData with completeSystemQ=True.";


Options[DSKiraPlan] = {ProgressReporting -> Automatic};


(* ::Section::Closed:: *)
(*既定积分顺序与显式重排*)

DSReorderIntegrals::badlinear = "DSReorderIntegrals 需要 DSLinear 返回的 backend-neutral linearData。";
DSReorderIntegrals::badorder = "积分顺序必须是由现有 J 或积分 ID 组成的非空列表。";


DSReorderIntegrals[linearData_Association, order_List] := Module[{result},
   If[Lookup[linearData, "dSIBPStatus", "failed"] =!= "generated" ||
     ! ListQ[Lookup[linearData, "integralList", Missing["integralList"]]],
    Message[DSReorderIntegrals::badlinear];
    Return[<|"status" -> "failed", "reason" -> "notLinearData"|>]
    ];
   If[Lookup[Lookup[linearData, "activeBasis", <||>], "status", "disabled"] === "configured",
    Message[DSReorderIntegrals::badorder];
    Return[<|"status" -> "failed", "reason" -> "reorderMustPrecedeUserMI"|>]
    ];
   If[order === {},
    Message[DSReorderIntegrals::badorder];
    Return[<|"status" -> "failed", "reason" -> "emptyIntegralOrder"|>]
    ];
   result = reorderLinearSystemIntegrals[linearData, order];
   Join[result, <|
     "integralOrderAuthority" -> "linearData.integralList",
     "integralOrderDigest" -> dsKiraExpressionDigest[result["integralList"]]
     |>]
   ];


DSReorderIntegrals[_, _] := (Message[DSReorderIntegrals::badorder]; <|"status" -> "failed", "reason" -> "invalidIntegralOrder"|>);

dsKiraPlanIntegralFromItem[item_, linearData_Association] := Which[
   Head[item] === J && MemberQ[linearData["integralList"], item], item,
   IntegerQ[item] && 1 <= item <= Length[linearData["integralList"]], linearData["integralList"][[item]],
   True, Missing["UnknownIntegral", item]
   ];


dsKiraPlanIntegralList[items_List, linearData_Association] := DeleteDuplicates@DeleteMissing[
   dsKiraPlanIntegralFromItem[#, linearData] & /@ items
   ];


dsKiraPlanCertificate[activeData_Association] := Module[{payload},
   payload = HoldComplete[
     Lookup[activeData, "activeExpressions", {}],
     Lookup[activeData, "derivativeVariables", {}],
     Lookup[activeData, "rawDerivatives", {}],
     Lookup[activeData, "derivativeTargetIntegrals", {}]
     ];
   <|
    "status" -> "frozenBeforeNumericalRules",
    "hashAlgorithm" -> "SHA256",
    "hash" -> IntegerString[Hash[payload, "SHA256"], 16, 64],
    "activeCount" -> Lookup[activeData, "activeCount", 0],
    "derivativeVariableCount" -> Length[Lookup[activeData, "derivativeVariables", {}]],
    "derivativeTargetCount" -> Length[Lookup[activeData, "derivativeTargetIntegrals", {}]]
    |>
   ];


(* ::Section:: *)
(*公开 Kira 两阶段计划*)

DSKiraPlan[linearData_Association, spec_Association, OptionsPattern[]] := Module[
   {stage, preferred, order, ordered, candidates, activeSetting, preview, activeData, targets,
    numericStage, coefficientRules, outputDirectory, jobOptions, certificate, progress},
   If[Lookup[linearData, "dSIBPStatus", "failed"] =!= "generated" || ! KeyExistsQ[linearData, "linearEquations"],
    Message[DSKiraPlan::badlinear]; Return[<|"status" -> "failed", "reason" -> "notLinearData"|>]
    ];
   progress = OptionValue[ProgressReporting];
   stage = Lookup[spec, "stage", Missing["stage"]];
   If[! MemberQ[{"preReduction", "formal"}, stage],
    Message[DSKiraPlan::badstage, stage]; Return[<|"status" -> "failed", "reason" -> "invalidStage"|>]
    ];
   If[stage === "formal" && ! TrueQ[Lookup[linearData, "completeSystemQ", False]],
    Message[DSKiraPlan::incomplete];
    Return[<|"status" -> "failed", "reason" -> "incompleteSystemForFormalReduction",
      "completeSystemQ" -> Lookup[linearData, "completeSystemQ", False]|>]
    ];
   preferred = dsKiraPlanIntegralList[Lookup[spec, "preferredIntegrals", {}], linearData];
   order = linearData["integralList"];
   ordered = linearData;
   coefficientRules = Lookup[spec, "coefficientRules", {}];
   outputDirectory = Lookup[spec, "outputDirectory", None];
   jobOptions = Lookup[spec, "jobOptions", Automatic];
   If[! MemberQ[{None, Automatic}, outputDirectory] && ! StringQ[outputDirectory],
    Message[DSKiraPlan::badspec, "outputDirectory must be a string or None"];
    Return[<|"status" -> "failed", "reason" -> "invalidOutputDirectory"|>]
    ];
   If[stage === "preReduction",
    candidates = Replace[Lookup[spec, "candidateIntegrals", Automatic], Automatic :> order];
    If[! ListQ[candidates], candidates = {candidates}];
    candidates = dsKiraPlanIntegralList[candidates, ordered];
    If[candidates === {},
     Message[DSKiraPlan::badspec, "empty candidateIntegrals"];
     Return[<|"status" -> "failed", "reason" -> "emptyCandidateIntegrals"|>]
     ];
    Return[<|
      "status" -> "planned", "kiraPlanQ" -> True, "stage" -> stage,
      "caseName" -> Lookup[linearData, "caseName", Missing["caseName"]],
      "linearData" -> linearData, "integralOrder" -> order,
      "orderingConvention" -> "linearDataIntegralList",
      "preferredIntegrals" -> preferred, "targetIntegrals" -> candidates,
      "targetCount" -> Length[candidates], "activeBasis" -> None,
      "numericStage" -> "symbolic", "coefficientRules" -> coefficientRules,
      "outputDirectory" -> outputDirectory, "jobOptions" -> jobOptions,
      "phaseIsolation" -> <|"stage" -> stage, "requiresSeparateOutputDirectoryQ" -> True|>
      |>]
    ];
   activeSetting = Lookup[spec, "activeBasis", Automatic];
   If[! AssociationQ[activeSetting] && ! (
       activeSetting === Automatic &&
        Lookup[Lookup[ordered, "activeBasis", <||>], "status", "disabled"] === "configured"
       ),
    Message[DSKiraPlan::badbasis, "missing activeBasis Association"];
    Return[<|"status" -> "failed", "reason" -> "missingActiveBasis"|>]
    ];
   preview = dsStageRun[
     "构造解析 active-basis 一阶导数与 target closure / Building analytic active-basis first derivatives and target closure",
     dsKiraAttachActiveBasis[ordered, activeSetting],
     progress
     ];
   If[Lookup[preview, "status", "missing"] =!= "generated" ||
     Lookup[Lookup[preview, "activeBasis", <||>], "status", "failed"] =!= "configured",
    Message[DSKiraPlan::badbasis, Lookup[preview, "reason", "closure failed"]];
    Return[<|"status" -> "failed", "reason" -> "activeBasisClosureFailed", "preview" -> preview|>]
    ];
   activeData = preview["activeBasis"];
   targets = activeData["targetIntegralIDs"];
   numericStage = Lookup[spec, "numericStage", "symbolic"];
   If[! MemberQ[{"symbolic", "postDerivative"}, numericStage],
    Message[DSKiraPlan::badspec, "numericStage"];
    Return[<|"status" -> "failed", "reason" -> "invalidNumericStage"|>]
    ];
   certificate = dsKiraPlanCertificate[activeData];
   <|
    "status" -> "planned", "kiraPlanQ" -> True, "stage" -> stage,
    "caseName" -> Lookup[linearData, "caseName", Missing["caseName"]],
    "linearData" -> linearData, "integralOrder" -> order,
    "orderingConvention" -> "linearDataIntegralList",
    "preferredIntegrals" -> preferred, "activeBasis" -> activeSetting,
    "activeBasisPreview" -> activeData,
    "preparedLinearData" -> preview,
    "analyticDerivativeCertificate" -> certificate,
    "targetIntegralIDsPreview" -> targets,
    "minimalTargetsQ" -> True,
    "numericStage" -> numericStage, "coefficientRules" -> coefficientRules,
    "outputDirectory" -> outputDirectory, "jobOptions" -> jobOptions,
    "phaseIsolation" -> <|"stage" -> stage, "requiresSeparateOutputDirectoryQ" -> True|>
    |>
   ];


dsKiraPlanQ[plan_Association] := TrueQ[Lookup[plan, "kiraPlanQ", False]] &&
   Lookup[plan, "status", "failed"] === "planned";


(* ::Section:: *)
(*计划驱动的 Kira 导出*)

DSKiraExport[plan_Association] /; dsKiraPlanQ[plan] := Module[
   {result, manifest, path, prepared, expectedCertificate, actualCertificate},
   prepared = If[
     Lookup[plan, "stage", "formal"] === "preReduction",
     Lookup[plan, "linearData", Missing["linearData"]],
     Lookup[plan, "preparedLinearData", Missing["preparedLinearData"]]
     ];
   If[! AssociationQ[prepared],
    Return[<|"status" -> "failed", "reason" -> "missingPreparedLinearData"|>]
    ];
   If[Lookup[plan, "stage", "formal"] === "formal",
    expectedCertificate = Lookup[plan, "analyticDerivativeCertificate", <||>];
    actualCertificate = dsKiraPlanCertificate[Lookup[prepared, "activeBasis", <||>]];
    If[Lookup[expectedCertificate, "hash", None] =!= Lookup[actualCertificate, "hash", Missing["hash"]],
     Return[<|"status" -> "failed", "reason" -> "preparedActiveBasisDigestMismatch",
       "expectedCertificate" -> expectedCertificate, "actualCertificate" -> actualCertificate|>]
     ]
    ];
   result = DSKiraExport[
     prepared,
     KiraActiveBasis -> If[Lookup[plan, "stage", "formal"] === "formal", Automatic, None],
     KiraTargetIntegrals -> Lookup[plan, "targetIntegrals", Automatic],
     KiraCoefficientRules -> plan["coefficientRules"],
     KiraJobOptions -> plan["jobOptions"],
     KiraNumericStage -> plan["numericStage"],
     KiraRequireCompleteSystem -> TrueQ[Lookup[plan, "stage", "formal"] === "formal"],
     OutputDirectory -> plan["outputDirectory"]
     ];
   If[Lookup[result, "status", "failed"] =!= "ready", Return[result]];
   manifest = Join[Lookup[result, "dSIBPExportManifest", <||>], <|
      "kiraPlan" -> KeyDrop[plan, {"linearData", "activeBasisPreview", "preparedLinearData"}]
      |>];
   path = Lookup[result, "dSIBPExportManifestPath", Missing["NotWritten"]];
   If[StringQ[path], Quiet[Check[Put[manifest, path], Null]]];
   Join[result, <|"dSIBPExportManifest" -> manifest|>]
   ];
