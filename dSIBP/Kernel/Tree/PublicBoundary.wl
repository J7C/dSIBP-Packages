(* ::Package:: *)
(* 本模块把论文 vertex-basis 公式实现限制在 Private 内部，并为 tree 公开接口提供
   J[sectorKey,timeShifts,stateBits] 适配。它不改变 massive-only 递推或 dlog 公式，也不绕过
   massless quotient 的 PendingRederivation 门禁。 *)

(* ::Chapter:: *)
(*统一 J 与私有 vertex basis 的双向映射*)

(* 公式内核按顶点保存 massive 状态；先构造既有 line-packed producer 对象，再通过
   中央 registry 唯一转换为公开 time-only 表示。 *)
dsTreeFormulaIntegralToPublic018[int : J[_List], sectorKey_String, context_Association] := Module[
   {familyContext, family, topo, packs, aList, baseline, linePacks, line, states,
    vertexIndex, legIndex, internalIntegral},
   familyContext = dsTreeFamilyContext[context];
   family = dsTreeFamilyBySector[sectorKey, familyContext];
   If[Head[family] === Missing || ! treeIntegralQ[int, family], Return[$Failed]];
   topo = family["topology"];
   packs = First[int];
   aList = packs[[All, 1]];
   baseline = Lookup[family, "loopBaselinePowers", <||>];
   linePacks = Table[
     line = topo["lines"][[e]];
     Switch[
      Lookup[line, "state", "full"] === "shrunk" || Lookup[line, "packType", ""] === "shrunk",
      True,
      If[lineIndexedPowerQ[line], {Lookup[baseline, line["id"], 0]}, {fixedLineSentinel018[]}],
      False,
      states = Table[
        vertexIndex = FirstPosition[
          family["vertexOrder"],
          line["endpoints"][[endpointSlot]],
          Missing["NoVertex"]
          ];
        If[Head[vertexIndex] === Missing, Return[$Failed]];
        legIndex = FirstPosition[
          Lookup[family["vertices"][[First[vertexIndex], "massiveLegs"]], "id", {}],
          {line["id"], endpointSlot},
          Missing["NoLeg"]
          ];
        If[Head[legIndex] === Missing, Return[$Failed]];
        packs[[First[vertexIndex], 1 + First[legIndex]]],
        {endpointSlot, 2}
        ];
      If[lineIndexedPowerQ[line],
       Prepend[states, Lookup[baseline, line["id"], 0]],
       Prepend[states, fixedLineSentinel018[]]
       ]
      ],
     {e, topo["nE"]}
     ];
   internalIntegral = J[aList, linePacks, {}];
   dsTimeOnlyInternalIntegralToPublic020[internalIntegral, context]
   ];


dsTreePublicIntegralToFormula018[
   int : J[_, _, _], sectorKey_String, context_Association
   ] := Module[
   {internalIntegral, metadata, familyContext, family, topo, aList, linePacks, lineIds, vertexPacks,
    legId, linePosition, endpointSlot, statePosition, states},
   internalIntegral = dsTimeOnlyPublicIntegralToInternal020[int, context];
   If[internalIntegral === $Failed, Return[$Failed]];
   metadata = SelectFirst[
     context["sectors"],
     Lookup[#, "sectorKey", None] === sectorKey && integralMatchesSectorMetadataQ[internalIntegral, #] &,
     Missing["NoSector"]
     ];
   If[Head[metadata] === Missing, Return[$Failed]];
   familyContext = dsTreeFamilyContext[context];
   family = dsTreeFamilyBySector[sectorKey, familyContext];
   If[Head[family] === Missing, Return[$Failed]];
   topo = family["topology"];
   aList = internalIntegral[[1]];
   linePacks = internalIntegral[[2]];
   lineIds = Lookup[topo["lines"], "id", {}];
   If[Length[aList] =!= Length[family["vertexOrder"]], Return[$Failed]];

   (* massiveLegs 的 id 已保存 root line 与端点槽；公开 full pack 的后两槽正是这两个端点态。
      fixed sentinel 只占第一槽，因此 cycle/fixed line 共用同一读取规则。 *)
   vertexPacks = Table[
     states = Map[
       Function[leg,
        legId = Lookup[leg, "id", Missing["NoLegId"]];
        If[! MatchQ[legId, {_, 1 | 2}], Return[$Failed]];
        {linePosition, endpointSlot} = {
          FirstPosition[lineIds, legId[[1]], Missing["NoLine"]],
          legId[[2]]
          };
        If[Head[linePosition] === Missing, Return[$Failed]];
        linePosition = First[linePosition];
        statePosition = 1 + endpointSlot;
        If[Length[linePacks[[linePosition]]] < statePosition, Return[$Failed]];
        linePacks[[linePosition, statePosition]]
        ],
       family["vertices"][[vertexIndex, "massiveLegs"]]
       ];
     Prepend[states, aList[[vertexIndex]]],
     {vertexIndex, Length[family["vertices"]]}
     ];
   If[FreeQ[vertexPacks, $Failed], J[vertexPacks], $Failed]
   ];


dsTreeTaggedTermToPublic018[term_Association, context_Association] := Module[
   {sectorKey, formulaIntegral, publicIntegral},
   sectorKey = Lookup[term, "sectorKey", Missing["NoSector"]];
   formulaIntegral = Lookup[term, "integral", Missing["NoIntegral"]];
   If[! StringQ[sectorKey] || ! MatchQ[formulaIntegral, J[_List]], Return[$Failed]];
   publicIntegral = dsTreeFormulaIntegralToPublic018[formulaIntegral, sectorKey, context];
   If[publicIntegral === $Failed, Return[$Failed]];
   Join[KeyDrop[term, {"integral"}], <|"integral" -> publicIntegral|>]
   ];


dsTreeTaggedTermToFormula018[term_Association, context_Association] := Module[
   {sectorKey, publicIntegral, formulaIntegral},
   sectorKey = Lookup[term, "sectorKey", Missing["NoSector"]];
   publicIntegral = Lookup[term, "integral", Missing["NoIntegral"]];
   If[! StringQ[sectorKey] || ! MatchQ[publicIntegral, J[_, _, _]], Return[$Failed]];
   formulaIntegral = dsTreePublicIntegralToFormula018[publicIntegral, sectorKey, context];
   If[formulaIntegral === $Failed, Return[$Failed]];
   Join[KeyDrop[term, {"integral"}], <|"integral" -> formulaIntegral|>]
   ];


(* ::Chapter:: *)
(*Seed 与 tagged linearData 公开化*)

dsTreeLinearDataToPublic018[data_Association, context_Association] := Module[{terms},
   terms = dsTreeTaggedTermToPublic018[#, context] & /@ Lookup[data, "terms", {}];
   If[MemberQ[terms, $Failed], Return[$Failed]];
   Join[
    KeyDrop[data, {"terms", "expression"}],
    <|
     "terms" -> terms,
     "expression" -> Total[Lookup[#, "coefficient"] Lookup[#, "integral"] & /@ terms],
     "representation" -> "sectorTaggedJ[sectorKey,timeShifts,stateBits]"
     |>
    ]
   ];


dsTreeLinearDataToFormula018[data_Association, context_Association] := Module[{terms},
   terms = dsTreeTaggedTermToFormula018[#, context] & /@ Lookup[data, "terms", {}];
   If[MemberQ[terms, $Failed], Return[$Failed]];
   Join[
    KeyDrop[data, {"terms", "expression"}],
    <|
     "terms" -> terms,
     "expression" -> Total[Lookup[#, "coefficient"] Lookup[#, "integral"] & /@ terms]
     |>
    ]
   ];


dsTreeSeedRecordToPublic018[record_Association, context_Association, referenceInt_J] := Module[
   {linearData, publicLinearData},
   linearData = Lookup[record, "treeLinearData", Missing["NoTreeLinearData"]];
   If[! AssociationQ[linearData], Return[$Failed]];
   publicLinearData = dsTreeLinearDataToPublic018[linearData, context];
   If[publicLinearData === $Failed, Return[$Failed]];
   publicLinearData = Join[
     KeyDrop[publicLinearData, {"referenceTreeIntegral", "referenceLoopIntegral"}],
     <|"referenceTreeIntegral" -> referenceInt|>
     ];
   Join[
    KeyDrop[record, {"treeSeed", "treeIntegral", "treeLinearData"}],
    <|
     "treeSeed" -> publicLinearData["expression"],
     "treeIntegral" -> referenceInt,
     "treeLinearData" -> publicLinearData,
     "representation" -> "J[sectorKey,timeShifts,stateBits]"
     |>
    ]
   ];


(* ::Chapter:: *)
(*批量 direct pure-time 模板*)

(* DSSeeds 的 timeOnly 路线与单积分 DSTreeSeeds 共用同一个原子生成器。
   连续 a 指标保持符号，端点态完整遍历 0/1；sector 与离散规则只作 provenance，
   不另建积分 Head，也不把 Private vertex basis 暴露给用户。 *)
dsPureTimeDirectTemplateRecord018[
   raw_Association,
   family_Association,
   vertexId_,
   formulaIntegral_J,
   context_Association
   ] := Module[
   {sectorKey, publicIntegral, internalIntegral, publicRecord, baseIntegral, discreteVariables,
    discretePositions, discreteRules, metadata},
   sectorKey = family["sector"];
   publicIntegral = dsTreeFormulaIntegralToPublic018[formulaIntegral, sectorKey, context];
   If[publicIntegral === $Failed, Return[$Failed]];
   internalIntegral = dsTimeOnlyPublicIntegralToInternal020[publicIntegral, context];
   If[internalIntegral === $Failed, Return[$Failed]];
   publicRecord = dsTreeSeedRecordToPublic018[raw, context, publicIntegral];
   If[publicRecord === $Failed, Return[$Failed]];
   baseIntegral = makeBaseIntegral[family["topology"]];
   discreteVariables = DeleteDuplicates[Cases[baseIntegral, _n, Infinity]];
   discretePositions = Position[baseIntegral, _n, Infinity];
   discreteRules = Thread[Extract[baseIntegral, discretePositions] -> Extract[internalIntegral, discretePositions]];
   metadata = SelectFirst[
     context["sectors"],
     Lookup[#, "sectorKey", None] === sectorKey &,
     Missing["NoSectorMetadata"]
     ];
   If[Head[metadata] === Missing, Return[$Failed]];
   <|
    "source" -> {"directPureTime", sectorKey},
    "generationRoute" -> "directPureTime",
    "sectorKey" -> sectorKey,
    "sectorShrunkLines" -> Lookup[family["topology"], "sectorShrunkLines", {}],
    "sectorMetadata" -> metadata,
    "generator" -> dtau[vertexId],
    "ibpClass" -> "tIBP",
    "continuousIndices" -> (a /@ family["vertexOrder"]),
    "sourceIntegral" -> publicIntegral,
    "discreteVariables" -> discreteVariables,
    "discreteRules" -> discreteRules,
    "discreteStateCountExpected" -> 2^Length[discreteVariables],
    "equation" -> publicRecord["treeSeed"],
    "forbiddenNData" -> {},
    "eomCanonicalQ" -> True,
    "representation" -> "J[sectorKey,timeShifts,stateBits]",
    "ibpMode" -> "timeOnly",
    "directSeedProvenance" -> KeyTake[
      publicRecord,
      {"generationRoute", "sectorKey", "generator", "treeIntegral", "representation"}
      ]
    |>
   ];


dsPureTimeDirectTemplateData018[context_Association] /; dsContextQ[context] := Module[
   {familyContext, families, records},
   If[treeFormulaMasslessPendingQ018[context],
    Return[treeFormulaPendingRederivation018["DSSeeds", context]]
    ];
   familyContext = dsTreeFamilyContext[context];
   If[familyContext === $Failed,
    Return[<|"status" -> "failed", "reason" -> "treeFamilyInitializationFailed"|>]
    ];
   families = familyContext["families"];
   records = Flatten@Table[
      With[{formulaIntegrals = treeMasterList[family, a /@ family["vertexOrder"]]},
       If[formulaIntegrals === $Failed,
        {$Failed},
        Table[
         With[{raw = dsDirectTreeSeedRecord[vertexId, formulaIntegral, family, familyContext]},
          If[! AssociationQ[raw] || Lookup[raw, "status", "failed"] =!= "generated",
           $Failed,
           dsPureTimeDirectTemplateRecord018[
            raw, family, vertexId, formulaIntegral, context
            ]
           ]
          ],
         {vertexId, family["vertexOrder"]},
         {formulaIntegral, formulaIntegrals}
         ]
        ]
       ],
      {family, families}
      ];
   If[MemberQ[records, $Failed],
    Return[<|"status" -> "failed", "reason" -> "directPureTimeTemplateFailed"|>]
    ];
   <|
    "status" -> "generated",
    "generationRoute" -> "directPureTime",
    "sectorCount" -> Length[families],
    "templateCount" -> Length[records],
    "allSeeds" -> records
    |>
   ];


(* ::Chapter:: *)
(*公开 DSTreeSeeds 与 repIterative*)

DSTreeSeeds[vertex_, int : J[_, _, _], context_Association] /; dsContextQ[context] := Module[
   {internalIntegral, metadata, sectorKey, familyContext, family, formulaIntegral, raw},
   If[treeFormulaMasslessPendingQ018[context],
    Return[treeFormulaPendingRederivation018["DSTreeSeeds", context]]
    ];
   internalIntegral = dsTimeOnlyPublicIntegralToInternal020[int, context];
   If[internalIntegral === $Failed,
    Return[<|"status" -> "failed", "reason" -> "invalidPublicTimeOnlyIntegral"|>]
    ];
   metadata = integralSectorMetadata018[context["topology"], internalIntegral];
   If[Head[metadata] === Missing, Return[<|"status" -> "failed", "reason" -> "unknownTreeSector"|>]];
   sectorKey = metadata["sectorKey"];
   familyContext = dsTreeFamilyContext[context];
   family = dsTreeFamilyBySector[sectorKey, familyContext];
   formulaIntegral = dsTreePublicIntegralToFormula018[int, sectorKey, context];
   If[Head[family] === Missing || formulaIntegral === $Failed,
    Return[<|"status" -> "failed", "reason" -> "treePublicToFormulaFailed"|>]
    ];
   raw = dsDirectTreeSeedRecord[vertex, formulaIntegral, family, familyContext];
   If[! AssociationQ[raw], $Failed, dsTreeSeedRecordToPublic018[raw, context, int]]
   ];


DSTreeSeeds[int : J[_, _, _], context_Association] /; dsContextQ[context] :=
  If[treeFormulaMasslessPendingQ018[context],
   treeFormulaPendingRederivation018["DSTreeSeeds", context],
   DSTreeSeeds[#, int, context] & /@ Lookup[
     Select[makeIBPGenerators[context["topology"]], #["type"] === "time" &],
     "vertex"
     ]
   ];


dsTreePublicExpressionData018[expr_, context_Association] := Module[
   {terms, parsed},
   terms = dsExpandedTerms[Expand[expr]];
   parsed = Map[
     Function[term,
      Module[{integrals, int, internalIntegral, coefficient, metadata},
       (* 裸 J 本身位于 level 0；若只从 level 1 搜索，单积分输入会被误判为空。 *)
       integrals = DeleteDuplicates[Cases[term, item : J[_, _, _] :> item, {0, Infinity}]];
       If[Length[integrals] =!= 1, Return[$Failed]];
       int = First[integrals];
       internalIntegral = dsTimeOnlyPublicIntegralToInternal020[int, context];
       If[internalIntegral === $Failed, Return[$Failed]];
       coefficient = Cancel[term/int];
       metadata = integralSectorMetadata018[context["topology"], internalIntegral];
       If[Head[metadata] === Missing, Return[$Failed]];
       <|"sectorKey" -> metadata["sectorKey"], "integral" -> int, "coefficient" -> coefficient|>
       ]
      ],
     terms
     ];
   If[MemberQ[parsed, $Failed], $Failed,
    <|"status" -> "generated", "terms" -> parsed,
      "expression" -> Total[Lookup[#, "coefficient"] Lookup[#, "integral"] & /@ parsed]|>]
   ];


repIterative[data_Association, end_: Automatic, context_Association, opts : OptionsPattern[]] /;
   dsContextQ[context] := Module[{formulaData, result},
   If[treeFormulaMasslessPendingQ018[context],
    Return[treeFormulaPendingRederivation018["repIterative", context]]
    ];
   formulaData = dsTreeLinearDataToFormula018[data, context];
   If[formulaData === $Failed, Return[<|"status" -> "error", "reason" -> "invalidUnifiedTreeData"|>]];
   result = dsRepIterativeTreeLinearData[formulaData, end, context, opts];
   If[! AssociationQ[result] || ! ListQ[Lookup[result, "terms", None]], result,
    dsTreePublicResult018[result, context]
    ]
   ];


repIterative[expr_, end_: Automatic, context_Association, opts : OptionsPattern[]] /;
   dsContextQ[context] := Module[{data},
   If[treeFormulaMasslessPendingQ018[context],
    Return[treeFormulaPendingRederivation018["repIterative", context]]
    ];
   data = dsTreePublicExpressionData018[expr, context];
   If[data === $Failed,
    <|"status" -> "error", "reason" -> "invalidUnifiedTreeExpression"|>,
    repIterative[data, end, context, opts]
    ]
   ];


(* ::Chapter:: *)
(*Naive 与 dlog 结果的公开结构化映射*)

dsTreeMasterToFormula018[master_Association, context_Association] := Module[{term},
   term = dsTreeTaggedTermToFormula018[
     KeyTake[master, {"sectorKey", "integral", "coefficient"}],
     context
     ];
   If[term === $Failed, $Failed,
    Join[KeyDrop[master, {"integral"}], <|"integral" -> term["integral"]|>]
    ]
   ];


(* 公开 naive-IBP 结果中的 lhs/rhsTerms 已经是统一 time-only J；本适配器只恢复
   Private 公式内核所需的 sector-tagged vertex basis，并保留全部系数与 sectorKey。 *)
dsTreeNaiveReductionToFormula018[record_Association, context_Association] := Module[
   {lhs, rhsTerms},
   lhs = dsTreeTaggedTermToFormula018[Lookup[record, "lhs", <||>], context];
   rhsTerms = dsTreeTaggedTermToFormula018[#, context] & /@ Lookup[record, "rhsTerms", {}];
   If[lhs === $Failed || MemberQ[rhsTerms, $Failed], Return[$Failed]];
   Join[
    KeyDrop[record, {"lhs", "rhsTerms", "rhsExpression"}],
    <|
     "lhs" -> lhs,
     "rhsTerms" -> rhsTerms,
     "rhsExpression" -> Total[Lookup[#, "coefficient"] Lookup[#, "integral"] & /@ rhsTerms]
     |>
    ]
   ];


(* 输入必须是同一 context 产生的公开 solved 数据；失败时不猜测 sector 或积分形状。 *)
dsTreeNaiveIBPDataToFormula018[data_Association, context_Association] := Module[
   {masters, reductions},
   masters = dsTreeMasterToFormula018[#, context] & /@ Lookup[data, "masters", {}];
   reductions = dsTreeNaiveReductionToFormula018[#, context] & /@ Lookup[data, "reductions", {}];
   If[masters === {} || reductions === {} || MemberQ[Join[masters, reductions], $Failed],
    Return[$Failed]
    ];
   Join[
    KeyDrop[data, {"masters", "reductions", "bareMasters", "publicRepresentationAudit"}],
    <|"masters" -> masters, "reductions" -> reductions, "representation" -> "J[vertexPacks]"|>
    ]
   ];


dsTreePublicizeData018[value_, context_Association, sectorHint_: None] := Module[
   {hint = sectorHint, mapped, publicIntegral},
   Which[
    MatchQ[value, dsTreeToken[_String, J[_List]]],
    publicIntegral = dsTreeFormulaIntegralToPublic018[value[[2]], value[[1]], context];
    If[publicIntegral === $Failed, value, publicIntegral],

    MatchQ[value, J[_List]] && StringQ[hint],
    publicIntegral = dsTreeFormulaIntegralToPublic018[value, hint, context];
    If[publicIntegral === $Failed, value, publicIntegral],

    AssociationQ[value],
    hint = Lookup[value, "sectorKey", Lookup[value, "sector", hint]];
    mapped = Association@KeyValueMap[
       Function[{key, item}, key -> dsTreePublicizeData018[item, context, hint]],
       value
       ];
    If[AssociationQ[Lookup[mapped, "treeLinearData", None]],
     mapped = Join[mapped, <|
        "treeSeed" -> mapped["treeLinearData", "expression"],
        "treeIntegral" -> Lookup[
          mapped["treeLinearData"],
          "referenceLoopIntegral",
          Lookup[mapped, "treeIntegral", Missing["NoReferenceIntegral"]]
          ]
        |>]
     ];
    If[KeyExistsQ[mapped, "terms"] && ListQ[mapped["terms"]],
     mapped = Join[mapped, <|
        "expression" -> Total[
          Lookup[#, "coefficient", 1] Lookup[#, "integral", 0] & /@ mapped["terms"]
          ]
        |>]
     ];
    If[KeyExistsQ[mapped, "rhsTerms"] && ListQ[mapped["rhsTerms"]],
     mapped = Join[mapped, <|
        "rhsExpression" -> Total[
          Lookup[#, "coefficient", 1] Lookup[#, "integral", 0] & /@ mapped["rhsTerms"]
          ]
        |>]
     ];
    mapped,

    ListQ[value], dsTreePublicizeData018[#, context, hint] & /@ value,
    Head[value] === Rule,
    Rule @@ (dsTreePublicizeData018[#, context, hint] & /@ List @@ value),
    True, value
    ]
   ];


dsTreePublicResult018[result_Association, context_Association, auditLevel_: "standard"] := Module[
   {public, masters, legacyIntegrals, oldPublicLeakQ, privateTokens, representationAudit},
   public = dsTreePublicizeData018[result, context];
   If[ListQ[Lookup[public, "masters", None]],
    masters = public["masters"];
    If[And @@ (AssociationQ /@ masters),
     public = Join[public, <|"bareMasters" -> Lookup[masters, "integral", {}]|>]
     ]
    ];
   If[ListQ[Lookup[public, "contactMaps", None]],
    public = Join[public, <|"sourceEquations" -> Lookup[public["contactMaps"], "rowsByVertex", {}]|>]
    ];
   representationAudit = If[
      auditLevel === "full",
      legacyIntegrals = DeleteDuplicates[Cases[public, item : J[_List] :> item, Infinity]];
      oldPublicLeakQ = dsTimeOnlyOldIntegralLeakQ020[public];
      privateTokens = DeleteDuplicates[Cases[public, _dsTreeToken, Infinity]];
      <|"status" -> "audited", "passQ" -> TrueQ[legacyIntegrals === {} && ! oldPublicLeakQ && privateTokens === {}],
        "legacyVertexIntegrals" -> legacyIntegrals, "oldPublicIntegralLeakQ" -> oldPublicLeakQ,
        "privateTreeTokens" -> privateTokens|>,
      <|"status" -> "producerGuaranteed", "passQ" -> True,
        "legacyVertexIntegrals" -> {}, "oldPublicIntegralLeakQ" -> False, "privateTreeTokens" -> {},
       "reason" -> "unified public conversion completed by the same producer"|>
     ];
   Join[public, <|
     "representation" -> "J[sectorKey,timeShifts,stateBits]",
     "auditLevel" -> auditLevel,
     "publicRepresentationAudit" -> representationAudit
     |>]
   ];


DSTreeNaiveIBP[context_Association, masters_: Automatic, OptionsPattern[]] /; dsContextQ[context] := Module[
   {formulaMasters, raw},
   If[treeFormulaMasslessPendingQ018[context],
    Return[treeFormulaPendingRederivation018["DSTreeNaiveIBP", context]]
    ];
   formulaMasters = If[masters === Automatic,
     Automatic,
     dsTreeMasterToFormula018[#, context] & /@ masters
     ];
   If[ListQ[formulaMasters] && MemberQ[formulaMasters, $Failed],
    Return[<|"status" -> "failed", "reason" -> "invalidUnifiedMasters"|>]
    ];
   raw = dsTreeNaiveIBPRaw018[
     context,
     formulaMasters,
     AuditLevel -> OptionValue[AuditLevel],
     ProgressReporting -> OptionValue[ProgressReporting]
     ];
   If[AssociationQ[raw], dsTreePublicResult018[raw, context, OptionValue[AuditLevel]], raw]
   ];


DSTreeNaiveDE[context_Association, variables_: Automatic, masters_: Automatic, OptionsPattern[]] /;
   dsContextQ[context] := Module[{formulaMasters, raw},
   If[treeFormulaMasslessPendingQ018[context],
    Return[treeFormulaPendingRederivation018["DSTreeNaiveDE", context]]
    ];
   formulaMasters = If[masters === Automatic,
     Automatic,
     dsTreeMasterToFormula018[#, context] & /@ masters
     ];
   If[ListQ[formulaMasters] && MemberQ[formulaMasters, $Failed],
    Return[<|"status" -> "failed", "reason" -> "invalidUnifiedMasters"|>]
    ];
    raw = dsTreeNaiveDERaw018[
     context,
     variables,
      formulaMasters,
      AuditLevel -> OptionValue[AuditLevel],
      ProgressReporting -> OptionValue[ProgressReporting]
     ];
   If[AssociationQ[raw], dsTreePublicResult018[raw, context, OptionValue[AuditLevel]], raw]
   ];


(* 支持把公开 DSTreeNaiveIBP 返回值直接传回 DSTreeNaiveDE；只在表示审计通过时跨越边界。 *)
DSTreeNaiveDE[ibpData_Association, variables_: Automatic, OptionsPattern[]] /;
   Lookup[ibpData, "status", "failed"] === "solved" &&
     Lookup[ibpData, "representation", None] === "J[sectorKey,timeShifts,stateBits]" := Module[
   {context, formulaData, raw},
   context = Lookup[ibpData, "context", Missing["context"]];
   If[! dsContextQ[context] ||
     ! TrueQ[Lookup[Lookup[ibpData, "publicRepresentationAudit", <||>], "passQ", False]],
    Return[<|"status" -> "failed", "reason" -> "invalidPublicNaiveIBPData"|>]
    ];
   formulaData = dsTreeNaiveIBPDataToFormula018[ibpData, context];
   If[formulaData === $Failed,
    Return[<|"status" -> "failed", "reason" -> "publicNaiveIBPConversionFailed"|>]
    ];
   raw = dsTreeNaiveDEFromIBP[
     formulaData,
     variables,
     ProgressReporting -> OptionValue[ProgressReporting]
     ];
   If[AssociationQ[raw], dsTreePublicResult018[raw, context, OptionValue[AuditLevel]], raw]
   ];


Options[DSTreeDLogDE] = {AuditLevel -> "standard"};


DSTreeDLogDE[context_Association, opts : OptionsPattern[]] /; dsContextQ[context] := Module[{raw},
   If[treeFormulaMasslessPendingQ018[context],
    Return[treeFormulaPendingRederivation018["DSTreeDLogDE", context]]
    ];
   raw = dsTreeMultiSectorDLog[context, Automatic, OptionValue[AuditLevel]];
   If[AssociationQ[raw], dsTreePublicResult018[raw, context, OptionValue[AuditLevel]], raw]
   ];


DSTreeDLogDE[context_Association, seedData_, opts : OptionsPattern[]] /; dsContextQ[context] := Module[{raw},
   If[treeFormulaMasslessPendingQ018[context],
    Return[treeFormulaPendingRederivation018["DSTreeDLogDE", context]]
    ];
   raw = dsTreeMultiSectorDLog[context, seedData, OptionValue[AuditLevel]];
   If[AssociationQ[raw], dsTreePublicResult018[raw, context, OptionValue[AuditLevel]], raw]
   ];
