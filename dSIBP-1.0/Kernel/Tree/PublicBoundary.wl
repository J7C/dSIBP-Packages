(* ::Package:: *)
(* 本模块把论文 vertex-basis 公式实现限制在 Private 内部，并为 naive tree 公开接口提供
   J[sectorKey,timeShifts,stateBits] 适配。它不改变 massive-only 递推或 dlog 公式，也不绕过
   massless quotient 的 PendingRederivation 门禁。 *)

(* ::Chapter:: *)
(*统一 J 与私有 vertex basis 的双向映射*)

(* 公式内核按顶点保存 massive 状态；先构造既有 line-packed producer 对象，再通过
   中央 registry 唯一转换为公开 time-only 表示。 *)
dsTreeFormulaIntegralToPublic[int : J[_List], sectorKey_String, context_Association] := Module[
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
      Lookup[line, "state", "full"] === "pinched" || Lookup[line, "packType", ""] === "pinched",
      True,
      If[lineIndexedPowerQ[line], {Lookup[baseline, line["id"], 0]}, {fixedLineSentinel[]}],
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
       Prepend[states, fixedLineSentinel[]]
       ]
      ],
     {e, topo["nE"]}
     ];
   internalIntegral = J[aList, linePacks, {}];
   dsTimeOnlyInternalIntegralToPublic[internalIntegral, context]
   ];


dsTreePublicIntegralToFormula[
   int : J[_, _, _], sectorKey_String, context_Association
   ] := Module[
   {internalIntegral, metadata, familyContext, family, topo, aList, linePacks, lineIds, vertexPacks,
    legId, linePosition, endpointSlot, statePosition, states},
   internalIntegral = dsTimeOnlyPublicIntegralToInternal[int, context];
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


dsTreeTaggedTermToFormula[term_Association, context_Association] := Module[
   {sectorKey, publicIntegral, formulaIntegral},
   sectorKey = Lookup[term, "sectorKey", Missing["NoSector"]];
   publicIntegral = Lookup[term, "integral", Missing["NoIntegral"]];
   If[! StringQ[sectorKey] || ! MatchQ[publicIntegral, J[_, _, _]], Return[$Failed]];
   formulaIntegral = dsTreePublicIntegralToFormula[publicIntegral, sectorKey, context];
   If[formulaIntegral === $Failed, Return[$Failed]];
   Join[KeyDrop[term, {"integral"}], <|"integral" -> formulaIntegral|>]
   ];


(* ::Chapter:: *)
(*Naive 结果的公开结构化映射*)

dsTreeMasterToFormula[master_Association, context_Association] := Module[{term},
   term = dsTreeTaggedTermToFormula[
     KeyTake[master, {"sectorKey", "integral", "coefficient"}],
     context
     ];
   If[term === $Failed, $Failed,
    Join[KeyDrop[master, {"integral"}], <|"integral" -> term["integral"]|>]
    ]
   ];


(* 公开 naive-IBP 结果中的 lhs/rhsTerms 已经是统一 time-only J；本适配器只恢复
   Private 公式内核所需的 sector-tagged vertex basis，并保留全部系数与 sectorKey。 *)
dsTreeNaiveReductionToFormula[record_Association, context_Association] := Module[
   {lhs, rhsTerms},
   lhs = dsTreeTaggedTermToFormula[Lookup[record, "lhs", <||>], context];
   rhsTerms = dsTreeTaggedTermToFormula[#, context] & /@ Lookup[record, "rhsTerms", {}];
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
dsTreeNaiveIBPDataToFormula[data_Association, context_Association] := Module[
   {masters, reductions},
   masters = dsTreeMasterToFormula[#, context] & /@ Lookup[data, "masters", {}];
   reductions = dsTreeNaiveReductionToFormula[#, context] & /@ Lookup[data, "reductions", {}];
   If[masters === {} || reductions === {} || MemberQ[Join[masters, reductions], $Failed],
    Return[$Failed]
    ];
   Join[
    KeyDrop[data, {"masters", "reductions", "bareMasters", "publicRepresentationAudit"}],
    <|"masters" -> masters, "reductions" -> reductions, "representation" -> "J[vertexPacks]"|>
    ]
   ];


dsTreePublicizeData[value_, context_Association, sectorHint_: None] := Module[
   {hint = sectorHint, mapped, publicIntegral},
   Which[
    MatchQ[value, dsTreeToken[_String, J[_List]]],
    publicIntegral = dsTreeFormulaIntegralToPublic[value[[2]], value[[1]], context];
    If[publicIntegral === $Failed, value, publicIntegral],

    MatchQ[value, J[_List]] && StringQ[hint],
    publicIntegral = dsTreeFormulaIntegralToPublic[value, hint, context];
    If[publicIntegral === $Failed, value, publicIntegral],

    AssociationQ[value],
    hint = Lookup[value, "sectorKey", Lookup[value, "sector", hint]];
    mapped = Association@KeyValueMap[
       Function[{key, item}, key -> dsTreePublicizeData[item, context, hint]],
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

    ListQ[value], dsTreePublicizeData[#, context, hint] & /@ value,
    Head[value] === Rule,
    Rule @@ (dsTreePublicizeData[#, context, hint] & /@ List @@ value),
    True, value
    ]
   ];


dsTreePublicResult[result_Association, context_Association, auditLevel_: "standard"] := Module[
   {public, masters, legacyIntegrals, oldPublicLeakQ, privateTokens, representationAudit},
   public = dsTreePublicizeData[result, context];
   If[ListQ[Lookup[public, "masters", None]],
    masters = public["masters"];
    If[And @@ (AssociationQ /@ masters),
     public = Join[public, <|"bareMasters" -> Lookup[masters, "integral", {}]|>]
     ]
    ];
   If[ListQ[Lookup[public, "pinchedMaps", None]],
    public = Join[public, <|"sourceEquations" -> Lookup[public["pinchedMaps"], "rowsByVertex", {}]|>]
    ];
   representationAudit = If[
      auditLevel === "full",
      legacyIntegrals = DeleteDuplicates[Cases[public, item : J[_List] :> item, Infinity]];
      oldPublicLeakQ = dsTimeOnlyOldIntegralLeakQ[public];
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
   If[treeFormulaMasslessPendingQ[context],
    Return[treeFormulaPendingRederivation["DSTreeNaiveIBP", context]]
    ];
   formulaMasters = If[masters === Automatic,
     Automatic,
     dsTreeMasterToFormula[#, context] & /@ masters
     ];
   If[ListQ[formulaMasters] && MemberQ[formulaMasters, $Failed],
    Return[<|"status" -> "failed", "reason" -> "invalidUnifiedMasters"|>]
    ];
   raw = dsTreeNaiveIBPRaw[
     context,
     formulaMasters,
     AuditLevel -> OptionValue[AuditLevel],
     ProgressReporting -> OptionValue[ProgressReporting]
     ];
   If[AssociationQ[raw], dsTreePublicResult[raw, context, OptionValue[AuditLevel]], raw]
   ];


DSTreeNaiveDE[context_Association, variables_: Automatic, masters_: Automatic, OptionsPattern[]] /;
   dsContextQ[context] := Module[{formulaMasters, raw},
   If[treeFormulaMasslessPendingQ[context],
    Return[treeFormulaPendingRederivation["DSTreeNaiveDE", context]]
    ];
   formulaMasters = If[masters === Automatic,
     Automatic,
     dsTreeMasterToFormula[#, context] & /@ masters
     ];
   If[ListQ[formulaMasters] && MemberQ[formulaMasters, $Failed],
    Return[<|"status" -> "failed", "reason" -> "invalidUnifiedMasters"|>]
    ];
    raw = dsTreeNaiveDERaw[
     context,
     variables,
      formulaMasters,
      AuditLevel -> OptionValue[AuditLevel],
      ProgressReporting -> OptionValue[ProgressReporting]
     ];
   If[AssociationQ[raw], dsTreePublicResult[raw, context, OptionValue[AuditLevel]], raw]
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
   formulaData = dsTreeNaiveIBPDataToFormula[ibpData, context];
   If[formulaData === $Failed,
    Return[<|"status" -> "failed", "reason" -> "publicNaiveIBPConversionFailed"|>]
    ];
   raw = dsTreeNaiveDEFromIBP[
     formulaData,
     variables,
     ProgressReporting -> OptionValue[ProgressReporting]
     ];
   If[AssociationQ[raw], dsTreePublicResult[raw, context, OptionValue[AuditLevel]], raw]
   ];
