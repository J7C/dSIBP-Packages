(* ::Package:: *)

(* ::Chapter:: *)
(*tree sector-tagged 数据边界*)

(* 裸 J[vertexPacks] 继续作为公开公式表示；进入跨 sector linearData 时必须由后续 token 同时携带 sectorKey。
   loop-to-tree 投影系数始终复用冻结核心的完整物理幂次公式，不能从 tree pack 反推零点。 *)

treeTaggedIntegral[sectorKey_, integral_J, coefficient_: 1] := <|
   "sectorKey" -> sectorKey,
   "integral" -> integral,
   "coefficient" -> coefficient
   |>;

dsExpandedTerms[expr_] := If[Head[Expand[expr]] === Plus, List @@ Expand[expr], {Expand[expr]}];

dsLoopTermTreeTag[term_, topo_Association, referenceInt_J] := Module[
   {loopIntegrals, loopIntegral, loopCoefficient, targetTopology, referenceTopology, projected, treeIntegrals,
    treeIntegral, projectionCoefficient, targetVertices, referenceVertices, targetAInteger, referenceAInteger,
    targetAZeroPoint, referenceAZeroPoint, targetAPhysical, referenceAPhysical, targetBInteger,
    referenceBInteger, targetBZeroPoint, referenceBZeroPoint, targetBPhysical, referenceBPhysical,
     deltaTimePower, deltaLinePowers, explicitMomentumMagnitudePowers, lineMomentumMagnitudes, auditedProjectionCoefficient},
   loopIntegrals = DeleteDuplicates[Cases[term, J[_, _, _], {0, Infinity}]];
   If[Length[loopIntegrals] =!= 1, Return[$Failed]];
   loopIntegral = First[loopIntegrals];
   loopCoefficient = term /. loopIntegral -> 1;
   targetTopology = loopTreeTargetTopology[loopIntegral, topo];
   referenceTopology = loopTreeTargetTopology[referenceInt, topo];
   projected = projectLoopIntegralToTree[loopIntegral, topo, referenceInt];
   If[MemberQ[{targetTopology, referenceTopology, projected}, $Failed], Return[$Failed]];
   treeIntegrals = DeleteDuplicates[Cases[projected, J[_List], {0, Infinity}]];
   If[Length[treeIntegrals] =!= 1, Return[$Failed]];
   treeIntegral = First[treeIntegrals];
   projectionCoefficient = projected /. treeIntegral -> 1;
   targetVertices = activeAVertexIds[targetTopology];
   referenceVertices = activeAVertexIds[referenceTopology];
   targetAInteger = First[loopIntegral];
   referenceAInteger = First[referenceInt];
   targetAZeroPoint = vertexZeroPoint[targetTopology, #] & /@ targetVertices;
   referenceAZeroPoint = vertexZeroPoint[referenceTopology, #] & /@ referenceVertices;
   targetAPhysical = MapThread[Plus, {targetAInteger, targetAZeroPoint}];
   referenceAPhysical = MapThread[Plus, {referenceAInteger, referenceAZeroPoint}];
   targetBInteger = Table[lineIntegerPowerIndex[targetTopology, loopIntegral, e], {e, targetTopology["nE"]}];
   referenceBInteger = Table[lineIntegerPowerIndex[referenceTopology, referenceInt, e], {e, referenceTopology["nE"]}];
   targetBZeroPoint = Table[
     If[actualLinePackType[targetTopology, e, loopIntegral[[2, e]]] === "shrunk",
      lineBSZeroPoint[targetTopology, e],
      lineBZeroPoint[targetTopology, e]
      ],
     {e, targetTopology["nE"]}
     ];
   referenceBZeroPoint = Table[
     If[actualLinePackType[referenceTopology, e, referenceInt[[2, e]]] === "shrunk",
      lineBSZeroPoint[referenceTopology, e],
      lineBZeroPoint[referenceTopology, e]
      ],
     {e, referenceTopology["nE"]}
     ];
   targetBPhysical = MapThread[Plus, {targetBInteger, targetBZeroPoint}];
   referenceBPhysical = MapThread[Plus, {referenceBInteger, referenceBZeroPoint}];
   deltaTimePower = Expand[Total[targetAPhysical] - Total[referenceAPhysical]];
   deltaLinePowers = MapThread[Expand[#1 - #2] &, {targetBPhysical, referenceBPhysical}];
   explicitMomentumMagnitudePowers = -deltaLinePowers;
   lineMomentumMagnitudes = Table[
     lineMomentumMagnitude[targetTopology, e],
     {e, targetTopology["nE"]}
     ];
   (* 两侧积分都使用 (-tau)^A；deltaTimePower 只作幂次审计，不产生重复的 (-1)^DeltaA。 *)
   auditedProjectionCoefficient = Times @@ MapThread[Power, {lineMomentumMagnitudes, explicitMomentumMagnitudePowers}];
   <|
    "sectorKey" -> sectorKeyFromShrunkLines[targetTopology, Lookup[targetTopology, "sectorShrunkLines", {}]],
    "integral" -> treeIntegral,
    "coefficient" -> Expand[loopCoefficient projectionCoefficient],
    "sourceLoopIntegral" -> loopIntegral,
    "projectionCoefficient" -> projectionCoefficient,
    "physicalPowerAudit" -> <|
      "target" -> <|
        "sectorKey" -> sectorKeyFromShrunkLines[targetTopology, Lookup[targetTopology, "sectorShrunkLines", {}]],
        "aInteger" -> targetAInteger,
        "aZeroPoint" -> targetAZeroPoint,
        "aPhysical" -> targetAPhysical,
        "bInteger" -> targetBInteger,
        "bZeroPoint" -> targetBZeroPoint,
        "bPhysical" -> targetBPhysical,
        "treeNu0" -> targetAZeroPoint
        |>,
      "reference" -> <|
        "sectorKey" -> sectorKeyFromShrunkLines[referenceTopology, Lookup[referenceTopology, "sectorShrunkLines", {}]],
        "aInteger" -> referenceAInteger,
        "aZeroPoint" -> referenceAZeroPoint,
        "aPhysical" -> referenceAPhysical,
        "bInteger" -> referenceBInteger,
        "bZeroPoint" -> referenceBZeroPoint,
        "bPhysical" -> referenceBPhysical
        |>,
      "deltaTimePower" -> deltaTimePower,
      "deltaLineIntegerPowers" -> MapThread[Expand[#1 - #2] &, {targetBInteger, referenceBInteger}],
      "deltaLineZeroPointPowers" -> MapThread[Expand[#1 - #2] &, {targetBZeroPoint, referenceBZeroPoint}],
      "deltaLinePhysicalPowers" -> deltaLinePowers,
      "explicitMomentumMagnitudePowers" -> explicitMomentumMagnitudePowers,
      "lineMomentumMagnitudes" -> lineMomentumMagnitudes,
      "auditedProjectionCoefficient" -> auditedProjectionCoefficient,
      "projectionCoefficientMatchesAudit" -> TrueQ[projectionCoefficient === auditedProjectionCoefficient],
      "zeroPointRules" -> targetTopology["zeroPointRules"],
      "unsafePowerExpand" -> False
      |>
    |>
   ];

dsCombineTreeTaggedTerms[terms_List] := Map[
   Function[group,
    Module[{optionalKeys, optionalData},
     optionalKeys = Select[
       {"sourceLoopIntegral", "projectionCoefficient", "physicalPowerAudit"},
       Function[key, AnyTrue[group, Function[item, KeyExistsQ[item, key]]]]
       ];
     optionalData = Association@Map[
        Switch[#,
          "sourceLoopIntegral", "sourceLoopIntegrals",
          "projectionCoefficient", "projectionCoefficients",
          "physicalPowerAudit", "physicalPowerAudits"
          ] -> Lookup[group, #] &,
        optionalKeys
        ];
     Join[
      First[group],
      <|
       "coefficient" -> Total[Lookup[group, "coefficient"]],
       "contributions" -> (KeyTake[#, Prepend[optionalKeys, "coefficient"]] & /@ group)
       |>,
      optionalData
      ]
     ]
    ],
   GatherBy[terms, {Lookup[#, "sectorKey"], Lookup[#, "integral"]} &]
   ];

dsTreeLinearData[record_Association, topo_Association, referenceInt_J] := Module[{taggedTerms, combined},
   taggedTerms = dsLoopTermTreeTag[#, topo, referenceInt] & /@ dsExpandedTerms[record["loopSeed"]];
   If[MemberQ[taggedTerms, $Failed], Return[<|"status" -> "failed", "reason" -> "treeTagProjectionFailed"|>]];
   combined = dsCombineTreeTaggedTerms[taggedTerms];
   <|
    "status" -> "generated",
    "terms" -> combined,
    "termCount" -> Length[combined],
    "sectorKeys" -> DeleteDuplicates[Lookup[combined, "sectorKey"]],
    "expression" -> Total[(Lookup[#, "coefficient"] Lookup[#, "integral"]) & /@ combined],
    "referenceLoopIntegral" -> referenceInt,
    "coefficientConvention" -> "complete physical powers: a+a0 and b+b0 or bS+bS0"
    |>
   ];

dsEnrichTreeSeedRecord[record_Association, topo_Association, referenceInt_J] := Join[
   record,
   <|"treeLinearData" -> dsTreeLinearData[record, topo, referenceInt]|>
   ];


(* ::Chapter:: *)
(*Sector-tagged tree 迭代*)

(* dsTreeToken 只在 Private context 中承载 sector 身份；对外结果仍序列化为 Association 与统一 Head J。 *)

dsTreeLinearDataQ[data_] := AssociationQ[data] &&
   Lookup[data, "status", Missing["status"]] === "generated" &&
   ListQ[Lookup[data, "terms", Missing["terms"]]] &&
   And @@ (AssociationQ[#] && StringQ[Lookup[#, "sectorKey", None]] && MatchQ[Lookup[#, "integral", None], _J] & /@ data["terms"]);


dsTreeFamilyBySector[sectorKey_String, context_Association] := Module[{matches},
   matches = Select[Lookup[context, "families", {}], Lookup[#, "sector", None] === sectorKey &];
   If[Length[matches] === 1, First[matches], Missing["TreeSectorFamily", sectorKey, Length[matches]]]
   ];


dsTreeTokenExpression[data_Association] /; dsTreeLinearDataQ[data] := Total[
   Lookup[#, "coefficient"] dsTreeToken[Lookup[#, "sectorKey"], Lookup[#, "integral"]] & /@ data["terms"]
   ];


dsTreeTokenTerms[expr_] := Module[{rawTerms, tagged},
   rawTerms = dsExpandedTerms[expr];
   tagged = Map[
     Function[term,
      With[{tokens = DeleteDuplicates[Cases[term, token : dsTreeToken[_String, _J] :> token, {0, Infinity}]]},
       If[Length[tokens] =!= 1,
        $Failed,
        treeTaggedIntegral[tokens[[1, 1]], tokens[[1, 2]], Expand[term /. tokens[[1]] -> 1]]
        ]
       ]
     ],
     rawTerms
     ];
   If[MemberQ[tagged, $Failed], Return[$Failed]];
   Map[
    Function[group, Join[First[group], <|"coefficient" -> Total[Lookup[group, "coefficient"]]|>]],
    GatherBy[tagged, {Lookup[#, "sectorKey"], Lookup[#, "integral"]} &]
    ]
   ];


dsTreeRecordTokenEquation[record_Association] := Module[{linearData = Lookup[record, "treeLinearData", Missing["treeLinearData"]]},
   If[dsTreeLinearDataQ[linearData], dsTreeTokenExpression[linearData], $Failed]
   ];


dsMakeTreeTaggedTimeReductionRules[records_List, data_Association] := Module[
   {good, grouped, minusRules = {}, plusRules = {}, sourceQ = False, sectorKey, vertexId, vertexIndex,
    vertex, states, ordered, currentInts, minusInts, currentTokens, minusTokens, equations,
    m1, m0, invM1, invM0t, tp, tpInv, regular, source, lowerRhs, upperRhs, singularQ = False},
   good = Select[Flatten[records], AssociationQ[#] && Lookup[#, "status", "error"] === "generated" &&
       dsTreeLinearDataQ[Lookup[#, "treeLinearData", <||>]] &];
   grouped = GroupBy[good, treeSeedRuleGroupKey[#, data] &];
   sectorKey = data["sector"];
   KeyValueMap[
    Function[{key, group},
     If[Head[key] === Missing, Return[]];
     vertexId = key[[1]];
     vertexIndex = First@FirstPosition[data["vertexOrder"], vertexId];
     vertex = data["vertices"][[vertexIndex]];
     states = treeBinaryStates[vertex["p"]];
     ordered = SortBy[group, Function[record,
        treeStateIndex[Rest[First[First[Cases[record["treeIntegral"], _J, {0, Infinity}]]][[vertexIndex]]]]
        ]];
     If[Length[ordered] =!= Length[states],
      Message[makeTreeTimeReductionRules::incomplete, dsTreeTimeGroupSentence[<|"key" -> key, "expected" -> Length[states], "actual" -> Length[ordered]|>]];
      Return[]
      ];
     currentInts = First[Cases[#"treeIntegral", _J, {0, Infinity}]] & /@ ordered;
     minusInts = MapThread[
       Function[{int, state}, J[ReplacePart[First[int], vertexIndex -> Prepend[state, First[int][[vertexIndex, 1]] - 1]]]],
       {currentInts, states}
       ];
     currentTokens = dsTreeToken[sectorKey, #] & /@ currentInts;
     minusTokens = dsTreeToken[sectorKey, #] & /@ minusInts;
     equations = dsTreeRecordTokenEquation /@ ordered;
     If[MemberQ[equations, $Failed], Return[]];
     m1 = treeM1[vertex, vertex["nu0"] + First[currentInts[[1]]][[vertexIndex, 1]]];
     m0 = treeM0[vertex];
     invM1 = treeDiagonalInverse[m1];
     invM0t = treeDiagonalInverse[treeM0Tilde[vertex]];
     If[MemberQ[{invM1, invM0t}, $Failed], singularQ = True; Return[]];
     tp = treeTp[vertex];
     tpInv = treeTpInverse[vertex];
     regular = m1.minusTokens + m0.currentTokens;
     source = Expand[equations - regular];
     If[! TrueQ[source === ConstantArray[0, Length[source]]], sourceQ = True];
     lowerRhs = Expand[treeAminusMatrix[vertex, First[currentInts[[1]]][[vertexIndex, 1]]] . currentTokens - invM1.source];
     upperRhs = Expand[treeAplusMatrix[vertex, First[currentInts[[1]]][[vertexIndex, 1]] - 1] . minusTokens - tpInv.invM0t.tp.source];
     minusRules = Join[minusRules, Thread[minusTokens -> lowerRhs]];
     plusRules = Join[plusRules, Thread[currentTokens -> upperRhs]];
     ],
    grouped
    ];
   <|
    "status" -> Which[singularQ, "singular", Length[grouped] === 0, "empty", True, "generated"],
    "minus" -> DeleteDuplicates[minusRules],
    "plus" -> DeleteDuplicates[plusRules],
    "sourceQ" -> sourceQ,
    "groupCount" -> Length[grouped]
    |>
   ];


dsTreeSeedRecordFromSector[vertex_, int : J[_, _, _], family_Association] := Module[
   {sectorTopology = family["topology"], rootTopology, loopSeed, treeSeed, treeIntegral, record},
   rootTopology = Lookup[family, "rootTopology", sectorTopology];
   loopSeed = dtau[vertex, int, sectorTopology];
   If[loopSeed === $Failed, Return[$Failed]];
   treeSeed = projectLoopTimeEquationToTree[loopSeed, rootTopology, int];
   treeIntegral = projectLoopIntegralToTree[int, rootTopology, int];
   If[MemberQ[{treeSeed, treeIntegral}, $Failed], Return[$Failed]];
   record = <|
     "status" -> "generated",
     "generator" -> dtau[vertex],
     "loopSeed" -> loopSeed,
     "treeSeed" -> treeSeed,
     "treeIntegral" -> treeIntegral,
     "referenceA" -> First[int],
     "sectorKey" -> family["sector"]
     |>;
   dsEnrichTreeSeedRecord[record, rootTopology, int]
   ];


dsTreeTaggedSourceAwareStep[
   token : dsTreeToken[sectorKey_String, int_J], vertexIndex_Integer, endpoint_Integer,
   family_Association, familyContext_Association
   ] := Module[
   {packs = First[int], current, seedA, states, vertexId, localIntegrals, records, ruleData, rules, result},
   current = packs[[vertexIndex, 1]];
   If[current === endpoint, Return[token]];
   seedA = If[current < endpoint, current + 1, current];
   states = treeBinaryStates[family["vertices"][[vertexIndex, "p"]]];
   vertexId = family["vertexOrder"][[vertexIndex]];
   localIntegrals = J[ReplacePart[packs, vertexIndex -> Prepend[#, seedA]]] & /@ states;
   records = dsDirectTreeSeedRecord[vertexId, #, family, familyContext] & /@ localIntegrals;
   If[! FreeQ[records, $Failed], Return[$Failed]];
   ruleData = dsMakeTreeTaggedTimeReductionRules[records, family];
   If[Lookup[ruleData, "status", "failed"] === "singular",
    Return[Failure["TreeRecurrenceSingular", <|"sector" -> sectorKey, "integral" -> int,
       "vertex" -> vertexId,
       "MessageTemplate" ->
        "tree 时间递推在所列 sector、积分和顶点处遇到奇异层：该层的两项递推矩阵有对角元为零，无法求逆，因此闭式递推不能继续。请把该顶点的时间幂移开这一层，或对它加入解析正规化后重新生成；也可以改用 general time IBP seed 路线，不依赖这条闭式递推。 The tree time recurrence hit a singular layer at the listed sector, integral and vertex: a diagonal entry of the two-term recurrence matrix vanishes there, so it cannot be inverted and the closed recurrence cannot continue. Move the time power of that vertex off this layer, or add an analytic regulator to it and regenerate; alternatively use the general time IBP seed route, which does not rely on this closed recurrence."|>]]
    ];
   rules = Lookup[ruleData, If[current < endpoint, "minus", "plus"], {}];
   result = Replace[token, rules];
   If[result === token,
    Message[makeTreeFamilyData::badinput, dsTreeIssueSentence[{<|"code" -> "missingSectorTaggedSourceStep", "sector" -> sectorKey,
       "integral" -> int, "vertex" -> vertexId|>}]];
    $Failed,
    result
    ]
   ];


dsTreeTaggedSingleStep[
   token : dsTreeToken[sectorKey_String, int_J], vertexIndex_Integer, endpoint_Integer,
   family_Association, familyContext_Association
   ] := Module[{bareResult},
   If[TrueQ[Lookup[family, "requiresSourceRules", False]] && AssociationQ[Lookup[family, "topology", Missing["NoLoopTopology"]]],
    Return[dsTreeTaggedSourceAwareStep[token, vertexIndex, endpoint, family, familyContext]]
    ];
   bareResult = treeSingleStepIntegral[int, vertexIndex, endpoint, family];
   If[bareResult === $Failed, $Failed, bareResult /. item_J :> dsTreeToken[sectorKey, item]]
   ];


dsTreeTokenDistanceData[token : dsTreeToken[sectorKey_String, int_J], end_, familyContext_Association] := Module[
   {family, endpoints},
   family = dsTreeFamilyBySector[sectorKey, familyContext];
   If[Head[family] === Missing || ! treeIntegralQ[int, family] ||
     ! And @@ (IntegerQ /@ First[int][[All, 1]]), Return[$Failed]];
   endpoints = normalizeTreeEndpoints[treeSectorEndpoints[end, family, familyContext], family];
   If[endpoints === $Failed, Return[$Failed]];
   <|"token" -> token, "sectorKey" -> sectorKey, "endpoints" -> endpoints,
    "remainingThetaLines" -> Length[thetaFullLineIndices[family["topology"]]],
    "distance" -> treeEndpointDistance[int, endpoints],
    "progressKey" -> {Length[thetaFullLineIndices[family["topology"]]], treeEndpointDistance[int, endpoints]}|>
   ];


(* tagged linearData 允许递推跨 contact sector；每条跨 sector 依赖仍须严格降低各自 endpoint 距离。 *)
dsTreeTaggedStepProgress[source_dsTreeToken, replacement_, end_, familyContext_Association] := Module[
   {sourceData, targets, targetData, invalidTargets},
   sourceData = dsTreeTokenDistanceData[source, end, familyContext];
   targets = DeleteDuplicates[Cases[replacement, token : dsTreeToken[_String, _J] :> token, {0, Infinity}]];
   targetData = dsTreeTokenDistanceData[#, end, familyContext] & /@ targets;
   invalidTargets = Pick[targets, (# === $Failed & /@ targetData)];
   <|
    "passQ" -> TrueQ[sourceData =!= $Failed && invalidTargets === {} &&
       And @@ (treeProgressKeyLessQ[Lookup[#, "progressKey", {Infinity, Infinity}], sourceData["progressKey"]] & /@ targetData)],
    "source" -> sourceData,
    "targets" -> targetData,
    "invalidTargets" -> invalidTargets
    |>
   ];


Options[dsRepIterativeTreeLinearData] = Options[repIterativeData];


dsRepIterativeTreeLinearData[data_Association, end_: Automatic, context_Association, OptionsPattern[]] := Module[
   {familyContext, result, steps = 0, tokens, token, sectorKey, int, family, endpoints, vertexIndex,
    reducedTerms, replacement, progressData, seenStates = <||>, stateHash},
   If[! dsTreeLinearDataQ[data], Return[<|"status" -> "error", "reason" -> "invalidTreeLinearData"|>]];
   familyContext = dsTreeFamilyContext[context];
   result = Expand[dsTreeTokenExpression[data]];
   AssociateTo[seenStates, treeRecurrenceStateHash[result] -> True];
   While[True,
    tokens = DeleteDuplicates[Cases[result, token : dsTreeToken[_String, _J] :> token, {0, Infinity}]];
    token = SelectFirst[
      tokens,
      Function[item,
       family = dsTreeFamilyBySector[item[[1]], familyContext];
       If[Head[family] === Missing, True,
        endpoints = normalizeTreeEndpoints[treeSectorEndpoints[end, family, familyContext], family];
        endpoints === $Failed || Or @@ MapThread[Unequal, {First[item[[2]]][[All, 1]], endpoints}]
        ]
       ],
      Missing["AllAtEndpoints"]
      ];
    If[Head[token] === Missing, Break[]];
    sectorKey = token[[1]];
    int = token[[2]];
    family = dsTreeFamilyBySector[sectorKey, familyContext];
    If[Head[family] === Missing,
     Message[repIterativeData::nosector, dsTreeNoSectorSentence[<|"sectorKey" -> sectorKey, "integral" -> int|>]];
     Return[<|"status" -> "error", "reason" -> "unknownSector", "steps" -> steps|>]
     ];
    If[! treeIntegralQ[int, family] || ! And @@ (IntegerQ /@ First[int][[All, 1]]),
     Message[repIterativeData::badindex, dsTreeBadIndexSentence[<|"sectorKey" -> sectorKey, "integral" -> int|>]];
     Return[<|"status" -> "error", "reason" -> "invalidTaggedIntegral", "steps" -> steps|>]
     ];
    endpoints = normalizeTreeEndpoints[treeSectorEndpoints[end, family, familyContext], family];
    If[endpoints === $Failed, Return[<|"status" -> "error", "reason" -> "invalidEndpoint", "steps" -> steps|>]];
    vertexIndex = SelectFirst[Range[Length[endpoints]], First[int][[#, 1]] =!= endpoints[[#]] &];
    replacement = dsTreeTaggedSingleStep[token, vertexIndex, endpoints[[vertexIndex]], family, familyContext];
    If[MatchQ[replacement, _Failure] && replacement[[1]] === "TreeRecurrenceSingular",
     Return[<|"status" -> "singular", "reason" -> "singularRecurrence", "steps" -> steps,
       "details" -> replacement[[2]]|>]
     ];
    If[replacement === $Failed, Return[<|"status" -> "error", "reason" -> "taggedStepFailed", "steps" -> steps|>]];
    progressData = dsTreeTaggedStepProgress[token, replacement, end, familyContext];
    If[! TrueQ[progressData["passQ"]],
     Message[repIterativeData::noprogress, dsTreeIssueSentence[<|"code" -> "noprogress", "integral" -> int, "steps" -> steps|>]];
     Return[<|"status" -> "error", "reason" -> "nonDecreasingRecurrence", "steps" -> steps,
       "progressData" -> progressData|>]
     ];
    result = Expand[result /. token -> replacement];
    stateHash = treeRecurrenceStateHash[result];
    If[KeyExistsQ[seenStates, stateHash],
     Message[repIterativeData::cycle, dsTreeIssueSentence[<|"code" -> "cycle", "steps" -> steps|>]];
     Return[<|"status" -> "error", "reason" -> "recurrenceCycle", "steps" -> steps,
       "stateHash" -> stateHash|>]
     ];
    AssociateTo[seenStates, stateHash -> True];
    steps++;
    ];
   reducedTerms = dsTreeTokenTerms[result];
   If[reducedTerms === $Failed, Return[<|"status" -> "error", "reason" -> "taggedTermExtractionFailed", "steps" -> steps|>]];
   <|
    "status" -> "reduced",
    "terms" -> reducedTerms,
    "termCount" -> Length[reducedTerms],
    "sectorKeys" -> DeleteDuplicates[Lookup[reducedTerms, "sectorKey"]],
    "expression" -> Total[Lookup[#, "coefficient"] Lookup[#, "integral"] & /@ reducedTerms],
    "steps" -> steps,
    "endpoints" -> end,
    "inputTreeLinearData" -> data,
    "coefficientConvention" -> Lookup[data, "coefficientConvention", Missing["coefficientConvention"]]
    |>
   ];


(* ::Chapter:: *)
(*Naive tree time-IBP 线性系统*)

(* 这一路只把投影后的 dtau 当作线性方程求解，不调用 A-/A+ 递推或 dlog 公式。
   master record 必须携带 sectorKey 与 normalization，避免同 shape lower sector 混淆。 *)

Options[DSTreeNaiveIBP] = {AuditLevel -> "standard", ProgressReporting -> Automatic};

DSTreeNaiveIBP::badmasters = "tree naive IBP 需要非空、无重复且可唯一匹配 sector 的 tagged master 列表。";
DSTreeNaiveIBP::nonsquare = "tree naive IBP 方程数 `1` 与待约化对象数 `2` 不相等。";
DSTreeNaiveIBP::solvefailed = "tree naive IBP 线性系统求解失败。";


dsTreeTaggedMasterQ[record_] := AssociationQ[record] &&
   StringQ[Lookup[record, "sectorKey", None]] &&
   MatchQ[Lookup[record, "integral", None], _J] &&
   ! TrueQ[Lookup[record, "coefficient", 0] === 0];


dsTreeResolveNaiveMasters[Automatic, context_Association] := Module[{dlog = dsTreeMultiSectorDLog[context]},
   If[Lookup[dlog, "status", "error"] === "generated", Lookup[dlog, "masters", $Failed], $Failed]
   ];
dsTreeResolveNaiveMasters[masters_List, _Association] := masters;
dsTreeResolveNaiveMasters[_, _Association] := $Failed;


dsTreeMasterFamilyRecords[masters_List, familyContext_Association] := Map[
   Function[master,
    With[{family = dsTreeFamilyBySector[master["sectorKey"], familyContext]},
     If[Head[family] === Missing || ! treeIntegralQ[master["integral"], family],
      $Failed,
      <|"master" -> master, "family" -> family|>
      ]
     ]
    ],
   masters
   ];


dsTreeNaiveSeedRecords[masterFamilyRecords_List, familyContext_Association, progress_] := Flatten@dsProgressMap[
    "正在生成 naive tree time-IBP / Building naive-tree time-IBP relations",
    masterFamilyRecords,
    Function[item,
     With[{master = item["master"], family = item["family"]},
      MapIndexed[
       Function[{vertexId, position},
        With[{seedIntegral = J[ReplacePart[
             First[master["integral"]],
             First[position] -> ReplacePart[First[master["integral"]][[First[position]]], 1 -> 1]
             ]]},
         dsDirectTreeSeedRecord[vertexId, seedIntegral, family, familyContext]
         ]
        ],
       family["vertexOrder"]
       ]
      ]
     ],
    progress
    ];


dsTreePublicReductionRecord[rule_Rule] := Module[{lhs, rhsTerms},
   lhs = First[rule];
   rhsTerms = dsTreeTokenTerms[Last[rule]];
   If[! MatchQ[lhs, dsTreeToken[_String, _J]] || rhsTerms === $Failed,
    $Failed,
    <|
     "lhs" -> treeTaggedIntegral[lhs[[1]], lhs[[2]]],
     "rhsTerms" -> rhsTerms,
     "rhsExpression" -> Total[Lookup[#, "coefficient"] Lookup[#, "integral"] & /@ rhsTerms]
     |>
    ]
   ];


dsTreeInternalReductionRules[data_Association] := Map[
   Function[record,
    dsTreeToken[record["lhs", "sectorKey"], record["lhs", "integral"]] ->
     Total[Lookup[#, "coefficient"] dsTreeToken[Lookup[#, "sectorKey"], Lookup[#, "integral"]] & /@ record["rhsTerms"]]
    ],
   Lookup[data, "reductions", {}]
   ];


dsTreeNaiveIBPRaw018[context_Association, masters_: Automatic, OptionsPattern[DSTreeNaiveIBP]] /; dsContextQ[context] := Module[
   {resolvedMasters, familyContext, masterFamilyRecords, seedRecords, equations, masterTokens, allTokens,
    unknownTokens, matrix, constantVector, solution, tokenRules, publicRules, solveResiduals, status},
   If[treeFormulaMasslessPendingQ018[context],
    Return[treeFormulaPendingRederivation018["DSTreeNaiveIBP", context]]
    ];
   resolvedMasters = dsTreeResolveNaiveMasters[masters, context];
   familyContext = dsTreeFamilyContext[context];
   If[! ListQ[resolvedMasters] || resolvedMasters === {} || ! And @@ (dsTreeTaggedMasterQ /@ resolvedMasters) ||
     DuplicateFreeQ[{Lookup[#, "sectorKey"], Lookup[#, "integral"]} & /@ resolvedMasters] =!= True || familyContext === $Failed,
    Message[DSTreeNaiveIBP::badmasters];
    Return[<|"status" -> "failed", "reason" -> "invalidMasters"|>]
    ];
   masterFamilyRecords = dsTreeMasterFamilyRecords[resolvedMasters, familyContext];
   If[MemberQ[masterFamilyRecords, $Failed],
    Message[DSTreeNaiveIBP::badmasters];
    Return[<|"status" -> "failed", "reason" -> "masterSectorMismatch", "masters" -> resolvedMasters|>]
    ];
   seedRecords = dsTreeNaiveSeedRecords[masterFamilyRecords, familyContext, OptionValue[ProgressReporting]];
   If[MemberQ[seedRecords, $Failed] || AnyTrue[seedRecords, Lookup[#, "status", "error"] =!= "generated" &],
    Return[<|"status" -> "failed", "reason" -> "seedGenerationFailed", "seedRecords" -> seedRecords|>]
    ];
   equations = DeleteDuplicates[dsTreeRecordTokenEquation /@ seedRecords];
   If[MemberQ[equations, $Failed], Return[<|"status" -> "failed", "reason" -> "taggedEquationFailed"|>]];
   masterTokens = dsTreeToken[Lookup[#, "sectorKey"], Lookup[#, "integral"]] & /@ resolvedMasters;
   allTokens = DeleteDuplicates[Cases[equations, token : dsTreeToken[_String, _J] :> token, Infinity]];
   unknownTokens = Select[allTokens, ! MemberQ[masterTokens, #] &];
   If[Length[equations] =!= Length[unknownTokens],
    Message[DSTreeNaiveIBP::nonsquare, Length[equations], Length[unknownTokens]];
    Return[<|"status" -> "failed", "reason" -> "nonSquareSystem", "equationCount" -> Length[equations],
      "unknownCount" -> Length[unknownTokens], "masters" -> resolvedMasters, "seedRecords" -> seedRecords|>]
    ];
   matrix = Table[Coefficient[equations[[row]], unknownTokens[[column]]],
     {row, Length[equations]}, {column, Length[unknownTokens]}];
   constantVector = Expand[equations /. Thread[unknownTokens -> 0]];
   solution = Quiet[Check[LinearSolve[matrix, -constantVector], $Failed]];
   If[solution === $Failed,
    Message[DSTreeNaiveIBP::solvefailed];
    Return[<|"status" -> "failed", "reason" -> "linearSolveFailed", "equationCount" -> Length[equations],
      "unknownCount" -> Length[unknownTokens], "masters" -> resolvedMasters, "seedRecords" -> seedRecords|>]
    ];
   tokenRules = Thread[unknownTokens -> solution];
   solveResiduals = If[
     OptionValue[AuditLevel] === "full",
     Together /@ Expand[equations /. tokenRules],
     Missing["NotRunAtStandardAuditLevel"]
     ];
   publicRules = dsTreePublicReductionRecord /@ tokenRules;
   status = If[
     FreeQ[publicRules, $Failed] &&
      (Head[solveResiduals] === Missing || And @@ (TrueQ[# === 0] & /@ solveResiduals)),
     "solved",
     "failed"
     ];
   <|
    "status" -> status,
    "reason" -> If[status === "solved", None, "nonzeroSolveResidual"],
    "masters" -> resolvedMasters,
    "masterCount" -> Length[resolvedMasters],
    "masterSource" -> If[masters === Automatic, "DSTreeDLogDE", "explicit"],
    "equationCount" -> Length[equations],
    "unknownCount" -> Length[unknownTokens],
    "seedRecords" -> seedRecords,
    "reductions" -> publicRules,
    "solveResiduals" -> solveResiduals,
    "auditLevel" -> OptionValue[AuditLevel],
    "context" -> context,
    "equationConvention" -> "projected loop dtau == 0; non-master tagged tree integrals solved in the supplied master basis",
    "formulaRecurrenceUsedQ" -> False
    |>
   ];


DSTreeNaiveIBP[_, ___] := (Message[DSTreeNaiveIBP::badmasters]; <|"status" -> "failed", "reason" -> "invalidContextOrMasters"|>);
