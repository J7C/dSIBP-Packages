(* ::Package:: *)
(* 本模块定义唯一积分表示、sector 身份与 massless seed quotient。所有 root
   line 槽位永久保留；full line 使用三槽，shrunk line 使用单槽，fixed line 用短
   字符串 "F" 标记无整数动量幂。masslessFull 的公开槽仍是 {n1,n2}，但 source
   seed 只枚举 n2->0 的两个代数代表，导数输出再由统一 canonical 映回该代表。 *)

(* ::Chapter:: *)
(*统一 line-pack schema*)

fixedLineSentinel[] := "F";


linePackNPositions[line_Association, packType_String] := Switch[packType,
   "massiveFull" | "massiveCross" | "masslessFull", {2, 3},
   _, {}
   ];


linePackBPosition[line_Association] := If[lineIndexedPowerQ[line], 1, Missing["FixedLinePower"]];


makeLinePack[line_Association] := Module[
   {id = line["id"], indexedQ = lineIndexedPowerQ[line], firstSlot},
   firstSlot = If[indexedQ, b[id], fixedLineSentinel[]];
   Switch[line["packType"],
    "massiveFull" | "massiveCross" | "masslessFull",
    {firstSlot, n[id, 1], n[id, 2]},
    "masslessCross",
    {firstSlot, 0, 0},
    "shrunk",
    If[indexedQ, {bS[id]}, {fixedLineSentinel[]}],
    _, Message[makeLinePack::badtype, line["packType"], id]; $Failed
    ]
   ];


actualLinePackType[topo_Association, e_Integer, pack_List] := Module[
   {line = topo["lines"][[e]], declared},
   declared = line["packType"];
   Which[
    declared === "shrunk", "shrunk",
    Length[pack] === 1 && lineIndexedPowerQ[line] && pack[[1]] =!= fixedLineSentinel[], "shrunk",
    Length[pack] === 1 && ! lineIndexedPowerQ[line] && pack[[1]] === fixedLineSentinel[], "shrunk",
    True, declared
    ]
   ];


lineIntegerPowerIndex[topo_Association, J[_, linePacks_, _], e_Integer] := If[
   lineIndexedPowerQ[topo["lines"][[e]]],
   linePacks[[e, 1]],
   0
   ];


discreteVarsForLine[line_Association] := If[
   Lookup[line, "state", "full"] === "shrunk" || Lookup[line, "packType", ""] === "shrunk",
   {},
   Switch[line["packType"],
    "massiveFull" | "massiveCross" | "masslessFull", {n[line["id"], 1], n[line["id"], 2]},
    _, {}
    ]
   ];


(* 输入输出仍保留两个端点槽；这里只缩减 DSSeeds 的 source representatives。
   massive 需要四态，masslessFull 只生成 00/10，cross 与 shrunk 不生成离散态。
   显式 parity 约束触及的 banana 束成员只生成 00——奇偶位折到束首线，与
   parityToFirstLine canonical 一致；无约束或未触及的束不折叠。 *)
discreteSeedRuleSetsForLine[line_Association, parityFoldedQ_: False] := Module[
   {id = line["id"], variables = discreteVarsForLine[line]},
   Switch[Lookup[line, "packType", ""],
    "masslessFull" /; Lookup[line, "state", "full"] =!= "shrunk",
    If[parityFoldedQ,
     {{n[id, 1] -> 0, n[id, 2] -> 0}},
     {
      {n[id, 1] -> 0, n[id, 2] -> 0},
      {n[id, 1] -> 1, n[id, 2] -> 0}
      }],
    _, ruleSetsForVars[variables]
    ]
   ];


rawBinaryDiscreteStateCountForLine[line_Association] :=
  2^Length[discreteVarsForLine[line]];


discreteStateCountForLine[line_Association] :=
  Length[discreteSeedRuleSetsForLine[line]];


(* enumerateDiscreteStates 是所有普通 time/momentum seed 的共同入口。返回值同时保存
   raw 双端点空间与 quotient representative 空间，供 release/audit 核对缩减边界；
   parity 触及的 banana 束只枚举首线奇偶位（parityFoldedLineIndices）。 *)
enumerateDiscreteStates[expr_, topo_Association] := Module[
   {parityFoldedLines, perLineRuleSets, allRuleSets, rawStateCount, result},
   parityFoldedLines = parityFoldedBundleMemberLineIndices[topo];
   perLineRuleSets = MapIndexed[
     discreteSeedRuleSetsForLine[#1, MemberQ[parityFoldedLines, First[#2]]] &,
     topo["lines"]];
   allRuleSets = Flatten[#, 1] & /@ Tuples[perLineRuleSets];
   rawStateCount = Times @@ (rawBinaryDiscreteStateCountForLine /@ topo["lines"]);
   result = <|
    "rules" -> allRuleSets,
    "integrals" -> (expr /. # & /@ allRuleSets),
    "mode" -> "masslessQuotientRepresentatives",
    "canonicalDirection" -> "n2ToZero",
    "rawBinaryStateCount" -> rawStateCount,
    "representativeStateCount" -> Length[allRuleSets],
    "eliminatedAlgebraicStateCount" -> rawStateCount - Length[allRuleSets]
    |>;
   If[parityFoldedLines =!= {},
    result = Join[result, <|"parityFoldedLineIndices" -> parityFoldedLines|>]];
   result
   ];


selectedDiscreteSeedRules[topo_Association, OptionsPattern[]] := Module[
   {data = enumerateDiscreteStates[makeBaseIntegral[topo], topo]},
   Join[
    <|"status" -> "generated", "ruleCount" -> Length[data["rules"]]|>,
    KeyTake[data, {
      "mode", "canonicalDirection", "rawBinaryStateCount",
      "representativeStateCount", "eliminatedAlgebraicStateCount", "rules"
      }]
    ]
   ];


publicExpectedPackLength[_Association, "shrunk"] := 1;
publicExpectedPackLength[_Association, "massiveFull" | "massiveCross" | "masslessFull" | "masslessCross"] := 3;
publicExpectedPackLength[_Association, packType_String] := Missing["UnknownPackType", packType];


makeBaseIntegral[topo_Association] := Module[
   {active = Lookup[topo, "activeVertexIds", topo["vertexIds"]],
    fixedA = Lookup[topo, "fixedAVertexValues", <||>], ispList},
   ispList = If[
     Lookup[topo, "ibpMode", "full"] === "timeOnly",
     {},
     Table[ispN[j], {j, Length[topo["ispData"]]}]
     ];
   J[
    Table[If[AssociationQ[fixedA] && KeyExistsQ[fixedA, v], fixedA[v], a[v]], {v, active}],
    makeLinePacks[topo],
    ispList
    ]
   ];


(* ::Chapter:: *)
(*Sector pattern 与唯一解析*)

sectorLinePattern[line_Association, pack_List] := <|
   "powerKind" -> If[lineIndexedPowerQ[line], "cycle", "fixed"],
   "state" -> If[Length[pack] === 1, "shrunk", "full"]
   |>;


sectorPattern[topo_Association, linePacks_List] := MapThread[
   sectorLinePattern,
   {topo["lines"], linePacks}
   ];


sectorShrunkLinesFromPattern[pattern_List] := Flatten@Position[Lookup[pattern, "state", {}], "shrunk"];


sectorKeyFromPattern[topo_Association, pattern_List] := sectorKeyFromShrunkLines[
   topo,
   sectorShrunkLinesFromPattern[pattern]
   ];


linePackMatchesSlotQ[pack_List, slot_Association] := Module[
   {template = Lookup[slot, "packTemplate", {}], firstTemplate, fixedQ},
   If[Length[pack] =!= Length[template], Return[False]];
   If[template === {}, Return[pack === {}]];
   firstTemplate = First[template];
   fixedQ = SameQ[firstTemplate, fixedLineSentinel[]];
   TrueQ[
    If[fixedQ, SameQ[First[pack], fixedLineSentinel[]], First[pack] =!= fixedLineSentinel[]]
    ]
   ];


(* ::Chapter:: *)
(*结构化 fixed-line sector prefactor*)

zeroPointRuleValueFromRules[rules_List, symbol_, default_: 0] := Module[{hits},
   hits = Cases[rules, (Rule | RuleDelayed)[lhs_, rhs_] /; lhs === symbol :> rhs];
   If[hits === {}, default, Last[hits]]
   ];


explicitZeroPointRuleValue[rules_List, symbol_] := Module[{hits},
   hits = Cases[rules, (Rule | RuleDelayed)[lhs_, rhs_] /; lhs === symbol :> rhs];
   If[hits === {}, Missing["NoExplicitZeroPoint", symbol], Last[hits]]
   ];


lineDerivedShrinkZeroPoint[topo_Association, e_Integer] := Expand[
   zeroPointRuleValueFromRules[
     Lookup[topo, "rootZeroPointRules", topo["zeroPointRules"]],
     b0[topo["lines"][[e, "id"]]],
     0
     ] + lineShrinkZeroPointShift[topo["lines"][[e]]]
   ];


lineTargetShrinkZeroPoint[topo_Association, e_Integer] := Module[
   {rootRules, id, explicit},
   rootRules = Lookup[topo, "rootZeroPointRules", topo["zeroPointRules"]];
   id = topo["lines"][[e, "id"]];
   explicit = explicitZeroPointRuleValue[rootRules, bS0[id]];
   If[Head[explicit] === Missing, lineDerivedShrinkZeroPoint[topo, e], explicit]
   ];


lineShrinkZeroPointSource[topo_Association, e_Integer] := Module[
   {rootRules, id},
   rootRules = Lookup[topo, "rootZeroPointRules", topo["zeroPointRules"]];
   id = topo["lines"][[e, "id"]];
   If[Head[explicitZeroPointRuleValue[rootRules, bS0[id]]] === Missing,
    "compiledDefault",
    "userBS0Override"
    ]
   ];


(* child sector 继承 root 的显式 bS0 override；未给出时仍严格使用函数系统编译的 z shift。 *)
sectorZeroPointRules[
   topo_Association,
   shrunkLines_List,
   repMap_Association,
   activeVertices_List
   ] := Module[
   {lineRules, vertexRules, classVertices, classShrunkLines, shiftSum, rep},
   lineRules = Table[
     If[MemberQ[shrunkLines, e],
      bS0[topo["lines"][[e, "id"]]] -> lineTargetShrinkZeroPoint[topo, e],
      Switch[actualLinePackType[topo, e, makeLinePack[topo["lines"][[e]]]],
       "shrunk", bS0[topo["lines"][[e, "id"]]] -> lineBSZeroPoint[topo, e],
       _, b0[topo["lines"][[e, "id"]]] -> lineBZeroPoint[topo, e]
       ]
      ],
     {e, topo["nE"]}
     ];
   vertexRules = Table[
     rep = activeVertices[[i]];
     classVertices = Select[topo["vertexIds"], Lookup[repMap, #] === rep &];
     classShrunkLines = Select[
       shrunkLines,
       Function[e, And @@ (Function[v, Lookup[repMap, v] === rep] /@ topo["lines"][[e, "endpoints"]])]
       ];
     shiftSum = Total[lineShrinkZeroPointShift[topo["lines"][[#]]] & /@ classShrunkLines];
     a0[rep] -> Total[vertexZeroPoint[topo, #] & /@ classVertices] - shiftSum,
     {i, Length[activeVertices]}
     ];
   DeleteDuplicates@Join[vertexRules, lineRules]
   ];


independentExternalMagnitudeExpressions[topo_Association] := Module[
   {count, bindingData},
   count = Length[Lookup[topo, "independentExternalMomenta", {}]];
   bindingData = externalLegMagnitudeBindingData[topo];
   Lookup[Take[bindingData, UpTo[count]], "userMagnitudeExpression", {}]
   ];


fixedLineKEIndex[topo_Association, e_Integer] := Module[
   {momentum, declared, matches},
   momentum = Expand[topo["lines"][[e, "momentum"]]];
   declared = Expand /@ Lookup[topo, "independentExternalMomenta", {}];
   matches = Flatten@Position[
      (TrueQ[Expand[# - momentum] === 0] || TrueQ[Expand[# + momentum] === 0]) & /@ declared,
      True
      ];
   If[matches === {}, Missing["NotIndependentExternalMagnitude"], First[matches]]
   ];


lineShrinkNormalizationFactor[topo_Association, e_Integer] := Module[
   {line, terms, coefficients, nonzero, contourSign},
   line = topo["lines"][[e]];
   If[Lookup[line, "massType", "massive"] =!= "massive", Return[1]];
   terms = lineCompiledShrinkTerms[line];
   coefficients = Expand[Lookup[#, "coefficient", 0]] & /@ terms;
   nonzero = Select[coefficients, ! TrueQ[# === 0] &];
   contourSign = sameBranchSKSign[line];
   (* shrinkTerms 的 coefficient 是 -W_T；sector normalization 必须吸收物理 W_T。 *)
   If[nonzero === {}, 1, -contourSign First[nonzero]]
   ];


fixedLinePrefactorRecord[topo_Association, e_Integer, state_: Automatic] := Module[
   {line = topo["lines"][[e]], id, shrunkQ, rootRules, sourcePower, targetPower,
     integerShift, zeroPointShift, targetMinusSource, residualPower, parameter,
     derivedTarget, overrideSource, kEIndex},
   id = line["id"];
   shrunkQ = If[
     state === Automatic,
     Lookup[line, "state", "full"] === "shrunk" || Lookup[line, "packType", ""] === "shrunk",
     TrueQ[state === "shrunk"]
     ];
   rootRules = Lookup[topo, "rootZeroPointRules", topo["zeroPointRules"]];
   sourcePower = zeroPointRuleValueFromRules[rootRules, b0[id], 0];
   (* 本函数也会以 root topology 和目标 J 重建 prefactor；此时 root rules 通常没有 bS0。
      必须沿现行 compiled-default/显式 override 合同求 target，而不能把缺失 bS0 当成零。 *)
   targetPower = If[shrunkQ, lineTargetShrinkZeroPoint[topo, e], lineBZeroPoint[topo, e]];
   integerShift = If[shrunkQ, lineShrinkBShift[line], 0];
   zeroPointShift = If[shrunkQ, lineShrinkZeroPointShift[line], 0];
   derivedTarget = If[shrunkQ, lineDerivedShrinkZeroPoint[topo, e], sourcePower];
   overrideSource = If[shrunkQ, lineShrinkZeroPointSource[topo, e], "sourceSector"];
   targetMinusSource = Expand[targetPower - sourcePower];
   residualPower = Expand[integerShift + zeroPointShift - targetMinusSource];
   parameter = fixedLineMomentumMagnitude[topo, e];
   kEIndex = fixedLineKEIndex[topo, e];
   <|
    "lineIndex" -> e,
    "lineId" -> id,
    "parameterKey" -> ("fixedLine:" <> ToString[id, InputForm]),
     "parameter" -> parameter,
     "kEIndex" -> kEIndex,
    "state" -> If[shrunkQ, "shrunk", "full"],
    "sourceZeroPoint" -> sourcePower,
    "targetZeroPoint" -> targetPower,
    "derivedTargetZeroPoint" -> derivedTarget,
    "zeroPointOverrideSource" -> overrideSource,
    "integerShift" -> integerShift,
    "compiledIntegerShift" -> integerShift,
    "zeroPointShift" -> zeroPointShift,
    "compiledZeroPointShift" -> zeroPointShift,
    "targetMinusSource" -> targetMinusSource,
    "normalizedResidualPower" -> residualPower,
    "prefactorPower" -> Expand[-targetPower]
    |>
   ];


sectorPrefactorData[topo_Association, shrunkLines_: Automatic] := Module[
   {effectiveShrunkLines, records, kEIndices, kEMomenta, kEExpressions, kEPowers,
    residualRecords, contractionRecords},
   effectiveShrunkLines = If[
     shrunkLines === Automatic,
     Select[
       Range[topo["nE"]],
       Lookup[topo["lines"][[#]], "state", "full"] === "shrunk" ||
         Lookup[topo["lines"][[#]], "packType", ""] === "shrunk" &
       ],
     Sort@DeleteDuplicates[shrunkLines]
     ];
   records = fixedLinePrefactorRecord[topo, #] & /@ Select[
      Range[topo["nE"]],
      ! lineIndexedPowerQ[topo["lines"][[#]]] &
      ];
   records = Map[
     fixedLinePrefactorRecord[
       topo,
       Lookup[#, "lineIndex"],
       If[MemberQ[effectiveShrunkLines, Lookup[#, "lineIndex"]], "shrunk", "full"]
       ] &,
     records
     ];
   kEIndices = Range[Length[Lookup[topo, "independentExternalMomenta", {}]]];
   kEMomenta = Lookup[topo, "independentExternalMomenta", {}];
   kEExpressions = independentExternalMagnitudeExpressions[topo];
   kEPowers = Table[
     Total[Lookup[Select[records, SameQ[Lookup[#, "kEIndex", None], i] &], "prefactorPower", {}]],
     {i, Length[kEIndices]}
     ];
   residualRecords = Select[records, Head[Lookup[#, "kEIndex", Missing[]]] === Missing &];
   contractionRecords = Map[
     <|
       "lineIndex" -> #,
       "lineId" -> topo["lines"][[#, "id"]],
       "normalizationFactor" -> lineShrinkNormalizationFactor[topo, #]
       |> &,
     effectiveShrunkLines
     ];
   <|
    "lineIndices" -> Lookup[records, "lineIndex", {}],
    "parameterKeys" -> Lookup[records, "parameterKey", {}],
    "parameterList" -> Lookup[records, "parameter", {}],
    "powerList" -> Lookup[records, "prefactorPower", {}],
     "powerParts" -> records,
     "kEIndices" -> kEIndices,
     "kEMomenta" -> kEMomenta,
     "kEParameterExpressions" -> kEExpressions,
     "kEPower" -> Apply[kEpower, kEPowers],
     "residualPowerParts" -> residualRecords,
     "contractionParts" -> contractionRecords,
     "constantPrefactor" -> Times @@ Lookup[contractionRecords, "normalizationFactor", {}],
     "normalizationConvention" -> "structuralKEPowerAndContact-v1"
     |>
   ];


materializeSectorPrefactor[data_Association] := Module[
    {powerHead, powers, expressions, residualParts},
    powerHead = Lookup[data, "kEPower", Missing["NoStructuralKEPower"]];
    If[Head[powerHead] === Missing,
     Return[Failure["MissingStructuralKEPower", <|"requiredKey" -> "kEPower",
        "MessageTemplate" ->
         "无法构造该 sector 的正规化前置因子，因为传入的 sector 数据缺少结构性 kEPower 记录，只靠表达式无法还原幂向量。请使用当前 context 的 sector 生成入口重新产生这份数据，不要手工拼装或跨 context 复用。 The normalized prefactor of this sector cannot be built because the supplied sector data has no structural kEPower record, so the power vector cannot be recovered from expressions alone. Regenerate the data with the sector producer of the current context instead of assembling it by hand or reusing it across contexts."|>]]
     ];
   powers = List @@ powerHead;
   expressions = Lookup[data, "kEParameterExpressions", {}];
   residualParts = Lookup[data, "residualPowerParts", {}];
   If[Length[powers] =!= Length[expressions], Return[$Failed]];
   Expand[
    Lookup[data, "constantPrefactor", 1]
     Times @@ MapThread[Power, {expressions, powers}]
     Times @@ (Power[Lookup[#, "parameter"], Lookup[#, "prefactorPower"]] & /@ residualParts)
    ]
   ];


sectorPrefactorDataForIntegral[topo_Association, int : J[_, linePacks_List, _]] := Module[
   {metadata, frozenData, pattern, shrunkLines},
   (* time-only 的 sector 身份已经由公开 key 和冻结 metadata 唯一确定；逆转换后的
      旧 pack 只服务 producer，不能再次推断 normalization，否则会丢失 child zero point 幂。 *)
   If[Lookup[topo, "ibpMode", "full"] === "timeOnly",
    metadata = integralSectorMetadata[topo, int];
    frozenData = If[
      AssociationQ[metadata],
      Lookup[metadata, "sectorPrefactorData", Missing["NoSectorPrefactorData"]],
      Missing["NoSectorMetadata"]
      ];
    If[AssociationQ[frozenData], Return[frozenData]]
    ];
   pattern = sectorPattern[topo, linePacks];
   shrunkLines = sectorShrunkLinesFromPattern[pattern];
   sectorPrefactorData[topo, shrunkLines]
   ];


normalizeContactTerm[topo_Association, sourceInt_J, rawTerm_] := Module[
   {targets, targetInt, rawCoefficient, sourceData, targetData, ratio},
   targets = DeleteDuplicates[Cases[rawTerm, _J, {0, Infinity}]];
   If[Length[targets] =!= 1, Return[$Failed]];
   targetInt = First[targets];
   rawCoefficient = Cancel[rawTerm/targetInt];
   sourceData = sectorPrefactorDataForIntegral[topo, sourceInt];
   targetData = sectorPrefactorDataForIntegral[topo, targetInt];
   ratio = Cancel[
     materializeSectorPrefactor[sourceData]/materializeSectorPrefactor[targetData]
     ];
   Expand[rawCoefficient ratio targetInt]
   ];


sectorPrefactorRatio[source_Association, target_Association] := Cancel[
   materializeSectorPrefactor[Lookup[source, "sectorPrefactorData", <||>]]/
    materializeSectorPrefactor[Lookup[target, "sectorPrefactorData", <||>]]
   ];


sectorPrefactorLogDerivative[metadata_Association, variable_] := Cancel[
   D[materializeSectorPrefactor[Lookup[metadata, "sectorPrefactorData", <||>]], variable]/
    materializeSectorPrefactor[Lookup[metadata, "sectorPrefactorData", <||>]]
   ];


(* 按 sectorKey 从 DSInit 冻结的 metadata 读取物理 prefactor。tree dlog、naive DE 与
   普通 ds 必须共用这一份数据，不能从 contact coefficient 或旧 pack 反推连续零点幂。 *)
sectorPrefactorRecordForKey[context_Association, sectorKey_] := Module[
   {metadata, data, prefactor},
   metadata = SelectFirst[
     Lookup[context, "sectors", {}],
     Lookup[#, "sectorKey", None] === sectorKey &,
     Missing["NoSectorMetadata", sectorKey]
     ];
   If[Head[metadata] === Missing, Return[metadata]];
   data = Lookup[metadata, "sectorPrefactorData", Missing["NoSectorPrefactorData", sectorKey]];
   If[Head[data] === Missing, Return[data]];
   prefactor = materializeSectorPrefactor[data];
   If[prefactor === $Failed, Return[Missing["InvalidSectorPrefactorData", sectorKey]]];
   <|
    "sectorKey" -> sectorKey,
    "physicalSectorPrefactorData" -> data,
    "physicalSectorPrefactor" -> prefactor
    |>
   ];


(* makeSectorMetadataBase 保存 vertex/line slot 的基础 metadata；本包装层只在其上追加
   sector 不变量，不复制既有构造。 *)
makeSectorMetadata[topo_Association] := Module[
   {base, pattern, shrunkLines, lineSlots, ispSlots, parityData, timeOnlyStateSlots},
   base = makeSectorMetadataBase[topo];
   pattern = sectorPattern[topo, Lookup[Lookup[base, "lineSlots", {}], "packTemplate", {}]];
   shrunkLines = sectorShrunkLinesFromPattern[pattern];
   lineSlots = MapThread[
     Join[#1, <|
        "rootLinePosition" -> #2,
        "powerSlotKind" -> Lookup[pattern[[#2]], "powerKind"],
        "shrunkQ" -> (Lookup[pattern[[#2]], "state"] === "shrunk")
        |>] &,
     {Lookup[base, "lineSlots", {}], Range[Length[pattern]]}
     ];
   ispSlots = If[Lookup[topo, "ibpMode", "full"] === "timeOnly", {}, Lookup[base, "ispSlots", {}]];
   parityData = If[Length[DownValues[parityMetadataForSector]] > 0,
     parityMetadataForSector[topo],
     <|"status" -> "notLoaded"|>
     ];
    timeOnlyStateSlots = If[
      Lookup[topo, "ibpMode", "full"] === "timeOnly" &&
       Length[DownValues[dsTimeOnlyStateSlots]] > 0,
      dsTimeOnlyStateSlots[topo, Join[base, <|"lineSlots" -> lineSlots|>]],
      {}
      ];
    Join[base, <|
      "rootLineCount" -> Length[pattern],
      "rootLineOrder" -> Lookup[topo["lines"], "id", Range[Length[pattern]]],
      "ibpMode" -> Lookup[topo, "ibpMode", "full"],
      "sectorShrunkLines" -> shrunkLines,
      "sectorPattern" -> pattern,
      "sectorKey" -> sectorKeyFromPattern[topo, pattern],
      "sectorBits" -> If[
        Lookup[topo, "ibpMode", "full"] === "timeOnly",
        Characters[sectorKeyFromPattern[topo, pattern]],
        Missing["NotTimeOnlyBitString"]
        ],
      "sectorKeySchema" -> sectorKeySchemaFromTopology[topo],
      "lineSlots" -> lineSlots,
      "ispSlots" -> ispSlots,
      "timeOnlyStateSlots" -> timeOnlyStateSlots,
      "timeOnlyStateCount" -> Length[timeOnlyStateSlots],
      "publicIntegralRepresentation" -> If[
        Lookup[topo, "ibpMode", "full"] === "timeOnly",
        "J[sectorKey,timeShifts,stateBits]",
        "J[aList,linePacks,ispList]"
        ],
     "sectorPrefactorData" -> sectorPrefactorData[topo],
     "builtInRelationData" -> If[Length[DownValues[masslessBuiltInRelationData]] > 0,
       Join[masslessBuiltInRelationData[topo], masslessParallelFoldRelationData[topo]],
       {}
       ],
     "parityData" -> parityData,
      "representation" -> If[
        Lookup[topo, "ibpMode", "full"] === "timeOnly",
        "J[sectorKey,timeShifts,stateBits]",
        "J[aList,linePacks,ispList]"
        ]
      |>]
   ];


integralSectorMetadata[topo_Association, int_J] := Module[{metadataList, matches},
   metadataList = Lookup[topo, "sectorMetadataList", {makeSectorMetadata[topo]}];
   matches = Select[metadataList, integralMatchesSectorMetadataQ[int, #] &];
   Which[
    Length[matches] === 1, First[matches],
    Length[matches] === 0, Missing["NoMatchingSector"],
    True, Missing["AmbiguousSector", Lookup[matches, "sectorKey", {}]]
    ]
   ];


sectorTopologyForMetadata[rootTopo_Association, metadata_Association] := Module[{shrunk},
   shrunk = Lookup[metadata, "sectorShrunkLines", {}];
   If[shrunk === {}, rootTopo, shrinkSectorTopology[rootTopo, shrunk]]
   ];


sectorTopologyForIntegral[rootTopo_Association, int_J] := Module[{metadata},
   metadata = integralSectorMetadata[rootTopo, int];
   If[Head[metadata] === Missing, metadata, sectorTopologyForMetadata[rootTopo, metadata]]
   ];


(* ::Chapter:: *)
(*Sector-aware 公开 shape gate*)

publicIntegralShapeIssues[topo_Association, int : J[aList_, linePacks_, ispList_]] := Module[
   {metadata, issues = {}, lineSlots, expectedISP},
   metadata = integralSectorMetadata[topo, int];
   If[Head[metadata] === Missing,
    Return[{<|"slot" -> "sector", "reason" -> metadata|>}]
    ];
   lineSlots = Lookup[metadata, "lineSlots", {}];
   expectedISP = Length[Lookup[metadata, "ispSlots", {}]];
   If[Length[aList] =!= Length[Lookup[metadata, "compactASlots", {}]],
    AppendTo[issues, <|"slot" -> "aList", "expected" -> Length[Lookup[metadata, "compactASlots", {}]], "actual" -> Length[aList]|>]
    ];
   If[Length[linePacks] =!= Length[lineSlots],
    AppendTo[issues, <|"slot" -> "linePacks", "expected" -> Length[lineSlots], "actual" -> Length[linePacks]|>],
    Do[
     If[! ListQ[linePacks[[e]]] || ! linePackMatchesSlotQ[linePacks[[e]], lineSlots[[e]]],
      AppendTo[issues, <|"slot" -> "linePack", "lineIndex" -> e,
        "expected" -> Lookup[lineSlots[[e]], "packTemplate", {}], "actual" -> linePacks[[e]]|>]
      ],
     {e, Length[linePacks]}
     ]
    ];
   If[Length[ispList] =!= expectedISP,
    AppendTo[issues, <|"slot" -> "ispList", "expected" -> expectedISP, "actual" -> Length[ispList]|>]
    ];
   issues
   ];


publicResolvedDiscreteStateIssues[topo_Association, int : J[_, linePacks_, _]] := Module[
   {metadata, issues = {}, lineSlots, positions},
   metadata = integralSectorMetadata[topo, int];
   If[Head[metadata] === Missing, Return[{<|"slot" -> "sector", "reason" -> metadata|>}]];
   lineSlots = Lookup[metadata, "lineSlots", {}];
   Do[
    positions = Lookup[lineSlots[[e]], "nPositions", {}];
    Do[
     If[! MemberQ[{0, 1}, linePacks[[e, p]]],
      AppendTo[issues, <|"lineIndex" -> e, "packPosition" -> p, "value" -> linePacks[[e, p]]|>]
      ],
     {p, positions}
     ],
    {e, Length[lineSlots]}
    ];
   issues
   ];


(* time-only 公开对象先按冻结的 sector/state registry 还原为 producer 表示；
   shape 与离散态门禁仍由同一套既有检查负责，避免维护第二套物理规则。 *)
publicIntegralShapeIssues[topo_Association, int : J[_String, _List, _List]] /;
   Lookup[topo, "ibpMode", "full"] === "timeOnly" := Module[{internal},
   internal = dsTimeOnlyPublicIntegralToInternal[int, topo];
   If[
    internal === $Failed,
    {<|"slot" -> "timeOnlyIntegral", "reason" -> "invalidPublicTimeOnlyIntegral", "integral" -> int|>},
    publicIntegralShapeIssues[topo, internal]
    ]
   ];


publicResolvedDiscreteStateIssues[topo_Association, int : J[_String, _List, _List]] /;
   Lookup[topo, "ibpMode", "full"] === "timeOnly" := Module[{internal},
   internal = dsTimeOnlyPublicIntegralToInternal[int, topo];
   If[
    internal === $Failed,
    {<|"slot" -> "timeOnlyIntegral", "reason" -> "invalidPublicTimeOnlyIntegral", "integral" -> int|>},
    publicResolvedDiscreteStateIssues[topo, internal]
    ]
   ];


validatePublicExpression[expr_, topo_Association, requireDiscreteQ_: False] := Module[
   {integrals = publicExpressionIntegrals[expr], shapeIssues, stateIssues},
   shapeIssues = Flatten[publicIntegralShapeIssues[topo, #] & /@ integrals];
   If[shapeIssues =!= {}, Message[dSIBPPublicAPI::badshape, dsPublicShapeIssueSentence[shapeIssues]]; Return[False]];
   If[TrueQ[requireDiscreteQ],
    stateIssues = Flatten[publicResolvedDiscreteStateIssues[topo, #] & /@ integrals];
    If[stateIssues =!= {}, Message[dSIBPPublicAPI::badstate, dsPublicStateIssueSentence[stateIssues]]; Return[False]]
    ];
   True
   ];
