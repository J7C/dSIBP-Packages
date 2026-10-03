(* ::Package:: *)
(* 本模块保留 direct pure-time 公式作为 massive-only 交叉检查；生产路径统一使用
   J[aList,linePacks,{}]。regular 项直接由
   vertex-family 的 M1/M0 构造，contact 项直接消费 compiled WT/shrinkTerms 与共同 theta。
   loop 三槽表示仅由独立交叉验证调用，不是本模块的生产路径。 *)

(* ::Chapter:: *)
(*Pure-time 表示能力与状态访问*)

(* 同分支 massless 内线的 regular 导数在 n=0/1 间翻转；当前用户约定的 tree pack
   只保存 massive h 状态，因此这类线必须显式拒绝，不能当作无状态相位。 *)
dsPureTimeUnsupportedLines[family_Association] := Select[
   Range[Lookup[family["topology"], "nE", 0]],
   With[{line = family["topology", "lines"][[#]]},
     Lookup[line, "state", "full"] =!= "shrunk" && Lookup[line, "packType", ""] === "masslessFull"
     ] &
   ];


dsPureTimeFamilyUsableQ[family_Association] := dsPureTimeUnsupportedLines[family] === {};


dsTreeVertexIndex[family_Association, vertexId_] := Module[{position},
   position = FirstPosition[family["vertexOrder"], vertexId, Missing["NoVertex"]];
   If[Head[position] === Missing, position, First[position]]
   ];


dsTreeLegState[int_J, family_Association, lineId_, endpointSlot_Integer] := Module[
   {matches, vertexIndex, legIndex},
   matches = Reap[
       Do[
        legIndex = FirstPosition[
          Lookup[family["vertices"][[vertexIndex, "massiveLegs"]], "id", {}],
          {lineId, endpointSlot},
          Missing["NoLeg"]
          ];
        If[Head[legIndex] =!= Missing, Sow[{vertexIndex, First[legIndex]}]],
        {vertexIndex, Length[family["vertices"]]}
        ]
       ][[2]];
   matches = If[matches === {}, {}, First[matches]];
   If[Length[matches] =!= 1,
    Missing["TreeLegState", lineId, endpointSlot, Length[matches]],
    {vertexIndex, legIndex} = First[matches];
    int[[1, vertexIndex, 1 + legIndex]]
    ]
   ];


(* ::Chapter:: *)
(*Regular M1/M0 time seed*)

dsTreeStateVectorIntegrals[int_J, family_Association, vertexIndex_Integer, aValue_] := Module[
   {packs = First[int], states},
   states = treeBinaryStates[family["vertices"][[vertexIndex, "p"]]];
   J[ReplacePart[packs, vertexIndex -> Prepend[#, aValue]]] & /@ states
   ];


dsDirectTreeRegularSeed[vertexId_, int_J, family_Association] := Module[
   {vertexIndex, vertex, pack, currentA, stateRow, currentIntegrals, lowerIntegrals, equationVector},
   If[! treeIntegralQ[int, family], Return[$Failed]];
   vertexIndex = dsTreeVertexIndex[family, vertexId];
   If[Head[vertexIndex] === Missing, Return[$Failed]];
   vertex = family["vertices"][[vertexIndex]];
   pack = First[int][[vertexIndex]];
   currentA = First[pack];
   stateRow = treeStateIndex[Rest[pack]];
   currentIntegrals = dsTreeStateVectorIntegrals[int, family, vertexIndex, currentA];
   lowerIntegrals = dsTreeStateVectorIntegrals[int, family, vertexIndex, currentA - 1];
   equationVector = Expand[
     treeM1[vertex, vertex["nu0"] + currentA] . lowerIntegrals +
      treeM0[vertex] . currentIntegrals
     ];
   Expand[equationVector[[stateRow]]]
   ];


(* ::Chapter:: *)
(*Direct 共同-theta contact*)

dsDirectTreeAtomicContactChoices[
   vertexId_, int_J, family_Association, lineIndex_Integer
   ] := Module[
   {topo = family["topology"], line, endpointSlots, endpointSlot, n1, n2, coefficient, shrinkTerms},
   line = topo["lines"][[lineIndex]];
   endpointSlots = lineEndpointSlotsAtVertex[line, vertexId];
   If[Length[endpointSlots] =!= 1 || Lookup[line, "packType", ""] =!= "massiveFull", Return[{}]];
   endpointSlot = First[endpointSlots];
   n1 = dsTreeLegState[int, family, line["id"], 1];
   n2 = dsTreeLegState[int, family, line["id"], 2];
   If[MemberQ[{n1, n2}, _Missing], Return[{}]];
    (* shrinkTerms 保存 -W_T；离散 prefactor 的负号把乘积恢复为物理 theta contact。 *)
    coefficient = -sameBranchSKSign[line] KroneckerDelta[n1 + n2, 1]
      (-1)^If[endpointSlot === 1, n1, n2];
   shrinkTerms = lineCompiledShrinkTerms[line];
   Map[
    <|
      "lineIndex" -> lineIndex,
      "coefficient" -> coefficient Lookup[#, "coefficient", 0],
      "bShift" -> Lookup[#, "bShift", 1],
      "zeroPointShift" -> Lookup[#, "zeroPointShift", lineShrinkZeroPointShift[line]],
      "aShift" -> Lookup[#, "bShift", 1]
      |> &,
    shrinkTerms
    ]
   ];


dsTreeContactTargetIntegral[
   int_J,
   sourceFamily_Association,
   targetFamily_Association,
   choices_List
   ] := Module[
   {sourcePacks = First[int], sourceTopo = sourceFamily["topology"], targetTopo = targetFamily["topology"],
    sourceVertices, targetVertices, targetRepMap, selectedLines, targetPacks, sourceClass, aValue, legStates},
   sourceVertices = sourceFamily["vertexOrder"];
   targetVertices = targetFamily["vertexOrder"];
   targetRepMap = Lookup[
     targetTopo,
     "sectorVertexRepresentativeMap",
     AssociationThread[targetTopo["vertexIds"] -> targetTopo["vertexIds"]]
     ];
   selectedLines = Lookup[choices, "lineIndex"];
   targetPacks = Table[
     sourceClass = Select[sourceVertices, Lookup[targetRepMap, #, #] === targetVertices[[targetIndex]] &];
     aValue = Total[
        First[sourcePacks[[dsTreeVertexIndex[sourceFamily, #]]]] & /@ sourceClass
        ] - Total[
        Map[
         Function[choice,
          If[
           And @@ (Lookup[targetRepMap, #, #] === targetVertices[[targetIndex]] & /@
              sourceTopo["lines"][[choice["lineIndex"], "endpoints"]]),
           choice["aShift"],
           0
           ]
          ],
         choices
         ]
        ];
     legStates = Map[
       Function[leg,
        dsTreeLegState[int, sourceFamily, leg["id"][[1]], leg["id"][[2]]]
        ],
       targetFamily["vertices"][[targetIndex, "massiveLegs"]]
       ];
     If[MemberQ[legStates, _Missing], Return[$Failed]];
     Prepend[legStates, aValue],
     {targetIndex, Length[targetVertices]}
     ];
   J[targetPacks]
   ];


dsTreeLineZeroPoint[topo_Association, e_Integer] := If[
   MemberQ[Lookup[topo, "sectorShrunkLines", {}], e] || Lookup[topo["lines"][[e]], "state", "full"] === "shrunk",
   lineBSZeroPoint[topo, e],
   lineBZeroPoint[topo, e]
   ];


dsDirectTreeContactCoefficient[
   int_J,
   targetInt_J,
   sourceFamily_Association,
   targetFamily_Association,
   choices_List,
   thetaBundleCoefficient_
   ] := Module[
    {sourceTopo = sourceFamily["topology"], targetTopo = targetFamily["topology"],
     shiftByLine, deltaLinePowers, lineMomentumMagnitudes, atomicCoefficient},
   (* source 与 target 都按 (-tau)^A 定义，合并时间幂不产生额外连续幂相位。
      vertex basis 已除去逐 sector prefactor；统一 J 中的 residual power 再乘 N_t/N_s 后，
      每条 contact line 只剩 compiled s+z，因而不随用户 bS0 重基改变。 *)
   shiftByLine = Association@Map[
      Lookup[#, "lineIndex"] -> (Lookup[#, "bShift", 0] + Lookup[#, "zeroPointShift", 0]) &,
      choices
      ];
   deltaLinePowers = Table[
     Expand[Lookup[shiftByLine, e, 0]],
     {e, sourceTopo["nE"]}
     ];
   lineMomentumMagnitudes = Table[lineMomentumMagnitude[targetTopo, e], {e, targetTopo["nE"]}];
   atomicCoefficient = Times @@ Lookup[choices, "coefficient"];
   Expand[
    thetaBundleCoefficient atomicCoefficient
      Times @@ MapThread[Power, {lineMomentumMagnitudes, -deltaLinePowers}]
      sectorPrefactorRatio018[
       makeSectorMetadata[sourceTopo],
       makeSectorMetadata[targetTopo]
       ]
    ]
   ];


dsDirectTreeContactTerms[
   vertexId_, int_J, sourceFamily_Association, familyContext_Association
   ] := Module[
   {topo = sourceFamily["topology"], connectedLines, eligibleLines, bundles, sourceShrunk,
    oddSubsets, atomicChoices, totalShrunk, targetSector, targetFamily, targetInt, coefficient, terms = {}},
   connectedLines = DeleteDuplicates@Flatten[Cases[
       Lookup[topo, "vertexLines", {}][[dsTreeVertexIndex[sourceFamily, vertexId]]],
       {e_Integer, _} :> e
       ]];
   eligibleLines = Select[
     connectedLines,
     Lookup[topo["lines"][[#]], "packType", ""] === "massiveFull" &&
       Lookup[topo["lines"][[#]], "state", "full"] =!= "shrunk" &&
       Length[lineEndpointSlotsAtVertex[topo["lines"][[#]], vertexId]] === 1 &
     ];
   bundles = GatherBy[eligibleLines, thetaBundleKey[topo, #] &];
   sourceShrunk = Lookup[topo, "sectorShrunkLines", {}];
   Do[
    oddSubsets = Select[Rest[Subsets[bundle]], OddQ[Length[#]] &];
    Do[
     atomicChoices = dsDirectTreeAtomicContactChoices[vertexId, int, sourceFamily, #] & /@ selected;
     If[AnyTrue[atomicChoices, # === {} &], Continue[]];
     Do[
      totalShrunk = Sort@Union[sourceShrunk, Lookup[choice, "lineIndex"]];
       targetSector = sectorKeyFromShrunkLines[topo, totalShrunk];
      targetFamily = dsTreeFamilyBySector[targetSector, familyContext];
      If[Head[targetFamily] === Missing, Return[$Failed]];
      targetInt = dsTreeContactTargetIntegral[int, sourceFamily, targetFamily, choice];
      If[targetInt === $Failed, Return[$Failed]];
      coefficient = dsDirectTreeContactCoefficient[
        int, targetInt, sourceFamily, targetFamily, choice, 2^(1 - Length[selected])
        ];
      If[! TrueQ[coefficient === 0],
       AppendTo[terms, <|
         "sectorKey" -> targetSector,
         "integral" -> targetInt,
         "coefficient" -> coefficient,
         "selectedLines" -> Lookup[choice, "lineIndex"],
         "route" -> "directCompiledTheta"
         |>]
       ],
      {choice, Tuples[atomicChoices]}
      ],
     {selected, oddSubsets}
     ],
    {bundle, bundles}
    ];
   terms
   ];


(* ::Chapter:: *)
(*公开 direct seed record*)

dsTreeExpressionTerms[expr_, sectorKey_String] := Module[{terms, records},
   terms = If[Head[Expand[expr]] === Plus, List @@ Expand[expr], {Expand[expr]}];
   records = Map[
     Function[term,
      With[{integrals = DeleteDuplicates[Cases[term, int_J :> int, {0, Infinity}]]},
       If[Length[integrals] =!= 1,
        $Failed,
        <|"sectorKey" -> sectorKey, "integral" -> First[integrals],
          "coefficient" -> Expand[term /. First[integrals] -> 1], "route" -> "directM1M0"|>
        ]
       ]
      ],
     terms
     ];
   If[MemberQ[records, $Failed], $Failed, records]
   ];


dsDirectTreeSeedRecord[
   vertexId_, int_J, sourceFamily_Association, familyContext_Association
   ] := Module[{regular, regularTerms, contactTerms, combined},
   If[! dsPureTimeFamilyUsableQ[sourceFamily],
    Return[<|"status" -> "failed", "reason" -> "masslessFullNeedsTreeState",
      "lineIndices" -> dsPureTimeUnsupportedLines[sourceFamily]|>]
    ];
   regular = dsDirectTreeRegularSeed[vertexId, int, sourceFamily];
   If[regular === $Failed, Return[<|"status" -> "failed", "reason" -> "regularSeedFailed"|>]];
   regularTerms = dsTreeExpressionTerms[regular, sourceFamily["sector"]];
   contactTerms = dsDirectTreeContactTerms[vertexId, int, sourceFamily, familyContext];
   If[MemberQ[{regularTerms, contactTerms}, $Failed],
    Return[<|"status" -> "failed", "reason" -> "contactSeedFailed"|>]
    ];
   combined = dsCombineTreeTaggedTerms[Join[regularTerms, contactTerms]];
   <|
    "status" -> "generated",
    "generator" -> dtau[vertexId],
    "treeSeed" -> Total[Lookup[#, "coefficient"] Lookup[#, "integral"] & /@ combined],
    "treeIntegral" -> int,
    "sectorKey" -> sourceFamily["sector"],
    "generationRoute" -> "directPureTime",
    "loopSeed" -> Missing["NotUsed"],
    "treeLinearData" -> <|
      "status" -> "generated",
      "terms" -> combined,
      "termCount" -> Length[combined],
      "sectorKeys" -> DeleteDuplicates[Lookup[combined, "sectorKey"]],
      "expression" -> Total[Lookup[#, "coefficient"] Lookup[#, "integral"] & /@ combined],
      "referenceTreeIntegral" -> int,
      "coefficientConvention" -> "direct tree physical powers from M1/M0 and compiled WT"
      |>
   |>
   ];


(* ::Chapter:: *)
(*Pure-time backend-neutral linearData*)

dsPureTimeTaggedIntegral[term_Association] := dsTreeToken[
   Lookup[term, "sectorKey", Missing["SectorKeyRequired"]],
   Lookup[term, "integral", Missing["TreeIntegral"]]
   ];


dsPureTimeLinearEquation[record_Association, integralIndex_Association] := Module[
   {terms, rules},
   terms = Lookup[Lookup[record, "treeLinearData", <||>], "terms", {}];
   rules = Merge[
     Rule[
        integralIndex[dsPureTimeTaggedIntegral[#]],
        Lookup[#, "coefficient", 0]
        ] & /@ terms,
     Total
     ];
   <|
    "source" -> "directPureTime",
    "generator" -> Lookup[record, "generator", Missing["generator"]],
     "sectorKey" -> Lookup[record, "sectorKey", Missing["SectorKeyRequired"]],
    "referenceTreeIntegral" -> Lookup[record, "treeIntegral", Missing["TreeIntegral"]],
    "coefficientRules" -> Normal[rules],
    "constantTerm" -> 0,
    "nonlinearTerms" -> {},
    "linearQ" -> True
    |>
   ];


makePureTimeLinearSystemData[batch_Association, context_Association] /; dsContextQ[context] := Module[
   {records, terms, tokens, integralIndex, linearEquations, publicIntegrals, coefficientDiagnostics},
   If[Lookup[batch, "status", "missing"] =!= "generated" ||
     Lookup[batch, "representation", None] =!= "J[vertexPacks]",
    Return[<|"status" -> "notGenerated", "reason" -> "notPureTimeSeedBatch"|>]
    ];
   records = Lookup[batch, "seedRecords", {}];
   terms = Flatten[Lookup[Lookup[records, "treeLinearData", <||>], "terms", {}], 1];
   If[AnyTrue[terms, ! StringQ[Lookup[#, "sectorKey", None]] || ! MatchQ[Lookup[#, "integral", None], J[_List]] &],
    Return[<|"status" -> "notReady", "reason" -> "invalidSectorTaggedTreeTerms"|>]
    ];
   tokens = DeleteDuplicates[dsPureTimeTaggedIntegral /@ terms];
   integralIndex = makeIntegralIndex[tokens];
   linearEquations = dsPureTimeLinearEquation[#, integralIndex] & /@ records;
   publicIntegrals = Map[
     <|"id" -> integralIndex[#], "sectorKey" -> #[[1]], "integral" -> #[[2]]|> &,
     tokens
     ];
   coefficientDiagnostics = linearCoefficientDiagnostics[linearEquations];
   <|
    "status" -> "generated",
    "caseName" -> Lookup[batch, "caseName", context["topology", "name"]],
    "ibpMode" -> "timeOnly",
    "representation" -> "sectorTaggedJ[vertexPacks]",
    "topology" -> context["topology"],
    "integralCount" -> Length[tokens],
    "equationCount" -> Length[linearEquations],
    "integralList" -> tokens,
    "publicIntegralList" -> publicIntegrals,
    "integralRules" -> Normal[integralIndex],
    "sectorMetadataList" -> Lookup[batch, "sectorMetadataList", {}],
    "topologyValidationReport" -> context["validationReport"],
    "seedCoverageReport" -> <|
      "passQ" -> True,
      "mode" -> "timeOnly",
      "completeTimeGenerationQ" -> True,
      "completeMomentumGenerationQ" -> False
      |>,
    "linearEquations" -> linearEquations,
    "linearQ" -> And @@ Lookup[linearEquations, "linearQ", False],
    "nonlinearEquationCount" -> Count[Lookup[linearEquations, "linearQ", False], False],
    "numericCoefficientSystemQ" -> coefficientDiagnostics["numericCoefficientSystemQ"],
    "coefficientVariables" -> coefficientDiagnostics["coefficientVariables"],
    "sourceSeedBatch" -> KeyDrop[batch, {"equations", "seedRecords"}]
    |>
   ];
