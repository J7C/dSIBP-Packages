(* ::Package:: *)
(* 本模块统一 massless 双端点 relation、sector-aware canonical/求导和 parity seed 域。
   parity 只筛选待作用生成元的 seed 点；生成后证书只报告错误，绝不把积分替换为零。 *)

(* ::Chapter:: *)
(*Massless 双端点 quotient*)

masslessFullLineQ018[line_Association] := Lookup[line, "packType", ""] === "masslessFull" &&
   Lookup[line, "state", "full"] =!= "shrunk";


canonicalizeMasslessIntegral018[
   topo_Association,
   int : J[aList_, linePacks_, ispList_]
   ] := Module[{newPacks = linePacks, coefficient = 1, positions, n1, n2},
   Do[
    If[actualLinePackType[topo, e, newPacks[[e]]] =!= "masslessFull", Continue[]];
    positions = linePackNPositions[topo["lines"][[e]], "masslessFull"];
    {n1, n2} = newPacks[[e, positions]];
    If[MemberQ[{0, 1}, n1] && MemberQ[{0, 1}, n2] && n2 === 1,
     coefficient = -coefficient;
     newPacks[[e, positions]] = {1 - n1, 0}
     ],
    {e, Length[topo["lines"]]}
    ];
   coefficient J[aList, newPacks, ispList]
   ];


masslessCoincidentAntisymmetricIntegralQ[
   topo_Association,
   int : J[_, linePacks_, _]
   ] := Module[{repMap, originalEndpoints, targetEndpoints, positions, states},
   repMap = integralTargetVertexRepresentativeMap[topo, int];
   AnyTrue[
    Range[Length[topo["lines"]]],
    Function[e,
     If[actualLinePackType[topo, e, linePacks[[e]]] =!= "masslessFull",
      False,
      originalEndpoints = Lookup[topo["lines"][[e]], "originalEndpoints", topo["lines"][[e, "endpoints"]]];
      targetEndpoints = Lookup[repMap, originalEndpoints];
      positions = linePackNPositions[topo["lines"][[e]], "masslessFull"];
      states = linePacks[[e, positions]];
      TrueQ[SameQ @@ targetEndpoints && And @@ (IntegerQ /@ states) && OddQ[Total[states]]]
      ]
     ]
    ]
   ];


applyMasslessEndpointCanonical[expr_, topo_Association] := Module[{quotient},
   quotient = Expand[expr /. int_J :> canonicalizeMasslessIntegral018[topo, int]];
   Expand[quotient /. (int_J /; masslessCoincidentAntisymmetricIntegralQ[topo, int]) :> 0]
   ];


masslessBuiltInRelationData018[topo_Association] := Map[
   Function[e,
    With[{id = topo["lines"][[e, "id"]]},
     <|
      "lineIndex" -> e,
      "lineId" -> id,
      "relations" -> {
        HoldForm[F[id, 0, 1] + F[id, 1, 0] == 0],
        HoldForm[F[id, 1, 1] + F[id, 0, 0] == 0]
        },
      "canonicalDirection" -> "n2ToZero"
      |>
     ]
    ],
   Select[Range[topo["nE"]], masslessFullLineQ018[topo["lines"][[#]]] &]
   ];


forbiddenNDataForIntegral[topo_Association, J[_, linePacks_, _]] := Module[
   {issues = {}, packType, positions, values},
   Do[
    packType = actualLinePackType[topo, e, linePacks[[e]]];
    positions = linePackNPositions[topo["lines"][[e]], packType];
    values = If[positions === {}, {}, linePacks[[e, positions]]];
    Switch[packType,
     "massiveFull" | "massiveCross",
     Do[
      If[IntegerQ[values[[endpointSlot]]] && values[[endpointSlot]] >= 2,
       AppendTo[issues, <|"lineIndex" -> e, "packType" -> packType,
         "endpointSlot" -> endpointSlot, "nValue" -> values[[endpointSlot]]|>]
       ],
      {endpointSlot, Length[values]}
      ],
     "masslessFull",
     Do[
      If[IntegerQ[values[[endpointSlot]]] && ! MemberQ[{0, 1}, values[[endpointSlot]]],
       AppendTo[issues, <|"lineIndex" -> e, "packType" -> packType,
         "endpointSlot" -> endpointSlot, "nValue" -> values[[endpointSlot]]|>]
       ],
      {endpointSlot, Length[values]}
      ],
     _, Null
     ],
    {e, Length[topo["lines"]]}
    ];
   issues
   ];


(* ::Chapter:: *)
(*Massless regular 与 contact 原子*)

toggleMasslessEndpointState018[
   J[aList_, linePacks_, ispList_], e_Integer, endpointSlot_Integer
   ] := Module[{newPacks = linePacks, position},
   position = 1 + endpointSlot;
   newPacks[[e, position]] = 1 - newPacks[[e, position]];
   J[aList, newPacks, ispList]
   ];


toggleMasslessLineState[topo_Association, int_J, e_Integer] :=
  toggleMasslessEndpointState018[int, e, 1];


timeMasslessEndpointDerivativeTerms[
   topo_Association,
   int : J[_, linePacks_, _],
   vertexId_
   ] := Module[{pos, connectedLines, line, endpointSlots, sigma, shiftedIntegral},
   pos = vertexPosition[topo, vertexId];
   If[Head[pos] === Missing, Return[0]];
   connectedLines = topo["vertexLines"][[pos]][[All, 1]];
   Total@Table[
     line = topo["lines"][[e]];
     endpointSlots = lineEndpointSlotsAtVertex[line, vertexId];
     Switch[actualLinePackType[topo, e, linePacks[[e]]],
      "masslessFull",
      sigma = sameBranchSKSign[line];
      Total@Table[
        shiftedIntegral = shiftLinePower[
          topo, toggleMasslessEndpointState018[int, e, endpointSlot], e, -1
          ];
        I sigma shiftedIntegral,
        {endpointSlot, endpointSlots}
        ],
      "masslessCross",
      Total@Table[
        I skEndpointPhaseSign[line, endpointSlot] shiftLinePower[topo, int, e, -1],
        {endpointSlot, endpointSlots}
        ],
      _, 0
      ],
     {e, connectedLines}
     ]
   ];


momentumBuildingBlockDerivativeTerms[
   topo_Association, int_J, gen_Association, repSP2ZRules_List
   ] := Module[
   {dLoop = gen["dLoop"], vector = gen["vector"], lineMomenta, line, loopCoeff,
    vDotQ, packType, sigma, shiftedInt},
   lineMomenta = Lookup[topo["lines"], "momentum"];
   Total@Table[
     line = topo["lines"][[e]];
     loopCoeff = Coefficient[lineMomenta[[e]], topo["loopMomenta"][[dLoop]]];
     If[zeroQ[loopCoeff],
      0,
      vDotQ = Expand[expandDotExpr[vector, lineMomenta[[e]], topo] /. repSP2ZRules];
      packType = actualLinePackType[topo, e, int[[2, e]]];
      Switch[packType,
       "massiveFull" | "massiveCross",
       Total@Table[
         loopCoeff compiledMomentumEndpointDerivativeTerms[topo, int, e, endpointSlot, vDotQ],
         {endpointSlot, 2}
         ],
       "masslessFull",
       sigma = sameBranchSKSign[line];
       loopCoeff Total@Table[
        (* bridge 的 fixed-line 幂移会产生显式模长系数；先对裸 J 做顶点移位，
           再移线幂，避免把 coefficient J 整体送入只接受裸积分的 helper。 *)
        shiftedInt = shiftVertexA[
          toggleMasslessEndpointState018[int, e, endpointSlot],
          topo, line["endpoints"][[endpointSlot]], 1
          ];
        -I sigma absorbLinearFactor[
          vDotQ,
          shiftLinePower[topo, shiftedInt, e, 1],
          topo
          ],
         {endpointSlot, 2}
         ],
       "masslessCross",
       loopCoeff Total@Table[
         shiftedInt = shiftVertexA[int, topo, line["endpoints"][[endpointSlot]], 1];
         -I skEndpointPhaseSign[line, endpointSlot]
          absorbLinearFactor[vDotQ, shiftLinePower[topo, shiftedInt, e, 1], topo],
         {endpointSlot, 2}
         ],
       _, 0
       ]
      ],
     {e, topo["nE"]}
     ]
   ];


externalVectorBuildingBlockDerivativeTerms[
   topo_Association, int_J, gen_Association, repSP2ZRules_List
   ] := Module[
   {dExternal = gen["dExternal"], vector = gen["vector"], lineMomenta, line, extCoeff,
    vDotQ, packType, sigma, shiftedInt},
   lineMomenta = Lookup[topo["lines"], "momentum"];
   Total@Table[
     line = topo["lines"][[e]];
     extCoeff = Coefficient[lineMomenta[[e]], topo["effectiveLoopExternalMomenta"][[dExternal]]];
     If[zeroQ[extCoeff] || externalLegCoordinateLineQ[lineMomenta[[e]], topo],
      0,
      vDotQ = Expand[expandDotExpr[vector, lineMomenta[[e]], topo] /. repSP2ZRules];
      packType = actualLinePackType[topo, e, int[[2, e]]];
      Switch[packType,
       "massiveFull" | "massiveCross",
       Total@Table[
         extCoeff compiledMomentumEndpointDerivativeTerms[topo, int, e, endpointSlot, vDotQ],
         {endpointSlot, 2}
         ],
       "masslessFull",
       sigma = sameBranchSKSign[line];
       extCoeff Total@Table[
        (* 与 loop-momentum 路线使用同一顺序：裸 J 顶点移位先于 fixed-line 系数生成。 *)
        shiftedInt = shiftVertexA[
          toggleMasslessEndpointState018[int, e, endpointSlot],
          topo, line["endpoints"][[endpointSlot]], 1
          ];
        -I sigma absorbLinearFactor[
          vDotQ,
          shiftLinePower[topo, shiftedInt, e, 1],
          topo
          ],
         {endpointSlot, 2}
         ],
       "masslessCross",
       extCoeff Total@Table[
         shiftedInt = shiftVertexA[int, topo, line["endpoints"][[endpointSlot]], 1];
         -I skEndpointPhaseSign[line, endpointSlot]
          absorbLinearFactor[vDotQ, shiftLinePower[topo, shiftedInt, e, 1], topo],
         {endpointSlot, 2}
         ],
       _, 0
       ]
      ],
     {e, topo["nE"]}
     ]
   ];


scalarMomentumBuildingBlockDerivativeTerms[topo_Association, int_J, e_Integer] := Module[
   {line = topo["lines"][[e]], packType, sigma},
   packType = actualLinePackType[topo, e, int[[2, e]]];
   Switch[packType,
    "massiveFull" | "massiveCross",
    Total[compiledScalarMomentumEndpointDerivativeTerms[topo, int, e, #] & /@ {1, 2}],
    "masslessFull",
     sigma = sameBranchSKSign[line];
    Total@Table[
      -I sigma shiftVertexA[
        toggleMasslessEndpointState018[int, e, endpointSlot],
        topo, line["endpoints"][[endpointSlot]], 1
        ],
      {endpointSlot, 2}
      ],
    "masslessCross",
    Total@Table[
      -I skEndpointPhaseSign[line, endpointSlot] shiftVertexA[
        int, topo, line["endpoints"][[endpointSlot]], 1
        ],
      {endpointSlot, 2}
      ],
    _, 0
    ]
   ];


thetaBoundaryAtomicTerms[
   topo_Association,
   J[aList_, linePacks_, ispList_],
   e_Integer,
   vertexId_
   ] := Module[
   {line = topo["lines"][[e]], endpointSlots, endpointSlot, endpointOrientation,
    pack = linePacks[[e]], packType, positions, coeff, shrinkTerms},
   endpointSlots = lineEndpointSlotsAtVertex[line, vertexId];
   If[Length[endpointSlots] =!= 1, Return[{}]];
   endpointSlot = First[endpointSlots];
   endpointOrientation = If[endpointSlot === 1, 1, -1];
   packType = actualLinePackType[topo, e, pack];
   Switch[packType,
    "massiveFull",
    positions = linePackNPositions[line, packType];
      coeff = -sameBranchSKSign[line] KroneckerDelta[Total[pack[[positions]]], 1]
        (-1)^pack[[positions[[endpointSlot]]]];
    shrinkTerms = lineCompiledShrinkTerms[line];
    (<|
        "lineIndex" -> e,
        "coefficient" -> coeff Lookup[#, "coefficient", 0],
        "bShift" -> Lookup[#, "bShift", 1],
        "zeroPointShift" -> Lookup[#, "zeroPointShift", lineShrinkZeroPointShift[line]],
        "aShift" -> Lookup[#, "bShift", 1]
        |> &) /@ shrinkTerms,
    "masslessFull",
    positions = linePackNPositions[line, packType];
    {<|
      "lineIndex" -> e,
      "coefficient" -> -2 endpointOrientation
        KroneckerDelta[Total[pack[[positions]]], 1] (-1)^pack[[positions[[2]]]],
      "bShift" -> 0,
      "zeroPointShift" -> 0,
      "aShift" -> 0
      |>},
    _, {}
    ]
   ];


(* ::Chapter:: *)
(*Fixed-line shrink 与 normalized coefficient*)

fixedLineShrinkResidualPower018[
   topo_Association,
   e_Integer,
   integerShift_,
   zeroPointShift_
   ] := Module[{rootRules, id, sourceZero, targetZero},
   rootRules = Lookup[topo, "rootZeroPointRules", topo["zeroPointRules"]];
   id = topo["lines"][[e, "id"]];
   sourceZero = zeroPointRuleValue018[rootRules, b0[id], 0];
   targetZero = lineTargetShrinkZeroPoint018[topo, e];
   Expand[integerShift + zeroPointShift - (targetZero - sourceZero)]
   ];


shrinkLineIntegral[
   topo_Association, int : J[aList_, linePacks_, ispList_], e_Integer,
   bShift_: Automatic, aShift_: Automatic
   ] := Module[
   {line = topo["lines"][[e]], uSlot, vSlot, oldActive, newRepMap, newActive,
    newAList, newLinePacks = linePacks, mergedRep, oldSlotsForNewRep, slotValues,
    effectiveBShift, effectiveAShift, effectiveZeroPointShift, powerCoefficient},
   effectiveBShift = If[bShift === Automatic, lineShrinkBShift[line], bShift];
   effectiveZeroPointShift = lineShrinkZeroPointShift[line];
   effectiveAShift = If[aShift === Automatic,
     If[Lookup[line, "massType", "massive"] === "massless", 0, effectiveBShift],
     aShift
     ];
   uSlot = vertexASlot[topo, line["endpoints"][[1]]];
   vSlot = vertexASlot[topo, line["endpoints"][[2]]];
   If[Head[uSlot] === Missing || Head[vSlot] === Missing, Return[$Failed]];
   oldActive = activeAVertexIds[topo];
   newRepMap = vertexRepresentativeMap[topo["vertexIds"], Join[
      ({#, vertexRepresentative[topo, #]} & /@ topo["vertexIds"]),
      {line["endpoints"]}
      ]];
   newActive = DeleteDuplicates[Lookup[newRepMap, topo["vertexIds"]]];
   mergedRep = Lookup[newRepMap, line["endpoints"][[1]]];
   newAList = Table[
     oldSlotsForNewRep = Flatten[Position[Lookup[newRepMap, oldActive], newActive[[i]]]];
     slotValues = aList[[oldSlotsForNewRep]];
     If[newActive[[i]] === mergedRep, Total[slotValues] - effectiveAShift, Total[slotValues]],
     {i, Length[newActive]}
     ];
   powerCoefficient = If[
     lineIndexedPowerQ[line],
     1,
     fixedLineMomentumMagnitude[topo, e]^(-effectiveBShift - effectiveZeroPointShift)
     ];
   newLinePacks[[e]] = If[
     lineIndexedPowerQ[line],
     {lineIntegerPowerIndex[topo, int, e] + effectiveBShift},
     {fixedLineSentinel018[]}
     ];
   powerCoefficient J[newAList, newLinePacks, ispList]
   ];


shrinkLinesIntegral[
   topo_Association,
   int : J[aList_, linePacks_, ispList_],
   specs_List
   ] := Module[
   {selectedLines, oldActive, pairs, newRepMap, newActive, newAList,
    newLinePacks = linePacks, oldSlotsForNewRep, selectedShiftForRep,
    powerCoefficient = 1},
   selectedLines = Lookup[specs, "lineIndex"];
   oldActive = activeAVertexIds[topo];
   pairs = topo["lines"][[#, "endpoints"]] & /@ selectedLines;
   newRepMap = vertexRepresentativeMap[topo["vertexIds"], Join[
      ({#, vertexRepresentative[topo, #]} & /@ topo["vertexIds"]), pairs
      ]];
   newActive = DeleteDuplicates[Lookup[newRepMap, topo["vertexIds"]]];
   newAList = Table[
     oldSlotsForNewRep = Flatten[Position[Lookup[newRepMap, oldActive], newActive[[i]]]];
     selectedShiftForRep = Total[MapThread[
        If[SameQ @@ Lookup[newRepMap, #1] && First[Lookup[newRepMap, #1]] === newActive[[i]], #2, 0] &,
        {pairs, Lookup[specs, "aShift"]}
        ]];
     Total[aList[[oldSlotsForNewRep]]] - selectedShiftForRep,
     {i, Length[newActive]}
     ];
   Scan[
    Function[spec,
     If[lineIndexedPowerQ[topo["lines"][[spec["lineIndex"]]]],
      newLinePacks[[spec["lineIndex"]]] = {
        lineIntegerPowerIndex[topo, int, spec["lineIndex"]] + spec["bShift"]
        },
      powerCoefficient *= fixedLineMomentumMagnitude[topo, spec["lineIndex"]]^(
        -spec["bShift"] -
         Lookup[spec, "zeroPointShift", lineShrinkZeroPointShift[topo["lines"][[spec["lineIndex"]]]]]
        );
      newLinePacks[[spec["lineIndex"]]] = {fixedLineSentinel018[]}
      ]
     ],
    specs
    ];
   powerCoefficient J[newAList, newLinePacks, ispList]
   ];


(* ::Chapter:: *)
(*Sector-aware canonical、ds 与 integrand*)

(* fixed/non-loop line 的零点幂属于 N_s；裸积分核导数只微分有 b/bS 槽的
   cycle denominator。模式函数本身的动量导数仍由原 building-block 路线处理。 *)
lineBarePowerIndex018[topo_Association, int_J, e_Integer] := If[
   lineIndexedPowerQ[topo["lines"][[e]]],
   linePowerIndex[topo, int, e],
   0
   ];


externalLegMagnitudeOccurrenceLineDerivativeSeed[topo_Association, int_J, coordinate_Association] := Module[
   {momentum = Lookup[coordinate, "momentum", Missing["NoMomentum"]], matchingLines},
   If[Head[momentum] === Missing, Return[0]];
   matchingLines = Flatten@Position[
      Lookup[topo["lines"], "momentum", {}],
      lineMomentum_ /; SameQ[canonicalExternalLegMomentum[lineMomentum], canonicalExternalLegMomentum[momentum]],
      {1},
      Heads -> False
      ];
   Total@Table[
     -lineBarePowerIndex018[topo, int, e] shiftLinePower[topo, int, e, 1] +
      scalarMomentumBuildingBlockDerivativeTerms[topo, int, e],
     {e, matchingLines}
     ]
   ];


externalVectorPropagatorDerivativeTerms[topo_Association, int_J, gen_Association, repSP2ZRules_List] := Module[
   {dExternal, vector, lineMomenta, extCoeff, vDotQ, shiftedInt},
   dExternal = gen["dExternal"];
   vector = gen["vector"];
   lineMomenta = Lookup[topo["lines"], "momentum"];
   Total@Table[
     extCoeff = Coefficient[lineMomenta[[e]], topo["effectiveLoopExternalMomenta"][[dExternal]]];
     If[zeroQ[extCoeff] || externalLegCoordinateLineQ[lineMomenta[[e]], topo],
      0,
      vDotQ = Expand[expandDotExpr[vector, lineMomenta[[e]], topo] /. repSP2ZRules];
      shiftedInt = shiftLinePower[topo, int, e, 2];
      -extCoeff lineBarePowerIndex018[topo, int, e] absorbLinearFactor[vDotQ, shiftedInt, topo]
      ],
     {e, topo["nE"]}
     ]
   ];

sectorAwareCanonicalTerm018[term_, rootTopo_Association] := Module[
   {integrals, int, coefficient, sectorTopo},
   integrals = DeleteDuplicates[Cases[term, _J, {0, Infinity}]];
   Which[
    integrals === {}, term,
    Length[integrals] =!= 1, $Failed,
    True,
    int = First[integrals];
    coefficient = Cancel[term/int];
    sectorTopo = sectorTopologyForIntegral018[rootTopo, int];
    If[Head[sectorTopo] === Missing, $Failed, applySeedCanonical[coefficient int, sectorTopo]]
    ]
   ];


sectorAwareCanonical018[expr_, rootTopo_Association] := Module[{terms, result},
   terms = linearTerms[Expand[expr]];
   result = sectorAwareCanonicalTerm018[#, rootTopo] & /@ terms;
   If[MemberQ[result, $Failed], $Failed, Expand[Total[result]]]
   ];


(* 连续撒点会反复命中同一批 sector。缓存只复用已经初始化的 sector topology，
   积分仍须与该 sector 的 shape metadata 匹配，不能用缓存绕过表示门禁。 *)
sectorTopologyCache018[rootTopo_Association, metadataList_List] := Module[
   {keys, entries, sectorTopo, symmetryRules, tadpoleData, postSamplingCanonicalRequiredQ},
   keys = Lookup[metadataList, "sectorKey", Missing["NoSectorKey"]];
   If[MemberQ[keys, _Missing] || DuplicateFreeQ[keys] =!= True,
    Return[Missing["InvalidSectorMetadataKeys", keys]]
    ];
   entries = Map[
     Function[metadata,
      sectorTopo = sectorTopologyForMetadata018[rootTopo, metadata];
      symmetryRules = effectiveSymmetryRules0[sectorTopo];
      If[
       ! ListQ[symmetryRules] || ! And @@ (validDiscreteReplacementRuleQ /@ symmetryRules),
       Return[Missing["InvalidSectorSymmetryRules", Lookup[metadata, "sectorKey"]], Module]
       ];
      tadpoleData = Select[
        tadpoleLoopReversalData[sectorTopo],
        MemberQ[{"massiveFull", "masslessFull"}, Lookup[#, "packType", None]] &&
          TrueQ[Lookup[#, "exclusiveLoopQ", False]] &
        ];
      (* EOM 与端点 canonical 已在 general template 层完成。只有用户规则可能依赖
         撒点后的具体指标，或 tadpole odd-ISP 判定需要整数 ISP 幂次时，才逐点重跑。 *)
      postSamplingCanonicalRequiredQ =
       Lookup[sectorTopo, "symmetryRules", {}] =!= {} ||
        (tadpoleData =!= {} && Lookup[sectorTopo, "ispData", {}] =!= {});
      Lookup[metadata, "sectorKey"] -> <|
        "metadata" -> metadata,
        "topology" -> sectorTopo,
        "symmetryRules" -> symmetryRules,
        "postSamplingCanonicalRequiredQ" -> postSamplingCanonicalRequiredQ
        |>
      ],
     metadataList
     ];
   Association[entries]
   ];


sectorCachePostSamplingCanonicalRequiredQ018[cache_Association] := AnyTrue[
   Values[cache],
   TrueQ[Lookup[#, "postSamplingCanonicalRequiredQ", True]] &
   ];


sectorKeyForIntegral018[
   rootTopo_Association,
   J[_, linePacks_List, _]
   ] := If[
   Length[linePacks] =!= Length[Lookup[rootTopo, "lines", {}]],
   Missing["LinePackCountMismatch"],
   sectorKeyFromPattern018[rootTopo, sectorPattern018[rootTopo, linePacks]]
   ];


sectorTopologyForIntegral018[
   rootTopo_Association,
   int_J,
   cache_Association
   ] := Module[{key, entry, metadata},
   key = sectorKeyForIntegral018[rootTopo, int];
   If[Head[key] === Missing, Return[key]];
   entry = Lookup[cache, key, Missing["NoMatchingSector", key]];
   If[Head[entry] === Missing, Return[entry]];
   metadata = Lookup[entry, "metadata", Missing["NoSectorMetadata", key]];
   If[Head[metadata] === Missing || ! TrueQ[integralMatchesSectorMetadataQ[int, metadata]],
    Return[Missing["SectorShapeMismatch", key]]
    ];
   Lookup[entry, "topology", Missing["NoSectorTopology", key]]
   ];


sectorEntryForIntegral018[
   rootTopo_Association,
   int_J,
   cache_Association
   ] := Module[{key, entry, metadata},
   key = sectorKeyForIntegral018[rootTopo, int];
   If[Head[key] === Missing, Return[key]];
   entry = Lookup[cache, key, Missing["NoMatchingSector", key]];
   If[Head[entry] === Missing, Return[entry]];
   metadata = Lookup[entry, "metadata", Missing["NoSectorMetadata", key]];
   If[Head[metadata] === Missing || ! TrueQ[integralMatchesSectorMetadataQ[int, metadata]],
    Return[Missing["SectorShapeMismatch", key]]
    ];
   entry
   ];


applySeedCanonicalWithRules018[
   expr_,
   topo_Association,
   symmetryRules_List
   ] := Expand[
   applyMasslessEndpointCanonical[
      applyMassiveCoincidentCanonical[applyEOM[expr, topo], topo],
      topo
      ] /. symmetryRules
   ];


sectorAwareSymmetryTerm018[
   term_,
   rootTopo_Association,
   cache_Association
    ] := Module[{integrals, int, coefficient, sectorEntry, symmetryRules, symmetricInt},
   integrals = DeleteDuplicates[Cases[term, _J, {0, Infinity}]];
   Which[
    integrals === {}, term,
    Length[integrals] =!= 1, $Failed,
    True,
    int = First[integrals];
    coefficient = Cancel[term/int];
    sectorEntry = sectorEntryForIntegral018[rootTopo, int, cache];
    If[Head[sectorEntry] === Missing, Return[$Failed, Module]];
    symmetryRules = Lookup[sectorEntry, "symmetryRules", Missing["NoSectorSymmetryRules"]];
    If[
     Head[symmetryRules] === Missing,
     $Failed,
     symmetricInt = Expand[int /. symmetryRules];
     Expand[coefficient symmetricInt]
     ]
    ]
   ];


sectorAwareSymmetry018[
   expr_,
   rootTopo_Association,
   cache_Association
   ] := Module[{terms, result},
   terms = linearTerms[Expand[expr]];
   result = sectorAwareSymmetryTerm018[#, rootTopo, cache] & /@ terms;
   If[MemberQ[result, $Failed], $Failed, Expand[Total[result]]]
   ];


sectorAwareCanonicalTerm018[
   term_,
   rootTopo_Association,
   cache_Association
    ] := Module[{integrals, int, coefficient, sectorEntry, sectorTopo, symmetryRules, canonicalInt},
   integrals = DeleteDuplicates[Cases[term, _J, {0, Infinity}]];
   Which[
    integrals === {}, term,
    Length[integrals] =!= 1, $Failed,
    True,
    int = First[integrals];
    coefficient = Cancel[term/int];
    sectorEntry = sectorEntryForIntegral018[rootTopo, int, cache];
    If[Head[sectorEntry] === Missing, Return[$Failed, Module]];
    sectorTopo = Lookup[sectorEntry, "topology", Missing["NoSectorTopology"]];
    symmetryRules = Lookup[sectorEntry, "symmetryRules", Missing["NoSectorSymmetryRules"]];
     If[Head[sectorTopo] === Missing || Head[symmetryRules] === Missing, Return[$Failed, Module]];
     (* EOM 与 canonical helper 的合同是裸 J。显式动力学系数只在外层保存，
        这样 fixed-line 幂移不能把系数误带进后续的顶点或指标移位。 *)
     canonicalInt = applySeedCanonicalWithRules018[int, sectorTopo, symmetryRules];
     If[canonicalInt === $Failed, $Failed, Expand[coefficient canonicalInt]]
     ]
   ];


sectorAwareCanonical018[
   expr_,
   rootTopo_Association,
   cache_Association
   ] := Module[{terms, result},
   terms = linearTerms[Expand[expr]];
   result = sectorAwareCanonicalTerm018[#, rootTopo, cache] & /@ terms;
   If[MemberQ[result, $Failed], $Failed, Expand[Total[result]]]
   ];


ds[expr_, userVariable_, topoSpec_Association] := Module[
   {rootTopo, internalExpr, variableData, userVariables, userExpr, linearData, internalVariable,
     coefficientDerivative, integralDerivativeTerms, sectorTopo, prefactorData,
     prefactorLogDerivative, result},
   rootTopo = resolvePublicTopologyContext[topoSpec];
   If[rootTopo === $Failed, Return[$Failed]];
   If[! dsTopologyCapabilityQ[rootTopo, "derivativeUsableQ"],
    dsErrorPrint["当前参数声明不支持唯一 ds 微分算符。 The current parameter declaration does not define a unique ds operator."];
    Return[$Failed]
    ];
   userVariables = Lookup[publicIndependentVariableDerivativeData[rootTopo], "userVariable", {}];
   variableData = resolvePublicIndependentVariableDerivativeData[rootTopo, userVariable];
    If[Head[variableData] === Missing,
     Message[dSIBPPublicAPI::badvar, userVariable, userVariables]; Return[$Failed]
     ];
    internalExpr = If[
      Lookup[rootTopo, "ibpMode", "full"] === "timeOnly",
      dsTimeOnlyExpressionToInternal020[expr, rootTopo],
      expr
      ];
    If[internalExpr === $Failed, Return[$Failed]];
    userExpr = rep2outform[internalExpr, rootTopo];
   If[userExpr === $Failed || ! validatePublicExpression[userExpr, rootTopo, True], Return[$Failed]];
   linearData = publicLinearIntegralDecomposition[userExpr];
   If[Lookup[linearData, "status", "failed"] =!= "linear",
    Message[dSIBPPublicAPI::nonlinear, Lookup[linearData, "heldExpression", userExpr]];
    Return[$Failed]
    ];
   internalVariable = variableData["variable"];
   coefficientDerivative = Expand[
     D[linearData["heldExpression"], userVariable] /. linearData["backwardRules"]
     ];
   integralDerivativeTerms = MapThread[
     Function[{coefficient, int},
      sectorTopo = sectorTopologyForIntegral018[rootTopo, int];
       If[Head[sectorTopo] === Missing,
        $Failed,
        prefactorData = sectorPrefactorDataForIntegral018[rootTopo, int];
        prefactorLogDerivative = sectorPrefactorLogDerivative018[
          <|"sectorPrefactorData" -> prefactorData|>,
          userVariable
          ];
        With[{term = applyIndependentVariableDerivativeSeed[sectorTopo, int, internalVariable]},
         If[term === $Failed || prefactorLogDerivative === $Failed,
          $Failed,
          coefficient (term + prefactorLogDerivative int)
          ]
         ]
       ]
      ],
     {linearData["coefficients"], linearData["integrals"]}
     ];
   If[MemberQ[integralDerivativeTerms, $Failed],
    Message[dSIBPPublicAPI::derivativefailed, userVariable]; Return[$Failed]
    ];
   result = sectorAwareCanonical018[
     Expand[coefficientDerivative + Total[integralDerivativeTerms]],
     rootTopo
     ];
    If[result === $Failed, Return[$Failed]];
    result = rep2outform[result, rootTopo];
    If[result === $Failed, Return[$Failed]];
    dsTimeOnlyExpressionToPublic020[result, rootTopo]
    ];


integrandBuildingBlock[line_Association, pack_List, momentumMagnitude_] := Module[
   {lineId = line["id"], endpoints = line["endpoints"], packType = line["packType"], positions},
   positions = linePackNPositions[line, packType];
   Switch[packType,
    "massiveFull" | "massiveCross",
    Hh[MassiveBlock[
       lineFunctionPreset018[line], Lookup[line, "nu", nu], Lookup[line, "skType", "++"],
      lineId, endpoints, momentumMagnitude, pack[[positions[[1]]]], pack[[positions[[2]]]]
      ]],
    "masslessFull",
    Hh[MasslessBlock[
      Lookup[line, "skType", "++"], lineId, endpoints, momentumMagnitude,
      pack[[positions[[1]]]], pack[[positions[[2]]]]
      ]],
    "masslessCross",
    Hh[MasslessCrossBlock[Lookup[line, "skType", "+-"], lineId, endpoints, momentumMagnitude]],
    _, 1
    ]
   ];


integrandLineBareFactor018[topo_Association, int_J, e_Integer] := Module[
   {line, pack, packType, momentumMagnitude, denominator},
   line = topo["lines"][[e]];
   pack = int[[2, e]];
   packType = actualLinePackType[topo, e, pack];
   momentumMagnitude = lineMomentumMagnitude[topo, e];
   denominator = If[
     lineIndexedPowerQ[line],
     momentumMagnitude^(-linePowerIndex[topo, int, e]),
     1
     ];
   If[
    packType === "shrunk",
    denominator,
    denominator integrandBuildingBlock[Join[line, <|"packType" -> packType|>], pack, momentumMagnitude]
    ]
   ];


integralToBareInertIntegrand018[topo_Association, int_J] := Times[
   integrandVertexFactor[topo, int],
   Times @@ Table[integrandLineBareFactor018[topo, int, e], {e, topo["nE"]}],
   integrandISPFactor[topo, int]
   ];


rep2Integrand[expr_, topoSpec_Association] := Module[
   {rootTopo, internalExpr, result, sectorTopo, prefactorData, prefactor},
   rootTopo = resolvePublicTopologyContext[topoSpec];
   If[rootTopo === $Failed, Return[$Failed]];
   internalExpr = If[
     Lookup[rootTopo, "ibpMode", "full"] === "timeOnly",
     dsTimeOnlyExpressionToInternal020[expr, rootTopo],
     expr
     ];
   If[internalExpr === $Failed || ! validatePublicExpression[internalExpr, rootTopo], Return[$Failed]];
   result = Expand[internalExpr /. int_J :> (
         sectorTopo = sectorTopologyForIntegral018[rootTopo, int];
         prefactorData = sectorPrefactorDataForIntegral018[rootTopo, int];
         prefactor = materializeSectorPrefactor018[prefactorData];
         If[
          Head[sectorTopo] === Missing || prefactor === $Failed,
          $Failed,
          prefactor integralToBareInertIntegrand018[sectorTopo, int]
          ]
         )];
   If[! FreeQ[result, $Failed], Return[$Failed]];
   rep2outform[result, rootTopo]
   ];


(* ::Chapter:: *)
(*GF(2) parity metadata 与 sector transport*)

normalizeParityConstraints018[constraints_] := Map[
   Function[item,
    Which[
     MatchQ[item, _Rule | _RuleDelayed], <|"expression" -> First[item], "remainder" -> Last[item]|>,
     AssociationQ[item] && KeyExistsQ[item, "expression"],
     <|"expression" -> item["expression"], "remainder" -> Lookup[item, "remainder", 0]|>,
     True, <|"status" -> "invalid", "input" -> item|>
     ]
    ],
   If[ListQ[constraints], constraints, {constraints}]
   ];


lineFunctionPreset018[line_Association] := Lookup[
   lineCompiledFunctionSystem[line],
   "preset",
    Lookup[Lookup[lineCompiledFunctionSystem[line], "input", <||>], "preset", "custom"]
   ];


(* Parity transport only needs every line building block to have a proved GF(2)
   closure. Massive h/H and massless exponential lines both satisfy this contract;
   an unknown custom function system remains fail closed. *)
parityLineFunctionSystemUsableQ018[line_Association] := Switch[
   Lookup[line, "massType", "massive"],
   "massive", MemberQ[{"h", "H"}, lineFunctionPreset018[line]],
    "massless", True,
   _, False
   ];


parityFunctionSystemUsableQ018[topo_Association] := Module[{lines},
   lines = Lookup[topo, "lines", {}];
   lines =!= {} && And @@ (parityLineFunctionSystemUsableQ018 /@ lines)
   ];


parityLineVariables018[line_Association] := {
   b[line["id"]], n[line["id"], 1], n[line["id"], 2]
   };


transportParityConstraint018[topo_Association, constraint_Association] := Module[
   {expr, remainder, rootRules, issues = {}, line, vars, coefficients, cb, c1, c2,
    integerShift, zeroShift, sourceZero, targetZero, delta, shrunkQ},
   If[Lookup[constraint, "status", "valid"] === "invalid", Return[constraint]];
   expr = Expand[constraint["expression"]];
   remainder = constraint["remainder"];
   rootRules = Lookup[topo, "rootZeroPointRules", topo["zeroPointRules"]];
   Do[
    line = topo["lines"][[e]];
    vars = parityLineVariables018[line];
    coefficients = Coefficient[expr, #] & /@ vars;
    {cb, c1, c2} = coefficients;
    If[And @@ (zeroQ /@ coefficients), Continue[]];
    If[! lineIndexedPowerQ[line] && ! And @@ (zeroQ /@ coefficients),
     AppendTo[issues, <|"reason" -> "fixedLineInParity", "lineIndex" -> e, "coefficients" -> coefficients|>]
     ];
    shrunkQ = Lookup[line, "state", "full"] === "shrunk" || Lookup[line, "packType", ""] === "shrunk";
    If[shrunkQ && lineIndexedPowerQ[line],
     If[! TrueQ[Expand[c1 - c2] === 0],
      AppendTo[issues, <|"reason" -> "endpointWeightsNotTransportable", "lineIndex" -> e, "coefficients" -> coefficients|>],
      integerShift = lineShrinkBShift[line];
      zeroShift = lineShrinkZeroPointShift[line];
      sourceZero = zeroPointRuleValue018[rootRules, b0[line["id"]], 0];
      targetZero = lineBSZeroPoint[topo, e];
      delta = Simplify[Expand[(targetZero - sourceZero) - zeroShift]];
      If[!(TrueQ[delta === 0] || IntegerQ[delta]),
       AppendTo[issues, <|"reason" -> "nonIntegerZeroPointRebase", "lineIndex" -> e, "delta" -> delta|>],
       expr = Expand[
         expr - cb vars[[1]] - c1 vars[[2]] - c2 vars[[3]] +
          cb (bS[line["id"]] - integerShift + delta) + c1
         ]
       ]
      ]
     ],
    {e, topo["nE"]}
    ];
   If[issues =!= {},
    <|"status" -> "invalid", "input" -> constraint, "issues" -> issues|>,
    <|"status" -> "transported", "expression" -> Expand[expr], "remainder" -> remainder|>
    ]
   ];


parityMetadataForSector018[topo_Association] := Module[
   {raw, normalized, transported, invalid, usableFunctionQ},
   raw = Lookup[topo, "parityConstraints", {}];
   usableFunctionQ = parityFunctionSystemUsableQ018[topo];
   If[raw === {} || raw === None,
    Return[<|"status" -> "disabled", "reason" -> "noParityConstraints",
      "parityUsableQ" -> usableFunctionQ, "constraints" -> {}|>]
    ];
   If[! usableFunctionQ,
    Return[<|"status" -> "disabled", "reason" -> "unsupportedParityFunctionSystem",
      "parityUsableQ" -> False, "constraints" -> {}|>]
    ];
   normalized = normalizeParityConstraints018[raw];
   transported = transportParityConstraint018[topo, #] & /@ normalized;
   invalid = Select[transported, Lookup[#, "status", "invalid"] =!= "transported" &];
   If[invalid =!= {},
    <|"status" -> "disabled", "reason" -> "constraintTransportFailed",
      "parityUsableQ" -> False, "constraints" -> {}, "issues" -> invalid|>,
    <|"status" -> "enabled", "parityUsableQ" -> True,
      "constraints" -> transported,
      "masslessCycleFlipCount" -> Mod[Count[
         Select[Range[topo["nE"]], MemberQ[Lookup[topo, "sectorShrunkLines", {}], #] &],
         e_ /; lineIndexedPowerQ[topo["lines"][[e]]] && Lookup[topo["lines"][[e]], "massType", "massive"] === "massless"
         ], 2]
      |>
    ]
   ];


parityIntegralIndexRules018[metadata_Association, int_J] := Module[{lineRules, ispRules},
   lineRules = Flatten@MapThread[
      Thread[Lookup[#1, "packTemplate", {}] -> #2] &,
      {Lookup[metadata, "lineSlots", {}], int[[2]]}
      ];
   ispRules = If[Lookup[metadata, "ispSlots", {}] === {},
     {},
     MapThread[Lookup[#1, "indexSymbol"] -> #2 &, {metadata["ispSlots"], int[[3]]}]
     ];
   Join[lineRules, ispRules]
   ];


parityIntegralSignature018[metadata_Association, int_J] := Module[{data, rules},
   data = Lookup[metadata, "parityData", <|"status" -> "disabled"|>];
   If[Lookup[data, "status", "disabled"] =!= "enabled", Return[{}]];
   rules = parityIntegralIndexRules018[metadata, int];
   Mod[Expand[(#expression - #remainder) /. rules], 2] & /@ Lookup[data, "constraints", {}]
   ];


parityIntegralAllowedQ018[metadata_Association, int_J] := Module[{signature},
   signature = parityIntegralSignature018[metadata, int];
   signature === {} || And @@ (TrueQ[# === 0] & /@ signature)
   ];


dsParityFilteredPointRules018[entry_Association, pointRules_List, context_Association] := Module[
   {metadata, parityData, sourceIntegral},
   metadata = Lookup[entry, "sectorMetadata", Missing["NoSectorMetadata"]];
   If[! AssociationQ[metadata], Return[pointRules]];
   parityData = Lookup[metadata, "parityData", <|"status" -> "disabled"|>];
   If[Lookup[parityData, "status", "disabled"] =!= "enabled", Return[pointRules]];
   sourceIntegral = Lookup[entry, "sourceIntegral", Missing["NoSourceIntegral"]];
   If[Head[sourceIntegral] =!= J, Return[pointRules]];
   Select[pointRules, parityIntegralAllowedQ018[metadata, sourceIntegral /. #] &]
   ];


(* 用户直接传入表达式时没有 sourceIntegral/sector parity provenance，不能在撒点前猜测
   seed parity。此时保留用户点域，并由生成后的 parity certificate 检查实际积分。 *)
dsParityFilteredPointRules018[_, pointRules_List, _Association] := pointRules;


dsParityCertificate018[records_List, context_Association] := Module[
   {metadataList = context["sectors"], failures = {}, integrals, metadata, signature},
   Do[
    integrals = DeleteDuplicates[Cases[Lookup[records[[i]], "equation", 0], _J, {0, Infinity}]];
    Do[
     metadata = SelectFirst[metadataList, integralMatchesSectorMetadataQ[int, #] &, Missing["NoMatchingSector"]];
     If[Head[metadata] === Missing,
      AppendTo[failures, <|"equationIndex" -> i, "integral" -> int, "reason" -> metadata|>],
      signature = parityIntegralSignature018[metadata, int];
      If[signature =!= {} && ! And @@ (TrueQ[# === 0] & /@ signature),
       AppendTo[failures, <|"equationIndex" -> i, "sectorKey" -> metadata["sectorKey"],
         "integral" -> int, "signature" -> signature|>]
       ]
      ],
     {int, integrals}
     ],
    {i, Length[records]}
    ];
   <|"passQ" -> (failures === {}), "checkedEquationCount" -> Length[records],
    "failureCount" -> Length[failures], "failures" -> failures|>
   ];


(* ::Chapter:: *)
(*Massless tree formula capability*)

treeFormulaMasslessLines018[context_Association] := Lookup[
   Select[
    Lookup[context["topology"], "lines", {}],
    Lookup[#, "massType", "massive"] === "massless" &
    ],
   "id",
   {}
   ];


treeFormulaMasslessPendingQ018[context_Association] := treeFormulaMasslessLines018[context] =!= {};


treeFormulaPendingRederivation018[operation_String, context_Association] := Module[{lines},
   lines = treeFormulaMasslessLines018[context];
   dsErrorPrint[
    operation <> " 尚未在 massless 三槽 quotient basis 上重新推导，公式型 tree 路线已停止。" <>
     " " <> operation <> " has not been rederived on the massless three-slot quotient basis; the formula-based tree route was stopped."
    ];
   <|
    "status" -> "PendingRederivation",
    "reason" -> "masslessQuotientFormulaNotCertified",
    "operation" -> operation,
    "masslessLineIds" -> lines,
    "representation" -> "J[aList,linePacks,ispList]",
    "availableAlternative" ->
     "Use DSSeeds -> DSGenerateIBP -> DSLinear on the unified line-pack basis; construct DSDE only after importing an external reduction."
    |>
   ];


(* ::Chapter:: *)
(*template-only DSSeeds*)

DSSeeds[context_: Automatic, opts : OptionsPattern[]] := Module[
   {resolved, seedSkeleton, templateData, sealedTemplates, seedGroups, seedGroupMetadata,
    seedRangeMetadata, discoveredIndices, progress = OptionValue[ProgressReporting]},
   resolved = dsResolveContext[context];
   If[Head[resolved] === Missing,
    Message[DSSeeds::noinit];
    dsErrorPrint["请先成功调用 DSInit。 Run DSInit successfully first."];
    Return[<|"status" -> "failed", "reason" -> "missingContext"|>]
    ];
   If[! dsContextCapabilityQ[resolved, "timeIBPUsableQ"] ||
     (Lookup[resolved["topology"], "ibpMode", "full"] === "full" &&
       ! dsContextCapabilityQ[resolved, "momentumIBPUsableQ"]),
    Message[DSSeeds::capability, dsContextCapabilities[resolved]];
    Return[<|"status" -> "failed", "reason" -> "capabilityGate"|>]
    ];
    seedSkeleton = <|
      "status" -> "generated",
      "representation" -> If[
        Lookup[resolved["topology"], "ibpMode", "full"] === "timeOnly",
        "J[sectorKey,timeShifts,stateBits]",
        "J[aList,linePacks,ispList]"
        ],
     "ibpMode" -> Lookup[resolved["topology"], "ibpMode", "full"],
     "sectorMetadata" -> First[resolved["sectors"]],
     "sectorMetadataList" -> resolved["sectors"]
     |>;
   templateData = dsStageRun[
     "构造全部 reachable-sector 离散态 seed 模板 / Building all reachable-sector discrete-state seed templates",
      dsLoopSeedTemplateData[resolved, seedSkeleton],
      progress
      ];
   If[Lookup[templateData, "status", "failed"] =!= "generated",
    Message[DSSeeds::failed, Lookup[templateData, "reason", "templateGenerationFailed"]];
    Return[Join[seedSkeleton, <|"status" -> "failed", "templateData" -> templateData|>]]
    ];
    If[Lookup[resolved["topology"], "ibpMode", "full"] === "timeOnly",
     templateData = Join[templateData, <|
        "generationRoute" -> "internalThreeSlotConverted",
        "representation" -> "J[sectorKey,timeShifts,stateBits]",
        "allSeeds" -> dsTimeOnlyDataToPublic020[Lookup[templateData, "allSeeds", {}], resolved]
        |>];
     If[Lookup[templateData, "allSeeds", $Failed] === $Failed,
      Message[DSSeeds::failed, "timeOnlyPublicConversionFailed"];
      Return[Join[seedSkeleton, <|"status" -> "failed", "reason" -> "timeOnlyPublicConversionFailed"|>]]
      ]
     ];
    sealedTemplates = dsSealSeedTemplates[templateData["allSeeds"], resolved];
   seedGroups = dsDefaultSeedGroups[sealedTemplates];
   seedGroupMetadata = dsSeedGroupMetadataFromGroups[seedGroups];
   discoveredIndices = DeleteDuplicates[Flatten[dsEntrySeedVariables /@ sealedTemplates]];
   seedRangeMetadata = DSMetaSeedRange[seedGroups, discoveredIndices];
   $dSIBPLastSeedTemplates = sealedTemplates;
   $dSIBPLastSeedGroups = seedGroups;
   $dSIBPLastSeedGroupMetadata = seedGroupMetadata;
   Join[seedSkeleton, <|
     "completeCanonicalQ" -> True,
     "completeMomentumIBPQ" -> (Lookup[resolved["topology"], "ibpMode", "full"] === "full"),
     "completeTimeIBPQ" -> True,
     "pendingFeatures" -> {},
     "forbiddenNData" -> {},
     "equationCount" -> 0,
     "equations" -> {},
     "allSeeds" -> sealedTemplates,
     "seedGroups" -> seedGroups,
     "seedGroupMetadata" -> seedGroupMetadata,
     "seedRangeMetadata" -> seedRangeMetadata,
     "seedTemplateSummary" -> KeyDrop[templateData, "allSeeds"],
     "dSIBPStatus" -> "generated",
     "dSIBPContextSummary" -> dsContextSummary[resolved]
     |>]
   ];
