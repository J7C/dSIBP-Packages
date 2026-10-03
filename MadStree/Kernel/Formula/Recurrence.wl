(* ::Package:: *)

(***
文件：Recurrence.wl
用途：构造 Eq. (3.65) contact maps，并用 (-tau)^A 约定执行单步和完整迭代约化。
核心边界：contact 始终携带精确 target sector 与 time shift；每个待约化积分还必须通过初始化时
冻结的小 t Gamma 参数证书，不能跨入非正整数发散层。
***)

(* ::Chapter:: *)
(* Contact matrices *)

msLineById[context_?MSContextQ, lineId_] := SelectFirst[
  context["lines"], #["id"] === lineId &, Missing["UnknownLine", lineId]
];

msTargetSectorForEvent[event_Association, context_?MSContextQ] := msSectorByKey[
  context,
  event["targetSector"]
];

msBitsIndexAssociation[states_List] := AssociationThread[
  ToString[#, InputForm] & /@ states,
  Range[Length[states]]
];

msSpectatorTargetBits[sourceSector_Association, targetSector_Association, sourceBits_List] := Module[
  {sourcePositionByKey, targetKeys},
  sourcePositionByKey = AssociationThread[
    ToString[#, InputForm] & /@ Lookup[sourceSector["slots"], "key"],
    Range[sourceSector["slotCount"]]
  ];
  targetKeys = If[targetSector["slots"] === {}, {}, Lookup[targetSector["slots"], "key"]];
  sourceBits[[sourcePositionByKey[ToString[#, InputForm]]]] & /@ targetKeys
];

msCanonicalTargetBits[targetSector_Association, rawBits_List] := Module[{row, imageRow, column},
  row = FirstPosition[targetSector["rawStateOrder"], rawBits, Missing["UnknownRawState"]];
  If[Head[row] === Missing, Return[row]];
  imageRow = targetSector["stateEmbedding"][[First[row]]];
  column = FirstPosition[imageRow, 1, Missing["ZeroQuotientState"]];
  If[Head[column] === Missing, column, targetSector["stateOrder"][[First[column]]]]
];

msEventEndpointIndex[line_Association, sourceSector_Association, componentPosition_Integer] := Module[
  {positions = sourceSector["rootToComponent"] /@ line["endpoints"]},
  First@FirstPosition[positions, componentPosition]
];


msEventAtomicCoefficient[
  line_Association,
  sourceSector_Association,
  sourceBits_List,
  componentPosition_Integer
] := Module[{endpointIndex, positions, endpointBits, sharedPosition},
  endpointIndex = msEventEndpointIndex[line, sourceSector, componentPosition];
  Switch[line["type"],
    "massiveFull",
      positions = (First@msSlotPositionByKey[sourceSector, {line["id"], #}] &) /@ {1, 2};
      endpointBits = sourceBits[[positions]];
      If[Total[endpointBits] =!= 1,
        0,
        If[endpointIndex === 1, (-1)^endpointBits[[2]], (-1)^endpointBits[[1]]]
      ],
    "masslessFull",
      sharedPosition = First@msSlotPositionByKey[sourceSector, {line["id"], "shared"}];
      If[sourceBits[[sharedPosition]] =!= 1, 0, If[endpointIndex === 1, -2, 2]],
    _, 0
  ]
];


msNormalizedEventAbsorption[
  sourceSector_Association,
  targetSector_Association,
  selectedLines_List
] := Simplify[
  sourceSector["normalization"]/targetSector["normalization"] Times @@ Lookup[selectedLines, "pinchNormalization"]
];


(* 纯 massive contact 不再乘轮廓号；轮廓只进入有符号能量、链式因子和最终 SK branch 权重。
   massless quotient 相对 h/Wronskian 原子仍为每条收缩线一个负号。 *)
msContactEventSign[contactMap_Association] :=
  (-1)^Count[contactMap["lineTypes"], "masslessFull"];


msEventContactMatrix[
  sourceSector_Association,
  targetSector_Association,
  event_Association,
  componentPosition_Integer,
  context_?MSContextQ
] := Module[
  {sourceStates, targetStates, targetIndex, selectedLines, matrix, targetBits, canonicalBits,
   column, coefficient, absorption},
  sourceStates = sourceSector["stateOrder"];
  targetStates = targetSector["stateOrder"];
  targetIndex = msBitsIndexAssociation[targetStates];
  selectedLines = msLineById[context, #] & /@ event["selectedLineIds"];
  matrix = ConstantArray[0, {Length[sourceStates], Length[targetStates]}];
  absorption = msNormalizedEventAbsorption[sourceSector, targetSector, selectedLines];
  Do[
    coefficient = event["thetaBundleCoefficient"] Times @@ (
      msEventAtomicCoefficient[#, sourceSector, sourceStates[[row]], componentPosition] & /@ selectedLines
    );
    If[! TrueQ[coefficient === 0],
      targetBits = msSpectatorTargetBits[sourceSector, targetSector, sourceStates[[row]]];
      canonicalBits = msCanonicalTargetBits[targetSector, targetBits];
      If[Head[canonicalBits] =!= Missing,
        column = targetIndex[ToString[canonicalBits, InputForm]];
        matrix[[row, column]] += Simplify[absorption coefficient]
      ]
    ],
    {row, Length[sourceStates]}
  ];
  matrix
];

msContactTargetShiftEvent[
  sourceSector_Association,
  targetSector_Association,
  sourceShift_List
] := Module[{shift, sourceComponentsInTarget},
  Map[
    Function[targetComponent,
      sourceComponentsInTarget = MapIndexed[
        If[Intersection[#1, targetComponent] =!= {}, First[#2], Nothing] &,
        sourceSector["vertexComponents"]
      ];
      shift = Total[sourceShift[[sourceComponentsInTarget]]];
      shift - (Length[sourceComponentsInTarget] - 1)
    ],
    targetSector["vertexComponents"]
  ]
];

msContactMapForEvent[sourceSector_Association, event_Association, context_?MSContextQ] := Module[
  {targetSector, endpointComponents, matrices, zeroShift, r1TargetShifts, selectedLines},
  targetSector = msTargetSectorForEvent[event, context];
  endpointComponents = event["componentPair"];
  selectedLines = msLineById[context, #] & /@ event["selectedLineIds"];
  matrices = Association@Table[
    componentPosition -> msEventContactMatrix[sourceSector, targetSector, event, componentPosition, context],
    {componentPosition, endpointComponents}
  ];
  zeroShift = ConstantArray[0, Length[sourceSector["vertexComponents"]]];
  r1TargetShifts = Association@Map[
    Function[componentPosition,
      componentPosition -> msContactTargetShiftEvent[
        sourceSector,
        targetSector,
        ReplacePart[zeroShift, componentPosition -> 1]
      ]
    ],
    DeleteDuplicates[endpointComponents]
  ];
  <|
    "sourceSector" -> sourceSector["sectorKey"],
    "targetSector" -> targetSector["sectorKey"],
    "eventId" -> event["eventId"],
    "lineId" -> If[Length[event["selectedLineIds"]] === 1, First[event["selectedLineIds"]], Missing["MultiEdgeEvent"]],
    "selectedLineIds" -> event["selectedLineIds"],
    "bundleLineIds" -> event["bundleLineIds"],
    "thetaBundleCoefficient" -> event["thetaBundleCoefficient"],
    "lineTypes" -> Lookup[selectedLines, "type"],
    "endpointComponents" -> endpointComponents,
    "matricesByComponent" -> matrices,
    "R0TargetShift" -> msContactTargetShiftEvent[sourceSector, targetSector, zeroShift],
    "R1TargetShiftsByComponent" -> r1TargetShifts,
    "normalizationRatio" -> Simplify[sourceSector["normalization"]/targetSector["normalization"]],
    "pinchNormalization" -> Times @@ Lookup[selectedLines, "pinchNormalization"],
    "absorbedNormalization" -> msNormalizedEventAbsorption[sourceSector, targetSector, selectedLines]
  |>
];

msContactMapsForSector[sourceSector_Association, context_?MSContextQ] := Map[
  msContactMapForEvent[sourceSector, #, context] &,
  Select[context["contactTransitions"], #["sourceSector"] === sourceSector["sectorKey"] &]
];

MSContactMaps[context_?MSContextQ, All] := AssociationThread[
  context["sectorOrder"] -> (msContactMapsForSector[#, context] & /@ context["sectors"])
];

MSContactMaps[context_?MSContextQ, key_String] := Module[{sector = msSectorByKey[context, key]},
  If[Head[sector] === Missing,
    Message[MSContactMaps::nosector, key];
    msFailure["UnknownSector", <|"sectorKey" -> key|>],
    msContactMapsForSector[sector, context]
  ]
];

MSContactMaps[context_?MSContextQ] := MSContactMaps[context, All];


(* ::Chapter:: *)
(* Formula recurrence *)

msIntegralVector[sector_Association, shifts_List] := MSIntegral[
  sector["sectorKey"],
  shifts,
  #
] & /@ sector["stateOrder"];

msResolveComponentPosition[sector_Association, component_] := Which[
  IntegerQ[component] && 1 <= component <= Length[sector["vertexComponents"]], component,
  ListQ[component], First@FirstPosition[sector["vertexComponents"], component, Missing["UnknownComponent"]],
  True, Missing["UnknownComponent"]
];

msRemainingVector[
  sector_Association,
  sourceShift_List,
  componentPosition_Integer,
  context_?MSContextQ
] := Module[{result, maps, map, matrix, targetSector, targetShift, targetVector, contactSign},
  result = ConstantArray[0, sector["masterCount"]];
  maps = Select[
    msContactMapsForSector[sector, context],
    KeyExistsQ[#["matricesByComponent"], componentPosition] &
  ];
  Do[
    map = contactMap;
    matrix = map["matricesByComponent"][componentPosition];
    targetSector = msSectorByKey[context, map["targetSector"]];
    targetShift = msContactTargetShiftEvent[
      sector,
      targetSector,
      sourceShift
    ];
    targetVector = msIntegralVector[targetSector, targetShift];
    contactSign = msContactEventSign[map];
    (* M0/M1 递推式已有自身的 contact 符号；这里的 event sign 只处理 massless quotient。 *)
    result += -contactSign matrix.targetVector,
    {contactMap, maps}
  ];
  result
];

msRecurrenceTowardZero[int : MSIntegral[_, _, _], component_, context_?MSContextQ] := Module[
  {data, sector, shifts, bits, componentPosition, stateRow, currentShift, nearerShift,
    regularMatrix, contactMatrix, remainingVector, nearerVector, result, denominators,
    m0Inverse, timeIBPIssues},
  data = msIntegralData[int, context];
  If[Head[data] === Failure, Return[data]];
  sector = data["sector"];
  shifts = data["shifts"];
  bits = data["bits"];
  timeIBPIssues = msTimeIBPShiftIssues[sector, shifts, bits, context];
  If[timeIBPIssues =!= {},
    Return[msFailure[
      "TimeIBPShiftHitsDivergentLayer",
      <|"input" -> int, "issues" -> timeIBPIssues|>
    ]]
  ];
  componentPosition = msResolveComponentPosition[sector, component];
  If[Head[componentPosition] === Missing, Return[msFailure["UnknownComponent", <|"component" -> component|>]]];
  currentShift = shifts[[componentPosition]];
  If[currentShift === 0, Return[msFailure["ZeroShiftAtComponent", <|"component" -> componentPosition|>]]];
  stateRow = First@FirstPosition[sector["stateOrder"], bits];
  If[currentShift > 0,
    m0Inverse = msM0Inverse[sector, componentPosition];
    If[Head[m0Inverse] === Failure, Return[m0Inverse]];
    nearerShift = ReplacePart[shifts, componentPosition -> currentShift - 1];
    regularMatrix = Simplify[
      m0Inverse.msM1Matrix[sector, componentPosition, currentShift]
    ];
    remainingVector = msRemainingVector[sector, shifts, componentPosition, context];
    contactMatrix = -m0Inverse;
    denominators = msEnergyLetters[sector, componentPosition],
    nearerShift = ReplacePart[shifts, componentPosition -> currentShift + 1];
    regularMatrix = Simplify[
      msM1Inverse[sector, componentPosition, currentShift + 1].msM0Matrix[sector, componentPosition]
    ];
    remainingVector = msRemainingVector[sector, nearerShift, componentPosition, context];
    contactMatrix = msM1Inverse[sector, componentPosition, currentShift + 1];
    denominators = msM1SingularSurfaces[sector, componentPosition, currentShift + 1]
  ];
  If[Head[denominators] === Failure, Return[denominators]];
  nearerVector = msIntegralVector[sector, nearerShift];
  result = Expand[
    regularMatrix[[stateRow]].nearerVector + contactMatrix[[stateRow]].remainingVector
  ];
  <|
    "status" -> "reducedOneStep",
    "input" -> int,
    "componentPosition" -> componentPosition,
    "direction" -> If[currentShift > 0, "raiseFormulaTowardZero", "lowerFormulaTowardZero"],
    "nearerShift" -> nearerShift,
    "result" -> result,
    "singularSurfaces" -> denominators
  |>
];

MSRecurrenceStep[int : MSIntegral[_, _, _], component_, context_?MSContextQ] := msRecurrenceTowardZero[
  int,
  component,
  context
];

MSRecurrenceStep[int : MSIntegral[_, shifts_List, _], context_?MSContextQ] := Module[
  {position = SelectFirst[Range[Length[shifts]], shifts[[#]] =!= 0 &, Missing["AlreadyMaster"]]},
  If[Head[position] === Missing,
    <|"status" -> "alreadyMaster", "input" -> int, "result" -> int|>,
    msRecurrenceTowardZero[int, position, context]
  ]
];

MSRecurrenceStep[other_, ___] := (
  Message[MSRecurrenceStep::badint, HoldForm[other]];
  msFailure["InvalidRecurrenceInput", <|"input" -> HoldForm[other]|>]
);

Options[MSReduce] = {MasterBasis -> Automatic};

msRequestedMasterBasis[context_?MSContextQ, Automatic] := Lookup[context["masters"], "integral"];
msRequestedMasterBasis[context_?MSContextQ, records_List] := Module[
  {basis, expected = Lookup[context["masters"], "integral"]},
  basis = Replace[records, item_Association :> Lookup[item, "integral", Missing["Integral"]], {1}];
  If[
    MemberQ[basis, _Missing] || ! DuplicateFreeQ[basis] ||
      Sort[ToString[#, InputForm] & /@ basis] =!= Sort[ToString[#, InputForm] & /@ expected],
    msFailure["MasterBasisMustBeCompletePermutation", <|"expected" -> expected, "actual" -> basis|>],
    basis
  ]
];

MSReduce[expr_, context_?MSContextQ, OptionsPattern[]] := Module[
  {memo = <||>, reduceIntegral, reduceScalar, reduceTree, fieldTree, resultTree,
   masterBasis, leafRecords, failedRecords, incompleteRecords, status},
  masterBasis = msRequestedMasterBasis[context, OptionValue[MasterBasis]];
  If[Head[masterBasis] === Failure, Return[masterBasis]];

  reduceIntegral[int : MSIntegral[_, shifts_List, _], stack_List] := Module[
    {key, data, step, childIntegrals, childEntries, childFailure, reduced, layer, entry},
    key = ToString[HoldComplete[int], InputForm];
    If[KeyExistsQ[memo, key], Return[memo[key]]];
    If[MemberQ[stack, key],
      Message[MSReduce::cycle, int];
      Return[msFailure["RecurrenceCycle", <|"integral" -> int|>]]
    ];
    data = msIntegralData[int, context];
    If[Head[data] === Failure, Return[data]];
    If[And @@ (# === 0 & /@ shifts),
      entry = <|"result" -> int, "singularLayers" -> {}|>;
      AssociateTo[memo, key -> entry];
      Return[entry]
    ];
    step = MSRecurrenceStep[int, context];
    If[Head[step] === Failure, Return[step]];
    childIntegrals = DeleteDuplicates@Cases[{step["result"]}, _MSIntegral, Infinity];
    childEntries = reduceIntegral[#, Append[stack, key]] & /@ childIntegrals;
    childFailure = FirstCase[childEntries, _Failure, Missing["NoFailure"], Infinity];
    If[Head[childFailure] === Failure, Return[childFailure]];
    reduced = Expand[
      step["result"] /. Thread[childIntegrals -> Lookup[childEntries, "result"]]
    ];
    layer = <|
      "integral" -> int,
      "componentPosition" -> step["componentPosition"],
      "direction" -> step["direction"],
      "surfaces" -> step["singularSurfaces"]
    |>;
    entry = <|
      "result" -> reduced,
      "singularLayers" -> DeleteDuplicates@Join[
        {layer},
        Flatten[Lookup[childEntries, "singularLayers", {}], 1]
      ]
    |>;
    AssociateTo[memo, key -> entry];
    entry
  ];

  reduceScalar[scalar_, position_] := Module[
    {integrals, entries, failure, result, variables, replaced, coefficients, residual,
     remainingShifted, masterRules, scalarStatus, scalarLayers, record},
    integrals = DeleteDuplicates@Cases[{scalar}, _MSIntegral, Infinity];
    entries = reduceIntegral[#, {}] & /@ integrals;
    failure = FirstCase[entries, _Failure, Missing["NoFailure"], Infinity];
    If[Head[failure] === Failure,
      If[position === None,
        Return[failure],
        Return[<|
          "kind" -> "MSReduceElementRecord",
          "status" -> "failed",
          "position" -> position,
          "input" -> scalar,
          "failure" -> failure,
          "result" -> failure,
          "masterBasis" -> masterBasis,
          "masterRules" -> failure,
          "coefficientVector" -> failure,
          "nonMasterResidual" -> failure,
          "remainingShiftedIntegrals" -> failure,
          "singularLayers" -> {}
        |>]
      ]
    ];
    result = Expand[scalar /. Thread[integrals -> Lookup[entries, "result"]]];
    scalarLayers = DeleteDuplicates@Flatten[Lookup[entries, "singularLayers", {}], 1];
    variables = Array[Unique["msMasterCoefficient"] &, Length[masterBasis]];
    replaced = Expand[result /. Thread[masterBasis -> variables]];
    coefficients = Simplify[Coefficient[replaced, #] & /@ variables];
    residual = Simplify[
      Together[(replaced - coefficients.variables) /. Thread[variables -> masterBasis]]
    ];
    remainingShifted = DeleteDuplicates@Cases[
      {result},
      int : MSIntegral[_, shifted_List, _] /; AnyTrue[shifted, # =!= 0 &] :> int,
      Infinity
    ];
    masterRules = Thread[masterBasis -> coefficients];
    scalarStatus = If[
      remainingShifted === {} && TrueQ[residual === 0],
      "reduced",
      "partiallyReduced"
    ];
    record = <|
      "status" -> scalarStatus,
      "input" -> scalar,
      "result" -> result,
      "masterBasis" -> masterBasis,
      "masterRules" -> masterRules,
      "coefficientVector" -> coefficients,
      "nonMasterResidual" -> residual,
      "memoizedIntegralCount" -> Length[memo],
      "remainingShiftedIntegrals" -> remainingShifted,
      "singularLayers" -> scalarLayers
    |>;
    If[
      position === None,
      record,
      Join[<|"kind" -> "MSReduceElementRecord", "position" -> position|>, record]
    ]
  ];

  If[! ListQ[expr], Return[reduceScalar[expr, None]]];

  reduceTree[list_List, path_List] := MapIndexed[
    Function[{item, index},
      If[
        ListQ[item],
        reduceTree[item, Join[path, index]],
        reduceScalar[item, Join[path, index]]
      ]
    ],
    list
  ];
  fieldTree[list_List, field_String] := fieldTree[#, field] & /@ list;
  fieldTree[record_Association, field_String] := Lookup[
    record,
    field,
    Lookup[record, "failure", Missing["ReductionField", field]]
  ];

  resultTree = reduceTree[expr, {}];
  leafRecords = Cases[
    resultTree,
    record_Association /; Lookup[record, "kind", None] === "MSReduceElementRecord",
    Infinity
  ];
  failedRecords = Select[leafRecords, Lookup[#, "status", None] === "failed" &];
  incompleteRecords = Select[leafRecords, Lookup[#, "status", None] =!= "reduced" &];
  status = Which[
    incompleteRecords === {}, "reduced",
    Length[failedRecords] === Length[leafRecords], "failed",
    True, "partiallyReduced"
  ];
  <|
    "status" -> status,
    "input" -> expr,
    "result" -> fieldTree[resultTree, "result"],
    "masterBasis" -> masterBasis,
    "masterRules" -> fieldTree[resultTree, "masterRules"],
    "coefficientVector" -> fieldTree[resultTree, "coefficientVector"],
    "nonMasterResidual" -> fieldTree[resultTree, "nonMasterResidual"],
    "remainingShiftedIntegrals" -> fieldTree[resultTree, "remainingShiftedIntegrals"],
    "singularLayers" -> fieldTree[resultTree, "singularLayers"],
    "elementResults" -> resultTree,
    "failurePositions" -> Lookup[failedRecords, "position", {}],
    "incompletePositions" -> Lookup[incompleteRecords, "position", {}],
    "memoizedIntegralCount" -> Length[memo]
  |>
];
