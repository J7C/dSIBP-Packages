(* ::Package:: *)

(***
文件：Boundary.wl
用途：按 dlog 主积分顺序，从 2411.03088 的 k0->Infinity Frobenius 解生成生产边界。
范围：统一消费通过公式 dlog 与 strict-rank chart 认证的 tree/time-only context，不按图名分派。
正确性边界：不调用 NIntegrate；每个 component 的 Gamma 参数独立检查，非正整数拒绝，其余参数
使用 Gamma 的亚纯延拓。
***)

(* ::Chapter:: *)
(* Input, point selection and common options *)

Options[MSBoundaryData] = {
  BoundaryScale -> 8,
  BoundarySeriesOrder -> 24,
  RankOrder -> Automatic,
  WorkingPrecision -> 200
};

msRuleList[rules_List] := rules;
msRuleList[association_Association] := Normal[association];
msRuleList[_] := $Failed;

msBoundaryRequiredSymbols[context_?MSContextQ] := Complement[
  DeleteDuplicates@Cases[
    {
      Lookup[context["vertices"], "externalLegEnergy"],
      Lookup[context["vertices"], "timePower"],
      If[context["lines"] === {}, {}, Lookup[context["lines"], "momentum"]],
      If[context["lines"] === {}, {}, Lookup[context["lines"], "nu"]],
      context["rootNormalization"]
    },
    symbol_Symbol /; Context[symbol] =!= "System`",
    Infinity
  ],
  msAuxiliaryExternalLegEnergySymbols[context]
];

(* The preferred point only determines the strict rank; larger damping energies correspond to smaller typical |tau|. *)
msPreferredVertexOrder[context_?MSContextQ, targetRules_List] := Lookup[
  SortBy[
    context["vertices"],
    Function[vertex,
      With[{value = N[I vertex["vertexSign"] vertex["externalLegEnergy"] /. targetRules]},
        {-Re[value], -Abs[value], vertex["position"]}
      ]
    ]
  ],
  "id"
];

msRetainNonEnergyRules[targetRules_List, vertices_List] := DeleteCases[
  targetRules,
  Rule[left_, _] /; MemberQ[Lookup[vertices, "externalLegEnergy"], left]
];


(* ::Chapter:: *)
(* Nested blow-up, theta fixing and normal-crossing certificates *)

msBoundaryRankAssociation[rankOrder_List] := AssociationThread[rankOrder -> Range[Length[rankOrder]]];


msBoundaryChartCoordinates[count_Integer] := Table[MSBlowupCoordinate[index], {index, count}];


msBoundaryEnergyChartData[context_?MSContextQ, rankOrder_List] := Module[
  {rank, coordinates, dampingSymbols, dampingRules, energyRules, vertexById, count},
  count = Length[rankOrder];
  rank = msBoundaryRankAssociation[rankOrder];
  coordinates = msBoundaryChartCoordinates[count];
  vertexById = AssociationThread[Lookup[context["vertices"], "id"] -> context["vertices"]];
  dampingSymbols = Association@Table[id -> MSDampingEnergy[id], {id, rankOrder}];
  dampingRules = Table[
    dampingSymbols[id] -> 1/Times @@ coordinates[[rank[id] ;; count]],
    {id, rankOrder}
  ];
  energyRules = Table[
    vertexById[id]["externalLegEnergy"] -> -I vertexById[id]["vertexSign"] dampingSymbols[id],
    {id, rankOrder}
  ];
  <|
    "rank" -> rank,
    "coordinates" -> coordinates,
    "dampingSymbols" -> dampingSymbols,
    "dampingRules" -> dampingRules,
    "energyRules" -> energyRules
  |>
];


msAffineLetterChartRecord[letter_, context_?MSContextQ, chart_Association] := Module[
  {dampingVariables, affineExpression, coefficientRules, activeRanks, dominantRank,
   monomial, strictTransform, boundaryValue, affineQ, unitQ, coordinates},
  dampingVariables = Values[chart["dampingSymbols"]];
  coordinates = chart["coordinates"];
  affineExpression = Together[letter /. chart["energyRules"]];
  coefficientRules = If[PolynomialQ[affineExpression, dampingVariables],
    CoefficientRules[affineExpression, dampingVariables], {}
  ];
  affineQ = PolynomialQ[affineExpression, dampingVariables] &&
    Max[0, Sequence @@ (Total[First[#]] & /@ coefficientRules)] <= 1;
  If[! TrueQ[affineQ],
    Return[<|
      "letter" -> letter,
      "status" -> "unsupportedNonAffineBoundaryDivisor",
      "affineExpression" -> affineExpression,
      "normalCrossingQ" -> False
    |>]
  ];
  activeRanks = Select[
    Range[Length[dampingVariables]],
    ! TrueQ[Simplify[Coefficient[affineExpression, dampingVariables[[#]]]] === 0] &
  ];
  dominantRank = If[activeRanks === {}, Missing["NoDampingEnergy"], Min[activeRanks]];
  monomial = If[
    Head[dominantRank] === Missing,
    1,
    1/Times @@ coordinates[[dominantRank ;; Length[coordinates]]]
  ];
  strictTransform = Simplify[(affineExpression /. chart["dampingRules"])/monomial];
  boundaryValue = Simplify[strictTransform /. Thread[coordinates -> 0]];
  unitQ = ! TrueQ[boundaryValue === 0] &&
    FreeQ[boundaryValue, Indeterminate | ComplexInfinity | DirectedInfinity];
  <|
    "letter" -> letter,
    "status" -> If[unitQ, "monomialTimesUnit", "strictTransformVanishes"],
    "affineExpression" -> affineExpression,
    "dominantRank" -> dominantRank,
    "coordinateMonomial" -> monomial,
    "strictTransform" -> strictTransform,
    "boundaryValue" -> boundaryValue,
    "normalCrossingQ" -> unitQ
  |>
];


msSectorThetaFixing[sector_Association, rank_Association] := Module[{componentDominantRank},
  componentDominantRank[component_List] := Min[rank /@ component];
  Map[
    Function[line,
      With[{components = sector["rootToComponent"] /@ line["endpoints"]},
        If[SameQ @@ components,
          <|"lineId" -> line["id"], "status" -> "coincident", "thetaUV" -> "equalTime"|>,
          With[{ranks = componentDominantRank /@ sector["vertexComponents"][[components]]},
            <|
              "lineId" -> line["id"],
              "status" -> "fixed",
              "endpointComponentRanks" -> ranks,
              "thetaUV" -> If[ranks[[1]] < ranks[[2]], 1, 0],
              "thetaVU" -> If[ranks[[2]] < ranks[[1]], 1, 0]
            |>
          ]
        ]
      ]
    ],
    Select[sector["activeLines"], msFullLineQ]
  ]
];


msDLogShiftDenominators[de_Association] := DeleteDuplicates@Flatten@Cases[
  Lookup[de, "contactBlocks", {}],
  layer_Association /; KeyExistsQ[layer, "surfaces"] :> layer["surfaces"],
  Infinity
];


msNormalizationDivisors[context_?MSContextQ] := DeleteDuplicates@Flatten[
  FactorList[Denominator[Together[#["normalization"]]]][[All, 1]] & /@ context["sectors"]
];


msSingleBoundaryChartCertificate[
  context_?MSContextQ,
  de_Association,
  rankOrder_List
] := Module[
  {chart, letters, records, theta, energyVector, jacobian, jacobianDeterminant,
   quotientQ, passedQ},
  chart = msBoundaryEnergyChartData[context, rankOrder];
  letters = DeleteDuplicates@Join[
    Lookup[de, "letters", {}],
    msDLogShiftDenominators[de],
    DeleteCases[msNormalizationDivisors[context], 1]
  ];
  records = msAffineLetterChartRecord[#, context, chart] & /@ letters;
  theta = AssociationThread[
    context["sectorOrder"] -> (msSectorThetaFixing[#, chart["rank"]] & /@ context["sectors"])
  ];
  energyVector = Lookup[
    Association[Rule @@@ ({#["externalLegEnergy"], #["externalLegEnergy"] /. chart["energyRules"] /. chart["dampingRules"]} & /@
      context["vertices"])],
    Lookup[context["vertices"], "externalLegEnergy"]
  ];
  jacobian = Outer[D, energyVector, chart["coordinates"]];
  jacobianDeterminant = Simplify[Det[jacobian]];
  quotientQ = And @@ Lookup[context["sectors"], "quotientCertifiedQ", False];
  passedQ = And @@ Lookup[records, "normalCrossingQ", True] &&
    ! TrueQ[jacobianDeterminant === 0] && quotientQ;
  <|
    "status" -> If[passedQ, "certifiedNormalCrossingChart", "normalCrossingCertificationFailed"],
    "normalCrossingQ" -> passedQ,
    "rankOrder" -> rankOrder,
    "coordinates" -> chart["coordinates"],
    "dampingDefinition" -> (I #["vertexSign"] #["externalLegEnergy"] & /@ context["vertices"]),
    "energyRules" -> chart["energyRules"],
    "dampingRules" -> chart["dampingRules"],
    "coordinateJacobian" -> jacobian,
    "coordinateJacobianDeterminant" -> jacobianDeterminant,
    "divisorRecords" -> records,
    "thetaFixingBySector" -> theta,
    "sectorQuotientCertifiedQ" -> quotientQ,
    "proofScope" -> "affine K-dependent dlog, normalization and shifted-recurrence divisors"
  |>
];


MSBoundaryChartCertificate[
  context_?MSContextQ,
  targetRulesInput_,
  OptionsPattern[{RankOrder -> Automatic}]
] := Module[{targetRules, vertexIds, requestedOrders, de, certificates},
  targetRules = msRuleList[targetRulesInput];
  If[targetRules === $Failed, Return[msFailure["NumericalRulesRequired", <||>]]];
  targetRules = Join[targetRules, msAuxiliaryExternalLegEnergyRules[context]];
  vertexIds = Lookup[context["vertices"], "id"];
  requestedOrders = Replace[
    OptionValue[RankOrder],
    {
      Automatic :> {msPreferredVertexOrder[context, targetRules]},
      All :> Permutations[vertexIds],
      order_List :> {order}
    }
  ];
  If[AnyTrue[requestedOrders, Sort[#] =!= Sort[vertexIds] || ! DuplicateFreeQ[#] &],
    Return[msFailure["InvalidRankOrder", <|"expected" -> vertexIds, "actual" -> requestedOrders|>]]
  ];
  de = MSDLogDE[context];
  certificates = msSingleBoundaryChartCertificate[context, de, #] & /@ requestedOrders;
  If[Length[certificates] === 1, First[certificates], <|"status" -> "generated", "charts" -> certificates|>]
];


(* ::Chapter:: *)
(* 2411 Sec. 3.3 single-vertex Frobenius coefficients *)

msPaperC[nu_] := 2^(-nu) Gamma[-nu]/(I Pi);

(* Returns the two endpoint branch coefficients of h=z^(+/-nu) H_nu; bit=0/1 corresponds to the MadStree state. *)
msVertexEndpointCoefficient["Negative", 1, 0, nu_] := Exp[-I Pi nu] msPaperC[nu];
msVertexEndpointCoefficient["Negative", 1, 1, nu_] := -msPaperC[-nu - 1];
msVertexEndpointCoefficient["Negative", 2, 0, nu_] := -Exp[I Pi nu] msPaperC[nu];
msVertexEndpointCoefficient["Negative", 2, 1, nu_] := msPaperC[-nu - 1];
msVertexEndpointCoefficient["Positive", 1, 0, nu_] := msPaperC[-nu];
msVertexEndpointCoefficient["Positive", 1, 1, nu_] := 2 nu Exp[-I Pi nu] msPaperC[nu];
msVertexEndpointCoefficient["Positive", 2, 0, nu_] := -msPaperC[-nu];
msVertexEndpointCoefficient["Positive", 2, 1, nu_] := -2 nu Exp[I Pi nu] msPaperC[nu];


(* ::Chapter:: *)
(* Frobenius leading system of an arbitrary sector DAG *)

(* The Hankel branches of a full propagator are determined only by the SK type and the strict time rank. *)
msBoundaryFullBranches["++", firstDominantQ_] := If[firstDominantQ, {1, 2}, {2, 1}];
msBoundaryFullBranches["--", firstDominantQ_] := If[firstDominantQ, {2, 1}, {1, 2}];


msBoundaryComponentRank[sector_Association, componentPosition_Integer, rank_Association] := Min[
  rank /@ sector["vertexComponents"][[componentPosition]]
];


(* Given an endpoint slot, returns the Hankel branch fixed at initialization or by the strict rank. *)
msBoundaryHankelBranch[
  slot_Association,
  sector_Association,
  rank_Association,
  context_?MSContextQ
] := Module[{line, endpointComponents, firstDominantQ, branches},
  line = msLineById[context, slot["lineId"]];
  Switch[line["type"],
    "massiveExternal",
      First[line["hankelBranches"]],
    "massiveCross",
      line["hankelBranches"][[slot["endpointIndex"]]],
    "massiveFull" | "masslessFull",
      endpointComponents = sector["rootToComponent"] /@ line["endpoints"];
      firstDominantQ = If[
        SameQ @@ endpointComponents,
        True,
        msBoundaryComponentRank[sector, First[endpointComponents], rank] <
          msBoundaryComponentRank[sector, Last[endpointComponents], rank]
      ];
      branches = msBoundaryFullBranches[line["contourType"], firstDominantQ];
      branches[[slot["endpointIndex"]]],
    _,
      Missing["NoHankelBranch", slot["key"]]
  ]
];


(* A quotient massless slot is a fixed exponent inside a strict chamber; bit=1 only leaves the direction sign. *)
msBoundaryMasslessSharedFactor[
  slot_Association,
  bit_Integer,
  sector_Association,
  rank_Association,
  context_?MSContextQ
] := Module[{line, endpointComponents, firstDominantQ},
  line = msLineById[context, slot["lineId"]];
  endpointComponents = sector["rootToComponent"] /@ line["endpoints"];
  If[SameQ @@ endpointComponents, Return[If[bit === 0, 1, 0]]];
  firstDominantQ = msBoundaryComponentRank[sector, First[endpointComponents], rank] <
    msBoundaryComponentRank[sector, Last[endpointComponents], rank];
  If[bit === 0, 1, If[firstDominantQ, -1, 1]]
];


(* All blow-up coordinates scale along the same t ray; t=1 is the finite matching point. *)
msGenericBoundaryCurveData[
  context_?MSContextQ,
  de_Association,
  targetRules_List,
  scale_,
  rankOrder_List,
  workingPrecision_Integer
] := Module[
  {rank, vertexCount, vertexById, momentumSizes, dampingBase, parameter, curveRules,
   retainedRules, anchorRules, omegaAlongCurve, connection, residue, convergenceRatio},
  rank = msBoundaryRankAssociation[rankOrder];
  vertexCount = Length[rankOrder];
  vertexById = AssociationThread[Lookup[context["vertices"], "id"] -> context["vertices"]];
  momentumSizes = Abs[N[Lookup[context["lines"], "momentum", {}] /. targetRules, workingPrecision]];
  dampingBase = Max[2, Ceiling[N[scale (1 + Total[momentumSizes]), 20]]];
  parameter = Unique["msBoundaryT"];
  curveRules = Table[
    With[{vertex = vertexById[id], weight = vertexCount - rank[id] + 1},
      vertex["externalLegEnergy"] -> -I vertex["vertexSign"] dampingBase^weight parameter^-weight
    ],
    {id, rankOrder}
  ];
  retainedRules = msRetainNonEnergyRules[targetRules, context["vertices"]];
  anchorRules = Join[curveRules /. parameter -> 1, retainedRules];
  omegaAlongCurve = de["omegaPotential"] /. curveRules /. retainedRules;
  connection = Map[Cancel[Together[D[#, parameter]]] &, omegaAlongCurve, {2}];
  If[! FreeQ[connection, _Real],
    Return[msFailure[
      "ExactSingularConnectionRequired",
      <|"reason" -> "generic sector pullback must remain exact Q(i)(t)"|>
    ]]
  ];
  residue = Map[Together[Limit[parameter #, parameter -> 0]] &, connection, {2}];
  If[! FreeQ[residue, Indeterminate | ComplexInfinity | DirectedInfinity],
    Return[msFailure["RegularSingularPullbackRequired", <|"residue" -> residue|>]]
  ];
  convergenceRatio = If[momentumSizes === {}, 1/dampingBase, (1 + Total[momentumSizes])/dampingBase];
  <|
    "rank" -> rank,
    "dampingBase" -> dampingBase,
    "parameter" -> parameter,
    "curveRules" -> curveRules,
    "retainedRules" -> retainedRules,
    "anchorRules" -> anchorRules,
    "connection" -> connection,
    "residue" -> residue,
    "convergenceRatio" -> convergenceRatio
  |>
];


msBoundaryAncestorSectorKeys[context_?MSContextQ, sectorKey_String] := FixedPoint[
  Function[keys,
    Union[
      keys,
      Lookup[
        Select[context["contactTransitions"], MemberQ[keys, #["targetSector"]] &],
        "sourceSector",
        {}
      ]
    ]
  ],
  {sectorKey}
];


(* Fixes the unit component of the root sector, sets non-ancestor sectors to zero, and solves all ancestor components linearly. *)
msBoundaryRootedIndicialVector[
  context_?MSContextQ,
  residue_?MatrixQ,
  exponent_,
  sectorKey_String,
  rootIndex_Integer
] := Module[
  {dimension, identity, matrix, ancestorKeys, allowedIndices, rootIndices, fixedIndices,
   fixedValues, constraints, solution, residual},
  dimension = Length[residue];
  identity = IdentityMatrix[dimension];
  matrix = exponent identity - residue;
  ancestorKeys = msBoundaryAncestorSectorKeys[context, sectorKey];
  allowedIndices = Lookup[
    Select[context["masters"], MemberQ[ancestorKeys, #["sectorKey"]] &],
    "globalIndex"
  ];
  rootIndices = Lookup[Select[context["masters"], #["sectorKey"] === sectorKey &], "globalIndex"];
  fixedIndices = Union[Complement[Range[dimension], allowedIndices], rootIndices];
  fixedValues = If[# === rootIndex, 1, 0] & /@ fixedIndices;
  constraints = UnitVector[dimension, #] & /@ fixedIndices;
  solution = Quiet@Check[
    LinearSolve[
      Join[matrix, constraints],
      Join[ConstantArray[0, dimension], fixedValues]
    ],
    $Failed
  ];
  If[solution === $Failed,
    Return[msFailure[
      "SectorLeadingSystemFailed",
      <|"sectorKey" -> sectorKey, "rootIndex" -> rootIndex, "frobeniusExponent" -> exponent|>
    ]]
  ];
  residual = Together[matrix.solution];
  If[
    ! TrueQ[And @@ (PossibleZeroQ /@ residual)] ||
      ! TrueQ[And @@ MapThread[PossibleZeroQ[#1 - #2] &, {solution[[fixedIndices]], fixedValues}]],
    Return[msFailure[
      "SectorLeadingSystemResidual",
      <|"sectorKey" -> sectorKey, "rootIndex" -> rootIndex, "residual" -> residual|>
    ]]
  ];
  Together[solution]
];


(* Leading coefficient of the defining integral of a single sector master; all data come from sector slots/components. *)
msBoundaryGammaPoleQ[value_] := msProvablyNonPositiveIntegerQ[value];


msGenericSectorLeadingRecord[
  context_?MSContextQ,
  sector_Association,
  stateBits_List,
  globalIndex_Integer,
  targetRules_List,
  curveData_Association,
  workingPrecision_Integer
] := Module[
  {rank, parameter, vertexCount, componentRecords, hSlots, timePower, alpha,
   componentRank, weight, energyExpression, energyConstant, endpointFactor, slotBit,
     line, momentum, physicalNu, formulaNu, branch, sharedFactor, coefficient, exponent,
     componentFailure},
  rank = curveData["rank"];
  parameter = curveData["parameter"];
  vertexCount = Length[context["vertices"]];
  componentRecords = Table[
    hSlots = msHankelSlotsAtComponent[sector, componentPosition];
    timePower = sector["baseTimePowers"][[componentPosition]] /. targetRules;
    alpha = Simplify[
      timePower + 1 - Total[
        Function[slot,
          slotBit = stateBits[[slot["slotPosition"]]];
          slotBit (2 (slot["formulaNu"] /. targetRules) + 1)
        ] /@ hSlots
      ]
    ];
    componentFailure = If[msBoundaryGammaPoleQ[alpha],
      msFailure[
        "BoundaryGammaPole",
        <|
          "sectorKey" -> sector["sectorKey"],
          "stateBits" -> stateBits,
          "componentPosition" -> componentPosition,
          "gammaArgument" -> alpha
        |>
      ],
      None
    ];
    If[Head[componentFailure] === Failure,
      componentFailure,
      componentRank = msBoundaryComponentRank[sector, componentPosition, rank];
      weight = vertexCount - componentRank + 1;
      energyExpression = Together[
        sector["componentEnergies"][[componentPosition]] /.
          curveData["curveRules"] /. curveData["retainedRules"]
      ];
      energyConstant = Together[Limit[parameter^weight energyExpression, parameter -> 0]];
      If[TrueQ[PossibleZeroQ[energyConstant]],
        msFailure[
          "VanishingComponentBoundaryEnergy",
          <|"sectorKey" -> sector["sectorKey"], "componentPosition" -> componentPosition|>
        ],
        endpointFactor = Times @@ Map[
          Function[slot,
            slotBit = stateBits[[slot["slotPosition"]]];
            line = msLineById[context, slot["lineId"]];
            momentum = line["momentum"] /. targetRules;
            physicalNu = line["nu"] /. targetRules;
            formulaNu = line["formulaNu"] /. targetRules;
            branch = msBoundaryHankelBranch[slot, sector, rank, context];
            (-I momentum)^(-slotBit (2 formulaNu + 1))
              msVertexEndpointCoefficient[
                context["convention", "nuConvention"], branch, slotBit, physicalNu
              ]
          ],
          hSlots
        ];
        <|
          "componentPosition" -> componentPosition,
          "timePower" -> timePower,
          "frobeniusPower" -> alpha,
          "gammaArgument" -> alpha,
          "rankWeight" -> weight,
          "energyLeadingConstant" -> energyConstant,
          "coefficient" -> Simplify[
            (-I)^(timePower + 1) Gamma[alpha] endpointFactor
          ],
          "anchorFactor" -> Simplify[energyConstant^-alpha]
        |>
      ]
    ],
    {componentPosition, Length[sector["vertexComponents"]]}
  ];
  If[Head[componentRecords] === Failure || MemberQ[componentRecords, _Failure],
    Return[FirstCase[componentRecords, _Failure, componentRecords]]
  ];
  sharedFactor = Times @@ MapIndexed[
    Function[{slot, position},
      If[
        slot["kind"] === "masslessShared",
        msBoundaryMasslessSharedFactor[
          slot, stateBits[[First[position]]], sector, rank, context
        ],
        1
      ]
    ],
    sector["slots"]
  ];
  coefficient = Simplify[
    (sector["normalization"] /. targetRules)
      sharedFactor Times @@ Lookup[componentRecords, "coefficient"]
  ];
  exponent = Simplify[Total[
    Lookup[componentRecords, "rankWeight"] Lookup[componentRecords, "frobeniusPower"]
  ]];
  <|
    "kind" -> "genericSectorDefinitionBranch",
    "sectorKey" -> sector["sectorKey"],
    "binaryState" -> stateBits,
    "rootGlobalIndex" -> globalIndex,
    "componentData" -> componentRecords,
    "frobeniusExponent" -> exponent,
    "logPower" -> 0,
    "coefficient" -> coefficient,
    "physicalWeight" -> Simplify[
      coefficient Times @@ Lookup[componentRecords, "anchorFactor"]
    ]
  |>
];


msGenericSectorFrobeniusData[
  context_?MSContextQ,
  de_Association,
  targetRules_List,
  scale_,
  rankOrder_List,
  order_Integer,
  workingPrecision_Integer
] := Module[
  {curveData, sectorByKey, localRecords, failure, branches, normalizedVector, thetaFixing},
  curveData = msGenericBoundaryCurveData[
    context, de, targetRules, scale, rankOrder, workingPrecision
  ];
  If[Head[curveData] === Failure, Return[curveData]];
  sectorByKey = AssociationThread[context["sectorOrder"] -> context["sectors"]];
  localRecords = Map[
    Function[master,
      msGenericSectorLeadingRecord[
        context,
        sectorByKey[master["sectorKey"]],
        master["stateBits"],
        master["globalIndex"],
        targetRules,
        curveData,
        workingPrecision
      ]
    ],
    context["masters"]
  ];
  failure = FirstCase[localRecords, _Failure, Missing["NoFailure"]];
  If[Head[failure] === Failure, Return[failure]];
  localRecords = Select[localRecords, ! TrueQ[PossibleZeroQ[#["physicalWeight"]]] &];
  branches = Map[
    Function[record,
      normalizedVector = msBoundaryRootedIndicialVector[
        context,
        curveData["residue"],
        record["frobeniusExponent"],
        record["sectorKey"],
        record["rootGlobalIndex"]
      ];
      If[
        Head[normalizedVector] === Failure,
        normalizedVector,
        Join[record, <|"normalizedLeadingVector" -> normalizedVector|>]
      ]
    ],
    localRecords
  ];
  failure = FirstCase[branches, _Failure, Missing["NoFailure"]];
  If[Head[failure] === Failure, Return[failure]];
  thetaFixing = AssociationThread[
    context["sectorOrder"] ->
      (msSectorThetaFixing[#, curveData["rank"]] & /@ context["sectors"])
  ];
  <|
    "boundaryKind" -> "singularFrobenius",
    "anchorRules" -> curveData["anchorRules"],
    "convergenceRatio" -> curveData["convergenceRatio"],
    "leadingBranches" -> branches,
    "seriesOrder" -> order,
    "blowupVariables" -> msBoundaryEnergyChartData[context, rankOrder],
    "singularParameter" -> curveData["parameter"],
    "singularCurveRules" -> curveData["curveRules"],
    "singularConnection" -> curveData["connection"],
    "singularResidue" -> curveData["residue"],
    "singularStart" -> 0,
    "singularTarget" -> 1,
    "thetaFixing" -> thetaFixing,
    "sectorLeadingSystemQ" -> True,
    "genericSectorProducerQ" -> True
  |>
];


(* ::Chapter:: *)
(* Public production boundary *)

MSBoundaryData[
  context_?MSContextQ,
  targetRulesInput_,
  OptionsPattern[]
] := Module[
  {de, targetRules, scale, order, workingPrecision, rankOrder, vertexIds, requiredSymbols,
   unresolved, chartCertificate, seriesData, anchorLetters},
  de = MSDLogDE[context];
  If[Lookup[de, "dlogStatus", None] =!= "certifiedByFormulaChecks",
    Return[msFailure["CertifiedDLogRequired", <|"status" -> Lookup[de, "dlogStatus", Missing["Absent"]]|>]]
  ];
  targetRules = msRuleList[targetRulesInput];
  If[targetRules === $Failed, Return[msFailure["NumericalRulesRequired", <||>]]];
  targetRules = Join[targetRules, msAuxiliaryExternalLegEnergyRules[context]];
  requiredSymbols = msBoundaryRequiredSymbols[context];
  unresolved = Select[requiredSymbols, ! NumericQ[N[# /. targetRules]] &];
  If[unresolved =!= {}, Return[msFailure["IncompleteNumericalPoint", <|"symbols" -> unresolved|>]]];
  scale = OptionValue[BoundaryScale];
  order = OptionValue[BoundarySeriesOrder];
  workingPrecision = OptionValue[WorkingPrecision];
  If[! NumericQ[scale] || ! TrueQ[N[scale] > 1],
    Return[msFailure["BoundaryScaleMustExceedOne", <|"value" -> scale|>]]
  ];
  If[! IntegerQ[order] || order < 0,
    Return[msFailure["BoundarySeriesOrderMustBeNonNegative", <|"value" -> order|>]]
  ];
  vertexIds = Lookup[context["vertices"], "id"];
  rankOrder = Replace[OptionValue[RankOrder], Automatic :> msPreferredVertexOrder[context, targetRules]];
  If[Sort[rankOrder] =!= Sort[vertexIds] || ! DuplicateFreeQ[rankOrder],
    Return[msFailure["InvalidRankOrder", <|"expected" -> vertexIds, "actual" -> rankOrder|>]]
  ];
  chartCertificate = MSBoundaryChartCertificate[context, targetRules, RankOrder -> rankOrder];
  If[Head[chartCertificate] === Failure || ! TrueQ[chartCertificate["normalCrossingQ"]],
    Return[msFailure["BoundaryChartNotCertified", <|"certificate" -> chartCertificate|>]]
  ];
  (* Every boundary runs through the generic sector-DAG route. *)
  seriesData = msGenericSectorFrobeniusData[
    context, de, targetRules, scale, rankOrder, order, workingPrecision
  ];
  If[Head[seriesData] === Failure, Return[seriesData]];
  anchorLetters = N[de["letters"] /. seriesData["anchorRules"], workingPrecision];
  If[AnyTrue[anchorLetters, ! NumericQ[#] || TrueQ[PossibleZeroQ[#]] &],
    Return[msFailure[
      "SingularBoundaryAnchor",
      <|"letters" -> anchorLetters, "rules" -> seriesData["anchorRules"]|>
    ]]
  ];
  Join[
  <|
    "status" -> "generated",
    "method" -> If[
      TrueQ[Lookup[seriesData, "sectorLeadingSystemQ", False]],
      "2411GenericSectorLeadingSeries",
      "2411VertexFrobeniusSeries"
    ],
    "masters" -> de["masters"],
    "masterDigest" -> de["masterDigest"],
    "boundaryKind" -> Lookup[seriesData, "boundaryKind", "finiteFrobeniusSeries"],
    "rankOrder" -> rankOrder,
    "boundaryChartCertificate" -> chartCertificate,
    "boundaryScale" -> scale,
    "seriesOrder" -> seriesData["seriesOrder"],
    "convergenceRatio" -> seriesData["convergenceRatio"],
    "lateTimeExponents" -> Lookup[seriesData, "lateTimeExponents", Lookup[seriesData["leadingBranches"], "frobeniusExponent", {}]],
    "leadingBranches" -> seriesData["leadingBranches"],
    "blowupVariables" -> Lookup[seriesData, "blowupVariables", Missing["SingleVariable"]],
    "thetaFixing" -> Lookup[seriesData, "thetaFixing", {}],
    "sectorLeadingSystemQ" -> Lookup[seriesData, "sectorLeadingSystemQ", False],
    "anchorRules" -> seriesData["anchorRules"],
    "targetRules" -> targetRules,
    "anchorLetters" -> anchorLetters,
    "branchConvention" -> "HankelH[1|2] and h prefactor fixed by initialized context",
    "normalizationIncludedQ" -> True,
    "ordinaryAnchorQ" -> Lookup[seriesData, "boundaryKind", "finiteFrobeniusSeries"] === "finiteFrobeniusSeries",
    "formulaAuthority" -> If[
      TrueQ[Lookup[seriesData, "sectorLeadingSystemQ", False]],
      "2411.03088 Secs. 3.3 and 4.2; sector-DAG power-log recurrence delegated to FlintNDE",
      "2411.03088 Eqs. (3.44)-(3.46), with endpoint coefficients from Eqs. (3.51)-(3.54)"
    ]
  |>,
  KeyTake[
    seriesData,
    {
      "values", "errorEstimate", "singularParameter", "singularCurveRules",
      "singularConnection", "singularResidue", "singularStart", "singularTarget"
    }
  ]
  ]
];

MSBoundaryData[___] := msFailure["InitializedContextRequired", <|"function" -> "MSBoundaryData"|>];
