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
    Select[sector["activeLines"], msTimeOrderedLineQ]
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
(* 2411.03088 Sec. 3.3 single-vertex Frobenius coefficients *)

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

(* The Hankel branches of a (anti-)time-ordered propagator are determined only by the SK type and the strict time rank. *)
msBoundaryTimeOrderedBranches["++", firstDominantQ_] := If[firstDominantQ, {1, 2}, {2, 1}];
msBoundaryTimeOrderedBranches["--", firstDominantQ_] := If[firstDominantQ, {2, 1}, {1, 2}];


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
    "massiveW",
      line["hankelBranches"][[slot["endpointIndex"]]],
    "massiveT" | "masslessT",
      endpointComponents = sector["rootToComponent"] /@ line["endpoints"];
      firstDominantQ = If[
        SameQ @@ endpointComponents,
        True,
        msBoundaryComponentRank[sector, First[endpointComponents], rank] <
          msBoundaryComponentRank[sector, Last[endpointComponents], rank]
      ];
      branches = msBoundaryTimeOrderedBranches[line["contourType"], firstDominantQ];
      branches[[slot["endpointIndex"]]],
    _,
      Missing["NoHankelBranch", slot["key"]]
  ]
];


(* A (anti-)time-ordered line whose two endpoints already merged into one component carries no time ordering:
   the two Hankel arrangements are the two sides of the same equal-time fibre, so a boundary
   coefficient that keeps only one of them is not constant along the fibre. *)
msBoundaryCoincidentTimeOrderedLineQ[sector_Association, line_Association] :=
  line["type"] === "massiveT" &&
    SameQ @@ (sector["rootToComponent"] /@ line["endpoints"]);


(* Fibre-symmetric endpoint coefficient of one coincident (anti-)time-ordered line, m = (alpha + beta)/2.
   Swapping the two bits gives the same pair of products, so the slot order is irrelevant.
   Averaging rescales coefficients only; it never moves the Frobenius power. *)
msBoundaryCoincidentEndpointCoefficient[
  convention_String,
  {firstBit_Integer, secondBit_Integer},
  nu_
] := FullSimplify[
  (msVertexEndpointCoefficient[convention, 1, firstBit, nu]
     msVertexEndpointCoefficient[convention, 2, secondBit, nu] +
   msVertexEndpointCoefficient[convention, 2, firstBit, nu]
     msVertexEndpointCoefficient[convention, 1, secondBit, nu])/2
];


(* Endpoint factor of one component, built line by line. A strictly ordered line keeps the
   branches fixed by the component rank, while every coincident (anti-)time-ordered line contributes exactly
   one averaged factor; with several spectator lines the result is therefore the product of
   their m factors, which is what the odd-subset theta bundle algebra requires. *)
msBoundaryComponentEndpointFactor[
  context_?MSContextQ,
  sector_Association,
  rank_Association,
  hSlots_List,
  stateBits_List,
  targetRules_List
] := Module[{convention},
  convention = context["convention", "nuConvention"];
  Times @@ Map[
    Function[lineSlots,
      Module[{line, bits, momentum, physicalNu, formulaNu, powerFactor, branchFactor},
        line = msLineById[context, First[lineSlots]["lineId"]];
        momentum = line["momentum"] /. targetRules;
        physicalNu = line["nu"] /. targetRules;
        formulaNu = line["formulaNu"] /. targetRules;
        bits = stateBits[[#["slotPosition"]]] & /@ lineSlots;
        powerFactor = Times @@ Map[
          Function[bit, (-I momentum)^(-bit (1 - 2 formulaNu))],
          bits
        ];
        branchFactor = If[
          msBoundaryCoincidentTimeOrderedLineQ[sector, line],
          msBoundaryCoincidentEndpointCoefficient[convention, bits, physicalNu],
          Times @@ Map[
            Function[slot,
              msVertexEndpointCoefficient[
                convention,
                msBoundaryHankelBranch[slot, sector, rank, context],
                stateBits[[slot["slotPosition"]]],
                physicalNu
              ]
            ],
            lineSlots
          ]
        ];
        powerFactor branchFactor
      ]
    ],
    GatherBy[hSlots, #["lineId"] &]
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
   retainedRules, anchorRules, perLetterRecords, zeroMatrix, connection, residue, convergenceRatio},
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
  (* 留数矩阵不依赖外能量；曲线代入、求导和取留数只在各 letter 标量上进行。
     letterAlong 逐 letter 保留：序列化层据此下发 letter list，后端按 dlog 分解
     装配认证极点，代替对 Total 后单矩阵做根隔离。 *)
  perLetterRecords = Table[
    Module[{letterAlong, dlogAlong, scalarResidue, matrix},
      letterAlong = letter /. curveRules /. retainedRules;
      dlogAlong = Cancel[Together[D[Log[letterAlong], parameter]]];
      scalarResidue = Together[Limit[parameter dlogAlong, parameter -> 0]];
      matrix = de["letterMatrices"][letter] /. retainedRules;
      <|
        "connection" -> dlogAlong matrix,
        "residue" -> scalarResidue matrix,
        "letter" -> letter,
        "letterAlong" -> letterAlong
      |>
    ],
    {letter, de["letters"]}
  ];
  zeroMatrix = ConstantArray[0, {de["masterCount"], de["masterCount"]}];
  connection = If[perLetterRecords === {}, zeroMatrix, Total[Lookup[perLetterRecords, "connection"]]];
  If[! FreeQ[connection, _Real],
    Return[msFailure[
      "ExactSingularConnectionRequired",
      <|"reason" -> "generic sector pullback must remain exact Q(i)(t)"|>
    ]]
  ];
  residue = If[perLetterRecords === {}, zeroMatrix, Total[Lookup[perLetterRecords, "residue"]]];
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
    "letterPullbacks" -> Table[
      <|"letter" -> record["letter"], "along" -> record["letterAlong"]|>,
      {record, perLetterRecords}
    ],
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


(* 固定 root 分量并解完整约束系统；共振或直接回代失败时使用这条精确路径。 *)
msBoundaryConstrainedIndicialVector[
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


msBoundaryDirectIndicialVector[
  context_?MSContextQ,
  residue_?MatrixQ,
  exponent_,
  sectorKey_String,
  rootIndex_Integer
] := Module[
  {dimension = Length[residue], indicesBySector, vector, previous = {rootIndex},
   rows, rhs, denominators, residual, ancestorKeys, allowedIndices},
  (* 非共振时按块三角结构逆序回代；零分母或约束/残差失败都回退完整精确求解。 *)
  indicesBySector = Association@Table[
    key -> Lookup[Select[context["masters"], #["sectorKey"] === key &], "globalIndex"],
    {key, context["sectorOrder"]}
  ];
  vector = ConstantArray[0, dimension];
  vector[[rootIndex]] = 1;
  Catch[
    Do[
      If[key =!= sectorKey,
        rows = indicesBySector[key];
        rhs = residue[[rows, previous]] . vector[[previous]];
        denominators = exponent - Diagonal[residue[[rows, rows]]];
        If[AnyTrue[denominators, TrueQ[PossibleZeroQ[#]] &], Throw[$Failed]];
        vector[[rows]] = rhs / denominators;
        previous = Join[previous, rows]
      ],
      {key, Reverse[context["sectorOrder"]]}
    ];
    ancestorKeys = msBoundaryAncestorSectorKeys[context, sectorKey];
    allowedIndices = Lookup[
      Select[context["masters"], MemberQ[ancestorKeys, #["sectorKey"]] &],
      "globalIndex"
    ];
    If[AnyTrue[Complement[Range[dimension], allowedIndices],
      ! TrueQ[PossibleZeroQ[vector[[#]]]] &], Throw[$Failed]];
    residual = Together[exponent vector - residue . vector];
    If[! And @@ (TrueQ[PossibleZeroQ[#]] & /@ residual), Throw[$Failed]];
    Together[vector]
  ]
];


msBoundaryRootedIndicialVector[
  context_?MSContextQ,
  residue_?MatrixQ,
  exponent_,
  sectorKey_String,
  rootIndex_Integer
] := Module[{direct},
  direct = msBoundaryDirectIndicialVector[context, residue, exponent, sectorKey, rootIndex];
  If[direct =!= $Failed,
    direct,
    msBoundaryConstrainedIndicialVector[context, residue, exponent, sectorKey, rootIndex]
  ]
];


(* ::Section::Closed:: *)
(*已知 roots 的 power-log 证书*)

(* 结构支撑检查：首项向量自身坐标为 1、同 sector 兄弟坐标为 0、sector 祖先以外坐标为 0。
   与构造器内的精确残差检查合用，即得 81 个向量线性独立的三角结构论证。 *)
msBoundaryVectorSupportQ[context_, vector_, sectorKey_, rootIndex_Integer] := Module[
  {ancestorKeys, allowedIndices, rootIndices, forbiddenIndices},
  ancestorKeys = msBoundaryAncestorSectorKeys[context, sectorKey];
  allowedIndices = Lookup[
    Select[context["masters"], MemberQ[ancestorKeys, #["sectorKey"]] &],
    "globalIndex"
  ];
  rootIndices = Lookup[Select[context["masters"], #["sectorKey"] === sectorKey &], "globalIndex"];
  forbiddenIndices = Join[
    Complement[Range[Length[vector]], allowedIndices],
    Complement[rootIndices, {rootIndex}]
  ];
  TrueQ[PossibleZeroQ[vector[[rootIndex]] - 1]] &&
    And @@ (TrueQ[PossibleZeroQ[vector[[#]]]] & /@ forbiddenIndices)
];


(* 已知指数的整数根差：返回 <|{alpha,beta} -> d|>，d 为正整数。 *)
msBoundaryPositiveIntegerDifferences[roots_List] := Module[{pairs},
  pairs = Select[
    Flatten[Table[{roots[[i]], roots[[j]]}, {i, Length[roots]}, {j, Length[roots]}], 1],
    (#[[1]] =!= #[[2]] &&
       IntegerQ[#[[1]] - #[[2]]] && #[[1]] - #[[2]] > 0) &
  ];
  Association[Thread[Rule[pairs, #[[1]] - #[[2]] & /@ pairs]]]
];


(* 缺员指数的 Jordan 链分析：alpha 已知，只解零空间与广义原像，不求特征多项式。
   自顶向下补基：从最高层核空间取模低层独立的链顶，逐层下移，得到全部链。 *)
msBoundaryJordanChains[residue_, alpha_] := Module[
  {dimension, identity, nilpotent, kernels, level, chains, images, basis, lower, tops,
   candidate, chain, vector, extended},
  dimension = Length[residue];
  identity = IdentityMatrix[dimension];
  nilpotent = alpha identity - residue;
  kernels = {NullSpace[nilpotent]};
  level = 1;
  While[level <= dimension,
    With[{next = NullSpace[MatrixPower[nilpotent, level + 1]]},
      If[Length[next] === Length[Last[kernels]], Break[]];
      AppendTo[kernels, next]
    ];
    level++
  ];
  chains = {};
  images = {};
  Do[
    basis = kernels[[order]];
    lower = If[order > 1, kernels[[order - 1]], {}];
    tops = {};
    Do[
      candidate = basis[[position]];
      extended = Join[images, tops];
      If[
        MatrixRank[Join[extended, {candidate}]] > MatrixRank[extended] &&
          MatrixRank[Join[extended, lower, {candidate}]] > MatrixRank[Join[extended, lower]],
        AppendTo[tops, candidate]
      ],
      {position, Length[basis]}
    ];
    Do[
      vector = tops[[position]];
      chain = Table[MatrixPower[nilpotent, drop] . vector, {drop, 0, order - 1}];
      AppendTo[chains, <|"tail" -> First[chain], "head" -> Last[chain], "length" -> order|>],
      {position, Length[tops]}
    ];
    images = Map[nilpotent . # &, tops],
    {order, Length[kernels], 1, -1}
  ];
  <|
    "geometric" -> Length[First[kernels]],
    "chains" -> chains,
    "chainLengths" -> Lookup[chains, "length"],
    "kernels" -> kernels
  |>
];


(* 已知 roots 的 power-log 证书：指数与代数重数来自 master 公式；几何重数由结构独立的首项
   向量数给出（免零空间）；缺员指数只用已知 alpha 解零空间与 Jordan 链；正整数根差直接列出。
   不解 indicial 特征多项式。 *)
msBoundaryKnownRootsCertificate[context_, residue_, allExponents_List, successRecords_List, failureRecords_List] := Module[
  {roots, algebraic, successCounts, geometric, supportViolation, differences,
   resonanceCandidates, defects, chainData, exponent, multiplicity, chains, route,
   chainTotal},
  roots = Union[allExponents];
  algebraic = Counts[allExponents];
  successCounts = Counts[Lookup[successRecords, "frobeniusExponent"]];
  supportViolation = Select[
    successRecords,
    ! msBoundaryVectorSupportQ[
      context, #["normalizedLeadingVector"], #["sectorKey"], #["rootGlobalIndex"]
    ] &
  ];
  If[Length[supportViolation] > 0,
    Return[msFailure[
      "CertificateSupportViolation",
      <|"records" -> Lookup[supportViolation, "rootGlobalIndex"]|>
    ]]
  ];
  geometric = Association[Thread[Rule[roots, Lookup[successCounts, roots, 0]]]];
  chainData = Association[];
  defects = {};
  Do[
    exponent = roots[[position]];
    multiplicity = algebraic[exponent];
    If[geometric[exponent] < multiplicity,
      chains = msBoundaryJordanChains[residue, exponent];
      chainTotal = Total[chains["chainLengths"]];
      If[chainTotal =!= multiplicity,
        Return[msFailure[
          "CertificateChainMismatch",
          <|"root" -> exponent, "algebraic" -> multiplicity, "chainTotal" -> chainTotal|>
        ]]
      ];
      AppendTo[defects, <|"root" -> exponent, "algebraic" -> multiplicity, "geometric" -> chains["geometric"]|>];
      chainData[exponent] = chains
    ],
    {position, Length[roots]}
  ];
  differences = msBoundaryPositiveIntegerDifferences[roots];
  resonanceCandidates = Map[
    <|"root" -> #[[2]], "degree" -> #[[1]] - #[[2]], "absolutePower" -> #[[1]]|> &,
    Keys[differences]
  ];
  route = Which[
    Length[defects] > 0, "defectiveKnownRoots",
    Length[differences] > 0, "resonantKnownRoots",
    True, "purePowerKnownRoots"
  ];
  <|
    "route" -> route,
    "roots" -> roots,
    "algebraicMultiplicities" -> Lookup[algebraic, roots],
    "geometricMultiplicities" -> Lookup[geometric, roots],
    "jordanDefects" -> defects,
    "positiveIntegerDifferences" -> differences,
    "resonanceCandidates" -> resonanceCandidates,
    "chainData" -> chainData,
    "requiredLogs" -> Join[
      Flatten[Map[
        Function[pair,
          Map[
            <|
              "root" -> pair[[1]],
              "kind" -> "jordan",
              "leadingLogDegree" -> #["length"] - 1,
              "head" -> #["head"],
              "tail" -> #["tail"]
            |> &,
            Select[pair[[2]]["chains"], #["length"] > 1 &]
          ]
        ],
        Normal[chainData]
      ]],
      Map[
        <|"root" -> #["root"], "kind" -> "resonance", "degree" -> #["degree"]|> &,
        resonanceCandidates
      ]
    ]
  |>
];


(* 缺员 master 的广义首项：在最小核空间 K_h 中解带 rooting 约束的线性系统，
   log 幂次为 h-1，{a,b,C} 中的 C 是 log^{h-1} 的系数 head/(h-1)!。 *)
msBoundaryGeneralizedLeadingVector[context_, residue_, exponent_, sectorKey_, rootIndex_, chains_] := Module[
  {dimension, identity, nilpotent, ancestorKeys, allowedIndices, rootIndices, forbiddenIndices,
   constraintRows, constraintValues, kernels, basis, coefficients, vector, head},
  dimension = Length[residue];
  identity = IdentityMatrix[dimension];
  nilpotent = exponent identity - residue;
  ancestorKeys = msBoundaryAncestorSectorKeys[context, sectorKey];
  allowedIndices = Lookup[
    Select[context["masters"], MemberQ[ancestorKeys, #["sectorKey"]] &],
    "globalIndex"
  ];
  rootIndices = Lookup[Select[context["masters"], #["sectorKey"] === sectorKey &], "globalIndex"];
  forbiddenIndices = Join[
    Complement[Range[dimension], allowedIndices],
    Complement[rootIndices, {rootIndex}]
  ];
  constraintRows = UnitVector[dimension, #] & /@ Join[{rootIndex}, forbiddenIndices];
  constraintValues = Join[{1}, ConstantArray[0, Length[forbiddenIndices]]];
  kernels = chains["kernels"];
  Catch[
    Do[
      basis = kernels[[height]];
      coefficients = Quiet@Check[
        LinearSolve[
          constraintRows . Transpose[basis],
          constraintValues
        ],
        $Failed
      ];
      If[coefficients =!= $Failed && VectorQ[coefficients],
        vector = basis . coefficients;
        head = MatrixPower[nilpotent, height - 1] . vector;
        If[AnyTrue[head, ! TrueQ[PossibleZeroQ[#]] &],
          Throw[<|"height" -> height, "vector" -> vector, "head" -> head|>]
        ]
      ],
      {height, Length[kernels]}
    ];
    $Failed
  ]
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
     sharedFactor, coefficient, exponent, componentFailure},
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
          slotBit (1 - 2 (slot["formulaNu"] /. targetRules))
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
        endpointFactor = msBoundaryComponentEndpointFactor[
          context, sector, rank, hSlots, stateBits, targetRules
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
  {curveData, sectorByKey, localRecords, failure, branchAttempts, certificate, branches,
   generalized, thetaFixing},
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
  (* 全部 master 都构造首项向量；求解失败者交由证书做 Jordan 链分析，不在此终止。 *)
  branchAttempts = Map[
    Function[record,
      With[
        {vector = msBoundaryRootedIndicialVector[
          context,
          curveData["residue"],
          record["frobeniusExponent"],
          record["sectorKey"],
          record["rootGlobalIndex"]
        ]},
        If[Head[vector] === Failure,
          Join[record, <|"normalizedLeadingVector" -> Missing["IndicialVectorFailed"]|>],
          Join[record, <|"normalizedLeadingVector" -> vector|>]
        ]
      ]
    ],
    localRecords
  ];
  certificate = msBoundaryKnownRootsCertificate[
    context,
    curveData["residue"],
    Lookup[localRecords, "frobeniusExponent"],
    Select[branchAttempts, ! MissingQ[#["normalizedLeadingVector"]] &],
    Select[branchAttempts, MissingQ[#["normalizedLeadingVector"]] &]
  ];
  If[Head[certificate] === Failure, Return[certificate]];
  If[certificate["route"] === "defectiveKnownRoots",
    branches = Map[
      Function[record,
        If[MissingQ[record["normalizedLeadingVector"]],
          generalized = msBoundaryGeneralizedLeadingVector[
            context,
            curveData["residue"],
            record["frobeniusExponent"],
            record["sectorKey"],
            record["rootGlobalIndex"],
            certificate["chainData"][record["frobeniusExponent"]]
          ];
          If[generalized === $Failed,
            msFailure[
              "SectorGeneralizedLeadingFailed",
              <|
                "sectorKey" -> record["sectorKey"],
                "rootIndex" -> record["rootGlobalIndex"],
                "frobeniusExponent" -> record["frobeniusExponent"]
              |>
            ],
            Join[
              record,
              <|
                "logPower" -> generalized["height"] - 1,
                "normalizedLeadingVector" -> generalized["head"]/(generalized["height"] - 1)!,
                "generalizedChainTail" -> generalized["vector"]
              |>
            ]
          ],
          record
        ]
      ],
      branchAttempts
    ],
    branches = branchAttempts
  ];
  failure = FirstCase[branches, _Failure, Missing["NoFailure"]];
  If[Head[failure] === Failure, Return[failure]];
  branches = Select[branches, ! TrueQ[PossibleZeroQ[#["physicalWeight"]]] &];
  thetaFixing = AssociationThread[
    context["sectorOrder"] ->
      (msSectorThetaFixing[#, curveData["rank"]] & /@ context["sectors"])
  ];
  <|
    "boundaryKind" -> "singularFrobenius",
    "anchorRules" -> curveData["anchorRules"],
    "convergenceRatio" -> curveData["convergenceRatio"],
    "leadingBranches" -> branches,
    "powerLogCertificate" -> KeyTake[certificate, {
      "route", "roots", "algebraicMultiplicities", "geometricMultiplicities",
      "jordanDefects", "positiveIntegerDifferences", "resonanceCandidates", "requiredLogs"
    }],
    "seriesOrder" -> order,
    "blowupVariables" -> msBoundaryEnergyChartData[context, rankOrder],
    "singularParameter" -> curveData["parameter"],
    "singularCurveRules" -> curveData["curveRules"],
    "singularLetterPullbacks" -> curveData["letterPullbacks"],
    "singularRetainedRules" -> curveData["retainedRules"],
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
  If[AssociationQ[Lookup[seriesData, "powerLogCertificate", Missing["Absent"]]],
    KeyTake[seriesData, {"powerLogCertificate"}],
    <||>
  ],
  KeyTake[
    seriesData,
    {
      "values", "errorEstimate", "singularParameter", "singularCurveRules",
      "singularLetterPullbacks", "singularRetainedRules",
      "singularConnection", "singularResidue", "singularStart", "singularTarget"
    }
  ]
  ]
];

MSBoundaryData[___] := msFailure["InitializedContextRequired", <|"function" -> "MSBoundaryData"|>];
