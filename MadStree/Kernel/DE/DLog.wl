(* ::Package:: *)

(***
File: DLog.wl
Purpose: Directly assembles the ordered master integrals and the block-triangular dlog potential according to 2401.00129 Eqs. (3.55), (3.65) and (3.68).
Gate: dlog is certified only when R^(1) lies exactly in the child zero shift, normalization is absorbed, and the residue does not depend on kinematics.
***)

(* ::Chapter:: *)
(* Sector diagonal blocks *)

msOmegaEx[sector_Association] := DiagonalMatrix@Map[
  Function[bits,
    -Total@MapIndexed[
       If[#1["kind"] === "massiveEndpoint",
        bits[[First[#2]]] (1 - 2 #1["formulaNu"]) Log[#1["momentum"]],
        0
      ] &,
      sector["slots"]
    ]
  ],
  sector["stateOrder"]
];

msComponentLogKernel[sector_Association, componentPosition_Integer] := Module[
  {u, uInverse, letters, omegaTilde},
  u = msSectorDiagonalizer[sector];
  If[Head[u] === Failure, Return[u]];
  uInverse = msSectorDiagonalizerInverse[sector];
  If[Head[uInverse] === Failure, Return[uInverse]];
  letters = msEnergyLetters[sector, componentPosition];
  If[Head[letters] === Failure, Return[letters]];
  omegaTilde = DiagonalMatrix[-I Log /@ letters];
  Simplify[-I uInverse.omegaTilde.u]
];

msSectorDLogBlock[sector_Association] := Module[
  {componentCount, componentKernels, omega, normalizationGauge, letters},
  componentCount = Length[sector["vertexComponents"]];
  componentKernels = Table[
    msComponentLogKernel[sector, componentPosition],
    {componentPosition, componentCount}
  ];
  If[! FreeQ[componentKernels, _Failure], Return[FirstCase[componentKernels, _Failure]]];
  omega = msOmegaEx[sector] + Total@Table[
    componentKernels[[componentPosition]].msM1Matrix[sector, componentPosition, 1],
    {componentPosition, componentCount}
  ];
  normalizationGauge = If[
    TrueQ[sector["normalization"] === 1],
    ConstantArray[0, Dimensions[omega]],
    Log[sector["normalization"]] IdentityMatrix[sector["masterCount"]]
  ];
  omega = Simplify[omega + normalizationGauge];
  letters = DeleteDuplicates@Cases[omega, Log[letter_] :> letter, Infinity];
  <|
    "sectorKey" -> sector["sectorKey"],
    "dimension" -> sector["masterCount"],
    "omega" -> omega,
    "omegaEx" -> msOmegaEx[sector],
    "normalizationGauge" -> normalizationGauge,
    "componentLogKernels" -> componentKernels,
    "letters" -> letters,
    "masters" -> Select[MSMasterIntegrals[$msActiveDLogContext], #["sectorKey"] === sector["sectorKey"] &]
  |>
];


(* ::Section:: *)
(*Eq. (3.68) parent-to-subsector block*)

msTargetMasterEmbedding[targetSector_Association, context_?MSContextQ] := Module[
  {globalBasis, targetBasis, globalIndex, matrix},
  globalBasis = Lookup[context["masters"], "integral"];
  targetBasis = MSIntegral[targetSector["sectorKey"], ConstantArray[0, Length[targetSector["vertexComponents"]]], #] & /@
    targetSector["stateOrder"];
  globalIndex = AssociationThread[ToString[#, InputForm] & /@ globalBasis, Range[Length[globalBasis]]];
  matrix = ConstantArray[0, {Length[targetBasis], Length[globalBasis]}];
  Do[
    matrix[[row, globalIndex[ToString[targetBasis[[row]], InputForm]]]] = 1,
    {row, Length[targetBasis]}
  ];
  matrix
];


msShiftedTargetReduction[
  targetSector_Association,
  targetShift_List,
  context_?MSContextQ
] := Module[{integrals, reductions, closedQ},
  integrals = MSIntegral[targetSector["sectorKey"], targetShift, #] & /@ targetSector["stateOrder"];
  reductions = MSReduce[#, context] & /@ integrals;
  closedQ = And @@ Map[
    Lookup[#, "status", None] === "reduced" &&
      TrueQ[Lookup[#, "nonMasterResidual", Missing["Absent"]] === 0] &&
      Lookup[#, "remainingShiftedIntegrals", {Missing["Absent"]}] === {} &,
    reductions
  ];
  <|
    "status" -> If[closedQ, "reducedToGlobalMasters", "shiftReductionFailed"],
    "targetShift" -> targetShift,
    "integrals" -> integrals,
    "matrix" -> If[closedQ, Lookup[reductions, "coefficientVector"], Missing["UnclosedReduction"]],
    "residuals" -> Lookup[reductions, "nonMasterResidual", Missing["Absent"]],
    "remainingShiftedIntegrals" -> DeleteDuplicates@Flatten[
      Lookup[reductions, "remainingShiftedIntegrals", {}]
    ],
    "singularLayers" -> DeleteDuplicates@Flatten[Lookup[reductions, "singularLayers", {}], 1]
  |>
];

msContactDLogBlock[
  sourceSector_Association,
  contactMap_Association,
  sourceBlock_Association,
  context_?MSContextQ
] := Module[
  {targetSector, endpointComponents, matrices, r1Shifts, contributions, reductionRecords,
    targetShift, reduction, contactSign, localBlock,
   globalBlock, closedQ, sectorSlices,
   sectorStart, sectorDimension},
  targetSector = msSectorByKey[context, contactMap["targetSector"]];
  endpointComponents = DeleteDuplicates[contactMap["endpointComponents"]];
  matrices = contactMap["matricesByComponent"];
  r1Shifts = Values[contactMap["R1TargetShiftsByComponent"]];
  reductionRecords = Table[
    targetShift = contactMap["R1TargetShiftsByComponent"][componentPosition];
    reduction = If[
      And @@ (# === 0 & /@ targetShift),
      <|
        "status" -> "alreadyTargetMasters",
        "targetShift" -> targetShift,
        "matrix" -> msTargetMasterEmbedding[targetSector, context],
        "residuals" -> ConstantArray[0, targetSector["masterCount"]],
        "remainingShiftedIntegrals" -> {},
        "singularLayers" -> {}
      |>,
      msShiftedTargetReduction[targetSector, targetShift, context]
    ];
    contactSign = msContactEventSign[contactMap];
    (* -I T^-1 OmegaTilde0 T 已在 component kernel 中；这里只补 massless quotient 相对 h/Wronskian 原子的 (-1)^N0。 *)
    localBlock = Simplify[
      contactSign
        sourceBlock["componentLogKernels"][[componentPosition]].matrices[componentPosition]
    ];
    <|
      "componentPosition" -> componentPosition,
      "localContactMatrix" -> localBlock,
      "shiftReduction" -> reduction,
      "globalContribution" -> If[
        MemberQ[{"alreadyTargetMasters", "reducedToGlobalMasters"}, reduction["status"]],
        Simplify[localBlock.reduction["matrix"]],
        Missing["UnclosedShiftReduction"]
      ]
    |>,
    {componentPosition, endpointComponents}
  ];
  closedQ = FreeQ[Lookup[reductionRecords, "globalContribution"], _Missing];
  globalBlock = If[
    closedQ,
    Simplify[Total[Lookup[reductionRecords, "globalContribution"]]],
    Missing["UnclosedShiftReduction"]
  ];
  sectorSlices = If[
    closedQ,
    Association@Table[
      sectorStart = First@Lookup[
        Select[context["masters"], #["sectorKey"] === sector["sectorKey"] &],
        "globalIndex"
      ];
      sectorDimension = sector["masterCount"];
      sector["sectorKey"] -> globalBlock[[All, sectorStart ;; sectorStart + sectorDimension - 1]],
      {sector, context["sectors"]}
    ],
    <||>
  ];
  <|
    "sourceSector" -> sourceSector["sectorKey"],
    "targetSector" -> targetSector["sectorKey"],
    "eventId" -> contactMap["eventId"],
    "lineId" -> contactMap["lineId"],
    "selectedLineIds" -> contactMap["selectedLineIds"],
    "matrix" -> If[closedQ, sectorSlices[targetSector["sectorKey"]], Missing["UnclosedShiftReduction"]],
    "globalMatrix" -> globalBlock,
    "matricesByTargetSector" -> sectorSlices,
    "dimensions" -> If[closedQ, Dimensions[globalBlock], Missing["UnclosedShiftReduction"]],
    "R1TargetShifts" -> contactMap["R1TargetShiftsByComponent"],
    "R1HitsTargetMasterQ" -> And @@ (And @@ (# === 0 & /@ #) & /@ r1Shifts),
    "R1ReductionClosedQ" -> closedQ,
    "shiftReductionRecords" -> reductionRecords,
    "absorbedNormalization" -> contactMap["absorbedNormalization"]
  |>
];


(* ::Chapter:: *)
(* Full connection and certification *)

msKinematicSymbols[context_?MSContextQ] := DeleteDuplicates@Cases[
  Join[
    Lookup[context["vertices"], "externalLegEnergy"],
    If[context["lines"] === {}, {}, Lookup[context["lines"], "momentum"]]
  ],
  symbol_Symbol /; Context[symbol] =!= "System`",
  Infinity
];

msConstantResidueQ[matrix_, kinematicSymbols_List] := If[
  kinematicSymbols === {},
  True,
  FreeQ[matrix, Alternatives @@ kinematicSymbols]
];

MSDLogDE[context_?MSContextQ] := Module[
  {sectors, sectorOrder, dimensions, sectorPositions, blocks, contactMaps, contactBlocks,
   omegaBlocks, sourceIndex, targetIndex, omega, letters, letterMatrices, residual,
   kinematicSymbols, constantResidueQ, shiftQ, normalizationQ, dlogStatus, offsets},
  sectors = context["sectors"];
  sectorOrder = context["sectorOrder"];
  dimensions = Lookup[sectors, "masterCount"];
  sectorPositions = AssociationThread[sectorOrder -> Range[Length[sectorOrder]]];
  Block[{$msActiveDLogContext = context}, blocks = msSectorDLogBlock /@ sectors];
  If[AnyTrue[blocks, Head[#] === Failure &],
    Return[FirstCase[blocks, failure_Failure :> failure]]
  ];
  contactMaps = Flatten[msContactMapsForSector[#, context] & /@ sectors];
  contactBlocks = Map[
    Function[contactMap,
      sourceIndex = sectorPositions[contactMap["sourceSector"]];
      msContactDLogBlock[sectors[[sourceIndex]], contactMap, blocks[[sourceIndex]], context]
    ],
    contactMaps
  ];
  If[AnyTrue[contactBlocks, Head[#] === Failure &],
    Return[FirstCase[contactBlocks, failure_Failure :> failure]]
  ];
  omegaBlocks = Table[
    If[sourceIndex === targetIndex,
      blocks[[sourceIndex, "omega"]],
      ConstantArray[0, {dimensions[[sourceIndex]], dimensions[[targetIndex]]}]
    ],
    {sourceIndex, Length[sectors]},
    {targetIndex, Length[sectors]}
  ];
  Do[
    sourceIndex = sectorPositions[contactBlock["sourceSector"]];
    If[TrueQ[contactBlock["R1ReductionClosedQ"]],
      Do[
        targetIndex = sectorPositions[targetSectorKey];
        omegaBlocks[[sourceIndex, targetIndex]] = Simplify[
          omegaBlocks[[sourceIndex, targetIndex]] +
            contactBlock["matricesByTargetSector"][targetSectorKey]
        ],
        {targetSectorKey, sectorOrder}
      ]
    ],
    {contactBlock, contactBlocks}
  ];
  omega = ArrayFlatten[omegaBlocks];
  letters = DeleteDuplicates@Cases[omega, Log[letter_] :> letter, Infinity];
  letterMatrices = Association@Table[
    letter -> Map[Coefficient[#, Log[letter]] &, omega, {2}],
    {letter, letters}
  ];
  residual = Simplify[omega - Total[(Log[#] letterMatrices[#]) & /@ letters]];
  kinematicSymbols = msKinematicSymbols[context];
  constantResidueQ = And @@ (msConstantResidueQ[#, kinematicSymbols] & /@ Values[letterMatrices]);
  shiftQ = If[contactBlocks === {}, True, And @@ Lookup[contactBlocks, "R1ReductionClosedQ", False]];
  normalizationQ = If[
    contactMaps === {},
    True,
    And @@ (TrueQ[Simplify[# - 1] === 0] & /@ Lookup[contactMaps, "absorbedNormalization", 1])
  ];
  dlogStatus = Which[
    ! shiftQ, "contactShiftReductionFailed",
    ! normalizationQ, "requiresNormalizationGauge",
    ! TrueQ[residual === ConstantArray[0, Dimensions[omega]]], "nonDLogResidual",
    ! constantResidueQ, "requiresGaugeTransformation",
    True, "certifiedByFormulaChecks"
  ];
  offsets = AssociationThread[
    sectorOrder -> Most@FoldList[Plus, 1, dimensions]
  ];
  <|
    "status" -> "generated",
    "dlogStatus" -> dlogStatus,
    "convention" -> context["convention"],
    "sectorOrder" -> sectorOrder,
    "sectorDimensions" -> dimensions,
    "masterSectorOffsets" -> offsets,
    "masters" -> context["masters"],
    "bareMasters" -> Lookup[context["masters"], "integral"],
    "masterCount" -> Length[context["masters"]],
    "masterDigest" -> context["masterDigest"],
    "sectorBlocks" -> blocks,
    "contactBlocks" -> contactBlocks,
    "omegaBlocks" -> omegaBlocks,
    "omegaPotential" -> omega,
    "matrixDimension" -> Dimensions[omega],
    "letters" -> letters,
    "letterMatrices" -> letterMatrices,
    "dlogResidual" -> residual,
    "constantResidueQ" -> constantResidueQ,
    "contactR1MasterShiftQ" -> And @@ Lookup[contactBlocks, "R1HitsTargetMasterQ", True],
    "contactR1ReductionClosedQ" -> shiftQ,
    "normalizationAbsorbedQ" -> normalizationQ,
    "kinematicSymbols" -> kinematicSymbols,
    "formulaAuthority" -> "2401.00129 Eq. (3.55), (3.65), (3.68)"
  |>
];
