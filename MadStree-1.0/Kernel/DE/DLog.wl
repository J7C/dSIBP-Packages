(* ::Package:: *)

(***
File: DLog.wl
Purpose: Directly assembles the ordered master integrals and the block-triangular dlog potential according to 2401.00129 Eqs. (3.55), (3.65) and (3.68).
Representation: The connection is first kept as structured 2x2-tensor kernels. Dense matrices are
materialized only by the public MSDLogDE query; the FlintNDE path consumes the structured system.
Gate: dlog is certified only when R^(1) lies exactly in the subsector zero shift, normalization is absorbed, and the residue does not depend on kinematics.
***)

(* ::Chapter:: *)
(* Structured sector diagonal blocks *)

(* Omega_ex 在 state-bit basis 对角；这里只构造对角向量，不做任何 2^n 直乘。 *)
msTensorOmegaExDiagonal[sector_Association] := Map[
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


(* 一个 sector 的结构化 dlog 数据：
   - omegaExDiagonal：state 对角的 Omega_ex 向量；
   - normalizationLog：Log[N_s]（恒等归一化时为 0）；
   - components：每个 component 的 U-对角化字母 lambda 与 M1(1) 对角向量。
   U^{-1} diag(Log lambda) U 永远按局部张量作用，不展开 U。 *)
msSectorDLogStructured[sector_Association] := Module[
  {componentCount, omegaExDiagonal, normalizationLog, components, kappa, letters,
   indicators, m1Vector},
  componentCount = Length[sector["vertexComponents"]];
  omegaExDiagonal = msTensorOmegaExDiagonal[sector];
  normalizationLog = If[
    TrueQ[sector["normalization"] === 1],
    0,
    Log[sector["normalization"]]
  ];
  components = Table[
    kappa = msEnergyLetters[sector, componentPosition];
    If[
      Head[kappa] === Failure,
      kappa,
      letters = DeleteDuplicates[Simplify /@ kappa];
      indicators = Table[
        Map[If[TrueQ[Simplify[# - letter] === 0], 1, 0] &, kappa],
        {letter, letters}
      ];
      m1Vector = msM1DiagonalVector[sector, componentPosition, 1];
      If[
        Head[m1Vector] === Failure,
        m1Vector,
        <|
          "componentPosition" -> componentPosition,
          "letters" -> letters,
          "indicators" -> indicators,
          "m1Vector" -> m1Vector,
          "kappa" -> kappa
        |>
      ]
    ],
    {componentPosition, componentCount}
  ];
  If[! FreeQ[components, _Failure],
    Return[FirstCase[components, _Failure, Missing["Absent"], Infinity]]];
  <|
    "sectorKey" -> sector["sectorKey"],
    "dimension" -> sector["masterCount"],
    "omegaExDiagonal" -> omegaExDiagonal,
    "normalization" -> sector["normalization"],
    "normalizationLog" -> normalizationLog,
    "components" -> components
  |>
];


(* ::Chapter:: *)
(* Eq. (3.68) supersector-to-subsector block *)

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


(* 只做约化与字段收集，不 materialize localBlock；结构化 dlog 系统与公开密集输出共用这一份。 *)
msContactDLogData[
  sourceSector_Association,
  contactMap_Association,
  context_?MSContextQ
] := Module[
  {targetSector, endpointComponents, matrices, r1Shifts, reductionRecords, targetShift,
   reduction, contactSign, statuses},
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
    <|
      "componentPosition" -> componentPosition,
      "matrices" -> matrices[componentPosition],
      "shiftReduction" -> reduction
    |>,
    {componentPosition, endpointComponents}
  ];
  contactSign = msContactEventSign[contactMap];
  statuses = Lookup[
    Lookup[reductionRecords, "shiftReduction"],
    "status",
    Missing["Absent"]
  ];
  <|
    "sourceSector" -> sourceSector["sectorKey"],
    "targetSector" -> targetSector["sectorKey"],
    "eventId" -> contactMap["eventId"],
    "lineId" -> contactMap["lineId"],
    "selectedLineIds" -> contactMap["selectedLineIds"],
    "endpointComponents" -> endpointComponents,
    "contactSign" -> contactSign,
    "reductionRecords" -> reductionRecords,
    "closedQ" -> And @@ (MemberQ[
      {"alreadyTargetMasters", "reducedToGlobalMasters"}, #] & /@ statuses),
    "R1TargetShifts" -> contactMap["R1TargetShiftsByComponent"],
    "R1HitsTargetMasterQ" -> And @@ (And @@ (# === 0 & /@ #) & /@ r1Shifts),
    "absorbedNormalization" -> contactMap["absorbedNormalization"]
  |>
];


(* ::Chapter:: *)
(* Structured connection system *)

(* sector-local 贡献只占 global 向量的一个 sector 行块；contact 贡献通过 middleMatrix 从
   global 列映回 source sector。两类都先算成 source 维向量，再嵌回 global 行块。 *)
msDLogEmbedSectorVector[system_Association, sectorKey_, localVector_List] := Module[
  {offset = system["masterSectorOffsets"][sectorKey], dimension, local},
  dimension = Length[localVector];
  Join[
    ConstantArray[0, offset - 1],
    localVector,
    ConstantArray[0, system["masterCount"] - offset - dimension + 1]
  ]
];


msDLogContributionApply[system_Association, contribution_Association, vector_List] := Module[
  {sector, sectorKey, offset, dimension, middle, value, indicator},
  Switch[contribution["type"],
    "diagonal",
      If[contribution["scope"] === "sector",
        sectorKey = contribution["sectorKey"];
        offset = system["masterSectorOffsets"][sectorKey];
        dimension = Length[contribution["vector"]];
        msDLogEmbedSectorVector[
          system, sectorKey,
          contribution["vector"] vector[[offset ;; offset + dimension - 1]]
        ],
        contribution["vector"] vector
      ],
    "conjugated",
      sectorKey = contribution["sectorKey"];
      sector = system["sectorByKey"][sectorKey];
      offset = system["masterSectorOffsets"][sectorKey];
      dimension = sector["masterCount"];
      middle = contribution["middleMatrix"];
      value = If[
        contribution["scope"] === "sector",
        vector[[offset ;; offset + dimension - 1]],
        middle . vector
      ];
      If[Head[value] === Failure, Return[value]];
      If[contribution["preVector"] =!= None, value = contribution["preVector"] value];
      value = msTensorApplyU[sector, value];
      If[Head[value] === Failure, Return[value]];
      indicator = contribution["letterIndicator"];
      value = indicator value;
      value = msTensorApplyUInverse[sector, value];
      If[Head[value] === Failure, Return[value]];
      msDLogEmbedSectorVector[system, sectorKey, contribution["postFactor"] value],
    _,
      msFailure[
        "UnsupportedDLogContribution",
        <|"type" -> contribution["type"]|>
      ]
  ]
];


msDLogContributionMatrix[system_Association, contribution_Association] := Module[
  {dimension = system["masterCount"], columns},
  Switch[contribution["type"],
    "diagonal",
      If[contribution["scope"] === "sector",
        DiagonalMatrix[msDLogEmbedSectorVector[
          system, contribution["sectorKey"], contribution["vector"]
        ]],
        DiagonalMatrix[contribution["vector"]]
      ],
    "conjugated",
      columns = Table[
        msDLogContributionApply[system, contribution, UnitVector[dimension, column]],
        {column, dimension}
      ];
      If[! FreeQ[columns, _Failure], Return[FirstCase[columns, _Failure]]];
      Transpose[columns],
    _,
      msFailure[
        "UnsupportedDLogContribution",
        <|"type" -> contribution["type"]|>
      ]
  ]
];


msDLogLetterMatrix[system_Association, letter_] := Module[
  {contributions = Lookup[system["residues"], letter, {}], matrices},
  matrices = msDLogContributionMatrix[system, #] & /@ contributions;
  If[! FreeQ[matrices, _Failure], Return[FirstCase[matrices, _Failure]]];
  Total[matrices]
];


msDLogContributionConstantQ[contribution_Association, kinematicSymbols_List] := Switch[
  contribution["type"],
  "diagonal",
    msConstantResidueQ[contribution["vector"], kinematicSymbols],
  "conjugated",
    msConstantResidueQ[contribution["preVector"], kinematicSymbols] &&
      msConstantResidueQ[contribution["middleMatrix"], kinematicSymbols] &&
      msConstantResidueQ[contribution["letterIndicator"], kinematicSymbols] &&
      msConstantResidueQ[contribution["postFactor"], kinematicSymbols],
  _,
    False
];


msDLogAddContribution[state_Association, letter_, contribution_] := Module[
  {existing, key,
   keyQ = MemberQ[state["letterOrder"], letter]},
  existing = Lookup[state["residues"], letter, Missing["Absent"]];
  <|
    "letterOrder" -> If[TrueQ[keyQ], state["letterOrder"], Append[state["letterOrder"], letter]],
    "residues" -> If[
      Head[existing] === Missing,
      Join[state["residues"], <|letter -> {contribution}|>],
      Join[state["residues"], <|letter -> Append[existing, contribution]|>]
    ]
  |>
];


msDLogTensorSystem[context_?MSContextQ] := Module[
  {sectors, sectorOrder, dimensions, sectorPositions, structuredSectors, sectorByKey,
   sectorStructuredByKey, contactMaps, contactData, state, kinematicSymbols,
   shiftQ, normalizationQ, structuralConstantQ, constantResidueQ, dlogStatus,
   contactRecords, data, rec, structured, component, contribution, contributions,
   slot, slotIndex, letters, offsets, letterMatrices, certificateSystem},
  sectors = context["sectors"];
  sectorOrder = context["sectorOrder"];
  dimensions = Lookup[sectors, "masterCount"];
  sectorPositions = AssociationThread[sectorOrder -> Range[Length[sectorOrder]]];
  structuredSectors = msSectorDLogStructured /@ sectors;
  If[! FreeQ[structuredSectors, _Failure],
    Return[FirstCase[structuredSectors, _Failure, Missing["Absent"], Infinity]]];
  sectorByKey = AssociationThread[sectorOrder -> sectors];
  sectorStructuredByKey = AssociationThread[sectorOrder -> structuredSectors];
  contactMaps = Flatten[msContactMapsForSector[#, context] & /@ sectors];
  contactData = msContactDLogData[
    sectors[[sectorPositions[#["sourceSector"]]]], #, context
  ] & /@ contactMaps;
  If[! FreeQ[contactData, _Failure],
    Return[FirstCase[contactData, _Failure, Missing["Absent"], Infinity]]];

  state = <|"residues" -> <||>, "letterOrder" -> {}|>;

  (* Omega_ex：massive slot 的 Log[q] 对角项；每个 slot 只贡献自身 bit，
     同一 momentum 的多个 slot 通过同一 letter 相加，恰好还原对角向量。 *)
  Do[
    structured = sectorStructuredByKey[sector["sectorKey"]];
    Do[
      slot = sector["slots"][[slotIndex]];
      If[slot["kind"] === "massiveEndpoint",
        With[{letter = slot["momentum"]},
          contribution = <|
            "type" -> "diagonal",
            "scope" -> "sector",
            "sectorKey" -> sector["sectorKey"],
            "vector" -> Map[
              Function[bits, -bits[[slotIndex]] (1 - 2 slot["formulaNu"])],
              sector["stateOrder"]
            ]
          |>;
          If[! TrueQ[msTensorExactZeroVectorQ[contribution["vector"]]],
            state = msDLogAddContribution[state, letter, contribution]
          ]
        ]
      ],
      {slotIndex, Length[sector["slots"]]}
    ];
    (* 归一化规范项 Log[N_s] Identity；对 letter N_s 的留数就是单位系数。 *)
    If[TrueQ[structured["normalizationLog"] =!= 0],
      state = msDLogAddContribution[
        state,
        structured["normalization"],
        <|
          "type" -> "diagonal",
          "scope" -> "sector",
          "sectorKey" -> sector["sectorKey"],
          "vector" -> ConstantArray[1, structured["dimension"]]
        |>
      ]
    ];
    Do[
      contributions = Table[
        <|
          "type" -> "conjugated",
          "scope" -> "sector",
          "sectorKey" -> sector["sectorKey"],
          "preVector" -> component["m1Vector"],
          "middleMatrix" -> None,
          "letterIndicator" -> component["indicators"][[index]],
          "postFactor" -> -1
        |>,
        {index, Length[component["letters"]]}
      ];
      Do[
        state = msDLogAddContribution[
          state, component["letters"][[index]], contributions[[index]]
        ],
        {index, Length[component["letters"]]}
      ],
      {component, structured["components"]}
    ],
    {sector, sectors}
  ];

  (* contact block：contactSign * K_c * C_c * R；U 的局部张量作用保持不变。 *)
  Do[
    data = contactRecord;
    If[TrueQ[data["closedQ"]],
      structured = sectorStructuredByKey[data["sourceSector"]];
      Do[
        rec = reductionRecord;
        If[MemberQ[{"alreadyTargetMasters", "reducedToGlobalMasters"},
            rec["shiftReduction"]["status"]],
          component = structured["components"][[rec["componentPosition"]]];
          With[{middle = Simplify[rec["matrices"] . rec["shiftReduction"]["matrix"]]},
            Do[
              state = msDLogAddContribution[
                state,
                component["letters"][[index]],
                <|
                  "type" -> "conjugated",
                  "scope" -> "global",
                  "sectorKey" -> data["sourceSector"],
                  "preVector" -> None,
                  "middleMatrix" -> middle,
                  "letterIndicator" -> component["indicators"][[index]],
                  "postFactor" -> -data["contactSign"]
                |>
              ],
              {index, Length[component["letters"]]}
            ]
          ]
        ],
        {reductionRecord, data["reductionRecords"]}
      ]
    ],
    {contactRecord, contactData}
  ];

  letters = state["letterOrder"];
  kinematicSymbols = msKinematicSymbols[context];
  shiftQ = If[contactData === {}, True, And @@ Lookup[contactData, "closedQ", False]];
  normalizationQ = If[
    contactMaps === {},
    True,
    And @@ (TrueQ[Simplify[# - 1] === 0] & /@ Lookup[contactMaps, "absorbedNormalization", 1])
  ];
  structuralConstantQ = And @@ Flatten@Table[
    msDLogContributionConstantQ[#, kinematicSymbols] & /@
      Lookup[state["residues"], letter, {}],
    {letter, letters}
  ];
  If[structuralConstantQ,
    constantResidueQ = True,
    (* 结构化系数含运动学符号时，只在证书需要时 materialize residue 判断是否真正抵消。 *)
    certificateSystem = <|
      "masterCount" -> Length[context["masters"]],
      "masterSectorOffsets" -> AssociationThread[sectorOrder -> Most@FoldList[Plus, 1, dimensions]],
      "residues" -> state["residues"],
      "sectorByKey" -> sectorByKey
    |>;
    letterMatrices = Association@Table[
      letter -> msDLogLetterMatrix[certificateSystem, letter],
      {letter, letters}
    ];
    If[! FreeQ[letterMatrices, _Failure], Return[FirstCase[letterMatrices, _Failure]]];
    constantResidueQ = And @@ (msConstantResidueQ[#, kinematicSymbols] & /@ Values[letterMatrices])
  ];
  dlogStatus = Which[
    ! shiftQ, "contactShiftReductionFailed",
    ! normalizationQ, "requiresNormalizationGauge",
    ! constantResidueQ, "requiresGaugeTransformation",
    True, "certifiedByFormulaChecks"
  ];
  offsets = AssociationThread[sectorOrder -> Most@FoldList[Plus, 1, dimensions]];
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
    "sectorStructured" -> sectorStructuredByKey,
    "sectorByKey" -> sectorByKey,
    "contactData" -> contactData,
    "contactR1ReductionClosedQ" -> shiftQ,
    "normalizationAbsorbedQ" -> normalizationQ,
    "constantResidueQ" -> constantResidueQ,
    "kinematicSymbols" -> kinematicSymbols,
    "letters" -> letters,
    "residues" -> state["residues"],
    "formulaAuthority" -> "2401.00129 Eq. (3.55), (3.65), (3.68)"
  |>
];


(* ::Chapter:: *)
(* Dense materialization at the public boundary *)

msComponentLogKernelMaterialized[sector_Association, component_Association] := Module[
  {dimension = sector["masterCount"], logKappa},
  logKappa = Map[Log, component["kappa"]];
  Transpose@Table[
    Simplify[-msTensorApplyUInverse[
      sector,
      logKappa msTensorApplyU[sector, UnitVector[dimension, column]]
    ]],
    {column, dimension}
  ]
];


msOmegaExDense[sector_Association, structured_Association] :=
  DiagonalMatrix[structured["omegaExDiagonal"]];


msSectorDLogBlock[
  sector_Association,
  structured_Association,
  contextMasters_List
] := Module[
  {componentCount, componentKernels, omega, normalizationGauge, letters},
  componentCount = Length[sector["vertexComponents"]];
  componentKernels = Map[
    msComponentLogKernelMaterialized[sector, #] &,
    structured["components"]
  ];
  If[! FreeQ[componentKernels, _Failure], Return[FirstCase[componentKernels, _Failure]]];
  omega = msOmegaExDense[sector, structured] + Total@Table[
    componentKernels[[componentPosition]] .
      DiagonalMatrix[structured["components"][[componentPosition, "m1Vector"]]],
    {componentPosition, componentCount}
  ];
  normalizationGauge = If[
    TrueQ[structured["normalizationLog"] === 0],
    ConstantArray[0, {sector["masterCount"], sector["masterCount"]}],
    structured["normalizationLog"] IdentityMatrix[sector["masterCount"]]
  ];
  omega = Simplify[omega + normalizationGauge];
  letters = DeleteDuplicates@Cases[omega, Log[letter_] :> letter, Infinity];
  <|
    "sectorKey" -> sector["sectorKey"],
    "dimension" -> sector["masterCount"],
    "omega" -> omega,
    "omegaEx" -> msOmegaExDense[sector, structured],
    "normalizationGauge" -> normalizationGauge,
    "componentLogKernels" -> componentKernels,
    "letters" -> letters,
    "masters" -> Select[contextMasters, #["sectorKey"] === sector["sectorKey"] &]
  |>
];


msContactDLogBlock[
  sourceSector_Association,
  data_Association,
  sourceBlock_Association,
  context_?MSContextQ
] := Module[
  {targetSector, reductionRecords, localBlocks, globalContributions, closedQ, sectorSlices,
   sectorStart, sectorDimension, records},
  targetSector = msSectorByKey[context, data["targetSector"]];
  reductionRecords = data["reductionRecords"];
  localBlocks = Table[
    Simplify[
      data["contactSign"] *
        sourceBlock["componentLogKernels"][[record["componentPosition"]]] .
        record["matrices"]
    ],
    {record, reductionRecords}
  ];
  globalContributions = Table[
    If[MemberQ[{"alreadyTargetMasters", "reducedToGlobalMasters"},
        reductionRecords[[index, "shiftReduction", "status"]]],
      Simplify[localBlocks[[index]] . reductionRecords[[index, "shiftReduction", "matrix"]]],
      Missing["UnclosedShiftReduction"]
    ],
    {index, Length[reductionRecords]}
  ];
  closedQ = FreeQ[globalContributions, _Missing];
  records = Table[
    <|
      "componentPosition" -> reductionRecords[[index, "componentPosition"]],
      "localContactMatrix" -> localBlocks[[index]],
      "shiftReduction" -> reductionRecords[[index, "shiftReduction"]],
      "globalContribution" -> globalContributions[[index]]
    |>,
    {index, Length[reductionRecords]}
  ];
  sectorSlices = If[
    closedQ,
    Association@Table[
      sectorStart = First@Lookup[
        Select[context["masters"], #["sectorKey"] === sector["sectorKey"] &],
        "globalIndex"
      ];
      sectorDimension = sector["masterCount"];
      sector["sectorKey"] -> Total[globalContributions][[All, sectorStart ;; sectorStart + sectorDimension - 1]],
      {sector, context["sectors"]}
    ],
    <||>
  ];
  <|
    "sourceSector" -> sourceSector["sectorKey"],
    "targetSector" -> targetSector["sectorKey"],
    "eventId" -> data["eventId"],
    "lineId" -> data["lineId"],
    "selectedLineIds" -> data["selectedLineIds"],
    "matrix" -> If[closedQ, sectorSlices[targetSector["sectorKey"]], Missing["UnclosedShiftReduction"]],
    "globalMatrix" -> If[closedQ, Total[globalContributions], Missing["UnclosedShiftReduction"]],
    "matricesByTargetSector" -> sectorSlices,
    "dimensions" -> If[closedQ, Dimensions[Total[globalContributions]], Missing["UnclosedShiftReduction"]],
    "R1TargetShifts" -> data["R1TargetShifts"],
    "R1HitsTargetMasterQ" -> data["R1HitsTargetMasterQ"],
    "R1ReductionClosedQ" -> closedQ,
    "shiftReductionRecords" -> records,
    "absorbedNormalization" -> data["absorbedNormalization"]
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
  {system, sectors, sectorOrder, dimensions, sectorPositions, blocks, contactBlocks,
   sourceIndex, targetIndex, omegaBlocks, omega, letters, letterMatrices, residual,
   kinematicSymbols, constantResidueQ, shiftQ, normalizationQ, dlogStatus, offsets,
   contextMasters, tensorSectors},
  system = msDLogTensorSystem[context];
  If[Head[system] === Failure, Return[system]];
  contextMasters = context["masters"];
  sectors = context["sectors"];
  sectorOrder = context["sectorOrder"];
  dimensions = system["sectorDimensions"];
  sectorPositions = AssociationThread[sectorOrder -> Range[Length[sectorOrder]]];
  blocks = Map[
    msSectorDLogBlock[#, system["sectorStructured"][#["sectorKey"]], contextMasters] &,
    sectors
  ];
  If[AnyTrue[blocks, Head[#] === Failure &],
    Return[FirstCase[blocks, failure_Failure :> failure]]
  ];
  contactBlocks = Map[
    Function[data,
      sourceIndex = sectorPositions[data["sourceSector"]];
      msContactDLogBlock[sectors[[sourceIndex]], data, blocks[[sourceIndex]], context]
    ],
    system["contactData"]
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
    system["contactData"] === {},
    True,
    And @@ (TrueQ[Simplify[# - 1] === 0] & /@ Lookup[system["contactData"], "absorbedNormalization", 1])
  ];
  dlogStatus = Which[
    ! shiftQ, "contactShiftReductionFailed",
    ! normalizationQ, "requiresNormalizationGauge",
    ! TrueQ[residual === ConstantArray[0, Dimensions[omega]]], "nonDLogResidual",
    ! constantResidueQ, "requiresGaugeTransformation",
    True, "certifiedByFormulaChecks"
  ];
  offsets = system["masterSectorOffsets"];
  tensorSectors = Association@Table[
    sectorKey -> msTensorSectorWireRecord[system["sectorByKey"][sectorKey], offsets[sectorKey]],
    {sectorKey, sectorOrder}
  ];
  If[! FreeQ[tensorSectors, _Failure],
    Return[FirstCase[tensorSectors, _Failure, Missing["Absent"], Infinity]]
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
    "tensorSectors" -> tensorSectors,
    "tensorResidues" -> system["residues"],
    "dlogResidual" -> residual,
    "constantResidueQ" -> constantResidueQ,
    "contactR1MasterShiftQ" -> And @@ Lookup[contactBlocks, "R1HitsTargetMasterQ", True],
    "contactR1ReductionClosedQ" -> shiftQ,
    "normalizationAbsorbedQ" -> normalizationQ,
    "kinematicSymbols" -> kinematicSymbols,
    "formulaAuthority" -> "2401.00129 Eq. (3.55), (3.65), (3.68)"
  |>
];
