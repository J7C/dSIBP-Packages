(* ::Package:: *)

(***
文件：TensorAtoms.wl
用途：从 sector slot registry 直接组装 Kronecker 原子、M1/M0、M0 局部张量对角基与能量 letters。
公式约定：包的基准公式即 +nu convention 形状，formulaNu = s_nu nu 是代入其中的带号 Hankel 阶（正 prefactor 给 +nu，负 prefactor 给 -nu）；arXiv:2401.00129 的印刷形状由 nu_ref=-nu 得到。M1 的 massive endpoint 原子取 (1-2 formulaNu)；massive
endpoint 只含 sigma2，massless shared slot 只含 sigma1；若 exact 张量变换不能对角化 M0，立即拒绝。
***)

(* ::Chapter:: *)
(* Kronecker tools *)

msKroneckerAll[{}] := {{1}};
msKroneckerAll[matrices_List] := Fold[KroneckerProduct, First[matrices], Rest[matrices]];

msEmbedMatrix[matrix_, position_Integer, slotCount_Integer] := msKroneckerAll@Table[
  If[index === position, matrix, msIdentity2],
  {index, slotCount}
];

msRawIdentityForSector[sector_Association] := IdentityMatrix[sector["rawStateCount"]];

msQuotientMatrix[sector_Association, rawMatrix_] := Simplify[
  sector["stateProjection"].rawMatrix.sector["stateEmbedding"]
];

msSlotPositionByKey[sector_Association, key_] := FirstPosition[
  Lookup[sector["slots"], "key"],
  key,
  Missing["UnknownSlot", key]
];


(* ::Section:: *)
(* Component and slot incidence *)

msComponentContainsRootQ[sector_Association, componentPosition_Integer, rootVertex_] := MemberQ[
  sector["vertexComponents"][[componentPosition]],
  rootVertex
];

msHankelSlotsAtComponent[sector_Association, componentPosition_Integer] := Select[
  MapIndexed[Append[#1, "slotPosition" -> First[#2]] &, sector["slots"]],
  #["kind"] === "massiveEndpoint" &&
    msComponentContainsRootQ[sector, componentPosition, #["rootVertex"]] &
];

msMasslessSharedIncidence[slot_Association, sector_Association, componentPosition_Integer] := Total@MapThread[
  If[MemberQ[sector["vertexComponents"][[componentPosition]], #1], #2, 0] &,
  {slot["endpoints"], {1, -1}}
];


(* ::Section:: *)
(* M1/M0 and simultaneous diagonalization *)

msM1Matrix[sector_Association, componentPosition_Integer, integerShift_: 0] := Module[
  {matrix, slots},
  matrix = (sector["baseTimePowers"][[componentPosition]] + integerShift) msRawIdentityForSector[sector];
  slots = msHankelSlotsAtComponent[sector, componentPosition];
  Do[
    matrix -= (1 - 2 slot["formulaNu"]) msEmbedMatrix[
      msProjector1,
      slot["slotPosition"],
      sector["slotCount"]
    ],
    {slot, slots}
  ];
  msQuotientMatrix[sector, matrix]
];

msM0Matrix[sector_Association, componentPosition_Integer] := Module[
  {matrix, slot, incidence},
  matrix = I sector["componentEnergies"][[componentPosition]] msRawIdentityForSector[sector];
  Do[
    slot = Append[sector["slots"][[slotPosition]], "slotPosition" -> slotPosition];
    Switch[
      slot["kind"],
      "massiveEndpoint",
        If[msComponentContainsRootQ[sector, componentPosition, slot["rootVertex"]],
          matrix += msEmbedMatrix[-I slot["momentum"] msSigma2, slotPosition, sector["slotCount"]]
        ],
      "masslessShared",
        incidence = msMasslessSharedIncidence[slot, sector, componentPosition];
        If[incidence =!= 0,
          matrix += msEmbedMatrix[
            I incidence slot["fullContourSign"] slot["momentum"] msSigma1,
            slotPosition,
            sector["slotCount"]
          ]
        ]
    ],
    {slotPosition, sector["slotCount"]}
  ];
  msQuotientMatrix[sector, matrix]
];

msLocalM0Diagonalizer[slot_Association, sector_Association] := Switch[
  slot["kind"],
  "massiveEndpoint", msPaperT,
  "masslessShared",
    If[MemberQ[sector["coincidentLineIds"], slot["lineId"]], msIdentity2, msHadamard],
  _, msFailure["UnsupportedTensorSlotKind", <|"slot" -> slot|>]
];


msLocalM0DiagonalizerInverse[slot_Association, sector_Association] := Switch[
  slot["kind"],
  "massiveEndpoint", msPaperTInverse,
  "masslessShared",
    If[MemberQ[sector["coincidentLineIds"], slot["lineId"]], msIdentity2, msHadamard],
  _, msFailure["UnsupportedTensorSlotKind", <|"slot" -> slot|>]
];


msSectorDiagonalizer[sector_Association] := Module[{localMatrices},
  localMatrices = msLocalM0Diagonalizer[#, sector] & /@ sector["slots"];
  If[! FreeQ[localMatrices, _Failure], Return[FirstCase[localMatrices, _Failure]]];
  msQuotientMatrix[sector, msKroneckerAll[localMatrices]]
];


msSectorDiagonalizerInverse[sector_Association] := Module[{localMatrices},
  localMatrices = msLocalM0DiagonalizerInverse[#, sector] & /@ sector["slots"];
  If[! FreeQ[localMatrices, _Failure], Return[FirstCase[localMatrices, _Failure]]];
  msQuotientMatrix[sector, msKroneckerAll[localMatrices]]
];


msExactZeroMatrixQ[matrix_] := MatrixQ[matrix] && TrueQ[
  And @@ Flatten[Map[PossibleZeroQ, Simplify[matrix], {2}]]
];


(* 只在 exact 检查通过后抽取对角元。这样未来同一 slot 若混入不能共同对角化的局部组分，
   不会被 Diagonal 静默丢弃。 *)
msSectorM0DiagonalizationData[sector_Association, componentPosition_Integer] := Module[
  {u, uInverse, roundTripResidual, transformed, offDiagonalResidual},
  u = msSectorDiagonalizer[sector];
  If[Head[u] === Failure, Return[u]];
  uInverse = msSectorDiagonalizerInverse[sector];
  If[Head[uInverse] === Failure, Return[uInverse]];
  roundTripResidual = Simplify[u.uInverse - IdentityMatrix[sector["masterCount"]]];
  transformed = Simplify[u.msM0Matrix[sector, componentPosition].uInverse];
  offDiagonalResidual = Simplify[transformed - DiagonalMatrix[Diagonal[transformed]]];
  If[! msExactZeroMatrixQ[roundTripResidual] || ! msExactZeroMatrixQ[offDiagonalResidual],
    Return[msFailure[
      "TensorM0DiagonalizationNotCertified",
      <|
        "sectorKey" -> sector["sectorKey"],
        "componentPosition" -> componentPosition,
        "roundTripResidual" -> roundTripResidual,
        "offDiagonalResidual" -> offDiagonalResidual
      |>
    ]]
  ];
  <|
    "U" -> u,
    "UInverse" -> uInverse,
    "diagonalMatrix" -> transformed,
    "diagonal" -> Diagonal[transformed]
  |>
];

(* M1 is strictly diagonal in the state-bit basis; here we only take reciprocals of the scalar eigenvalues. *)
msM1Inverse[sector_Association, componentPosition_Integer, integerShift_: 0] := Module[
  {diagonal = Diagonal[msM1Matrix[sector, componentPosition, integerShift]]},
  DiagonalMatrix[Simplify[1/diagonal]]
];

msM0Inverse[sector_Association, componentPosition_Integer] := Module[
  {data = msSectorM0DiagonalizationData[sector, componentPosition]},
  If[Head[data] === Failure,
    data,
    Simplify[data["UInverse"].DiagonalMatrix[1/data["diagonal"]].data["U"]]
  ]
];

msEnergyLetters[sector_Association, componentPosition_Integer] := Module[
  {data = msSectorM0DiagonalizationData[sector, componentPosition]},
  If[Head[data] === Failure, data, Simplify[data["diagonal"]/I]]
];

msM1SingularSurfaces[sector_Association, componentPosition_Integer, integerShift_: 0] := DeleteDuplicates[
  Simplify[Diagonal[msM1Matrix[sector, componentPosition, integerShift]]]
];


(* ::Section:: *)
(* Public matrix queries *)

msFormulaMatricesForSector[sector_Association] := Module[
  {componentCount, u, uInverse, components, diagonalizationData},
  componentCount = Length[sector["vertexComponents"]];
  u = msSectorDiagonalizer[sector];
  If[Head[u] === Failure, Return[u]];
  uInverse = msSectorDiagonalizerInverse[sector];
  If[Head[uInverse] === Failure, Return[uInverse]];
  components = Table[
    diagonalizationData = msSectorM0DiagonalizationData[sector, componentPosition];
    If[Head[diagonalizationData] === Failure, Return[diagonalizationData]];
    <|
      "componentPosition" -> componentPosition,
      "rootVertices" -> sector["vertexComponents"][[componentPosition]],
      "baseTimePower" -> sector["baseTimePowers"][[componentPosition]],
      "externalLegEnergy" -> sector["componentEnergies"][[componentPosition]],
      "M1" -> msM1Matrix[sector, componentPosition, 0],
      "M1PlusOne" -> msM1Matrix[sector, componentPosition, 1],
      "M0" -> msM0Matrix[sector, componentPosition],
      "M0Diagonal" -> diagonalizationData["diagonalMatrix"],
      "energyLetters" -> Simplify[diagonalizationData["diagonal"]/I],
      "M1SingularSurfaces" -> msM1SingularSurfaces[sector, componentPosition, 0]
    |>,
    {componentPosition, componentCount}
  ];
  <|
    "status" -> "generated",
    "sectorKey" -> sector["sectorKey"],
    "slotRegistry" -> sector["slots"],
      "stateOrder" -> sector["stateOrder"],
      "rawStateOrder" -> sector["rawStateOrder"],
      "stateEmbedding" -> sector["stateEmbedding"],
      "stateProjection" -> sector["stateProjection"],
    "dimension" -> sector["masterCount"],
    "U" -> u,
    "UInverse" -> uInverse,
    "components" -> components,
    "letters" -> DeleteDuplicates[Flatten[Lookup[components, "energyLetters"]]]
  |>
];

MSFormulaMatrices[context_?MSContextQ, All] := AssociationThread[
  context["sectorOrder"] -> (msFormulaMatricesForSector /@ context["sectors"])
];

MSFormulaMatrices[context_?MSContextQ, key_String] := Module[{sector = msSectorByKey[context, key]},
  If[Head[sector] === Missing,
    Message[MSFormulaMatrices::nosector, key];
    msFailure["UnknownSector", <|"sectorKey" -> key|>],
    msFormulaMatricesForSector[sector]
  ]
];

MSFormulaMatrices[context_?MSContextQ] := MSFormulaMatrices[context, All];
