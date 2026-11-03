(* ::Package:: *)

(***
文件：TensorAtoms.wl
用途：只保留公开公式查询所需的密集 materialization；M1/M0/U 的日常运算由 TensorProduct.wl 的逐 slot 张量表示承担。
约定：包的基准公式即 +nu convention 形状，formulaNu = s_nu nu；massive endpoint 只含 sigma2，massless shared slot 只含 sigma1。
***)

(* ::Chapter:: *)
(* Kronecker tools used only by the public dense query *)

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
(* Dense formulas materialized at the public boundary *)

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
            I incidence slot["timeOrderedContourSign"] slot["momentum"] msSigma1,
            slotPosition,
            sector["slotCount"]
          ]
        ]
    ],
    {slotPosition, sector["slotCount"]}
  ];
  msQuotientMatrix[sector, matrix]
];


(* U/U^{-1} 由已认证的局部对角化子 materialize；U M0 U^{-1}=diag(I kappa)，与 v0.18 的
   M0Diagonal 口径一致，energyLetters 再除以 I 还原成能级字母。 *)
msSectorM0DiagonalizationData[sector_Association, componentPosition_Integer] := Module[
  {kappa, u, uInverse},
  kappa = I msEnergyLetters[sector, componentPosition];
  If[Head[kappa] === Failure, Return[kappa]];
  u = msSectorDiagonalizer[sector];
  If[Head[u] === Failure, Return[u]];
  uInverse = msSectorDiagonalizerInverse[sector];
  If[Head[uInverse] === Failure, Return[uInverse]];
  <|
    "U" -> u,
    "UInverse" -> uInverse,
    "diagonalMatrix" -> DiagonalMatrix[kappa],
    "diagonal" -> kappa
  |>
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
