(* ::Package:: *)

(***
文件：TensorProduct.wl
用途：把 sector 的 M0/M1、M0 对角化基以及 U/U^{-1} 的作用全部保留为逐 slot 的 2x2 张量积，按局部矩阵乘法完成运算。
范围：本模块不展开 2^n 直乘；只有 MSFormulaMatrices 与 MSDLogDE 的公开密集输出在边界处 materialize。
约定：slot 1 是最高位，raw 状态编号与 msStateList 的 IntegerDigits 顺序一致；U=P U_raw S，U^{-1}=P U_raw^{-1} S。
***)

(* ::Chapter:: *)
(* Local 2x2 algebra *)

msTensorKroneckerAll[{}] := {{1}};
msTensorKroneckerAll[matrices_List] := Fold[KroneckerProduct, First[matrices], Rest[matrices]];


(* 沿一个 slot 轴作用局部 2x2 矩阵；其余轴保持原序。
   Transpose[t, p] 的第 k 个下标对应原数组的 p^{-1}(k)，所以「把 axis 移到最前」
   需要 p = Ordering[moveToFront]，而还原用 p^{-1} = moveToFront。 *)
msTensorApplyLocalToAxis[tensor_, matrix_, axis_Integer] := Module[
  {rank = TensorRank[tensor], moveToFront, permutation, inversePermutation, moved, applied},
  If[rank === 1, Return[matrix . tensor]];
  moveToFront = Prepend[DeleteCases[Range[rank], axis], axis];
  permutation = Ordering[moveToFront];
  inversePermutation = moveToFront;
  moved = Transpose[tensor, permutation];
  applied = matrix . ArrayReshape[moved, {2, 2^(rank - 1)}];
  Transpose[ArrayReshape[applied, ConstantArray[2, rank]], inversePermutation]
];


(* n 个 2x2 局部矩阵的直乘作用在长度 2^n 的向量上，复杂度 O(n 2^n)。 *)
msTensorApplyLocalChain[locals_List, vector_List] := Module[
  {rank = Length[locals], tensor},
  If[rank === 0, Return[vector]];
  If[Length[vector] =!= 2^rank, Return[$Failed]];
  tensor = ArrayReshape[vector, ConstantArray[2, rank]];
  Do[tensor = msTensorApplyLocalToAxis[tensor, locals[[axis]], axis], {axis, rank}];
  Flatten[tensor]
];


(* ::Chapter:: *)
(* Slot-local M0/M1 data *)

msTensorSlotM0Local[
  slot_Association,
  sector_Association,
  componentPosition_Integer
] := Switch[
  slot["kind"],
  "massiveEndpoint",
    If[MemberQ[sector["vertexComponents"][[componentPosition]], slot["rootVertex"]],
      -I slot["momentum"] msSigma2,
      None
    ],
  "masslessShared",
    With[{incidence = Total@MapThread[
        If[MemberQ[sector["vertexComponents"][[componentPosition]], #1], #2, 0] &,
        {slot["endpoints"], {1, -1}}
      ]},
      If[incidence =!= 0,
        I incidence slot["timeOrderedContourSign"] slot["momentum"] msSigma1,
        None
      ]
    ],
  _, $Failed
];


msTensorSlotM1Coefficient[slot_Association, sector_Association, componentPosition_Integer] :=
  Switch[
    slot["kind"],
    "massiveEndpoint",
      If[MemberQ[sector["vertexComponents"][[componentPosition]], slot["rootVertex"]],
        1 - 2 slot["formulaNu"],
        0
      ],
    _, 0
  ];


(* 局部对角化基的稳定标签。NDE 线格式只传标签：U 与 U^{-1} 共用的
   1/Sqrt[2] 因子在 U^{-1} D U 中整体抵消，后端因此可以只用整数高斯矩阵
   逐 slot 作用，完全不需要 2^n 维矩阵。 *)
msTensorLocalBaseTag[slot_Association, sector_Association] := Switch[
  slot["kind"],
  "massiveEndpoint", "massiveT",
  "masslessShared",
    If[MemberQ[sector["coincidentLineIds"], slot["lineId"]], "identity", "masslessH"],
  _, Missing["UnsupportedTensorSlotKind", slot]
];


msTensorLocalBaseMatrix[tag_String] := Switch[
  tag,
  "massiveT", msPaperT,
  "masslessH", msHadamard,
  "identity", msIdentity2,
  _, $Failed
];


msTensorLocalBaseMatrixInverse[tag_String] := Switch[
  tag,
  "massiveT", msPaperTInverse,
  "masslessH", msHadamard,
  "identity", msIdentity2,
  _, $Failed
];


msLocalM0Diagonalizer[slot_Association, sector_Association] := Module[
  {tag = msTensorLocalBaseTag[slot, sector]},
  If[
    Head[tag] === Missing,
    msFailure["UnsupportedTensorSlotKind", <|"slot" -> slot|>],
    msTensorLocalBaseMatrix[tag]
  ]
];


msLocalM0DiagonalizerInverse[slot_Association, sector_Association] := Module[
  {tag = msTensorLocalBaseTag[slot, sector]},
  If[
    Head[tag] === Missing,
    msFailure["UnsupportedTensorSlotKind", <|"slot" -> slot|>],
    msTensorLocalBaseMatrixInverse[tag]
  ]
];


(* ::Chapter:: *)
(* Exact checks *)

msExactZeroMatrixQ[matrix_] := MatrixQ[matrix] && TrueQ[
  And @@ Flatten[Map[PossibleZeroQ, Simplify[matrix], {2}]]
];

msTensorExactZeroVectorQ[vector_List] := TrueQ[
  And @@ (PossibleZeroQ[Simplify[#]] & /@ vector)
];


(* ::Chapter:: *)
(* Sector tensor data *)

msTensorRawPositionOrZero[vector_List] := With[
  {position = FirstPosition[vector, 1, Missing["None"]]},
  If[Head[position] === Missing, 0, First[position]]
];


(* U=P U_raw S 的四类稀疏作用。S 在 raw 行上按 canonical 列取值；P 取 representative raw 行；
   P^T 只点亮 representative；S^T 对同一 canonical 的全部 raw 原像求和。 *)
msTensorEmbedVector[data_Association, vector_List] := Map[
  If[# === 0, 0, vector[[#]]] &,
  data["canonicalOfRaw"]
];

msTensorProjectVector[data_Association, vector_List] := vector[[data["representativeRaw"]]];

msTensorProjectTransposeVector[data_Association, vector_List] := Module[
  {raw = ConstantArray[0, data["rawStateCount"]]},
  Do[raw[[data["representativeRaw"][[column]]]] = vector[[column]], {column, Length[vector]}];
  raw
];

msTensorEmbedTransposeVector[data_Association, vector_List] := Map[
  Total[vector[[#]]] &,
  data["rawPreimages"]
];


msTensorApplyUFromData[data_Association, vector_List] := msTensorProjectVector[
  data,
  msTensorApplyLocalChain[data["locals"], msTensorEmbedVector[data, vector]]
];

msTensorApplyUInverseFromData[data_Association, vector_List] := msTensorProjectVector[
  data,
  msTensorApplyLocalChain[data["localsInverse"], msTensorEmbedVector[data, vector]]
];

msTensorApplyUTransposeFromData[data_Association, vector_List] := msTensorEmbedTransposeVector[
  data,
  msTensorApplyLocalChain[Transpose /@ data["locals"],
    msTensorProjectTransposeVector[data, vector]]
];

msTensorApplyUInverseTransposeFromData[data_Association, vector_List] := msTensorEmbedTransposeVector[
  data,
  msTensorApplyLocalChain[Transpose /@ data["localsInverse"],
    msTensorProjectTransposeVector[data, vector]]
];


(* M0_raw 在 raw 张量上逐 slot 相加；最后按 P(...S) 映回 master 基。 *)
msTensorRawM0Apply[data_Association, componentPosition_Integer, rawVector_List] := Module[
  {tensor, accumulator, local},
  (* 无根传播子的 sector 只有 1 个 raw 状态：ArrayReshape 到 0 阶张量后 Flatten 会把整个向量
     塌成标量，下游 msTensorProjectVector 因此拿不到 List。这里直接退化为纯对角作用。 *)
  If[data["slotCount"] === 0,
    Return[I data["componentEnergies"][[componentPosition]] rawVector]
  ];
  tensor = ArrayReshape[rawVector, ConstantArray[2, data["slotCount"]]];
  accumulator = I data["componentEnergies"][[componentPosition]] tensor;
  Do[
    local = data["slotM0Locals"][[componentPosition]][[slotPosition]];
    If[local =!= None,
      accumulator = accumulator + msTensorApplyLocalToAxis[tensor, local, slotPosition]
    ],
    {slotPosition, data["slotCount"]}
  ];
  Flatten[accumulator]
];

msTensorM0ApplyFromData[data_Association, componentPosition_Integer, vector_List] :=
  msTensorProjectVector[
    data,
    msTensorRawM0Apply[data, componentPosition,
      msTensorEmbedVector[data, vector]]
  ];


(* M1 在 state-bit basis 严格对角：master 对角元只需读取 representative raw 状态。 *)
msTensorRawM1Diagonal[data_Association, componentPosition_Integer, integerShift_: 0] := Module[
  {vector = {data["baseTimePowers"][[componentPosition]] + integerShift}, coefficient},
  Do[
    coefficient = data["slotM1Coefficients"][[componentPosition]][[slotPosition]];
    vector = Join[vector, vector - coefficient],
    {slotPosition, data["slotCount"], 1, -1}
  ];
  vector
];

msTensorM1Diagonal[data_Association, componentPosition_Integer, integerShift_: 0] :=
  msTensorRawM1Diagonal[data, componentPosition, integerShift][[
    data["representativeRaw"]
  ]];


(* v0.18 约定是 U M0_c U^{-1} = DiagonalMatrix[kappa]；局部对角化因此按 t L t^{-1} 作用。 *)
msTensorRawComponentKappa[data_Association, componentPosition_Integer] := Module[
  {vector = {I data["componentEnergies"][[componentPosition]]}, pair},
  Do[
    pair = data["slotKappaPairs"][[componentPosition]][[slotPosition]];
    If[Head[pair] === Failure, Return[pair]];
    vector = Join[vector + pair[[1]], vector + pair[[2]]],
    {slotPosition, data["slotCount"], 1, -1}
  ];
  vector
];

msTensorMasterComponentKappa[data_Association, componentPosition_Integer] := Module[
  {raw = msTensorRawComponentKappa[data, componentPosition]},
  If[Head[raw] === Failure, raw, raw[[data["representativeRaw"]]]]
];

msTensorMasterComponentLetters[data_Association, componentPosition_Integer] := Module[
  {kappa = msTensorMasterComponentKappa[data, componentPosition]},
  If[Head[kappa] === Failure, kappa, Simplify[kappa/I]]
];


(* ::Chapter:: *)
(* Certificates *)

msTensorLocalCertificateQ[data_Association] := Module[
  {inverseResiduals, diagonalResiduals},
  inverseResiduals = MapThread[
    Simplify[#1 . #2 - msIdentity2] &,
    {data["locals"], data["localsInverse"]}
  ];
  If[! And @@ (msExactZeroMatrixQ /@ inverseResiduals), Return[False]];
  diagonalResiduals = Flatten@Table[
    With[{local = data["slotM0Locals"][[componentPosition, slotPosition]]},
      If[local === None, True,
        With[{transformed = Simplify[
            data["locals"][[slotPosition]] . local .
              data["localsInverse"][[slotPosition]]
          ]},
          msExactZeroMatrixQ[transformed - DiagonalMatrix[Diagonal[transformed]]]
        ]
      ]
    ],
    {componentPosition, Length[data["slotM0Locals"]]},
    {slotPosition, Length[data["slots"]]}
  ];
  And @@ diagonalResiduals
];


(* 非平凡 quotient 下逐列检查 U U^{-1}=I；恒等 quotient 由局部精确逆直接保证。 *)
msTensorQuotientCertificateQ[data_Association] := Module[
  {dimension = data["masterCount"], column, value},
  If[! data["quotientNontrivialQ"], Return[True]];
  Do[
    value = msTensorApplyUFromData[
      data,
      msTensorApplyUInverseFromData[data, UnitVector[dimension, column]]
    ];
    If[! msTensorExactZeroVectorQ[value - UnitVector[dimension, column]], Return[False]],
    {column, dimension}
  ];
  True
];


(* 对每个 component 检查 U M0_c U^{-1} = DiagonalMatrix[kappa_c]。
   逐列只做一次 U^{-1} e_c 与一次 M0 作用，U 的行由转置访问器单独给出；两侧合成矩阵后
   一次判零。这样避免把每个 U^{-1} e_c 再送进 U 的全列展开（后者按 4^n 增长）。
   注意：列向量列表按矩阵读取时是转置，必须显式 Transpose。
   平凡 quotient 的 sector 没有 P/S 混合，逐 slot 局部证书已经蕴含该恒等式。 *)
msTensorM0CertificateQ[data_Association, componentPosition_Integer] := Module[
  {dimension = data["masterCount"], kappa, columns, rows, residual},
  kappa = msTensorMasterComponentKappa[data, componentPosition];
  If[Head[kappa] === Failure, Return[kappa]];
  If[! data["quotientNontrivialQ"], Return[True]];
  columns = Table[
    msTensorM0ApplyFromData[data, componentPosition,
      msTensorApplyUInverseFromData[data, UnitVector[dimension, column]]],
    {column, dimension}
  ];
  rows = Table[
    msTensorApplyUTransposeFromData[data, UnitVector[dimension, row]],
    {row, dimension}
  ];
  residual = rows . Transpose[columns] - DiagonalMatrix[kappa];
  TrueQ[And @@ Flatten[Map[PossibleZeroQ, Map[Together, residual, {2}], {2}]]]
];


msTensorSectorData[sector_Association] := Module[
  {
    slots, slotCount, rawStateCount, masterCount, componentCount,
    locals, localsInverse, canonicalOfRaw, representativeRaw, rawPreimages,
    slotM0Locals, slotM1Coefficients, slotKappaPairs, data, certificate
  },
  slots = sector["slots"];
  slotCount = Length[slots];
  rawStateCount = sector["rawStateCount"];
  masterCount = sector["masterCount"];
  componentCount = Length[sector["vertexComponents"]];
  locals = msLocalM0Diagonalizer[#, sector] & /@ slots;
  If[! FreeQ[locals, _Failure], Return[FirstCase[locals, _Failure]]];
  localsInverse = msLocalM0DiagonalizerInverse[#, sector] & /@ slots;
  If[! FreeQ[localsInverse, _Failure], Return[FirstCase[localsInverse, _Failure]]];
  canonicalOfRaw = msTensorRawPositionOrZero /@ sector["stateEmbedding"];
  representativeRaw = First@FirstPosition[#, 1] & /@ sector["stateProjection"];
  rawPreimages = Table[
    Select[Range[rawStateCount], canonicalOfRaw[[#]] === column &],
    {column, masterCount}
  ];
  slotM0Locals = Table[
    Table[
      msTensorSlotM0Local[slots[[slotPosition]], sector, componentPosition],
      {slotPosition, slotCount}
    ],
    {componentPosition, componentCount}
  ];
  If[! FreeQ[slotM0Locals, _Failure],
    Return[FirstCase[slotM0Locals, _Failure, Missing["UnknownSlotKind"], Infinity]]
  ];
  slotM1Coefficients = Table[
    Table[
      msTensorSlotM1Coefficient[slots[[slotPosition]], sector, componentPosition],
      {slotPosition, slotCount}
    ],
    {componentPosition, componentCount}
  ];
  data = <|
    "sectorKey" -> sector["sectorKey"],
    "slots" -> slots,
    "slotCount" -> slotCount,
    "rawStateCount" -> rawStateCount,
    "masterCount" -> masterCount,
    "componentCount" -> componentCount,
    "vertexComponents" -> sector["vertexComponents"],
    "componentEnergies" -> sector["componentEnergies"],
    "baseTimePowers" -> sector["baseTimePowers"],
    "locals" -> locals,
    "localsInverse" -> localsInverse,
    "canonicalOfRaw" -> canonicalOfRaw,
    "representativeRaw" -> representativeRaw,
    "rawPreimages" -> rawPreimages,
    "slotM0Locals" -> slotM0Locals,
    "slotM1Coefficients" -> slotM1Coefficients,
    "quotientNontrivialQ" -> (sector["coincidentLineIds"] =!= {})
  |>;
  If[! msTensorLocalCertificateQ[data],
    Return[msFailure[
      "TensorM0DiagonalizationNotCertified",
      <|"sectorKey" -> sector["sectorKey"], "reason" -> "localM0"|>
    ]]
  ];
  slotKappaPairs = Table[
    Table[
      With[{local = slotM0Locals[[componentPosition, slotPosition]]},
        If[local === None,
          {0, 0},
          With[{transformed = Simplify[
              locals[[slotPosition]] . local . localsInverse[[slotPosition]]
            ]},
            If[! msExactZeroMatrixQ[transformed - DiagonalMatrix[Diagonal[transformed]]],
              $Failed,
              Simplify[Diagonal[transformed]]
            ]
          ]
        ]
      ],
      {slotPosition, slotCount}
    ],
    {componentPosition, componentCount}
  ];
  If[! FreeQ[slotKappaPairs, $Failed],
    Return[msFailure[
      "TensorM0DiagonalizationNotCertified",
      <|"sectorKey" -> sector["sectorKey"]|>
    ]]
  ];
  data = Join[data, <|"slotKappaPairs" -> slotKappaPairs|>];
  If[! msTensorQuotientCertificateQ[data],
    Return[msFailure[
      "TensorM0DiagonalizationNotCertified",
      <|"sectorKey" -> sector["sectorKey"], "reason" -> "quotientRoundTrip"|>
    ]]
  ];
  Do[
    certificate = msTensorM0CertificateQ[data, componentPosition];
    If[certificate === False || Head[certificate] === Failure, Break[]],
    {componentPosition, componentCount}
  ];
  If[certificate === False || Head[certificate] === Failure,
    Return[If[Head[certificate] === Failure, certificate,
      msFailure[
        "TensorM0DiagonalizationNotCertified",
        <|"sectorKey" -> sector["sectorKey"]|>
      ]
    ]]
  ];
  data
];


msAttachTensorData[sector_Association] := Module[{data = msTensorSectorData[sector]},
  If[Head[data] === Failure, data, Join[sector, <|"tensor" -> data|>]]
];


(* 缓存必须和当前 slots/sectorKey 一致；测试会替换 slots 构造非法 sector，不能复用旧张量数据。 *)
msTensorDataOf[sector_Association] := Module[{cached = Lookup[sector, "tensor", None]},
  If[
    AssociationQ[cached] && cached["sectorKey"] === sector["sectorKey"] &&
      cached["slots"] === sector["slots"],
    cached,
    msTensorSectorData[sector]
  ]
];


(* NDE 线格式的 sector 记录：局部基标签 + quotient 稀疏映射 + 全局偏移。
   不包含任何 2^n x 2^n 数据；后端按 slot 轴逐层作用重建 U。 *)
msTensorSectorWireRecord[sector_Association, offset_Integer] := Module[
  {data = msTensorDataOf[sector], tags},
  If[Head[data] === Failure, Return[data]];
  tags = msTensorLocalBaseTag[#, sector] & /@ sector["slots"];
  If[! FreeQ[tags, _Missing],
    Return[msFailure[
      "UnsupportedLocalBase",
      <|"sectorKey" -> sector["sectorKey"], "tags" -> tags|>
    ]]
  ];
  <|
    "sectorKey" -> sector["sectorKey"],
    "dimension" -> data["masterCount"],
    "offset" -> offset,
    "rawCount" -> data["rawStateCount"],
    "slotBases" -> tags,
    "rawToMasterColumn" -> data["canonicalOfRaw"],
    "masterRepresentativeRaw" -> data["representativeRaw"]
  |>
];


(* ::Chapter:: *)
(* Master-space wrappers *)

msTensorApplyU[sector_Association, vector_List] := Module[
  {data = msTensorDataOf[sector]},
  If[Head[data] === Failure, data, msTensorApplyUFromData[data, vector]]
];

msTensorApplyUInverse[sector_Association, vector_List] := Module[
  {data = msTensorDataOf[sector]},
  If[Head[data] === Failure, data, msTensorApplyUInverseFromData[data, vector]]
];

msTensorApplyUTranspose[sector_Association, vector_List] := Module[
  {data = msTensorDataOf[sector]},
  If[Head[data] === Failure, data, msTensorApplyUTransposeFromData[data, vector]]
];

msTensorApplyUInverseTranspose[sector_Association, vector_List] := Module[
  {data = msTensorDataOf[sector]},
  If[Head[data] === Failure, data, msTensorApplyUInverseTransposeFromData[data, vector]]
];

msTensorApplyM0[sector_Association, componentPosition_Integer, vector_List] := Module[
  {data = msTensorDataOf[sector]},
  If[Head[data] === Failure, data, msTensorM0ApplyFromData[data, componentPosition, vector]]
];

msEnergyLetters[sector_Association, componentPosition_Integer] := Module[
  {data = msTensorDataOf[sector]},
  If[Head[data] === Failure, data,
    msTensorMasterComponentLetters[data, componentPosition]]
];

msM1DiagonalVector[sector_Association, componentPosition_Integer, integerShift_: 0] := Module[
  {data = msTensorDataOf[sector]},
  If[Head[data] === Failure, data,
    msTensorM1Diagonal[data, componentPosition, integerShift]]
];

msM1SingularSurfaces[sector_Association, componentPosition_Integer, integerShift_: 0] := Module[
  {diagonal = msM1DiagonalVector[sector, componentPosition, integerShift]},
  If[Head[diagonal] === Failure, diagonal,
    DeleteDuplicates[Simplify[diagonal]]]
];

msTensorDenseLocalsProduct[sector_Association] := Module[
  {data = msTensorDataOf[sector]},
  If[Head[data] === Failure, data, msTensorKroneckerAll[data["locals"]]]
];

msSectorDiagonalizer[sector_Association] := Module[
  {data = msTensorDataOf[sector], product},
  If[Head[data] === Failure, Return[data]];
  product = msTensorKroneckerAll[data["locals"]];
  sector["stateProjection"] . product . sector["stateEmbedding"]
];

msSectorDiagonalizerInverse[sector_Association] := Module[
  {data = msTensorDataOf[sector], product},
  If[Head[data] === Failure, Return[data]];
  product = msTensorKroneckerAll[data["localsInverse"]];
  sector["stateProjection"] . product . sector["stateEmbedding"]
];
