(* ::Package:: *)

(***
文件：Representation.wl
用途：验证 MSIntegral、公开 normalized-to-bare 定义、执行 H/h 状态变换，并适配 dSIBP time-only J。
范围：所有表示共用 sectorKey、timeShifts、stateBits 三槽；massless shared quotient 在 H/h
变换中恒为单位作用，不重建未约四态表示或第二套指标约定。
***)

(* ::Chapter:: *)
(* Native integral validation *)

msIntegralData[int : MSIntegral[key_String, shifts_List, bits_List], context_?MSContextQ] := Module[
  {sector = msSectorByKey[context, key]},
  If[Head[sector] === Missing, Return[msFailure["UnknownSector", <|"sectorKey" -> key|>]]];
  If[Length[shifts] =!= Length[sector["vertexComponents"]] || ! And @@ (IntegerQ /@ shifts),
    Return[msFailure["BadTimeShifts", <|"expected" -> Length[sector["vertexComponents"]], "actual" -> shifts|>]]
  ];
  If[! MemberQ[sector["stateOrder"], bits],
    Return[msFailure[
      "BadStateBits",
      <|"expectedStateOrder" -> sector["stateOrder"], "actual" -> bits|>
    ]]
  ];
  <|"sector" -> sector, "shifts" -> shifts, "bits" -> bits|>
];

msIntegralData[other_, _] := msFailure["NotMSIntegral", <|"input" -> HoldForm[other]|>];


(* ::Chapter:: *)
(*Normalized master 与裸指标积分*)

(* 裸积分只复用 MSIntegral 的三个既有指标。normalization 始终读取当前 sector，
   definition 仅用于显示和审计，不参加递推、求导或数值输运。 *)
MSIntegralDefinition[
  int : MSIntegral[key_String, shifts_List, bits_List],
  context_?MSContextQ
] := Module[{data, bareIntegral, normalization, definition},
  data = msIntegralData[int, context];
  If[Head[data] === Failure, Return[data]];
  bareIntegral = MSBareIntegral[key, shifts, bits];
  normalization = data["sector", "normalization"];
  definition = If[
    TrueQ[normalization === 1],
    With[{normalized = int, bare = bareIntegral}, HoldForm[normalized == bare]],
    With[
      {normalized = int, bare = bareIntegral, coefficient = normalization},
      HoldForm[normalized == coefficient bare]
    ]
  ];
  <|
    "integral" -> int,
    "bareIntegral" -> bareIntegral,
    "normalization" -> normalization,
    "definition" -> definition
  |>
];


MSIntegralDefinition[other_, context_?MSContextQ] := msIntegralData[other, context];


(* 主积分列表继续以 context 中的认证顺序为 authority；新增字段全部来自上面的单对象入口。 *)
msMasterRecordWithDefinition[record_Association, context_?MSContextQ] := Module[{definition},
  definition = MSIntegralDefinition[record["integral"], context];
  If[Head[definition] === Failure, Return[definition]];
  Join[
    record,
    KeyTake[definition, {"bareIntegral", "normalization", "definition"}]
  ]
];


MSMasterIntegrals[context_?MSContextQ] := msMasterRecordWithDefinition[#, context] & /@ context["masters"];


(* ::Section:: *)
(* All-sector H/h state-vector transformation *)

(* 输入按 sector stateBits 排序。massive endpoint 使用 z=-k tau；massless shared quotient
   不是裸 Hankel endpoint pair，因此在 H/h 变换下恒为单位作用。 *)
msSectorBasisTransformMatrix[
  direction_Rule,
  sector_Association,
  componentTimes_List,
  nuConvention_String
] := Module[{localMatrices, componentPosition, zValue},
  If[! MemberQ[{"H" -> "h", "h" -> "H"}, direction],
    Return[msFailure["UnsupportedBasisDirection", <|"direction" -> direction|>]]
  ];
  localMatrices = Map[
    Function[slot,
      Switch[slot["kind"],
        "masslessShared", msIdentity2,
        "massiveEndpoint",
        componentPosition = sector["rootToComponent"][slot["rootVertex"]];
        zValue = -slot["momentum"] componentTimes[[componentPosition]];
        If[direction === ("H" -> "h"),
          msHTohMatrix[slot["nu"], zValue, nuConvention],
          msHToHMatrix[slot["nu"], zValue, nuConvention]
        ],
        _, msFailure["UnsupportedTensorSlotKind", <|"slot" -> slot|>]
      ]
    ],
    sector["slots"]
  ];
  If[! FreeQ[localMatrices, _Failure], Return[FirstCase[localMatrices, _Failure]]];
  msQuotientMatrix[sector, msKroneckerAll[localMatrices]]
];

(*
Returns the sector state vector in the same order as the input. This interface acts on integrand states and does not change the base time powers of MSIntegral; basis changes of integral objects require additional target-family metadata and remain rejected for now.
*)
MSConvertBasis[
  vector_List,
  direction_Rule,
  key_String,
  componentTimes_List,
  context_?MSContextQ
] := Module[{sector, matrix, nuConvention},
  sector = msSectorByKey[context, key];
  If[Head[sector] === Missing, Return[msFailure["UnknownSector", <|"sectorKey" -> key|>]]];
  If[Length[componentTimes] =!= Length[sector["vertexComponents"]],
    Return[msFailure[
      "ComponentTimeCount",
      <|"expected" -> Length[sector["vertexComponents"]], "actual" -> Length[componentTimes]|>
    ]]
  ];
  If[Length[vector] =!= sector["masterCount"],
    Return[msFailure[
      "StateVectorDimension",
      <|"expected" -> sector["masterCount"], "actual" -> Length[vector]|>
    ]]
  ];
  nuConvention = context["convention"]["nuConvention"];
  matrix = msSectorBasisTransformMatrix[direction, sector, componentTimes, nuConvention];
  If[Head[matrix] === Failure, matrix, Simplify[matrix.vector]]
];


(* ::Chapter:: *)
(*dSIBP time-only J[sectorKey,timeShifts,stateBits] adapter*)


(* Input must first pass validation against the current MadStree context; output stays lazy to avoid evaluation when dSIBP is not loaded. *)
MSToDSIBPJ[int : MSIntegral[_, _, _], context_?MSContextQ] := Module[
  {data},
  data = msIntegralData[int, context];
  If[Head[data] === Failure, Return[data]];
  Inactive[dSIBP`J][data["sector"]["sectorKey"], data["shifts"], data["bits"]]
];


msDSIBPActiveJToInactive[expr_] := Replace[
  Unevaluated[expr],
  dSIBP`J[key_String, shifts_List, bits_List] :> Inactive[dSIBP`J][key, shifts, bits]
];


(* The reverse entry validates sector, time-shift dimension and discrete states slot by slot; inconsistent schemas return a native Failure. *)
MSFromDSIBPJ[input_, context_?MSContextQ] := Module[
  {held, key, shifts, bits, candidate, validation},
  held = msDSIBPActiveJToInactive[input];
  If[! MatchQ[held, Inactive[dSIBP`J][_String, _List, _List]],
    Return[msFailure["UnsupportedDSIBPJ", <|"input" -> HoldForm[input]|>]]
  ];
  {key, shifts, bits} = List @@ held;
  candidate = MSIntegral[key, shifts, bits];
  validation = msIntegralData[candidate, context];
  If[Head[validation] === Failure, validation, candidate]
];


MSFromDSIBPExpression[expr_, context_?MSContextQ] := Module[
  {integrals, mapped, failures},
  integrals = DeleteDuplicates@Cases[expr, _dSIBP`J, Infinity];
  mapped = MSFromDSIBPJ[#, context] & /@ integrals;
  failures = Pick[mapped, Head /@ mapped, Failure];
  If[failures =!= {},
    Return[msFailure[
      "DSIBPExpressionConversionFailed",
      <|"integrals" -> integrals, "failures" -> failures|>
    ]]
  ];
  Expand[expr /. Thread[integrals -> mapped]]
];
