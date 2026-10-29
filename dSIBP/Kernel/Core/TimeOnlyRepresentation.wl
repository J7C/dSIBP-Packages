(* ::Package:: *)
(***
文件：TimeOnlyRepresentation.wl
用途：定义 time-only 公开积分 J[sectorKey,timeShifts,stateBits] 与内部
      J[aList,linePacks,{}] 的唯一双向转换，并提供表达式/数据容器级转换。
核心逻辑：sectorKey 只负责 sector 身份；timeShifts 保存当前 sector 的紧致时间幂指标；
          stateBits 按初始化 metadata 冻结的 root line/endpoint 顺序保存离散 n 态。
边界：full 模式不经过本模块；三槽 time-only 对象只允许存在于 Private 内部。
***)


(* ::Chapter:: *)
(*Time-only sector 与 slot metadata*)

(* 从 context 或 parsed topology 取得同一份 root topology；转换规则不得依赖调用方猜测。 *)
dsTimeOnlyTopology[data_Association] := If[
   KeyExistsQ[data, "topology"] && AssociationQ[data["topology"]],
   data["topology"],
   data
   ];


dsTimeOnlyMetadataList[data_Association] := Module[{topo = dsTimeOnlyTopology[data]},
   If[
    KeyExistsQ[data, "sectors"],
    Lookup[data, "sectors", {}],
    Lookup[topo, "sectorMetadataList", {makeSectorMetadata[topo]}]
    ]
   ];


dsTimeOnlyModeQ[data_Association] :=
  Lookup[dsTimeOnlyTopology[data], "ibpMode", "full"] === "timeOnly";


dsTimeOnlySectorMetadata[data_Association, sectorKey_String] := SelectFirst[
   dsTimeOnlyMetadataList[data],
   Lookup[#, "sectorKey", None] === sectorKey &,
   Missing["UnknownTimeOnlySector", sectorKey]
   ];


(* 每个 slot 直接保存旧 line pack 中的位置；双向转换因此只读取 metadata，不重复推断图。 *)
dsTimeOnlyStateSlots[topo_Association, metadata_Association] := Module[
   {lines = Lookup[topo, "lines", {}], lineSlots = Lookup[metadata, "lineSlots", {}]},
   Flatten@MapThread[
     Function[{line, lineSlot, linePosition},
      If[TrueQ[Lookup[lineSlot, "shrunkQ", False]],
       {},
       Switch[
        Lookup[line, "packType", ""],
        "massiveFull" | "massiveCross",
        Table[
         <|
          "lineId" -> Lookup[line, "id", linePosition],
          "rootLinePosition" -> linePosition,
          "kind" -> "massiveEndpoint",
          "endpointSlot" -> endpointSlot,
          "packPosition" -> 1 + endpointSlot
          |>,
         {endpointSlot, 2}
         ],
        "masslessFull",
        {<|
          "lineId" -> Lookup[line, "id", linePosition],
          "rootLinePosition" -> linePosition,
          "kind" -> "masslessShared",
          "endpointSlot" -> "shared",
          "packPosition" -> 2
          |>},
        _,
        {}
        ]
       ]
      ],
     {lines, lineSlots, Range[Length[lines]]}
     ]
   ];


dsTimeOnlyStateSlotsForMetadata[topo_Association, metadata_Association] := Lookup[
   metadata,
   "timeOnlyStateSlots",
   dsTimeOnlyStateSlots[topo, metadata]
   ];


(* ::Chapter:: *)
(*单积分双向转换*)

dsTimeOnlyInternalIntegralToPublic[
   int : J[aList_List, linePacks_List, {}],
   data_Association
   ] := Module[{topo, metadata, sectorKey, slots, stateBits},
   topo = dsTimeOnlyTopology[data];
   If[! dsTimeOnlyModeQ[data], Return[int]];
   metadata = integralSectorMetadata[topo, int];
   If[Head[metadata] === Missing, Return[$Failed]];
   sectorKey = Lookup[metadata, "sectorKey", Missing["NoSectorKey"]];
   slots = dsTimeOnlyStateSlotsForMetadata[topo, metadata];
   If[! StringQ[sectorKey] || Length[aList] =!= Length[Lookup[metadata, "compactASlots", {}]],
    Return[$Failed]
    ];
   stateBits = Extract[
      linePacks,
      ({#rootLinePosition, #packPosition} &) /@ slots
      ];
   If[! VectorQ[stateBits, MemberQ[{0, 1}, #] &], Return[$Failed]];
   J[sectorKey, aList, stateBits]
   ];


dsTimeOnlyInternalIntegralToPublic[int_J, data_Association] := If[
   dsTimeOnlyModeQ[data],
   $Failed,
   int
   ];


dsTimeOnlyPublicIntegralToInternal[
   int : J[sectorKey_String, timeShifts_List, stateBits_List],
   data_Association
   ] := Module[{topo, metadata, slots, linePacks, lines, lineSlots},
   topo = dsTimeOnlyTopology[data];
   If[! dsTimeOnlyModeQ[data], Return[$Failed]];
   metadata = dsTimeOnlySectorMetadata[data, sectorKey];
   If[Head[metadata] === Missing, Return[$Failed]];
   slots = dsTimeOnlyStateSlotsForMetadata[topo, metadata];
   If[
    Length[timeShifts] =!= Length[Lookup[metadata, "compactASlots", {}]] ||
     Length[stateBits] =!= Length[slots] ||
     ! VectorQ[stateBits, MemberQ[{0, 1}, #] &],
    Return[$Failed]
    ];
   lineSlots = Lookup[metadata, "lineSlots", {}];
   lines = Lookup[topo, "lines", {}];
   linePacks = Lookup[lineSlots, "packTemplate", {}];
   If[Length[linePacks] =!= Length[Lookup[topo, "lines", {}]], Return[$Failed]];
   Do[
    linePacks = ReplacePart[
      linePacks,
      {slots[[slotPosition, "rootLinePosition"]], slots[[slotPosition, "packPosition"]]} ->
       stateBits[[slotPosition]]
      ],
    {slotPosition, Length[slots]}
    ];
   (* quotient massless 的被消去端点槽固定为 0；cross/收缩占位由 sector packTemplate 恢复。 *)
   Do[
    If[
     ! TrueQ[Lookup[lineSlots[[linePosition]], "shrunkQ", False]] &&
      Lookup[lines[[linePosition]], "packType", ""] === "masslessFull",
     linePacks = ReplacePart[linePacks, {linePosition, 3} -> 0]
     ],
    {linePosition, Length[lines]}
    ];
   J[timeShifts, linePacks, {}]
   ];


dsTimeOnlyPublicIntegralToInternal[int_J, data_Association] := If[
   dsTimeOnlyModeQ[data],
   $Failed,
   int
   ];


(* ::Chapter:: *)
(*表达式与 Association 边界转换*)

dsTimeOnlyExpressionToPublic[expr_, data_Association] := Module[{result},
   If[! dsTimeOnlyModeQ[data], Return[expr]];
   result = expr /. int : J[_List, _List, {}] :>
      dsTimeOnlyInternalIntegralToPublic[int, data];
   If[FreeQ[result, $Failed], result, $Failed]
   ];


dsTimeOnlyExpressionToInternal[expr_, data_Association] := Module[{result, oldIntegrals},
   If[! dsTimeOnlyModeQ[data], Return[expr]];
   oldIntegrals = ! FreeQ[expr, J[_List, _List, _List]];
   If[TrueQ[oldIntegrals], Return[$Failed]];
   result = expr /. int : J[_String, _List, _List] :>
      dsTimeOnlyPublicIntegralToInternal[int, data];
   If[FreeQ[result, $Failed], result, $Failed]
   ];


dsTimeOnlyDataToPublic[value_Association, data_Association] := Association@KeyValueMap[
   dsTimeOnlyDataToPublic[#1, data] -> dsTimeOnlyDataToPublic[#2, data] &,
   value
   ];
dsTimeOnlyDataToPublic[value_List, data_Association] :=
  dsTimeOnlyDataToPublic[#, data] & /@ value;
dsTimeOnlyDataToPublic[value_, data_Association] :=
  dsTimeOnlyExpressionToPublic[value, data];


dsTimeOnlyDataToInternal[value_Association, data_Association] := Association@KeyValueMap[
   dsTimeOnlyDataToInternal[#1, data] -> dsTimeOnlyDataToInternal[#2, data] &,
   value
   ];
dsTimeOnlyDataToInternal[value_List, data_Association] :=
  dsTimeOnlyDataToInternal[#, data] & /@ value;
dsTimeOnlyDataToInternal[value_, data_Association] :=
  dsTimeOnlyExpressionToInternal[value, data];


dsTimeOnlyOldIntegralLeakQ[value_Association] := AnyTrue[
   Values[value],
   dsTimeOnlyOldIntegralLeakQ
   ];
dsTimeOnlyOldIntegralLeakQ[value_List] := AnyTrue[
   value,
   dsTimeOnlyOldIntegralLeakQ
   ];
dsTimeOnlyOldIntegralLeakQ[value_] := ! FreeQ[value, J[_List, _List, _List]];
