(* ::Package:: *)

(***
文件：07_zero_external_leg_energy.wl
用途：演示没有物理外腿的顶点可省略 externalLegEnergy 或显式输入 0；两种输入生成同一解析 dlog，数值接口不暴露私有辅助能量并自动输运到物理零点。
运行：可在 Mathematica 前端逐节执行，或用 wolframscript -file 运行整个文件。
***)

(* ::Chapter:: *)
(*加载与拓扑*)

packageRoot = DirectoryName[DirectoryName[$InputFileName]];
AppendTo[$Path, FileNameJoin[{packageRoot, "Kernel"}]];
Needs["MadStree`"];

commonLine = {
  <|"type" -> "massless", "endpoints" -> {1, 2}, "momentum" -> q|>
};
omittedSpec = <|
  "vertices" -> {
    <|"id" -> 1, "timePower" -> 1, "vertexType" -> "+"|>,
    <|"id" -> 2, "externalLegEnergy" -> k2, "timePower" -> 1,
      "vertexType" -> "+"|>
  },
  "lines" -> commonLine
|>;
zeroSpec = ReplacePart[
  omittedSpec,
  {"vertices", 1} -> Append[omittedSpec["vertices"][[1]], "externalLegEnergy" -> 0]
];

omittedContext = MSInitTree[omittedSpec];
zeroContext = MSInitTree[zeroSpec];


(* ::Chapter:: *)
(*只输入物理坐标*)

pointSequence = {{k2}, {3 I}};
result = MSEvaluatePath[
  omittedContext,
  pointSequence,
  ParameterRules -> {q -> 1},
  BoundaryScale -> 4,
  WorkingPrecision -> 40,
  TransportOrder -> 80,
  ReferenceTransportOrder -> 104,
  TargetRelativeError -> "1e-20",
  MessageLanguage -> "CN"
];


(* ::Chapter:: *)
(*验收*)

checks = <|
  "sameContext" -> omittedContext["masterDigest"] === zeroContext["masterDigest"],
  "sameAnalyticDLog" -> MSDLogDE[omittedContext] === MSDLogDE[zeroContext],
  "auxiliaryTargetIsZero" ->
    omittedContext["vertices"][[1, "externalLegEnergyTarget"]] === 0,
  "publicPointSequenceHasOnlyPhysicalCoordinate" -> First[pointSequence] === {k2},
  "computed" -> Lookup[result, "status", None] === "computed"
|>;

If[! And @@ Values[checks],
  Print["Example FAILED: ", InputForm[Select[checks, Not]]];
  Exit[1]
];
Print["Example PASSED: ", Count[Values[checks], True], "/", Length[checks],
  " checks; omitted and explicit-zero external-leg inputs are equivalent."];
