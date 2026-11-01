(* ::Package:: *)
(* 文件用途：演示 FlintNDE Wolfram 接口对不同固定 ep 的 NDE 任务做有界并行。
   方程 y'(x)=ep/(1+x)y(x)、y(0)=1 的闭式终值为 y(1)=2^ep。
   ParallelTaskCount 缺省为 12；实际并发自动取该值与 ep 数量的较小者。 *)

(* ::Chapter:: *)
(*程序包加载与 ep 任务*)

exampleDirectory = DirectoryName[$InputFileName];
(* 发布布局直接取上级目录；发展仓布局取 versions/ 下版本号最大的当前版本目录。 *)
packageRoot = ExpandFileName[FileNameJoin[{exampleDirectory, ".."}]];
If[! DirectoryQ[FileNameJoin[{packageRoot, "flintnde"}]],
  versionDirs =
    FileNames["FlintNDE-*", FileNameJoin[{packageRoot, "versions"}]];
  If[Length[versionDirs] == 0,
    Print["no versions/FlintNDE-* directory beside examples/"];
    Abort[]
  ];
  versionKey[name_] := Module[{parts},
    parts = StringPadLeft[StringSplit[
      StringReplace[StringDrop[name, StringLength["FlintNDE-"]],
        StartOfString ~~ "v" -> ""], "."], 4, "0"];
    StringRiffle[Join[parts,
      ConstantArray["0000", Max[0, 4 - Length[parts]]]], "."]
  ];
  packageRoot = ExpandFileName[Last[SortBy[versionDirs, versionKey]]]
];
versionDirectory = packageRoot;
runtimeDirectory = FileNameJoin[{exampleDirectory, "results_temp"}];
PrependTo[$Path, versionDirectory];
Needs["FlintNDE`"];

epValues = {1/5, 1/6, 1/7};
epParallelTaskCount = 12; (* 缺省值；例如改为 4 即最多并行四个 ep。 *)
jobs = Map[
  <|
    "ep" -> #,
    "system" -> FlintNDERationalSystem[{{#/(1 + x)}}, x],
    "start" -> 0,
    "points" -> {1},
    "initialVector" -> {1}
  |> &,
  epValues
];


(* ::Chapter:: *)
(*有界并行执行与闭式检查*)

batch = FlintNDEEvaluateEpBatch[
  jobs,
  FlintNDE`ParallelTaskCount -> epParallelTaskCount,
  FlintNDE`MessageLanguage -> "CN",
  "WorkDirectory" -> runtimeDirectory,
  "WorkingPrecisionDigits" -> 60,
  "OutputDigits" -> 40,
  "PrimaryOrder" -> 48,
  "ReferenceOrder" -> 64,
  "TargetRelativeError" -> "1e-30"
];

If[Head[batch] === Failure,
  Print["ep_parallel_mathematica: FAILED: ", InputForm[batch]];
  Exit[1]
];

values = First[#1["execution", "referenceFinalVector"]] & /@ batch["results"];
maximumDifference = Max[Abs[values - N[2^epValues, 40]]];
checks = <|
  "defaultParallelTaskCount" ->
    MemberQ[Options[FlintNDEEvaluateEpBatch], FlintNDE`ParallelTaskCount -> 12],
  "effectiveParallelTaskCount" ->
    batch["parallelTaskCountEffective"] === Min[epParallelTaskCount, Length[epValues]],
  "inputOrder" -> Lookup[batch["results"], "ep"] === epValues,
  "closedForm" -> TrueQ[maximumDifference < 10^-30]
|>;

Print["parallel requested/effective: ",
  batch["parallelTaskCountRequested"], "/", batch["parallelTaskCountEffective"]];
Print["maximum difference from 2^ep: ", maximumDifference];
Print["checks: ", Count[Values[checks], True], "/", Length[checks]];
If[! And @@ Values[checks], Exit[1]];
Print["ep_parallel_mathematica.wl: PASSED"];
