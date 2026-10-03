(* ::Package:: *)
(* 文件用途：演示 FlintNDE 标准 Wolfram Language 接口。
   功能范围：通过 Needs 加载程序包，构造单变量 exact 有理矩阵 DE，先规划路径再执行已有
   计划，并与闭式结果比较。运行期 bridge 文件只写入调用目录的 results_temp。 *)


(* ::Chapter:: *)
(*程序包加载与路径配置*)

exampleDirectory = DirectoryName[$InputFileName];
versionDirectory = ExpandFileName[FileNameJoin[{
  exampleDirectory, ".."
}]];
runtimeDirectory = FileNameJoin[{exampleDirectory, "results_temp"}];
resultDirectory = FileNameJoin[{
  exampleDirectory, "results", "mathematica_interface_example"
}];
$Path = Prepend[$Path, versionDirectory];
Needs["FlintNDE`"];


(* ::Chapter:: *)
(*单变量微分方程与路径规划*)

(* ::Section:: *)
(*Exact 有理矩阵系统*)

(* A(z)=diag(1/(z-1),-2/(z+3)); y(0)={1,1}. *)
system = FlintNDERationalSystem[
  {{1/(z - 1), 0}, {0, -2/(z + 3)}},
  z,
  "Name" -> "mathematica-interface-example"
];


(* ::Section:: *)
(*从原始用户点生成一次计划*)

plan = FlintNDEPlanPath[
  system,
  0,
  {1/4, 1/2},
  FlintNDE`SingularityMode -> "Avoid",
  FlintNDE`MessageLanguage -> "EN",
  "WorkDirectory" -> runtimeDirectory,
  "WorkingPrecisionDigits" -> 80,
  "OutputDigits" -> 40
];

If[! AssociationQ[plan] || Lookup[plan, "status", "error"] =!= "complete",
  Print["mathematica_interface_example: FAILED while planning: ", InputForm[plan]];
  Exit[1]
];


(* ::Chapter:: *)
(*已有计划执行与验证*)

(* ::Section:: *)
(*执行阶段不重新规划*)

result = FlintNDEExecutePath[
  system,
  {1, 1},
  plan,
  FlintNDE`MessageLanguage -> "EN",
  "WorkDirectory" -> runtimeDirectory,
  "WorkingPrecisionDigits" -> 80,
  "OutputDigits" -> 40,
  "PrimaryOrder" -> 40,
  "ReferenceOrder" -> 48,
  "TargetRelativeError" -> "1e-30"
];

If[! AssociationQ[result] || Lookup[result, "status", "error"] =!= "complete",
  Print["mathematica_interface_example: FAILED while executing: ", InputForm[result]];
  Exit[1]
];


(* ::Section:: *)
(*闭式比较*)

expectedFinal = N[{1/2, 36/49}, 40];
finalDifference = Max[Abs[N[result["primaryFinalVector"], 40] - expectedFinal]];
Print["certification mode: ", result["certificationMode"]];
Print["relative difference (primary vs reference): ", result["relativeDifferenceInf"]];
Print["final |delta| vs closed form: ", finalDifference];

If[! TrueQ[finalDifference < 10^-30],
  Print["mathematica_interface_example: FAILED numerical check"];
  Exit[1]
];


(* ::Section:: *)
(*用户选择的结果保存*)

If[! DirectoryQ[resultDirectory],
  CreateDirectory[resultDirectory, CreateIntermediateDirectories -> True]
];
Put[result, FileNameJoin[{resultDirectory, "mathematica_interface_result.m"}]];
Print["mathematica_interface_example: PASSED"];
