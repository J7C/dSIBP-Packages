(* ::Package:: *)
(* Massive function-system 示例：输入裸 Hankel P,Q,T,W，并通过公开工作流编译和生成 seed。 *)

(* ::Chapter:: *)
(*加载标准 package*)

exampleDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{exampleDir, "load_current_package.wl"}], CharacterEncoding -> "UTF-8"];


(* ::Chapter:: *)
(*函数系统与 topology 输入*)

prefactor = (4 I/Pi) Exp[Pi Im[nuM]];
hankelSystem = <|
   "variable" -> x,
   "P" -> 1/x,
   "Q" -> 1 - nuM^2/x^2,
   "T" -> IdentityMatrix[2],
   "W" -> -prefactor/x,
   "WT" -> Automatic,
   "shrinkBShift" -> 1,
   "shrinkZeroPointShift" -> 0
   |>;

case = <|
   "name" -> "bareHankelExample",
   "vertices" -> {
     <|"id" -> v1, "vertexType" -> "+", "externalLegEnergy" -> E1|>,
     <|"id" -> v2, "vertexType" -> "+", "externalLegEnergy" -> E2|>
     },
   "lines" -> {
     <|"id" -> e1, "endpoints" -> {v1, v2}, "momentum" -> q,
       "nu" -> nuM, "massType" -> "massive", "functionSystem" -> hankelSystem|>
     },
   (* 单边连接两顶点的结构圈数为零；q 在本例只是该线的固定模长。 *)
   "loopMomenta" -> {},
   "loopExternalMomenta" -> {},
   "independentExternalMomenta" -> {q},
   "ibpMode" -> "timeOnly",
   "zeroPointRules" -> {
     a0[v1] -> alpha1, a0[v2] -> alpha2, b0[e1] -> beta1
     },
   "symmetryRules" -> {}
   |>;


(* ::Chapter:: *)
(*公开工作流*)

(* 缺省不写初始化文件；本例枚举全部离散态但保持最小连续范围。 *)
context = DSInit[
   case,
   WriteInitializationFiles -> False,
   GenerateDerivativeMetadata -> False
   ];
seedData = DSSeeds[context];

fullInfo = DSInfo[context, "Full"];
summary = <|
   "initStatus" -> Lookup[context, "status", "missing"],
   "seedStatus" -> Lookup[seedData, "dSIBPStatus", "missing"],
   "sectorCount" -> Length[Lookup[context, "sectors", {}]],
   "functionSystemInputPreserved" -> ! FreeQ[fullInfo, hankelSystem]
   |>;

Print[summary];
If[! And[
    ToString[dSIBP`$dSIBPVersion] === currentVersion,
    summary["initStatus"] === "initialized",
    summary["seedStatus"] === "generated",
    TrueQ[summary["functionSystemInputPreserved"]]
    ], Exit[1]];
