(* ::Package:: *)
(* Mixed massive/massless 示例：用统一公开入口完成离散模板、连续撒点和 linearData。 *)

(* ::Chapter:: *)
(*加载标准 package*)

exampleDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{exampleDir, "load_current_package.wl"}], CharacterEncoding -> "UTF-8"];


(* ::Chapter:: *)
(*详细物理输入*)

case = <|
   "name" -> "mixedBubbleExample",
   "vertices" -> {
     <|"id" -> v1, "vertexType" -> "+", "externalLegEnergy" -> E1|>,
     <|"id" -> v2, "vertexType" -> "+", "externalLegEnergy" -> E2|>
     },
   "lines" -> {
     <|"id" -> e1, "endpoints" -> {v1, v2}, "momentum" -> q,
       "nu" -> nuM, "massType" -> "massive"|>,
     <|"id" -> e2, "endpoints" -> {v1, v2}, "momentum" -> q - k,
       "massType" -> "massless"|>
     },
   "loopMomenta" -> {q},
   "loopExternalMomenta" -> {k},
   "independentExternalMomenta" -> {},
   "ibpMode" -> "full",
   "ispData" -> {},
   "zeroPointRules" -> {
     a0[v1] -> alpha1, a0[v2] -> alpha2,
     b0[e1] -> beta1, b0[e2] -> beta2
     },
   "symmetryRules" -> {}
   |>;


(* ::Chapter:: *)
(*缺省选项*)

(* 缺省不写初始化文件、不生成导数 metadata；DSSeeds 始终完整枚举离散态。 *)
initOptions = {WriteInitializationFiles -> False, GenerateDerivativeMetadata -> False};
linearOptions = {LinearSystemMode -> "symbolic"};


(* ::Chapter:: *)
(*公开工作流*)

context = DSInit[case, Sequence @@ initOptions];
templateSeedData = DSSeeds[context];
allSeeds = DSAllSeeds[templateSeedData];

(* 统一目标包络覆盖当前 producer 的最宽 {-2,2} shift；更窄范围会令部分分组无 seed 点。 *)
seedData = DSGenerateIBP[allSeeds, {-2, 2}];

(* 精细范围：必须把 allSeeds 中实际出现的连续指标逐个完整覆盖。 *)
continuousIndices = DeleteDuplicates@Flatten[Lookup[allSeeds, "continuousIndices", {}], Infinity];
detailedRanges = ({#, -2, 2} & /@ continuousIndices);
seedDataDetailed = DSGenerateIBP[allSeeds, Sequence @@ detailedRanges];

linearData = DSLinear[seedData, context, Sequence @@ linearOptions];
(* 后端需要不同编号时只在 linearData 层显式重排一次；plan 不再自行改序。 *)
integralOrder = Lookup[linearData, "integralList", {}];
If[Length[integralOrder] >= 2,
 linearData = DSReorderIntegrals[linearData, Reverse@Take[integralOrder, 2]]
 ];
kiraPlan = DSKiraPlan[linearData, <|
    "stage" -> "preReduction",
    "candidateIntegrals" -> Take[Lookup[linearData, "integralList", {}], UpTo[3]]
    |>];

summary = <|
   "initStatus" -> Lookup[context, "status", "missing"],
   "sectorCount" -> Length[Lookup[context, "sectors", {}]],
   "templateStatus" -> Lookup[templateSeedData, "dSIBPStatus", "missing"],
   "templateCount" -> Length[allSeeds],
   "seedStatus" -> Lookup[seedData, "dSIBPStatus", "missing"],
   "equationCount" -> Lookup[seedData, "equationCount", Missing["NotAvailable"]],
   "uniformDetailedEqual" -> (Lookup[seedData, "equations", {}] === Lookup[seedDataDetailed, "equations", {}]),
   "linearStatus" -> Lookup[linearData, "dSIBPStatus", "missing"],
   "integralOrderAuthority" -> Lookup[linearData, "integralOrderAuthority", Missing["NotReordered"]],
   "kiraPlanStatus" -> Lookup[kiraPlan, "status", "missing"]
   |>;

Print[summary];
If[! And[
    ToString[dSIBP`$dSIBPVersion] === currentVersion,
    summary["initStatus"] === "initialized",
    summary["templateStatus"] === "generated",
    summary["seedStatus"] === "generated",
    TrueQ[summary["uniformDetailedEqual"]],
    summary["linearStatus"] === "generated",
    summary["kiraPlanStatus"] === "planned"
    ], Exit[1]];
