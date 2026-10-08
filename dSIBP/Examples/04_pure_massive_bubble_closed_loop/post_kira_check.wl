(* ::Package:: *)
(* 本脚本只验证已经完成的真实 Kira reduction，不重新生成 seed、linearData 或 Kira 输入。
   它重建与 main.wl 相同的等能量 family context，再依次检查导入、19 维 DE 和标度关系。 *)

(* ::Chapter:: *)
(*标准 package 与 family convention*)

exampleDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{exampleDir, "..", "load_current_package.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{exampleDir, "dlog_basis.wl"}]];
Get[FileNameJoin[{exampleDir, "reference_user_mi_basis.wl"}]];
Get[FileNameJoin[{exampleDir, "family_conventions.wl"}]];

externalWorkspaceRoot = Environment["DSIBP_KIRA_WORKSPACE"];
If[! StringQ[externalWorkspaceRoot] || StringTrim[externalWorkspaceRoot] === "",
  Print["Set DSIBP_KIRA_WORKSPACE to a directory outside this repository before running this example."];
  Exit[2]
];
externalWorkspaceRoot = ExpandFileName[externalWorkspaceRoot];
repositoryRoot = ExpandFileName@FileNameJoin[{exampleDir, "..", "..", ".."}];
If[
  externalWorkspaceRoot === repositoryRoot ||
    StringStartsQ[externalWorkspaceRoot, repositoryRoot <> $PathnameSeparator],
  Print["The Kira workspace must live outside the current Git repository; change DSIBP_KIRA_WORKSPACE."];
  Exit[2]
];
exampleWorkspace = FileNameJoin[{externalWorkspaceRoot, "pure_massive_bubble"}];


(* 与 export 端共享冻结的精确参数点，不在后处理阶段重新抽取或改值。 *)
parameterProbeSeed = 20260722;
parameterProbeRules = {
   dim -> 37/11, nu -> 7/13, etaNu -> 23/17,
   analyticRegulator -> 0
   };


(* ::Chapter:: *)
(*等能量 family 初始化*)

(* 顶点交换 symmetry 只在 reference P1=P2=-P0 时成立；独立 P1/P2 family 不得复用本配置。 *)
caseInput = <|
   "name" -> "phase2PureMassiveBubbleMinusMinus",
   "vertices" -> {
     <|"id" -> v1, "vertexType" -> "-", "externalLegEnergy" -> P0|>,
     <|"id" -> v2, "vertexType" -> "-", "externalLegEnergy" -> P0|>
     },
   "lines" -> {
     <|"id" -> 1, "endpoints" -> {v1, v2}, "momentum" -> q,
       "massType" -> "massive", "nu" -> nu|>,
     <|"id" -> 2, "endpoints" -> {v1, v2}, "momentum" -> q - k,
       "massType" -> "massive", "nu" -> nu|>
     },
   "loopMomenta" -> {q},
   "loopExternalMomenta" -> {k},
   "independentExternalMomenta" -> {},
   "ibpMode" -> "full",
   "ispData" -> {},
   "zeroPointRules" -> {
     a0[v1] -> 2 nu + analyticRegulator,
     a0[v2] -> 2 nu + analyticRegulator,
     b0[1] -> -2 nu, b0[2] -> -2 nu
     },
   "symmetryRules" -> exampleSymmetryRules0
   |>;

context = DSInit[
   caseInput,
   WriteInitializationFiles -> False,
   GenerateDerivativeMetadata -> False,
   OverwriteInitialization -> False
   ];
If[Lookup[context, "status", "failed"] =!= "initialized",
  Print["Pure massive bubble 初始化失败；请先修正上方报告的输入或正规化问题，Kira 结果不会导入。"];
  Exit[1]
];


(* ::Chapter:: *)
(*真实 Kira 结果取回与微分方程*)

kiraDir = FileNameJoin[{exampleWorkspace, "kira"}];
deDir = FileNameJoin[{exampleWorkspace, "results", "dlogDE"}];
savedDEFile = FileNameJoin[{deDir, "manifest.wl"}];
packageSourcePath = currentPackagePath;
summaryFile = FileNameJoin[{exampleWorkspace, "results", "post_kira_summary.wl"}];

reductionData = DSKiraImport[kiraDir, context];
(* 已有 DE 是本脚本上一轮从同一 manifest 重算的正式结果；删除该目录即可强制重新调用 DSDE。 *)
deData = If[
   FileExistsQ[savedDEFile] &&
    FileDate[savedDEFile] >= FileDate[FileNameJoin[{kiraDir, "results", "Tuserweight", "kira_list.m"}]] &&
    FileDate[savedDEFile] >= FileDate[packageSourcePath],
   Join[Get[savedDEFile], <|"context" -> context|>],
   DSDE[reductionData, {ss11, P0}, OutputDirectory -> deDir]
   ];
scaleData = DSScaleCheck[
   deData,
   <|
    "relation" -> "PureMassiveBubble",
    "variables" -> {ss11, P0},
    "weights" -> {1, 1},
    "degrees" -> (pureMassiveBubbleUserMIScalingDegrees /. parameterProbeRules)
    |>
   ];


(* ::Chapter:: *)
(*闭环验收*)

(* postDerivative 阶段只把参数规则用于 DE 系数；masters 保持 basis 原始定义，逐项同序比较。 *)
expectedMasters = pureMassiveBubbleUserMIExpressions[[pureMassiveBubbleUserMIActiveIndices]];
expectedMasterIDs = Range[19];
expectedMasterTokens = userMI /@ expectedMasterIDs;
expectedBackendMasterTokens = Tuserweight /@ expectedMasterIDs;
sourceManifest = Lookup[reductionData, "sourceManifest", <||>];
userMIData = DSUserMI[reductionData];
probeSymbols = First /@ parameterProbeRules;
dePublicData = HoldComplete[
   Lookup[deData, "masters", {}],
   Lookup[deData, "matrices", <||>],
   Lookup[deData, "sources", <||>]
   ];
fixedParameterResiduals = Select[probeSymbols, ! FreeQ[dePublicData, #] &];

(* sourceManifest/backendReductionRules 有意保留可逆后端映射；用户结果面不得出现这些原子。 *)
backendLeakSymbols = DeleteDuplicates @ Cases[
    HoldComplete[
     Lookup[reductionData, "reductionRules", {}],
     Lookup[deData, "masters", {}],
     Lookup[deData, "matrices", <||>],
     Lookup[deData, "sources", <||>]
     ],
    symbol_Symbol /; StringMatchQ[SymbolName[Unevaluated[symbol]], ("dsc" ~~ ___) | "dsii"],
    Infinity
    ];

validationChecks = Lookup[Lookup[reductionData, "validationReport", <||>], "checks", <||>];
checks = <|
   "importStatus" -> (Lookup[reductionData, "status", "missing"] === "imported"),
   "importValidation" -> (Lookup[Lookup[reductionData, "validationReport", <||>], "status", "missing"] === "passed"),
   "activeMasterCount" -> (Length[Lookup[reductionData, "masters", {}]] === 19),
   "activeMasterIDs" -> (Lookup[reductionData, "masterIDs", {}] === expectedMasterIDs),
   "activeMasterTokens" -> (Lookup[reductionData, "masterTokens", {}] === expectedMasterTokens),
   "backendMasterTokens" -> (Lookup[reductionData, "backendMasterTokens", {}] === expectedBackendMasterTokens),
   "activeMasterOrder" -> (Lookup[reductionData, "masters", {}] === expectedMasters),
   "userMIMapping" -> TrueQ[Lookup[userMIData, "reversibleQ", False]],
   (* 新 basis 十九项全部 active、没有辅助关系；改为检查 19 个 active 关系均被选为 master。 *)
   "backendMasterIDsCoverActive" -> (Select[Lookup[reductionData, "backendMasterIDs", {}], 1 <= # <= 19 &] === Range[19]),
   "completeTargetCoverage" -> TrueQ[Lookup[validationChecks, "completeTargetCoverage", False]],
   "rhsContainsOnlyMasters" -> TrueQ[Lookup[validationChecks, "rhsContainsOnlyMasters", False]],
   "deStatus" -> (Lookup[deData, "status", "missing"] === "generated"),
   "deMasterOrder" -> (Lookup[deData, "masters", {}] === expectedMasters),
   "deVariables" -> (Lookup[deData, "variables", {}] === {ss11, P0}),
   "deDimensions" -> (And @@ (Dimensions[#] === {19, 19} & /@ Values[Lookup[deData, "matrices", <||>]])),
   "noResidualJ" -> FreeQ[Lookup[deData, "residualIntegrals", <||>], _J],
   "noResidualBackendTokens" -> FreeQ[Lookup[deData, "residualBackendTokens", <||>], Tuserweight[_Integer]],
   "noBackendCoefficientSymbols" -> (backendLeakSymbols === {}),
   "manifestParameterRules" -> (Lookup[sourceManifest, "userCoefficientRulesApplied", Missing["coefficientRules"]] === parameterProbeRules),
   "noFixedParameterResidual" -> (fixedParameterResiduals === {}),
   "scalingStatus" -> (Lookup[scaleData, "status", "missing"] === "passed"),
   (* 当前 package 的标度证书按 certificateScope 区分键名；本例的定点参数不含 Euler 变量，
      证书范围为 symbolic，故读取 symbolicMatrixRelation/symbolicSourceRelation。 *)
   "scalingMatrixResidual" -> TrueQ[Lookup[Lookup[scaleData, "checks", <||>], "symbolicMatrixRelation", False]],
   "scalingSourceResidual" -> TrueQ[Lookup[Lookup[scaleData, "checks", <||>], "symbolicSourceRelation", False]]
   |>;

failedChecks = Keys @ Select[checks, ! TrueQ[#] &];
summary = <|
   "status" -> If[failedChecks === {}, "passed", "failed"],
   "passed" -> Count[Values[checks], True],
   "total" -> Length[checks],
   "failedChecks" -> failedChecks,
   "backendMasterIDs" -> Lookup[reductionData, "backendMasterIDs", {}],
   "activeMasterIDs" -> Lookup[reductionData, "masterIDs", {}],
   "deVariables" -> Lookup[deData, "variables", {}],
   "parameterProbeSeed" -> parameterProbeSeed,
   "parameterProbeRules" -> parameterProbeRules,
   "fixedParameterResiduals" -> fixedParameterResiduals,
   "scalingDegrees" -> Lookup[scaleData, "degrees", {}]
   |>;

Print[summary];
Quiet[CreateDirectory[DirectoryName[summaryFile], CreateIntermediateDirectories -> True]];
Put[summary, summaryFile];
If[failedChecks =!= {}, Exit[1]];
