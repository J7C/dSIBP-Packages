(* ::Package:: *)
(* pure massive bubble：固定 -- branch/parity，从 topology 到 formal Kira 输入。
   package 只生成和序列化关系，不启动 reduction；Kira 工作树必须位于包所在 Git 仓库或交付目录之外。 *)

(* ::Chapter:: *)
(*标准 package 加载*)

exampleDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{exampleDir, "..", "load_current_package.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{exampleDir, "dlog_basis.wl"}]];
Get[FileNameJoin[{exampleDir, "reference_user_mi_basis.wl"}]];
Get[FileNameJoin[{exampleDir, "family_conventions.wl"}]];
bundledInputSummary = Get[FileNameJoin[{exampleDir, "kira_input_summary.wl"}]];
bundledResultSummary = Get[FileNameJoin[{exampleDir, "kira_result_summary.wl"}]];

externalWorkspaceRoot = Environment["DSIBP_KIRA_WORKSPACE"];
If[! StringQ[externalWorkspaceRoot] || StringTrim[externalWorkspaceRoot] === "",
  Print["请先用 DSIBP_KIRA_WORKSPACE 指定仓库外的 Kira 工作目录，再运行本例。"];
  Exit[2]
];
externalWorkspaceRoot = ExpandFileName[externalWorkspaceRoot];
(* 包可放在任意位置：以上层最近的 .git 为仓库根；找不到 Git 仓库时以交付目录为界。 *)
workspaceGuardCandidates = Select[
   NestList[DirectoryName, exampleDir, 8],
   FileExistsQ[FileNameJoin[{#, ".git"}]] &
   ];
workspaceGuardRoot = ExpandFileName@If[workspaceGuardCandidates === {},
   FileNameJoin[{exampleDir, "..", ".."}],
   First[workspaceGuardCandidates]
   ];
(* ExpandFileName 保留尾分隔符，须先去掉再做前缀判定。 *)
workspaceGuardRoot = StringReplace[workspaceGuardRoot, RegularExpression["[\\\\/]+$"] -> ""];
If[$OperatingSystem === "Windows",
  workspaceGuardKey = ToLowerCase[workspaceGuardRoot];
  externalWorkspaceKey = ToLowerCase[externalWorkspaceRoot],
  workspaceGuardKey = workspaceGuardRoot;
  externalWorkspaceKey = externalWorkspaceRoot
  ];
If[
  externalWorkspaceKey === workspaceGuardKey ||
    StringStartsQ[externalWorkspaceKey, workspaceGuardKey <> $PathnameSeparator],
  Print["Kira 工作目录必须位于包所在 Git 仓库或交付目录之外；请修改 DSIBP_KIRA_WORKSPACE。"];
  Exit[2]
];
exampleWorkspace = FileNameJoin[{externalWorkspaceRoot, "pure_massive_bubble"}];
If[! DirectoryQ[exampleWorkspace],
  CreateDirectory[exampleWorkspace, CreateIntermediateDirectories -> True]
];


(* 固定随机种子只记录参数点的来源；规则本身冻结，Kira 与回读端必须逐项复用。 *)
parameterProbeSeed = 20260722;
parameterProbeRules = {
   dim -> 37/11, nu -> 7/13, etaNu -> 23/17,
   analyticRegulator -> 0
   };

(* ::Chapter:: *)
(*详细物理输入*)

(* 两个顶点均在 - branch；这与 reference 的 Vpm=0 convention 对齐。 *)
(* 两条 massive h 圈传播子依次取 q 与 q-k；effectiveLoopExternalMomenta 只含实际进入 loop 动量的独立向量 k。 *)
(* ss11=Sqrt[k.k]=ks 与 P0 是 ds 的独立变量；P_pkg=P0，reference P1=P2=-P0。 *)
(* J 只保存整数指标；a0=2 nu、b0=-2 nu 留在 metadata，并在 shrink/tree 投影时进入完整物理幂次。 *)
caseInput = <|
   "name" -> "pureMassiveBubbleMinusMinus",
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

(* ::Chapter:: *)
(*缺省选项与本例覆盖*)

(* 缺省：WriteInitializationFiles->False；GenerateDerivativeMetadata->False；OverwriteInitialization->False。 *)
initOptions = {
   WriteInitializationFiles -> True,
   InitializationDirectory -> FileNameJoin[{exampleWorkspace, "init"}],
   GenerateDerivativeMetadata -> True,
   (* 本例的 Kira 输入、结果摘要与 dlog 基准是 025 冻结基线；Negative 逐字节复现该口径。
      新计算请用缺省 Positive 记法。 *)
   NuConvention -> "Negative",
   (* 本闭环例会随输入配置刷新同目录 metadata；缺省仍为 False。 *)
   OverwriteInitialization -> True
   };

(* 缺省：LinearSystemMode->"symbolic"；CoefficientRules->Automatic；KiraOrdering->Automatic。 *)
linearOptions = {
   LinearSystemMode -> "symbolic",
   KiraOrdering -> Automatic
   };

(* 缺省：OutputDirectory->None；KiraJobOptions->Automatic；package 只写文件，不运行 Kira。 *)
kiraJobOptions = <|
     "RunInitiate" -> True,
     "RunFirefly" -> True,
     "WriteKira2MathJob" -> True,
     "WriteRunScript" -> False
     |>;

(* ::Chapter:: *)
(*初始化、IBP 与 Kira 输入*)

context = DSInit[caseInput, Sequence @@ initOptions];
If[Lookup[context, "status", "failed"] =!= "initialized",
  Print["Pure massive bubble 初始化失败；请先修正上方报告的输入或正规化问题，后续 Kira 步骤不会执行。"];
  Exit[1]
];
DSInfo[context]

seedData = DSSeeds[context];
allSeeds = DSAllSeeds[seedData];
seedGroups = DSSeedGroups[seedData];
seedGroupMetadata = DSSeedGroupMetadata[seedData];
seedRangeMetadata = DSMetaSeedRange[seedGroups, referenceSeedIndices];
(* 输入描述 top 目标积分包络；程序逐组反推 seed 点域、先筛 parity，再代入数值点。 *)
generatedIBP = DSGenerateIBP[allSeeds, Sequence @@ referenceTopTargetEnvelope];
linearData = DSLinear[
   generatedIBP,
   context,
   LinearSystemMode -> "numeric",
   CoefficientRules -> parameterProbeRules,
   KiraOrdering -> Automatic
   ];
linearDataWithUserMI = DSUserMI[
   linearData,
   pureMassiveBubbleUserMIExpressions,
   <|
    "names" -> pureMassiveBubbleUserMINames,
    "activeIndices" -> pureMassiveBubbleUserMIActiveIndices,
    "derivativeVariables" -> pureMassiveBubbleUserMIDerivativeVariables,
    "scalingDegrees" -> pureMassiveBubbleUserMIScalingDegrees
    |>
   ];
userMIData = DSUserMI[linearDataWithUserMI];
formalPlan = DSKiraPlan[
   linearDataWithUserMI,
   <|
    "stage" -> "formal",
    "activeBasis" -> Automatic,
     "numericStage" -> "postDerivative",
     "coefficientRules" -> parameterProbeRules,
     "outputDirectory" -> FileNameJoin[{exampleWorkspace, "kira"}],
    "jobOptions" -> kiraJobOptions
    |>
   ];
kiraExport = DSKiraExport[formalPlan];
exportSyntaxReport = Lookup[
   kiraExport,
   "backendCoefficientSyntaxReport",
   Lookup[Lookup[kiraExport, "kiraInput", <||>], "backendCoefficientSyntaxReport", <||>]
   ];

(* ::Chapter:: *)
(*完整 Kira 结果取回、DE 与标度关系*)

kiraDir = FileNameJoin[{exampleWorkspace, "kira"}];
requiredKiraResults = {
   FileNameJoin[{kiraDir, "kira.log"}],
   FileNameJoin[{kiraDir, "results", "Tuserweight", "kira_list.m"}],
   FileNameJoin[{kiraDir, "results", "Tuserweight", "masters"}]
   };
exportReadyQ = Lookup[kiraExport, "status", "missing"] === "ready";
reductionReadyQ = exportReadyQ && And @@ (FileExistsQ /@ requiredKiraResults);
reductionReadyQ = reductionReadyQ && Quiet[Check[
     FileDate[FileNameJoin[{kiraDir, "kira.log"}]] >= FileDate[FileNameJoin[{kiraDir, "dsibp-export-manifest.wl"}]],
     False
     ]];

closedLoopResult = If[
   reductionReadyQ,
   reductionData = DSKiraImport[kiraDir, context];
   deData = DSDE[
     reductionData,
     {ss11, P0},
      OutputDirectory -> FileNameJoin[{exampleWorkspace, "results", "dlogDE"}]
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
   <|"status" -> scaleData["status"], "reduction" -> reductionData, "de" -> deData, "scaling" -> scaleData|>,
   <|"status" -> If[exportReadyQ, "awaitingExternalKira", "exportNotReady"],
    "requiredFiles" -> requiredKiraResults|>
   ];

closedLoopSummary = If[
   reductionReadyQ,
   <|
    "status" -> Lookup[closedLoopResult, "status", "missing"],
    "exportStatus" -> Lookup[kiraExport, "status", "missing"],
    "exportReason" -> Lookup[kiraExport, "reason", None],
    "parameterProbeSeed" -> parameterProbeSeed,
     "parameterProbeRules" -> parameterProbeRules,
     "coefficientRulesApplied" -> Lookup[Lookup[kiraExport, "dSIBPExportManifest", <||>], "userCoefficientRulesApplied", {}],
    "exportSyntaxBadStrings" -> Take[DeleteDuplicates[Lookup[exportSyntaxReport, "badStrings", {}]], UpTo[8]],
    "reductionStatus" -> Lookup[reductionData, "status", "missing"],
    "masterCount" -> Length[Lookup[reductionData, "masters", {}]],
    "deStatus" -> Lookup[deData, "status", "missing"],
    "deResidualCounts" -> Map[Length, Lookup[deData, "residualIntegrals", <||>]],
    "scalingStatus" -> Lookup[scaleData, "status", "missing"],
    "scalingReason" -> Lookup[scaleData, "reason", None]
    |>,
   <|
    "status" -> If[exportReadyQ, "awaitingExternalKira", "exportNotReady"],
    "exportStatus" -> Lookup[kiraExport, "status", "missing"],
    "exportReason" -> Lookup[kiraExport, "reason", None],
    "parameterProbeSeed" -> parameterProbeSeed,
    "parameterProbeRules" -> parameterProbeRules,
    "exportSyntaxBadStrings" -> Take[DeleteDuplicates[Lookup[exportSyntaxReport, "badStrings", {}]], UpTo[8]],
    "requiredFiles" -> requiredKiraResults
    |>
   ];
Print["pure massive bubble closed-loop summary: ", closedLoopSummary];

If[! exportReadyQ, Exit[1]];

closedLoopResult
