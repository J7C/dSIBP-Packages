(* ::Package:: *)

(* ::Chapter:: *)
(*上下文、消息、进度与公开 API 清单*)

If[! ValueQ[$dSIBPMessagesEnabled], $dSIBPMessagesEnabled = True];
If[! ValueQ[$dSIBPCurrentContext], $dSIBPCurrentContext = Missing["NotInitialized"]];

DSMessagesOn[] := ($dSIBPMessagesEnabled = True);
DSMessagesOff[] := ($dSIBPMessagesEnabled = False);
DSMessagesQ[] := TrueQ[$dSIBPMessagesEnabled];

dsMessagesEnabledQ[setting_: Automatic] := If[setting === Automatic, DSMessagesQ[], TrueQ[setting]];

(* 所有调用点直接提供逐句中英文本；此层只统一非字符串对象的显示，不伪造占位翻译。 *)
dsBilingualRuntimeText[text_] := ToString[text];


dsInfoPrint[text_, setting_: Automatic] := If[
   dsMessagesEnabledQ[setting], Print["[dSIBP] ", dsBilingualRuntimeText[text]]
   ];

dsWarningPrint[text_, setting_: Automatic] := If[
   dsMessagesEnabledQ[setting],
   If[TrueQ[$Notebooks],
    Print[Style["警告 / Warning: " <> dsBilingualRuntimeText[text], Darker[Red]]],
    Print["[dSIBP 警告 / Warning] ", dsBilingualRuntimeText[text]]
    ]
   ];

(* fatal error 不读取全局开关；即使用户关闭可选提醒也必须可见。 *)
dsErrorPrint[text_] := If[
   TrueQ[$Notebooks],
   Print[Style["错误 / Error: " <> dsBilingualRuntimeText[text], Red]],
   Print["[dSIBP 错误 / Error] ", dsBilingualRuntimeText[text]]
   ];

dsStageRun[label_String, expression_, setting_: Automatic] := Module[{result, elapsed},
   dsInfoPrint["开始 / Start: " <> label, setting];
   {elapsed, result} = AbsoluteTiming[expression];
   dsInfoPrint["完成 / Completed: " <> label <> " (" <> ToString[Round[elapsed, 0.001], InputForm] <> " s)", setting];
   result
   ];

(* Notebook 由 Monitor 在同一输出区域刷新；headless 只报告 25% 里程碑，避免逐项刷屏。 *)
dsProgressMap[label_String, items_List, function_, setting_: Automatic] := Module[
   {index = 0, total = Length[items], marks, result},
   If[total == 0, Return[{}]];
   If[! dsMessagesEnabledQ[setting], Return[function /@ items]];
   If[TrueQ[$Notebooks],
    Monitor[
     Map[(index++; function[#]) &, items],
     Row[{label, " ", Dynamic[index], "/", total, "  ", ProgressIndicator[Dynamic[index], {0, total}]}]
     ],
    marks = If[total < 8, {total}, DeleteDuplicates[Ceiling[total {1/4, 1/2, 3/4, 1}]]];
    result = Map[
      Function[item,
       index++;
       If[MemberQ[marks, index], dsInfoPrint[label <> " " <> ToString[index] <> "/" <> ToString[total], setting]];
       function[item]
       ],
      items
      ];
    result
    ]
   ];

dsContextQ[context_] := AssociationQ[context] && Lookup[context, "status", Missing["status"]] === "initialized" &&
   parsedTopologyQ[Lookup[context, "topology", Missing["topology"]]];

(* 失败初始化不得携带可被下游误读的部分能力；原始 declaration/kinematics audit
   仍保留各自的局部状态，公开 context capability 一律关闭。 *)
dsDisabledCapabilities[] := AssociationMap[
   False &,
   {
    "initializationUsableQ", "timeIBPUsableQ", "momentumIBPUsableQ",
    "derivativeUsableQ", "inverseKinematicsUsableQ", "backendExportUsableQ", "parityUsableQ"
    }
   ];


dsTopologyWithDisabledCapabilities[topo_Association] := Join[
   topo,
   <|"capabilities" -> dsDisabledCapabilities[]|>
   ];


dsTopologyWithDisabledCapabilities[topo_] := topo;


dsFailedInitializationData[reason_String, data_Association : <||>] := Module[{result},
   result = Join[data, <|
      "status" -> "failed",
      "reason" -> reason,
      "capabilities" -> dsDisabledCapabilities[]
      |>];
   If[KeyExistsQ[result, "topologyData"],
    result = Join[result, <|
       "topologyData" -> dsTopologyWithDisabledCapabilities[result["topologyData"]]
       |>]
    ];
   result
   ];


dsResolveContext[Automatic] := If[dsContextQ[$dSIBPCurrentContext], $dSIBPCurrentContext, Missing["NotInitialized"]];
dsResolveContext[context_Association] := If[dsContextQ[context], context, Missing["InvalidContext"]];
dsResolveContext[_] := Missing["InvalidContext"];

dsContextCapabilities[context_Association] := Lookup[
   context,
   "capabilities",
   Lookup[Lookup[context, "topology", <||>], "capabilities", <||>]
   ];

dsTopologyCapabilityQ[topo_Association, capability_String] := TrueQ[
   Lookup[Lookup[topo, "capabilities", <||>], capability, False]
   ];

dsContextCapabilityQ[context_Association, capability_String] := TrueQ[
   Lookup[dsContextCapabilities[context], capability, False]
   ];

dsContextSummary[context_Association] := <|
   "packageVersion" -> Lookup[context, "packageVersion", Missing["packageVersion"]],
   "inputHash" -> Lookup[context, "inputHash", Missing["inputHash"]],
   "caseName" -> Lookup[context, "caseName", Missing["caseName"]],
   "sectorKeys" -> Lookup[Lookup[context, "sectors", {}], "sectorKey", {}],
   "capabilities" -> dsContextCapabilities[context],
   "loopTreeProjectionConvention" -> Lookup[context, "loopTreeProjectionConvention", <||>]
   |>;


(* 公开函数清单是手册汇总表和成品 example 覆盖检查的共同来源。
   函数按用户工作流分组；符号 Head 与 option 名不混入函数覆盖计数。 *)
dsPublicAPISections[] := <|
   "initialization" -> {
     "DSInit", "DSInfo", "DSKinematics", "DSParameterNotation", "DSRedefineParameters"
     },
   "messages" -> {"DSMessagesOn", "DSMessagesOff", "DSMessagesQ"},
   "atomicOperations" -> {
     "dtau", "dqq", "dqk", "ds", "rep2innerform", "rep2outform", "rep2Integrand",
     "symmetry", "repSymmetry0"
     },
    "loopWorkflow" -> {"DSSeeds", "DSAllSeeds", "DSSeedGroups", "DSSeedGroupMetadata", "DSMetaSeedRange", "DSGenerateIBP", "DSLinear", "DSReorderIntegrals", "DSUserMI", "DSKiraPlan", "DSKiraExport", "DSKiraImport", "DSDE", "DSScaleCheck"},
   "pureTimeWorkflow" -> {
     "DSTreeSeeds", "repIterative", "DSTreeNaiveIBP", "DSTreeNaiveDE", "DSTreeDLogDE"
     },
   "introspection" -> {"DSPublicAPI"}
   |>;


dsPublicAPIOptions[] := <|
   "DSInit" -> Options[DSInit],
   "DSSeeds" -> Options[DSSeeds],
   "DSLinear" -> Options[DSLinear],
   "DSGenerateIBP" -> Options[DSGenerateIBP],
   "DSKiraPlan" -> Options[DSKiraPlan],
   "DSKiraExport" -> Options[DSKiraExport],
   "DSKiraImport" -> Options[DSKiraImport],
   "DSDE" -> Options[DSDE],
   "DSScaleCheck" -> Options[DSScaleCheck],
   "DSTreeNaiveIBP" -> Options[DSTreeNaiveIBP],
   "DSTreeNaiveDE" -> Options[DSTreeNaiveDE],
   "DSTreeDLogDE" -> Options[DSTreeDLogDE],
   "repIterative" -> Options[repIterative]
   |>;


DSPublicAPI[] := Module[{sections = dsPublicAPISections[]},
   <|
    "version" -> $dSIBPVersion,
    "sections" -> sections,
    "functions" -> DeleteDuplicates@Flatten[Values[sections]],
    "options" -> dsPublicAPIOptions[]
    |>
   ];
