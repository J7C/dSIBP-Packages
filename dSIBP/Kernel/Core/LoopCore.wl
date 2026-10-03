(* ::Package:: *)
(* 本文件是 dS IBP package 的通用生成器骨架。
   目标是先把 topology-driven 的结构层做实：拓扑解析、传播子 metadata、统一 J 指标包、
   离散态枚举、完整圈动量 IBP 生成元列表、标量积/ISP 覆盖性验证。
   当前文件生成 momentum seed、time-core seed 与受保护的自动 shrink-sector seed；EOM、有方向的 per-line massless 单 n 规则与 massive/massless theta boundary shrink 已作为 seed 门禁接入，并提供 backend-neutral linear-system 与 Kira serializer。
   公开原子 API 为 dtau/dqq/dqk 与 rep2innerform/rep2outform/rep2Integrand；函数系统编译层提供
   P,Q,T,W，缺省 h 取 T=IdentityMatrix[2]、WT=W_h，Hankel 作为独立 preset。massive ++ Wronskian
   theta-boundary 的 Vpm 符号、按每个 J 的 shrunk packs 重建的目标顶点映射、跨 sector 的 massless
   coincident n=1 canonical、共同-theta odd-subset contact 与 contact 可达 sector、massless shrink
   幂次修正和 massive coincidence canonical 都在本模块内实现。pure time-IBP/tree 表示、loop-time
   投影、通用迭代矩阵与直接 dlog DE 由本模块与后面的 Tree 章节共同提供。
   性能原则：默认只定义函数和示例输入，不自动运行检查；验证必须是 seed/metadata 层或代数赋值后的小检查。 *)

(* 加载主线脚本不得清空用户 Global` 上下文；开发时需要彻底重载应使用新 kernel。 *)


(* ::Chapter:: *)
(*环境与通用工具*)

(* 本章只解析脚本目录并定义小工具函数。加载 package 不改变用户当前工作目录；输出路径由调用者显式决定。 *)

baseDir = If[StringQ[$InputFileName] && $InputFileName =!= "",
   DirectoryName[$InputFileName],
   If[TrueQ[$Notebooks],
    With[{nd = Quiet[NotebookDirectory[]]},
     If[StringQ[nd] && nd =!= $Failed, nd, Directory[]]
     ],
    Directory[]
    ]
   ];


(* 只做结构零判断；解析化简应放到单独的小规模 specialized check 中。 *)
zeroQ[expr_] := TrueQ[expr === 0];


(* 用户口统一用 sp 表示 scalar product；内部仍用 qq/qk/kk 编号坐标做线性代数。 *)
SetAttributes[sp, Orderless];
SetAttributes[qq, Orderless];
SetAttributes[kk, Orderless];


(* 对称标量积记号。qq/kk 只保留上三角索引，避免同一个标量积出现两种名字。 *)
qqSym[i_, j_] := If[i <= j, qq[i, j], qq[j, i]];
kkSym[i_, j_] := If[i <= j, kk[i, j], kk[j, i]];


(* ::Chapter:: *)
(*拓扑输入规范化*)

(* 本章把用户输入的 vertices/lines/extLegs/loopMomenta/ispData 规范化为 Association。
   顶点和传播子只接受 Association schema；额外键不进入内部 topology。 *)


(* 顶点相位的唯一公开 authority 是 vertexType；内部保留规范化后的同名字段。 *)
normalizeVertex[vertex_Association] := KeyTake[
   vertex,
   {"id", "vertexType", "externalLegEnergy"}
   ];


(* line 的 SK、pack 和 contact state 都由端点与 sector producer 派生，不能从用户输入读取。 *)
normalizeLine[line_Association] := Module[{massType, allowed},
   massType = Lookup[line, "massType", Missing["massType"]];
   allowed = If[
     massType === "massive",
     {"id", "massType", "endpoints", "momentum", "nu", "functionSystem"},
     {"id", "massType", "endpoints", "momentum"}
     ];
   KeyTake[line, allowed]
   ];


(* ISP 可用 {name, expr, range} 或 Association。其坐标由不可约 numerator 标量积定义，
   零点固定为 0；正幂是 numerator，用户显式选择的负幂作为额外 denominator 保留。 *)
normalizeISP[isp_Association] := Join[isp, <|"zeroPoint" -> 0|>];
normalizeISP[{name_, expr_, range_}] := <|
   "name" -> name,
   "expr" -> expr,
   "range" -> range,
   "zeroPoint" -> 0
   |>;


requiredCaseInputKeys[] := {"vertices", "lines", "loopMomenta"};


optionalCaseInputKeys[] := {
   "name", "extLegs", "ibpMode",
   "loopExternalMomenta", "independentExternalMomenta",
   "kinematicRules",
   "ispData", "zeroPointRules", "symmetryRules", "parityConstraints"
   };


vertexEntryAssociationMissingKeys[entry_Association] := Complement[
   {"id", "vertexType", "externalLegEnergy"},
   Keys[entry]
   ];
vertexEntryAssociationMissingKeys[_] := {};


validVertexEntryQ[entry_] := AssociationQ[entry];


validLineEntryShapeQ[entry_] := AssociationQ[entry];


lineEntryAssociationMissingKeys[entry_Association] := Module[{required},
   required = Join[
     {"id", "massType", "endpoints", "momentum"},
     If[Lookup[entry, "massType", Missing["massType"]] === "massive", {"nu"}, {}]
     ];
   Complement[required, Keys[entry]]
   ];
lineEntryAssociationMissingKeys[_] := {};


lineEntryEndpointValue[entry_Association] := Lookup[entry, "endpoints", Missing["NoEndpoints"]];
lineEntryEndpointValue[_] := Missing["NoEndpoints"];


validEndpointPairQ[endpoints_] := ListQ[endpoints] && Length[endpoints] == 2;


validISPEntryShapeQ[entry_] := AssociationQ[entry] || MatchQ[entry, {_, _, _}];


ispEntryAssociationMissingKeys[entry_Association] := Complement[{"name", "expr"}, Keys[entry]];
ispEntryAssociationMissingKeys[_] := {};


validIndexRangeSpecQ[spec_] := IntegerQ[spec] || (ListQ[spec] && Length[spec] > 0 && And @@ (IntegerQ /@ spec));


validDiscreteReplacementRuleQ[rule_] := MatchQ[Unevaluated[rule], _Rule | _RuleDelayed];


caseInputMalformedIssues[case_Association] := Module[
   {issues = {}, vertices, lineInput, loopMomenta, loopExternalMomenta,
    independentExternalMomenta, ispData, badVertexPositions, vertexMissingKeyData,
    badLineShapePositions, lineMissingKeyData, badEndpointData, badISPShapePositions, ispMissingKeyData,
    symmetryRules, badSymmetryRulePositions},
   If[KeyExistsQ[case, "vertices"],
    vertices = case["vertices"];
    If[! ListQ[vertices],
     AppendTo[issues, <|"severity" -> "error", "code" -> "malformedVertices", "reason" -> "vertices must be a list of Associations"|>],
     badVertexPositions = Flatten @ Position[vertices, entry_ /; ! validVertexEntryQ[entry], {1}, Heads -> False];
     If[badVertexPositions =!= {},
      AppendTo[issues, <|"severity" -> "error", "code" -> "malformedVertices", "badPositions" -> badVertexPositions|>]
      ];
     vertexMissingKeyData = DeleteCases[
       MapIndexed[
        If[AssociationQ[#1] && vertexEntryAssociationMissingKeys[#1] =!= {},
          <|"vertexPosition" -> First[#2], "missingKeys" -> vertexEntryAssociationMissingKeys[#1]|>,
          Nothing
          ] &,
        vertices
        ],
       Nothing
       ];
     If[vertexMissingKeyData =!= {},
      AppendTo[issues, <|"severity" -> "error", "code" -> "verticesMissingRequiredKeys", "vertices" -> vertexMissingKeyData|>]
      ]
     ]
    ];
   If[KeyExistsQ[case, "lines"],
    lineInput = case["lines"];
    If[! ListQ[lineInput],
     AppendTo[issues, <|"severity" -> "error", "code" -> "malformedLines", "reason" -> "lines must be a list of Associations"|>],
     badLineShapePositions = Flatten @ Position[lineInput, entry_ /; ! validLineEntryShapeQ[entry], {1}, Heads -> False];
     If[badLineShapePositions =!= {},
      AppendTo[issues, <|"severity" -> "error", "code" -> "malformedLines", "badPositions" -> badLineShapePositions|>]
      ];
     lineMissingKeyData = DeleteCases[
       MapIndexed[
        If[AssociationQ[#1] && lineEntryAssociationMissingKeys[#1] =!= {},
          <|"linePosition" -> First[#2], "missingKeys" -> lineEntryAssociationMissingKeys[#1]|>,
          Nothing
          ] &,
        lineInput
        ],
       Nothing
       ];
     If[lineMissingKeyData =!= {},
      AppendTo[issues, <|"severity" -> "error", "code" -> "linesMissingRequiredKeys", "lines" -> lineMissingKeyData|>]
      ];
     badEndpointData = DeleteCases[
       MapIndexed[
        If[validLineEntryShapeQ[#1] && lineEntryAssociationMissingKeys[#1] === {} && ! validEndpointPairQ[lineEntryEndpointValue[#1]],
          <|"linePosition" -> First[#2], "endpoints" -> lineEntryEndpointValue[#1]|>,
          Nothing
          ] &,
        lineInput
        ],
       Nothing
       ];
     If[badEndpointData =!= {},
      AppendTo[issues, <|"severity" -> "error", "code" -> "malformedLineEndpoints", "lines" -> badEndpointData|>]
      ]
     ]
    ];
   If[KeyExistsQ[case, "loopMomenta"],
    loopMomenta = case["loopMomenta"];
    If[! ListQ[loopMomenta],
     AppendTo[issues, <|"severity" -> "error", "code" -> "malformedLoopMomenta", "reason" -> "loopMomenta must be a list"|>]
     ]
    ];
   If[KeyExistsQ[case, "loopExternalMomenta"],
    loopExternalMomenta = case["loopExternalMomenta"];
    If[! ListQ[loopExternalMomenta],
     AppendTo[issues, <|"severity" -> "error", "code" -> "malformedLoopExternalMomenta", "reason" -> "loopExternalMomenta must be a list"|>]
     ]
    ];
   If[KeyExistsQ[case, "independentExternalMomenta"],
    independentExternalMomenta = case["independentExternalMomenta"];
    If[! ListQ[independentExternalMomenta],
     AppendTo[issues, <|"severity" -> "error", "code" -> "malformedIndependentExternalMomenta", "reason" -> "independentExternalMomenta must be a list"|>]
     ]
    ];
   If[KeyExistsQ[case, "ibpMode"] && ! MemberQ[{Automatic, "full", "timeOnly"}, case["ibpMode"]],
    AppendTo[issues, <|"severity" -> "error", "code" -> "malformedIBPMode", "value" -> case["ibpMode"], "allowed" -> {"full", "timeOnly"}|>]
    ];
   If[KeyExistsQ[case, "symmetryRules"],
    symmetryRules = case["symmetryRules"];
    If[! ListQ[symmetryRules],
     AppendTo[issues, <|"severity" -> "error", "code" -> "malformedSymmetryRules", "reason" -> "symmetryRules must be a list of Rule or RuleDelayed entries"|>],
     badSymmetryRulePositions = Flatten @ Position[
        symmetryRules,
        rule_ /; ! validDiscreteReplacementRuleQ[rule],
        {1},
        Heads -> False
        ];
     If[badSymmetryRulePositions =!= {},
      AppendTo[issues, <|"severity" -> "error", "code" -> "malformedSymmetryRules", "badPositions" -> badSymmetryRulePositions|>]
      ]
     ]
    ];
   If[KeyExistsQ[case, "ispData"],
    ispData = case["ispData"];
    If[! ListQ[ispData],
     AppendTo[issues, <|"severity" -> "error", "code" -> "malformedISPData", "reason" -> "ispData must be a list of ISP associations or {name,expr,range} entries"|>],
     badISPShapePositions = Flatten @ Position[ispData, entry_ /; ! validISPEntryShapeQ[entry], {1}, Heads -> False];
     If[badISPShapePositions =!= {},
      AppendTo[issues, <|"severity" -> "error", "code" -> "malformedISPData", "badPositions" -> badISPShapePositions|>]
      ];
     ispMissingKeyData = DeleteCases[
       MapIndexed[
        If[AssociationQ[#1] && ispEntryAssociationMissingKeys[#1] =!= {},
          <|"ispPosition" -> First[#2], "missingKeys" -> ispEntryAssociationMissingKeys[#1]|>,
          Nothing
          ] &,
        ispData
        ],
       Nothing
       ];
     If[ispMissingKeyData =!= {},
      AppendTo[issues, <|"severity" -> "error", "code" -> "ispDataMissingRequiredKeys", "isps" -> ispMissingKeyData|>]
      ]
     ]
     ];
    issues
    ];


(* raw case 入口的轻量 preflight；避免缺必需字段时 parser 先抛底层 Part 消息。 *)
caseInputRequirementReport[case_Association] := Module[
   {keys = Keys[case], required = requiredCaseInputKeys[], optional = optionalCaseInputKeys[],
    malformedIssues},
   malformedIssues = caseInputMalformedIssues[case];
   <|
    "providedKeys" -> keys,
    "requiredKeys" -> required,
    "optionalKeys" -> optional,
    "missingRequiredKeys" -> Complement[required, keys],
    "unknownKeys" -> Complement[keys, Join[required, optional]],
    "malformedInputIssues" -> malformedIssues,
    "completeRequiredKeysQ" -> TrueQ[Complement[required, keys] === {}],
    "inputPreflightPassQ" -> TrueQ[
      Complement[required, keys] === {} &&
       malformedIssues === {}
      ]
    |>
   ];


caseInputErrorReport[case_Association] := Module[{report = caseInputRequirementReport[case]},
   With[{issues = Join[
       If[report["missingRequiredKeys"] === {}, {}, {<|"severity" -> "error", "code" -> "missingRequiredCaseKeys", "missingRequiredKeys" -> report["missingRequiredKeys"]|>}],
        report["malformedInputIssues"]
       ]},
   <|
    "status" -> If[issues === {}, "ok", "issues"],
    "errorCount" -> Length[issues],
    "warningCount" -> 0,
    "pendingCount" -> 0,
    "pendingFeatures" -> {},
    "inputRequirementReport" -> report,
    "issues" -> issues
    |>
    ]
   ];


caseInputMissingRequiredKeysQ[case_Association] := ! TrueQ[caseInputRequirementReport[case]["completeRequiredKeysQ"]];


caseInputPreflightErrorQ[case_Association] := ! TrueQ[caseInputRequirementReport[case]["inputPreflightPassQ"]];


(* 由端点顶点的 +/- 标记推断 SK 类型。 *)
inferSKType[endpoints_, vertexSignAssoc_] := StringJoin[
   ToString[Lookup[vertexSignAssoc, endpoints[[1]], "+"]],
   ToString[Lookup[vertexSignAssoc, endpoints[[2]], "+"]]
   ];


(* packType 是后续统一 J 表示的核心分派键。massless 统一走双 theta 合并路线。 *)
inferPackType[massType_, skType_, state_] := Which[
   state === "shrunk", "shrunk",
   massType === "massive" && MemberQ[{"+-", "-+"}, skType], "massiveCross",
   massType === "massless" && MemberQ[{"++", "--"}, skType], "masslessFull",
   massType === "massless" && MemberQ[{"+-", "-+"}, skType], "masslessCross",
   True, "massiveFull"
   ];


(* ::Section:: *)
(*特殊函数系统编译*)

functionSystemZeroQ[expr_] := TrueQ[Quiet[FullSimplify[Together[expr] == 0]]];


functionSystemAdditiveTerms[expr_] := Module[{expanded = Expand[expr]},
   Which[
    functionSystemZeroQ[expanded], {},
    Head[expanded] === Plus, List @@ expanded,
    True, {expanded}
    ]
   ];


functionSystemMonomialData[term_, variable_] := Module[{factors, power, coefficient},
   factors = If[Head[term] === Times, List @@ term, {term}];
   power = Total[Replace[
      factors,
      {
       factor_ /; SameQ[factor, variable] :> 1,
       HoldPattern[Power[base_, exponent_]] /; SameQ[base, variable] :> exponent,
       _ -> 0
       },
      {1}
      ]];
   coefficient = Quiet[FullSimplify[term/variable^power]];
   If[FreeQ[coefficient, variable],
    <|"coefficient" -> coefficient, "xPower" -> power|>,
    Missing["NotFiniteLaurentMonomial", term]
    ]
   ];


compileDerivativeTerms[at_, variable_] := Catch[
   Module[{entries = {}, data, grouped},
    Do[
     Do[
      Do[
       data = functionSystemMonomialData[term, variable];
       If[Head[data] === Missing || ! IntegerQ[data["xPower"]], Throw[$Failed, derivativeCompileFailure]];
       AppendTo[entries, Join[
         <|"sourceState" -> i - 1, "targetState" -> j - 1|>,
         data
         ]],
       {term, functionSystemAdditiveTerms[at[[i, j]]]}
       ],
      {j, 2}
      ],
     {i, 2}
     ];
    grouped = GroupBy[entries, {#["sourceState"], #["targetState"], #["xPower"]} &];
    KeyValueMap[
     <|
       "sourceState" -> #1[[1]], "targetState" -> #1[[2]], "xPower" -> #1[[3]],
       "coefficient" -> Quiet[FullSimplify[Total[Lookup[#2, "coefficient"]]]]
       |> &,
     grouped
     ]
    ],
   derivativeCompileFailure
   ];


compileShrinkTerms[wt_, variable_, spec_Association] := Module[
   {rawTerms, data, powers, integerPowers, bShifts, zeroPointShifts, commonZeroPointShift},
   rawTerms = functionSystemAdditiveTerms[wt];
   If[rawTerms === {}, Return[$Failed]];
   data = functionSystemMonomialData[#, variable] & /@ rawTerms;
   If[AnyTrue[data, Head[#] === Missing &], Return[$Failed]];
   powers = Lookup[data, "xPower"];
   If[KeyExistsQ[spec, "shrinkBShift"] || KeyExistsQ[spec, "shrinkZeroPointShift"],
    If[Length[data] =!= 1 || ! KeyExistsQ[spec, "shrinkBShift"] || ! KeyExistsQ[spec, "shrinkZeroPointShift"], Return[$Failed]];
    bShifts = {spec["shrinkBShift"]};
    zeroPointShifts = {spec["shrinkZeroPointShift"]},
    integerPowers = Quiet[FullSimplify[# /. Thread[Variables[{#}] -> 0]]] & /@ powers;
    If[AnyTrue[integerPowers, ! IntegerQ[#] &], Return[$Failed]];
    bShifts = -integerPowers;
    zeroPointShifts = MapThread[Quiet[FullSimplify[-#1 - #2]] &, {powers, bShifts}]
    ];
   If[AnyTrue[bShifts, ! IntegerQ[#] &], Return[$Failed]];
   commonZeroPointShift = First[zeroPointShifts];
   If[AnyTrue[Rest[zeroPointShifts], ! functionSystemZeroQ[# - commonZeroPointShift] &], Return[$Failed]];
   <|
    "terms" -> MapThread[
      Join[#1, <|"coefficient" -> -#1["coefficient"], "bShift" -> #2, "zeroPointShift" -> #3|>] &,
      {data, bShifts, zeroPointShifts}
      ],
    "zeroPointShift" -> commonZeroPointShift
    |>
   ];


functionSystemPreset[preset_String, line_Association] := Module[
   {variable = x, nuValue = Lookup[line, "nu", nu], prefactor},
   prefactor = (4 I/Pi) Exp[Pi Im[nuValue]];
   Switch[preset,
    "h",
    <|
     "preset" -> "h", "variable" -> variable,
     "P" -> (2 nuValue + 1)/variable, "Q" -> 1,
     "T" -> IdentityMatrix[2],
     "W" -> -prefactor variable^(-2 nuValue - 1), "WT" -> Automatic,
     "shrinkBShift" -> 1, "shrinkZeroPointShift" -> 2 nuValue
     |>,
    "H",
    <|
     "preset" -> "H", "variable" -> variable,
     "P" -> 1/variable, "Q" -> 1 - nuValue^2/variable^2,
     "T" -> IdentityMatrix[2],
     "W" -> -prefactor/variable, "WT" -> Automatic,
     "shrinkBShift" -> 1, "shrinkZeroPointShift" -> 0
     |>,
    _,
    Missing["UnknownFunctionSystemPreset", preset]
    ]
   ];


normalizeLineFunctionSystem[line_Association] := Module[{explicit},
   explicit = Lookup[line, "functionSystem", Automatic];
   Which[
     AssociationQ[explicit], explicit,
     MemberQ[{"h", "H"}, explicit], functionSystemPreset[explicit, line],
     explicit =!= Automatic, Missing["MalformedFunctionSystem", explicit],
     True, functionSystemPreset["h", line]
     ]
   ];


compileFunctionSystem[spec_Association] := Module[
   {required = {"P", "Q", "W"}, missing, variable, p, q, t, w, explicitWT, a0, at, wt,
    detT, derivativeTerms, shrinkData, issues = {}},
   missing = Complement[required, Keys[spec]];
   If[missing =!= {}, Return[<|"status" -> "invalid", "issues" -> {<|"code" -> "missingFunctionSystemKeys", "missingKeys" -> missing|>}|>]];
   variable = Lookup[spec, "variable", x];
   p = spec["P"];
   q = spec["Q"];
   t = Replace[Lookup[spec, "T", Automatic], Automatic -> IdentityMatrix[2]];
   w = spec["W"];
   explicitWT = Lookup[spec, "WT", Automatic];
   If[! MatrixQ[t] || Dimensions[t] =!= {2, 2},
    Return[<|"status" -> "invalid", "issues" -> {<|"code" -> "functionSystemTNot2x2", "T" -> t|>}|>]
    ];
   detT = Quiet[FullSimplify[Det[t]]];
   If[functionSystemZeroQ[detT],
    Return[<|"status" -> "invalid", "issues" -> {<|"code" -> "functionSystemTSingular", "T" -> t|>}|>]
    ];
   a0 = {{0, 1}, {-q, -p}};
   at = Quiet[FullSimplify[D[t, variable].Inverse[t] + t.a0.Inverse[t]]];
   wt = Quiet[FullSimplify[detT w]];
   If[! functionSystemZeroQ[D[w, variable] + p w],
    AppendTo[issues, <|"code" -> "functionSystemWInconsistent", "residual" -> Quiet[FullSimplify[D[w, variable] + p w]]|>]
    ];
   If[explicitWT =!= Automatic && ! functionSystemZeroQ[explicitWT - wt],
    AppendTo[issues, <|"code" -> "functionSystemWTInconsistent", "expectedWT" -> wt, "providedWT" -> explicitWT|>]
    ];
   If[! functionSystemZeroQ[D[wt, variable]/wt - Tr[at]],
    AppendTo[issues, <|"code" -> "functionSystemTraceInconsistent", "residual" -> Quiet[FullSimplify[D[wt, variable]/wt - Tr[at]]]|>]
    ];
   derivativeTerms = compileDerivativeTerms[at, variable];
   If[derivativeTerms === $Failed,
    AppendTo[issues, <|"code" -> "functionSystemDerivativeNotFiniteLaurent", "AT" -> at|>]
    ];
   shrinkData = compileShrinkTerms[wt, variable, spec];
   If[shrinkData === $Failed,
    AppendTo[issues, <|"code" -> "functionSystemWTNotFiniteCompatibleLaurent", "WT" -> wt|>]
    ];
   If[issues =!= {}, Return[<|"status" -> "invalid", "issues" -> issues|>]];
   <|
    "status" -> "compiled", "input" -> Join[spec, <|"variable" -> variable, "T" -> t|>],
    "variable" -> variable, "P" -> p, "Q" -> q, "T" -> t, "W" -> w,
    "A0" -> a0, "AT" -> at, "WT" -> wt,
    "derivativeTerms" -> derivativeTerms,
    "shrinkTerms" -> shrinkData["terms"],
    "shrinkZeroPointShift" -> shrinkData["zeroPointShift"]
    |>
   ];


compileFunctionSystem[data_Association] /; KeyExistsQ[data, "endpoints"] := Module[{spec = normalizeLineFunctionSystem[data]},
   If[Head[spec] === Missing,
    <|"status" -> "invalid", "issues" -> {<|"code" -> "invalidLineFunctionSystem", "value" -> spec|>}|>,
    compileFunctionSystem[spec]
    ]
   ];


(* 补齐每条线的默认 metadata。thetaConvention 固定为 mergedTwoTheta。
   endpoints 是有序输入；对 masslessFull，第一端点定义单 n=1 的反对称方向。 *)
completeLineMetadata[line_, vertexSignAssoc_] := Module[
   {endpoints, referenceEndpoints, massType, skType, state, packType, metadata, functionSpec, compiled},
   endpoints = line["endpoints"];
   referenceEndpoints = endpoints;
   massType = line["massType"];
   skType = inferSKType[endpoints, vertexSignAssoc];
   state = "full";
   packType = inferPackType[massType, skType, state];
   metadata = Join[
     line,
    <|
     "massType" -> massType,
     "skType" -> skType,
     "state" -> state,
     "thetaConvention" -> "mergedTwoTheta",
     "packType" -> packType,
     "masslessN1ReferenceEndpoint" -> If[packType === "masslessFull", referenceEndpoints[[1]], Missing["NotApplicable"]],
     "masslessN1OppositeEndpoint" -> If[packType === "masslessFull", referenceEndpoints[[2]], Missing["NotApplicable"]]
      |>
     ];
   If[massType =!= "massive", Return[metadata]];
   functionSpec = normalizeLineFunctionSystem[metadata];
   compiled = If[Head[functionSpec] === Missing,
     <|"status" -> "invalid", "issues" -> {<|"code" -> "invalidLineFunctionSystem", "value" -> functionSpec|>}|>,
     compileFunctionSystem[functionSpec]
     ];
   Join[metadata, <|
     "functionSystem" -> Lookup[compiled, "input", Lookup[metadata, "functionSystem", Missing["InvalidFunctionSystem"]]],
     "compiledFunctionSystem" -> compiled
     |>]
   ];


(* 把用户 case 解析为通用拓扑对象。公开输入分别声明 loopExternalMomenta 与
   independentExternalMomenta；额外字段不会进入规范化 topology。 *)
parseTopology::missingkeys = "case 缺少必需字段：`1`。";
parseTopology::badinput = "case 输入 preflight 失败：`1`。";
parseTopology::badfunction = "massive line 的函数系统编译失败：`1`。";
parseTopology::badsectorkey = "sectorKey 与当前 ibpMode 不一致：`1`。";


parseTopology[case_Association] := Module[
   {vertices, vertexIds, vertexSignAssoc, rootExternalLegEnergies, rawLines, lines, badFunctionLines, loopMomenta,
   effectiveLoopExternalMomenta, declaredLoopExternalMomenta, loopKinematicRules, resolvedLoopKinematicRules,
    independentExternalMomenta, declaredIndependentExternalMomenta, loopDeclarationAudit,
    loopVectorAtoms, fixedExternalVectorAtoms, momentumDecompositionBasis,
    magnitudeKinematicRules, resolvedMagnitudeKinematicRules, kinematicRules, kinematicAudit,
    ispData, nV, nE, nL, nK, bMatrix, vertexLines,
    eMomenta, momentumBasis, momentumData, loopCoeffMatrix, externalCoeffMatrix, externalPartList,
    topologyMomentumAudit, graphAudit, routingAudit, declarationAudit, capabilities, ibpMode,
    rawLineMomenta, normalizedLineMomenta, referenceLoopMatrix, referenceLoopResiduals,
    loopMomentumInputToNormalizedRules,
      topoContext},
   If[caseInputPreflightErrorQ[case],
    If[caseInputMissingRequiredKeysQ[case],
     Message[parseTopology::missingkeys, caseInputRequirementReport[case]["missingRequiredKeys"]],
     Message[parseTopology::badinput, Lookup[caseInputErrorReport[case], "issues", {}]]
     ];
    Return[$Failed]
    ];
   vertices = normalizeVertex /@ case["vertices"];
   vertexIds = Lookup[vertices, "id"];
   vertexSignAssoc = AssociationThread[vertexIds -> Lookup[vertices, "vertexType"]];
   rootExternalLegEnergies = AssociationThread[vertexIds -> Lookup[vertices, "externalLegEnergy"]];
   rawLines = normalizeLine /@ case["lines"];
   lines = completeLineMetadata[#, vertexSignAssoc] & /@ rawLines;
   badFunctionLines = Select[
     MapIndexed[Join[<|"lineIndex" -> First[#2], "lineId" -> Lookup[#1, "id", Missing["NoLineId"]]|>, Lookup[#1, "compiledFunctionSystem", <||>]] &, lines],
     KeyExistsQ[#, "status"] && #["status"] =!= "compiled" &
     ];
   If[badFunctionLines =!= {},
    Message[parseTopology::badfunction, Lookup[badFunctionLines, {"lineIndex", "lineId", "issues"}]];
    Return[$Failed]
    ];
   loopMomenta = case["loopMomenta"];
    topologyMomentumAudit = ds016TopologyAndMomentumAudit[case, lines, vertexIds];
   graphAudit = Lookup[topologyMomentumAudit, "graph", <||>];
   routingAudit = Lookup[topologyMomentumAudit, "routing", <||>];
    declarationAudit = Lookup[topologyMomentumAudit, "declarations", <||>];
    capabilities = Lookup[topologyMomentumAudit, "capabilities", <||>];
    ibpMode = Lookup[routingAudit, "ibpMode", ds016ResolveIBPMode[case, Lookup[graphAudit, "graphLoopCount", 0]]];
   declaredLoopExternalMomenta = Lookup[declarationAudit, "loopExternalMomenta",
     Lookup[case, "loopExternalMomenta", {}]];
   declaredIndependentExternalMomenta = Lookup[declarationAudit, "independentExternalMomenta",
     Lookup[case, "independentExternalMomenta", {}]];
   loopDeclarationAudit = Lookup[declarationAudit, "loopExternalAudit", <||>];
   (* 过完备声明仍可生成 symbolic IBP，但不能把可吸收到圈变量平移中的额外方向
      虚增为 q.k 坐标。核心闭合使用 affine quotient 的必要基，原声明单独保留给审计。 *)
   effectiveLoopExternalMomenta = If[
     Lookup[loopDeclarationAudit, "status", "exact"] === "overcomplete",
     Lookup[loopDeclarationAudit, "requiredBasisDirections",
      Lookup[declarationAudit, "requiredLoopExternalDirections", declaredLoopExternalMomenta]],
     declaredLoopExternalMomenta
     ];
   independentExternalMomenta = declaredIndependentExternalMomenta;
   rawLineMomenta = Lookup[lines, "momentum", 0];
   referenceLoopMatrix = Lookup[routingAudit, "referenceLoopMatrix", {}];
   referenceLoopResiduals = Lookup[routingAudit, "referenceExternalResiduals", {}];
   normalizedLineMomenta = If[
     Lookup[routingAudit, "ibpMode", "invalid"] === "full" &&
      Lookup[routingAudit, "status", "invalid"] === "valid" && loopMomenta =!= {} &&
      MatrixQ[referenceLoopMatrix] && Length[referenceLoopMatrix] === Length[loopMomenta],
     Expand /@ (
       Lookup[routingAudit, "loopCoefficientMatrix", {}] . Inverse[referenceLoopMatrix] . loopMomenta +
        Lookup[routingAudit, "shiftInvariantLineResiduals", {}]
       ),
      rawLineMomenta
      ];
   (* 规范圈动量以参考传播子为新变量时，用户 ISP 仍定义在输入圈动量上。
      这里先建立 input -> normalized 的仿射反变换，再用同一变换改写 ISP；否则
      qIBP 会把 sp[ell,k] 错当成平移后的 sp[ell_new,k]。 *)
   loopMomentumInputToNormalizedRules = If[
     Lookup[routingAudit, "ibpMode", "invalid"] === "full" &&
      Lookup[routingAudit, "status", "invalid"] === "valid" && loopMomenta =!= {} &&
      MatrixQ[referenceLoopMatrix] && Length[referenceLoopMatrix] === Length[loopMomenta] &&
      Length[referenceLoopResiduals] === Length[loopMomenta],
     Thread[
      loopMomenta ->
       (Expand /@ (Inverse[referenceLoopMatrix] . (loopMomenta - referenceLoopResiduals)))
      ],
     {}
     ];
    (* timeOnly 是表示边界，不只是少生成一类算符。此模式下即使图有结构圈，
       所有传播子的动量幂也进入显式系数，tree pack 不得保留 b/bS。 *)
    lines = MapIndexed[
      Join[#1, <|
         "rawMomentum" -> rawLineMomenta[[First[#2]]],
         "momentum" -> normalizedLineMomenta[[First[#2]]],
         "loopLineQ" -> (ibpMode === "full" && MemberQ[Lookup[graphAudit, "cycleLineIndices", {}], First[#2]]),
         "bridgeQ" -> MemberQ[Lookup[graphAudit, "bridgeLineIndices", {}], First[#2]],
         "linePowerMode" -> If[
           ibpMode === "full" && MemberQ[Lookup[graphAudit, "cycleLineIndices", {}], First[#2]],
           "indexed",
           "fixedCoefficient"
           ]
         |>] &,
      lines
      ];
   ispData = Map[
     Function[entry,
      With[{normalizedEntry = normalizeISP[entry]},
       Join[
        normalizedEntry,
        <|"expr" -> Expand[Lookup[normalizedEntry, "expr", 0] /.
             loopMomentumInputToNormalizedRules]|>
        ]
       ]
      ],
     Lookup[case, "ispData", {}]
     ];
   nV = Length[vertexIds];
   nE = Length[lines];
   nL = Length[loopMomenta];
    nK = Length[effectiveLoopExternalMomenta];
   eMomenta = Lookup[lines, "momentum"];
   bMatrix = Table[
     Which[
      lines[[e]]["endpoints"][[1]] === vertexIds[[v]], 1,
      lines[[e]]["endpoints"][[2]] === vertexIds[[v]], -1,
      True, 0
      ],
     {v, nV}, {e, nE}
     ];
   vertexLines = Table[
     Select[Table[{e, bMatrix[[v, e]]}, {e, nE}], #[[2]] =!= 0 &],
     {v, nV}
     ];
   (* independentExternalMomenta 定义模长坐标，不是复合 momentum-IBP 的外向量基。
      这里只补入其原子向量，以便合法解析 bridge 上的 p1+p2 等固定动量。 *)
   loopVectorAtoms = ds016MomentumAtoms[Join[loopMomenta, effectiveLoopExternalMomenta]];
   fixedExternalVectorAtoms = Complement[
     ds016MomentumAtoms[independentExternalMomenta],
     loopVectorAtoms
     ];
   momentumDecompositionBasis = Join[loopMomenta, effectiveLoopExternalMomenta, fixedExternalVectorAtoms];
   momentumBasis = momentumDecompositionBasis;
   momentumData = linearMomentumExpressionData[#, momentumBasis] & /@ eMomenta;
   loopCoeffMatrix = If[nL === 0, ConstantArray[{}, nE],
     Take[Lookup[momentumData, "coefficients", {}], All, UpTo[nL]]
     ];
   externalCoeffMatrix = If[nK === 0, ConstantArray[{}, nE],
     Take[Lookup[momentumData, "coefficients", {}], All, {nL + 1, nL + nK}]
     ];
   externalPartList = Table[
     eMomenta[[e]] - Sum[loopCoeffMatrix[[e, l]] loopMomenta[[l]], {l, nL}],
     {e, nE}
     ];
   topoContext = <|"loopMomenta" -> loopMomenta, "effectiveLoopExternalMomenta" -> effectiveLoopExternalMomenta,
     "independentExternalMomenta" -> independentExternalMomenta, "lines" -> lines,
     "sectorExternalLegEnergyByVertex" -> rootExternalLegEnergies, "extLegs" -> Lookup[case, "extLegs", {}],
     "nL" -> nL, "nK" -> nK|>;
   kinematicRules = Lookup[case, "kinematicRules", Automatic];
   kinematicAudit = resolveKinematicRulesForCase[case, topoContext];
   (* 动量声明和坐标 Jacobian 是两层独立门禁；坐标过完备不得重新开启
      declaration audit 已关闭的能力，也不得让一般满秩混合坐标失去 ds。 *)
   capabilities = Join[capabilities, <|
      "derivativeUsableQ" -> TrueQ[
        Lookup[capabilities, "derivativeUsableQ", False] &&
         Lookup[kinematicAudit, "completeQ", False] &&
         ! Lookup[kinematicAudit, "overcompleteQ", False]
        ],
      "inverseKinematicsUsableQ" -> TrueQ[
        Lookup[capabilities, "inverseKinematicsUsableQ", False] &&
         Lookup[kinematicAudit, "inverseAvailableQ", False]
        ],
      "backendExportUsableQ" -> TrueQ[
        Lookup[capabilities, "backendExportUsableQ", False] &&
         Lookup[kinematicAudit, "completeQ", False] &&
         ! Lookup[kinematicAudit, "overcompleteQ", False]
        ]
      |>];
   loopKinematicRules = Lookup[
     kinematicAudit,
     "rawLoopRules",
     Automatic
     ];
   resolvedLoopKinematicRules = Lookup[
     kinematicAudit,
     "resolvedLoopRules",
     normalizeLoopKinematicRulesForTopology[loopKinematicRules, topoContext]
     ];
   magnitudeKinematicRules = Lookup[
     kinematicAudit,
     "rawExternalLegRules",
     Automatic
     ];
   resolvedMagnitudeKinematicRules = Lookup[
     kinematicAudit,
     "resolvedExternalLegRules",
     normalizeMagnitudeKinematicRulesForTopology[magnitudeKinematicRules, topoContext]
     ];
   topoContext = Join[topoContext, <|
      "resolvedLoopKinematicRules" -> resolvedLoopKinematicRules,
      "magnitudeKinematicRules" -> magnitudeKinematicRules,
      "resolvedMagnitudeKinematicRules" -> resolvedMagnitudeKinematicRules,
      "kinematicRules" -> kinematicRules,
      "kinematicCoordinateAudit" -> kinematicAudit
      |>];
   <|
    "name" -> Lookup[case, "name", "unnamed"],
    "vertices" -> vertices,
    "vertexIds" -> vertexIds,
    "vertexSignAssoc" -> vertexSignAssoc,
    "lines" -> lines,
    "extLegs" -> Lookup[case, "extLegs", {}],
    "sectorExternalLegEnergyByVertex" -> rootExternalLegEnergies,
    "activeVertexIds" -> vertexIds,
    "fixedAVertexValues" -> <||>,
    "loopMomenta" -> loopMomenta,
     "ibpMode" -> ibpMode,
    "graphLoopCount" -> Lookup[graphAudit, "graphLoopCount", Missing["graphLoopCount"]],
    "graphTopologyAudit" -> graphAudit,
     "loopMomentumRoutingAudit" -> routingAudit,
     "loopMomentumInputToNormalizedRules" -> loopMomentumInputToNormalizedRules,
    "normalizedLineMomenta" -> normalizedLineMomenta,
    "momentumDeclarationAudit" -> declarationAudit,
    "capabilities" -> capabilities,
    "cycleLineIndices" -> Lookup[graphAudit, "cycleLineIndices", {}],
    "bridgeLineIndices" -> Lookup[graphAudit, "bridgeLineIndices", {}],
    "selfLoopLineIndices" -> Lookup[graphAudit, "selfLoopLineIndices", {}],
    "loopExternalMomenta" -> declaredLoopExternalMomenta,
    "momentumDecompositionBasis" -> momentumDecompositionBasis,
    "fixedExternalVectorAtoms" -> fixedExternalVectorAtoms,
    "effectiveLoopExternalMomenta" -> effectiveLoopExternalMomenta,
    "independentExternalMomenta" -> independentExternalMomenta,
    "loopKinematicRules" -> loopKinematicRules,
    "resolvedLoopKinematicRules" -> resolvedLoopKinematicRules,
    "magnitudeKinematicRules" -> magnitudeKinematicRules,
    "resolvedMagnitudeKinematicRules" -> resolvedMagnitudeKinematicRules,
    "kinematicRules" -> kinematicRules,
    "kinematicCoordinateAudit" -> kinematicAudit,
    "ispData" -> ispData,
    "nV" -> nV,
    "nE" -> nE,
    "nL" -> nL,
    "nK" -> nK,
    "bMatrix" -> bMatrix,
    "vertexLines" -> vertexLines,
    "loopCoeffMatrix" -> loopCoeffMatrix,
    "externalCoeffMatrix" -> externalCoeffMatrix,
    "externalPartList" -> externalPartList,
    "zeroPointRules" -> Lookup[case, "zeroPointRules", {}],
    "rootZeroPointRules" -> Lookup[case, "zeroPointRules", {}],
    "symmetryRules" -> Lookup[case, "symmetryRules", {}],
    "parityConstraints" -> Lookup[case, "parityConstraints", {}],
    "sectorVertexRepresentativeMap" -> AssociationThread[vertexIds -> vertexIds]
    |>
   ];


(* ::Chapter:: *)
(*统一 J 指标包与离散态枚举*)

(* 本章实现 massive/massless 的自然指标打包。离散态按每条线 metadata 枚举，
   不再使用旧 bubble 原型中的 Tuples[{0,1}, 2 nE]。 *)


lineIndexedPowerQ[line_Association] := Lookup[line, "linePowerMode", "indexed"] === "indexed";


linePackNPositions[line_Association, packType_String] := Module[{offset = If[lineIndexedPowerQ[line], 1, 0]},
   Switch[packType,
    "massiveFull" | "massiveCross", offset + {1, 2},
    "masslessFull", {offset + 1},
    _, {}
    ]
   ];


(* 只有 cycle line 的 pack 第一槽是整数动量幂。bridge 的第一槽可能是 n，
   因而任何物理幂次、排序或投影代码都必须通过本访问器读取。 *)
lineIntegerPowerIndex[topo_Association, J[_, linePacks_, _], e_Integer] := If[
   lineIndexedPowerQ[topo["lines"][[e]]],
   linePacks[[e, 1]],
   0
   ];


linePackNPosition[topo_Association, e_Integer, endpointSlot_Integer] := Module[
   {positions = linePackNPositions[topo["lines"][[e]], actualLinePackType[topo, e, makeLinePack[topo["lines"][[e]]]]]},
   If[1 <= endpointSlot <= Length[positions], positions[[endpointSlot]], Missing["NoEndpointState"]]
   ];


makeLinePack[line_Association] := Module[{id = line["id"], indexedQ = lineIndexedPowerQ[line]},
   Switch[line["packType"],
    "massiveFull", If[indexedQ, {b[id], n[id, 1], n[id, 2]}, {n[id, 1], n[id, 2]}],
    "massiveCross", If[indexedQ, {b[id], n[id, 1], n[id, 2]}, {n[id, 1], n[id, 2]}],
    "masslessFull", If[indexedQ, {b[id], n[id]}, {n[id]}],
    "masslessCross", If[indexedQ, {b[id]}, {}],
    "shrunk", If[indexedQ, {bS[id]}, {}],
    _, Message[makeLinePack::badtype, line["packType"], id]; $Failed
    ]
   ];
makeLinePack::badtype = "未知 packType `1`，line id = `2`。";


makeLinePacks[topo_Association] := makeLinePack /@ topo["lines"];


makeSectorMetadataBase018[topo_Association] := Module[
   {fixedA = Lookup[topo, "fixedAVertexValues", <||>], active = Lookup[topo, "activeVertexIds", topo["vertexIds"]],
    lines = topo["lines"], vertexIds = topo["vertexIds"], repMap, originalSlotByVertex, compactSlotByVertex,
    activeOriginalGroups, vertexSlots, compactASlots, lineSlots},
   repMap = Lookup[topo, "sectorVertexRepresentativeMap", AssociationThread[topo["vertexIds"] -> topo["vertexIds"]]];
   originalSlotByVertex = AssociationThread[vertexIds -> Range[Length[vertexIds]]];
   activeOriginalGroups = Table[
     v -> Select[vertexIds, Lookup[repMap, #, #] === v &],
     {v, active}
     ];
   compactSlotByVertex = Association @ Flatten[
      Table[Thread[activeOriginalGroups[[i, 2]] -> i], {i, Length[activeOriginalGroups]}],
      1
      ];
   vertexSlots = Table[
     <|
      "slot" -> i,
      "vertexId" -> vertexIds[[i]],
      "representativeVertexId" -> Lookup[repMap, vertexIds[[i]], vertexIds[[i]]],
      "aSymbol" -> a[vertexIds[[i]]],
      "activeQ" -> MemberQ[active, vertexIds[[i]]],
      "fixedValue" -> If[AssociationQ[fixedA] && KeyExistsQ[fixedA, vertexIds[[i]]], fixedA[vertexIds[[i]]], None],
      "compactASlot" -> Lookup[compactSlotByVertex, vertexIds[[i]], None]
      |>,
     {i, Length[vertexIds]}
     ];
   compactASlots = Table[
     <|
      "compactSlot" -> i,
      "representativeVertexId" -> active[[i]],
      "originalVertexIds" -> activeOriginalGroups[[i, 2]],
      "originalSlots" -> Lookup[originalSlotByVertex, activeOriginalGroups[[i, 2]]],
      "aSymbol" -> a[active[[i]]]
      |>,
     {i, Length[active]}
     ];
   lineSlots = Table[
     With[{pack = makeLinePack[lines[[e]]], endpoints = lines[[e]]["endpoints"], originalEndpoints = Lookup[lines[[e]], "originalEndpoints", lines[[e]]["endpoints"]]},
      <|
       "slot" -> e,
       "lineId" -> lines[[e]]["id"],
       "packType" -> lines[[e]]["packType"],
       "massType" -> lines[[e]]["massType"],
       "state" -> lines[[e]]["state"],
       "endpoints" -> endpoints,
       "originalEndpoints" -> originalEndpoints,
       "masslessN1ReferenceEndpoint" -> Lookup[lines[[e]], "masslessN1ReferenceEndpoint", Missing["NotApplicable"]],
       "masslessN1OppositeEndpoint" -> Lookup[lines[[e]], "masslessN1OppositeEndpoint", Missing["NotApplicable"]],
       "endpointOriginalASlots" -> Lookup[originalSlotByVertex, originalEndpoints],
       "endpointCompactASlots" -> Lookup[compactSlotByVertex, endpoints, None],
       "linePowerMode" -> Lookup[lines[[e]], "linePowerMode", "indexed"],
       "bPosition" -> linePackBPosition[lines[[e]]],
       "bSymbol" -> If[lineIndexedPowerQ[lines[[e]]], First[pack], None],
       "nPositions" -> linePackNPositions[lines[[e]], lines[[e]]["packType"]],
       "packTemplate" -> pack
       |>
      ],
     {e, Length[lines]}
     ];
   <|
    "caseName" -> topo["name"],
    "sectorShrunkLines" -> Lookup[topo, "sectorShrunkLines", {}],
    "sectorKey" -> sectorKeyFromShrunkLines[topo, Lookup[topo, "sectorShrunkLines", {}]],
    "aSlotMode" -> "compactActiveSlots",
    "sectorVertexRepresentativeMap" -> repMap,
    "vertexIdToOriginalASlot" -> originalSlotByVertex,
    "vertexIdToCompactASlot" -> compactSlotByVertex,
    "vertexSlots" -> vertexSlots,
    "compactASlots" -> compactASlots,
    "activeASlots" -> Range[Length[active]],
    "lineSlots" -> lineSlots,
    "lineIdToSlot" -> AssociationThread[Lookup[lines, "id"] -> Range[Length[lines]]],
     "bSymbolToLineSlot" -> Association@Cases[
       MapIndexed[Lookup[#1, "bSymbol", None] -> First[#2] &, lineSlots],
       Rule[symbol_, slot_] /; symbol =!= None :> Rule[symbol, slot]
       ],
    "ispSlots" -> Table[
      <|"slot" -> j, "indexSymbol" -> ispN[j], "data" -> topo["ispData"][[j]]|>,
      {j, Length[topo["ispData"]]}
      ]
    |>
   ];


(* sector key 的位序严格等于 root propagator 顺序。结果必须保留为字符串，
   因为前导零编码了排在最前的传播子已经收缩。 *)
sectorKeyFromShrunkLines[rootLineCount_Integer?NonNegative, shrunkLines_List] := StringJoin[
   If[MemberQ[shrunkLines, #], "0", "1"] & /@ Range[rootLineCount]
   ];


sectorKeyFromShrunkLines[topo_Association, shrunkLines_List] := sectorKeyFromShrunkLines[
   Length[Lookup[topo, "lines", {}]],
   shrunkLines,
   Lookup[topo, "ibpMode", "full"]
   ];


sectorKeyFromShrunkLines[
   rootLineCount_Integer?NonNegative,
   shrunkLines_List,
   "timeOnly"
   ] := sectorKeyFromShrunkLines[rootLineCount, shrunkLines];


sectorKeyFromShrunkLines[
   _Integer?NonNegative,
   shrunkLines_List,
   _
   ] := If[shrunkLines === {}, "top", StringRiffle["e" <> ToString[#] & /@ shrunkLines, "_"]];


sectorKeyForModeQ[key_, "timeOnly", rootLineCount_Integer?NonNegative] :=
   StringQ[key] && StringLength[key] === rootLineCount &&
    And @@ (MemberQ[{"0", "1"}, #] & /@ Characters[key]);


sectorKeyForModeQ[key_, _String, _Integer?NonNegative] := StringQ[key] &&
   (key === "top" || StringMatchQ[key, "e" ~~ DigitCharacter .. ~~ ("_e" ~~ DigitCharacter ..) ...]);


sectorKeyTopQ[key_, "timeOnly", rootLineCount_Integer?NonNegative] :=
   sectorKeyForModeQ[key, "timeOnly", rootLineCount] &&
    And @@ (# === "1" & /@ Characters[key]);


sectorKeyTopQ[key_, _String, _Integer?NonNegative] := key === "top";


sectorKeySchemaFromTopology[topo_Association] := If[
   Lookup[topo, "ibpMode", "full"] === "timeOnly",
   <|
    "type" -> "rootPropagatorBitString",
    "rootLineOrder" -> Lookup[Lookup[topo, "lines", {}], "id", Range[Length[Lookup[topo, "lines", {}]]]],
    "width" -> Length[Lookup[topo, "lines", {}]],
    "contractedBit" -> "0",
    "uncontractedBit" -> "1",
    "storageType" -> "String"
    |>,
   <|"type" -> "legacyContractedLineList", "storageType" -> "String"|>
   ];


sectorMetadataKey[metadata_Association] := Lookup[
   metadata,
   "sectorKey",
   sectorKeyFromShrunkLines[
    Lookup[metadata, "rootLineCount", Length[Lookup[metadata, "lineSlots", {}]]],
    Lookup[metadata, "sectorShrunkLines", {}],
    Lookup[metadata, "ibpMode", "full"]
    ]
   ];


sectorMetadataLinePackLengths[metadata_Association] := Length /@ Lookup[Lookup[metadata, "lineSlots", {}], "packTemplate", {}];


parsedTopologyQ[spec_Association] := TrueQ[KeyExistsQ[spec, "lines"] && KeyExistsQ[spec, "nL"]];
parsedTopologyQ[_] := False;


normalizeTopologySpec[spec_Association] := If[parsedTopologyQ[spec], spec, parseTopology[spec]];
normalizeTopologySpec[spec_] := spec;


indexContainsLineSymbolQ[index_, symbol_] := ! FreeQ[index, symbol];


indexHasAnyLineSymbolQ[index_] := ! FreeQ[index, _b | _bS];


linePackMatchesSlotQ[pack_List, slot_Association] := Module[
   {template = Lookup[slot, "packTemplate", {}], bSymbol = Lookup[slot, "bSymbol", None],
    bPosition = Lookup[slot, "bPosition", Missing["FixedLinePower"]]},
   TrueQ[Length[pack] === Length[template]] &&
     TrueQ[
      bSymbol === None || Head[bPosition] === Missing ||
       indexContainsLineSymbolQ[pack[[bPosition]], bSymbol] ||
       ! indexHasAnyLineSymbolQ[pack[[bPosition]]]
      ]
   ];


integralMatchesSectorMetadataQ[J[aList_, linePacks_, ispList_], metadata_Association] := Module[
   {lineSlots = Lookup[metadata, "lineSlots", {}], compactASlots = Lookup[metadata, "compactASlots", {}], ispSlots = Lookup[metadata, "ispSlots", {}]},
   TrueQ[
    Length[aList] === Length[compactASlots] &&
     Length[linePacks] === Length[lineSlots] &&
     Length[ispList] === Length[ispSlots] &&
     And @@ MapThread[linePackMatchesSlotQ, {linePacks, lineSlots}]
    ]
   ];


integralSectorKey[J[aList_, linePacks_, ispList_], metadataList_List] := Module[
   {matches},
   matches = Select[metadataList, integralMatchesSectorMetadataQ[J[aList, linePacks, ispList], #] &];
   Which[
    Length[matches] == 1, sectorMetadataKey[First[matches]],
    Length[matches] == 0, Missing["NoMatchingSector"],
    True, Missing["AmbiguousSector", sectorMetadataKey /@ matches]
    ]
   ];


batchSectorMetadataList[batch_Association, topoSpec_: Automatic] := Module[
   {fromBatch, topo, topMetadata},
   fromBatch = Lookup[batch, "sectorMetadataList", Missing["NoSectorMetadataList"]];
   If[ListQ[fromBatch], Return[fromBatch]];
   If[AssociationQ[Lookup[batch, "sectorMetadata", Missing["NoSectorMetadata"]]], Return[{batch["sectorMetadata"]}]];
   topo = normalizeTopologySpec[topoSpec];
   topMetadata = If[parsedTopologyQ[topo], makeSectorMetadata[topo], Missing["NoSectorMetadata"]];
   If[Head[topMetadata] === Missing, {}, {topMetadata}]
   ];


normaliseKiraOrderingSpec[spec_] := Which[
   AssociationQ[spec], spec,
   True, <||>
   ];


allowedKiraOrderingKeys[] := {"IntegralOrder", "ManualIntegralOrder", "PreferredIntegrals", "PreferredPriority", "SectorRank", "SectorOrder"};


validateKiraOrderingSpec[Automatic] := <|"status" -> "ok", "issues" -> {}|>;
validateKiraOrderingSpec[spec_] /; ! AssociationQ[spec] := <|"status" -> "invalidKiraOrdering", "reason" -> "KiraOrdering must be Automatic or an Association", "kiraOrdering" -> spec|>;
validateKiraOrderingSpec[spec_Association] := Module[
   {unknownKeys, badValueData},
   unknownKeys = Complement[Keys[spec], allowedKiraOrderingKeys[]];
   badValueData = DeleteCases[
     {
      If[KeyExistsQ[spec, "IntegralOrder"] && ! ListQ[spec["IntegralOrder"]], <|"optionKey" -> "IntegralOrder", "optionValue" -> spec["IntegralOrder"]|>, Nothing],
      If[KeyExistsQ[spec, "ManualIntegralOrder"] && ! ListQ[spec["ManualIntegralOrder"]], <|"optionKey" -> "ManualIntegralOrder", "optionValue" -> spec["ManualIntegralOrder"]|>, Nothing],
      If[KeyExistsQ[spec, "PreferredIntegrals"] && ! ListQ[spec["PreferredIntegrals"]], <|"optionKey" -> "PreferredIntegrals", "optionValue" -> spec["PreferredIntegrals"]|>, Nothing],
      If[KeyExistsQ[spec, "PreferredPriority"] && ! MemberQ[{"BeforeB", "AfterB"}, spec["PreferredPriority"]], <|"optionKey" -> "PreferredPriority", "optionValue" -> spec["PreferredPriority"]|>, Nothing],
      If[KeyExistsQ[spec, "SectorRank"] && ! AssociationQ[spec["SectorRank"]], <|"optionKey" -> "SectorRank", "optionValue" -> spec["SectorRank"]|>, Nothing],
      If[KeyExistsQ[spec, "SectorOrder"] && ! ListQ[spec["SectorOrder"]], <|"optionKey" -> "SectorOrder", "optionValue" -> spec["SectorOrder"]|>, Nothing]
      },
     Nothing
     ];
   If[unknownKeys === {} && badValueData === {},
    <|"status" -> "ok", "issues" -> {}|>,
    <|"status" -> "invalidKiraOrdering", "reason" -> "KiraOrdering contains unknown keys or invalid values", "unknownKiraOrderingKeys" -> unknownKeys, "malformedKiraOrderingValues" -> badValueData, "allowedKiraOrderingKeys" -> allowedKiraOrderingKeys[]|>
    ]
   ];


resolveKiraOrderingSpec[_Association, _, optSpec_] :=
  normaliseKiraOrderingSpec[Replace[optSpec, Automatic -> <||>]];


preferredIntegralRank[int_, orderingSpec_Association] := Module[
   {manual = Lookup[orderingSpec, "IntegralOrder", Lookup[orderingSpec, "ManualIntegralOrder", {}]], preferred, pos},
   If[! ListQ[manual], manual = {}];
   preferred = DeleteDuplicates@Join[manual, Lookup[orderingSpec, "PreferredIntegrals", {}]];
   pos = FirstPosition[preferred, int, Missing["NotPreferred"], {1}, Heads -> False];
   If[Head[pos] === Missing, 10^9, First[pos] - 1]
   ];


automaticSectorRank[sectorKey_, metadataList_List] := Module[{metadata},
   metadata = SelectFirst[
     metadataList,
     AssociationQ[#] && sectorMetadataKey[#] === sectorKey &,
     Missing["NoMatchingSector"]
     ];
   If[
    AssociationQ[metadata],
    -Length[Lookup[metadata, "sectorShrunkLines", {}]],
    0
    ]
   ];


sectorRankFromOrdering[sectorKey_, orderingSpec_Association, metadataList_: {}] := Module[
   {rank = Lookup[orderingSpec, "SectorRank", <||>], order = Lookup[orderingSpec, "SectorOrder", {}], pos},
   Which[
    AssociationQ[rank] && KeyExistsQ[rank, sectorKey], rank[sectorKey],
    ListQ[order] && MemberQ[order, sectorKey], First[FirstPosition[order, sectorKey]] - 1,
    True, automaticSectorRank[sectorKey, metadataList]
    ]
   ];


numericIndexValue[x_Integer] := x;
numericIndexValue[x_Rational] := x;
numericIndexValue[x_Real] := x;
numericIndexValue[_] := 0;


integralSortKey[J[aList_, linePacks_, ispList_], orderingSpec_: <||>, metadataList_: {}] := Module[
    {int = J[aList, linePacks, ispList], spec = normaliseKiraOrderingSpec[orderingSpec], bValues, aValues, ispValues,
     continuousValues, negativeIndexCount, negativeGroupRank, effectiveContinuousLevel, nValues,
     preferredRank, sectorKey, sectorRank, baseKey, priority},
    (* 进入 linearData 的 line pack 已经数值化；连续幂位于长度为 1 或 3 的 pack 首槽。
       所有非负连续指标先组成完整前缀，负指标组内部再按 penalty level 排序；固定 F 和
       长度为 3 的离散 building-block 指标均不参与负指标判定。 *)
   bValues = Cases[
     linePacks,
     pack_List /; MemberQ[{1, 3}, Length[pack]] && NumberQ[First[pack]] :>
       numericIndexValue[First[pack]]
     ];
    aValues = numericIndexValue /@ aList;
    ispValues = numericIndexValue /@ ispList;
    continuousValues = Join[aValues, bValues, ispValues];
    negativeIndexCount = Count[continuousValues, value_ /; TrueQ[value < 0]];
     negativeGroupRank = Boole[negativeIndexCount > 0];
     effectiveContinuousLevel = Total[Abs /@ continuousValues] + If[
       negativeIndexCount > 0,
       3 + 2 negativeIndexCount,
       0
       ];
    nValues = Flatten@Cases[
      linePacks,
      pack_List /; Length[pack] === 3 :> (numericIndexValue /@ Rest[pack])
      ];
   preferredRank = preferredIntegralRank[int, spec];
   sectorKey = integralSectorKey[int, metadataList];
     sectorRank = sectorRankFromOrdering[sectorKey, spec, metadataList];
      baseKey = {
        negativeGroupRank,
        effectiveContinuousLevel,
        negativeIndexCount,
       Total[Abs /@ bValues],
       Total[Clip[bValues, {0, Infinity}]],
      sectorRank,
      Total[Abs /@ aValues],
      Total[Abs /@ ispValues],
      Total[nValues],
      ToString[InputForm[int]]
      };
   priority = Lookup[spec, "PreferredPriority", "AfterB"];
    If[priority === "BeforeB",
      Prepend[baseKey, preferredRank],
      Insert[baseKey, preferredRank, 6]
      ]
   ];


sortIntegralsForKira[integrals_List, orderingSpec_: <||>, metadataList_: {}] := SortBy[integrals, integralSortKey[#, orderingSpec, metadataList] &];


integralOrderItemToIntegral[item_, integrals_List] := Which[
   Head[item] === J && MemberQ[integrals, item], item,
   IntegerQ[item] && 1 <= item <= Length[integrals], integrals[[item]],
   True, Missing["UnknownIntegralOrderItem", item]
   ];


normaliseIntegralOrder[order_List, integrals_List] := DeleteDuplicates @ DeleteMissing[
    integralOrderItemToIntegral[#, integrals] & /@ order
    ];


missingIntegralOrderItems[order_List, integrals_List] := Cases[
   integralOrderItemToIntegral[#, integrals] & /@ order,
   Missing["UnknownIntegralOrderItem", item_] :> item
   ];


kiraOrderingIntegralRequestList[orderingSpec_Association] := Module[
   {manual = Lookup[orderingSpec, "IntegralOrder", Lookup[orderingSpec, "ManualIntegralOrder", {}]],
    preferred = Lookup[orderingSpec, "PreferredIntegrals", {}]},
   If[! ListQ[manual], manual = {}];
   If[! ListQ[preferred], preferred = {}];
   DeleteDuplicates@Join[manual, preferred]
   ];


kiraOrderingMatchReport[orderingSpec_Association, integrals_List] := Module[
   {requested = kiraOrderingIntegralRequestList[orderingSpec], matched, missing},
   matched = normaliseIntegralOrder[requested, integrals];
   missing = DeleteDuplicates @ missingIntegralOrderItems[requested, integrals];
   <|
    "requestedIntegrals" -> requested,
    "matchedIntegrals" -> matched,
    "missingIntegralOrderItems" -> missing,
    "allRequestedIntegralsMatchedQ" -> TrueQ[missing === {}]
    |>
   ];


integralMetadataList[integrals_List, metadataList_List, orderingSpec_: <||>] := Table[
   <|
    "id" -> i,
    "sectorKey" -> integralSectorKey[integrals[[i]], metadataList],
    "sortKey" -> integralSortKey[integrals[[i]], orderingSpec, metadataList]
    |>,
   {i, Length[integrals]}
   ];




ruleSetsForVars[vars_List] := If[Length[vars] == 0,
   {{}},
   Thread[vars -> #] & /@ Tuples[{0, 1}, Length[vars]]
   ];


(* 离散态计数只用于完整性证书；实际 seed 仍由 enumerateDiscreteStates 全量展开。 *)
discreteStateCount[topo_Association] := Times @@ (discreteStateCountForLine /@ topo["lines"]);


(* ::Chapter:: *)
(*EOM canonical 与 forbidden-n 扫描*)

(* 本章只处理已经出现在 J 指标中的 Hankel 二阶导数态。
   seed 生成时应先枚举 n=0/1；若生成元把 massiveFull 的端点 n 推到 2 或更高，
   这里立即用 EOM 递推回 n=0/1。massless 的 {b,n} 只有两个状态，不用 Hankel EOM。 *)

vertexPosition[topo_Association, vertexId_] := Module[{pos},
   pos = FirstPosition[topo["vertexIds"], vertexId, Missing["NotFound"]];
   If[Head[pos] === Missing, Missing["VertexNotFound", vertexId], First[pos]]
   ];


activeAVertexIds[topo_Association] := Lookup[topo, "activeVertexIds", topo["vertexIds"]];


vertexRepresentative[topo_Association, vertexId_] := Lookup[
   Lookup[topo, "sectorVertexRepresentativeMap", AssociationThread[topo["vertexIds"] -> topo["vertexIds"]]],
   vertexId,
   vertexId
   ];


vertexASlot[topo_Association, vertexId_] := Module[
   {rep = vertexRepresentative[topo, vertexId], pos},
   pos = FirstPosition[activeAVertexIds[topo], rep, Missing["NotActiveVertex"]];
   If[Head[pos] === Missing, Missing["ActiveVertexNotFound", vertexId, rep], First[pos]]
   ];


shiftVertexA[J[aList_, linePacks_, ispList_], topo_Association, vertexId_, delta_] := Module[
   {newAList = aList, pos},
   pos = vertexASlot[topo, vertexId];
   (* compact slot 必须同时属于当前 sector topology 和当前 J；任何失配都直接失败，
      不能让 Part/Set 消息留下一个看似可继续使用的积分。 *)
   If[
    Head[pos] === Missing || ! IntegerQ[pos] || pos < 1 || pos > Length[newAList],
    Return[$Failed]
    ];
   newAList[[pos]] = newAList[[pos]] + delta;
   J[newAList, linePacks, ispList]
   ];


shiftLinePackEntry[J[aList_, linePacks_, ispList_], e_Integer, packPos_Integer, delta_] := Module[
   {newLinePacks = linePacks},
   newLinePacks[[e, packPos]] = newLinePacks[[e, packPos]] + delta;
   J[aList, newLinePacks, ispList]
   ];


setLinePackEntry[J[aList_, linePacks_, ispList_], e_Integer, packPos_Integer, value_] := Module[
   {newLinePacks = linePacks},
   newLinePacks[[e, packPos]] = value;
   J[aList, newLinePacks, ispList]
   ];


actualLinePackType[topo_Association, e_Integer, pack_List] := Module[
   {line = topo["lines"][[e]], declared, fullLength, shrunkLength},
   declared = line["packType"];
   If[declared === "shrunk", Return["shrunk"]];
   fullLength = Length[makeLinePack[line]];
   shrunkLength = If[lineIndexedPowerQ[line], 1, 0];
   If[Length[pack] === shrunkLength && fullLength =!= shrunkLength,
    "shrunk",
    declared
    ]
   ];


lineCompiledFunctionSystem[line_Association] := Lookup[line, "compiledFunctionSystem", <||>];


lineDerivativeTerms[line_Association, sourceState_Integer] := Select[
   Lookup[lineCompiledFunctionSystem[line], "derivativeTerms", {}],
   Lookup[#, "sourceState", Missing["NoSourceState"]] === sourceState &
   ];


applyCompiledEOMTerm[topo_Association, int_J, e_Integer, endpointSlot_Integer, term_Association] := Module[
   {result, endpointVertex, nPosition, xPower = term["xPower"]},
   endpointVertex = topo["lines"][[e, "endpoints", endpointSlot]];
   nPosition = linePackNPositions[topo["lines"][[e]], actualLinePackType[topo, e, int[[2, e]]]][[endpointSlot]];
   result = setLinePackEntry[int, e, nPosition, term["targetState"]];
   result = shiftVertexA[result, topo, endpointVertex, xPower];
   (* bridge 的固定动量幂会由 shiftLinePower 提到显式系数，因此先完成只接受裸 J 的顶点移位。 *)
   result = shiftLinePower[topo, result, e, -xPower];
   term["coefficient"] result
   ];


massiveEOMTarget[topo_Association, J[aList_, linePacks_, ispList_]] := Module[
   {lines = topo["lines"], packType, pack, nValue, target = Missing["NoEOMTarget"]},
   Do[
    pack = linePacks[[e]];
    packType = actualLinePackType[topo, e, pack];
    If[MemberQ[{"massiveFull", "massiveCross"}, packType],
     Do[
      nValue = pack[[linePackNPositions[lines[[e]], packType][[endpointSlot]]]];
      If[IntegerQ[nValue] && nValue >= 2 && Head[target] === Missing,
       target = <|"lineIndex" -> e, "endpointSlot" -> endpointSlot, "nValue" -> nValue|>
       ],
      {endpointSlot, 2}
      ]
     ],
    {e, Length[lines]}
    ];
   target
   ];


eomReduceIntegralAt[topo_Association, int_J, target_Association] := Module[
   {e, endpointSlot, nValue, line, terms},
   e = target["lineIndex"];
   endpointSlot = target["endpointSlot"];
   nValue = target["nValue"];
   line = topo["lines"][[e]];
   If[nValue =!= 2, Return[int]];
   terms = lineDerivativeTerms[line, 1];
   If[terms === {}, Return[int]];
   Expand[Total[applyCompiledEOMTerm[topo, int, e, endpointSlot, #] & /@ terms]]
   ];


applyEOMToIntegral[topo_Association, int_J] := Module[{target},
   target = massiveEOMTarget[topo, int];
   If[Head[target] === Missing,
    int,
    applyEOM[eomReduceIntegralAt[topo, int, target], topo]
    ]
   ];


applyEOM[expr_, topo_Association] := Expand[expr /. int_J :> applyEOMToIntegral[topo, int]];


(* masslessFull 的单 n 只允许 0/1。它以 line endpoints 的第一端点定义方向，
   所有导数直接在 0/1 间翻转；不再用有歧义的临时 n=2 混写 {20} 与 {11}。
   若其它缩并使该线两端点重合，反对称态 n=1 在等时点恒为零。 *)
integralShrunkLineIndices[
   topo_Association,
   J[aList_, linePacks_, ispList_]
   ] := Select[
   Range[Length[topo["lines"]]],
   actualLinePackType[topo, #, linePacks[[#]]] === "shrunk" &
   ];


integralTargetVertexRepresentativeMap[
   topo_Association,
   int_J
   ] := Module[
   {shrunkLines, originalPairs},
   shrunkLines = integralShrunkLineIndices[topo, int];
   originalPairs = Lookup[
       topo["lines"][[#]],
       "originalEndpoints",
       topo["lines"][[#, "endpoints"]]
       ] & /@ shrunkLines;
   vertexRepresentativeMap[topo["vertexIds"], originalPairs]
   ];


(* J 中的 shrunk packs 决定目标 sector；不能只读取 source topology 的当前 endpoints。 *)
masslessCoincidentAntisymmetricIntegralQ[
   topo_Association,
   int : J[aList_, linePacks_, ispList_]
   ] := Module[
   {lines = topo["lines"], repMap},
   repMap = integralTargetVertexRepresentativeMap[topo, int];
   AnyTrue[
    Range[Length[lines]],
    Function[e,
     Module[{originalEndpoints, targetEndpoints},
      originalEndpoints = Lookup[
        lines[[e]],
        "originalEndpoints",
        lines[[e]]["endpoints"]
        ];
      targetEndpoints = Lookup[repMap, originalEndpoints];
      actualLinePackType[topo, e, linePacks[[e]]] === "masslessFull" &&
       SameQ @@ targetEndpoints &&
       linePacks[[e, First[linePackNPositions[lines[[e]], "masslessFull"]]]] === 1
      ]
     ]
    ]
   ];


applyMasslessEndpointCanonical[expr_, topo_Association] := Expand[
   expr /. (int_J /; masslessCoincidentAntisymmetricIntegralQ[topo, int]) :> 0
   ];


canonicalizeMassiveCoincidentIntegral[
   topo_Association,
   int : J[aList_, linePacks_, ispList_]
   ] := Module[{lines = topo["lines"], repMap, newPacks = linePacks, originalEndpoints, targetEndpoints},
   repMap = integralTargetVertexRepresentativeMap[topo, int];
   Do[
    originalEndpoints = Lookup[lines[[e]], "originalEndpoints", lines[[e]]["endpoints"]];
    targetEndpoints = Lookup[repMap, originalEndpoints];
    If[
     actualLinePackType[topo, e, linePacks[[e]]] === "massiveFull" &&
      SameQ @@ targetEndpoints &&
      linePacks[[e, linePackNPositions[lines[[e]], "massiveFull"]]] === {1, 0},
     newPacks[[e, linePackNPositions[lines[[e]], "massiveFull"]]] =
       Reverse[linePacks[[e, linePackNPositions[lines[[e]], "massiveFull"]]]]
     ],
    {e, Length[lines]}
    ];
   J[aList, newPacks, ispList]
   ];


applyMassiveCoincidentCanonical[expr_, topo_Association] := Expand[
   expr /. int_J :> canonicalizeMassiveCoincidentIntegral[topo, int]
   ];


applySeedCanonical[expr_, topo_Association] :=
  symmetry[
   applyMasslessEndpointCanonical[
    applyMassiveCoincidentCanonical[applyEOM[expr, topo], topo],
    topo
    ],
   topo
   ];


forbiddenNDataForIntegral[topo_Association, J[aList_, linePacks_, ispList_]] := Module[
   {lines = topo["lines"], issues = {}, packType, pack},
   Do[
    pack = linePacks[[e]];
    packType = actualLinePackType[topo, e, pack];
    Switch[packType,
     "massiveFull",
     Do[
      If[IntegerQ[pack[[linePackNPositions[lines[[e]], packType][[endpointSlot]]]]] &&
        pack[[linePackNPositions[lines[[e]], packType][[endpointSlot]]]] >= 2,
       AppendTo[issues, <|"lineIndex" -> e, "packType" -> packType, "endpointSlot" -> endpointSlot,
         "nValue" -> pack[[linePackNPositions[lines[[e]], packType][[endpointSlot]]]]|>]
       ],
     {endpointSlot, 2}
     ],
     "masslessFull",
     If[IntegerQ[pack[[First[linePackNPositions[lines[[e]], packType]]]]] &&
       ! MemberQ[{0, 1}, pack[[First[linePackNPositions[lines[[e]], packType]]]]],
      AppendTo[issues, <|"lineIndex" -> e, "packType" -> packType,
        "nValue" -> pack[[First[linePackNPositions[lines[[e]], packType]]]]|>]
      ],
     "massiveCross",
     Do[
      If[IntegerQ[pack[[linePackNPositions[lines[[e]], packType][[endpointSlot]]]]] &&
        pack[[linePackNPositions[lines[[e]], packType][[endpointSlot]]]] >= 2,
       AppendTo[issues, <|"lineIndex" -> e, "packType" -> packType, "endpointSlot" -> endpointSlot,
         "nValue" -> pack[[linePackNPositions[lines[[e]], packType][[endpointSlot]]]]|>]
       ],
      {endpointSlot, 2}
      ],
     _, Null
     ],
    {e, Length[lines]}
    ];
   issues
   ];


forbiddenNData[topo_Association, expr_] := DeleteCases[
   Flatten[forbiddenNDataForIntegral[topo, #] & /@ DeleteDuplicates[Cases[expr, _J, {0, Infinity}]]],
   Null
   ];


containsForbiddenNQ[topo_Association, expr_] := Length[forbiddenNData[topo, expr]] > 0;


assertNoForbiddenN::badn = "表达式仍含 forbidden n 指标：`1`。";
assertNoForbiddenN[expr_, topo_Association] := Module[{bad = forbiddenNData[topo, expr]},
   If[bad === {}, expr, Message[assertNoForbiddenN::badn, bad]; $Failed]
   ];

(* ::Chapter:: *)
(*用户输入的积分族对称性*)

(* sp 的 Orderless 只处理标量积交换性；本章只应用用户确认物理条件后的积分族规则。 *)
repSymmetry0[topo_Association] := Lookup[topo, "symmetryRules", {}];


symmetry::badrules = "symmetryRules 必须是 Rule/RuleDelayed 的列表。";


(* tadpole symmetry 只识别可由单个 loop momentum reversal 证明的 self-loop。*)
tadpoleLoopReversalData[topo_Association] := Module[
   {loops = Lookup[topo, "loopMomenta", {}], lines = Lookup[topo, "lines", {}], candidates, momentum, endpoints, loopIndex, exclusiveLoopQ, data = {}},
   Do[
    endpoints = Lookup[lines[[e]], "originalEndpoints", Lookup[lines[[e]], "endpoints", {}]];
    If[SameQ @@ endpoints,
     momentum = Expand[Lookup[lines[[e]], "momentum", 0]];
     candidates = Select[
       Range[Length[loops]],
       Function[l,
        With[{coefficient = Coefficient[momentum, loops[[l]]]},
         coefficient =!= 0 && Expand[momentum - coefficient loops[[l]]] === 0
         ]
        ]
       ];
     If[Length[candidates] === 1,
      loopIndex = First[candidates];
      exclusiveLoopQ = And @@ Table[
         otherLine === e || Coefficient[Expand[Lookup[lines[[otherLine]], "momentum", 0]], loops[[loopIndex]]] === 0,
         {otherLine, Length[lines]}
         ];
      AppendTo[data, <|
        "lineIndex" -> e,
        "lineId" -> Lookup[lines[[e]], "id", e],
        "packType" -> Lookup[lines[[e]], "packType", Missing["packType"]],
        "endpoints" -> endpoints,
        "loopIndex" -> loopIndex,
        "loopMomentum" -> loops[[loopIndex]],
        "exclusiveLoopQ" -> exclusiveLoopQ
        |>]
      ]
     ],
    {e, Length[lines]}
    ];
   data
   ];


reverseLoopScalarProductExpr[expr_, loopIndex_Integer] := Expand[
   expr /. {
     HoldPattern[qq[i_Integer, j_Integer]] :> (-1)^Count[{i, j}, loopIndex] qq[i, j],
     HoldPattern[qk[i_Integer, j_Integer]] :> (-1)^Count[{i}, loopIndex] qk[i, j]
     }
   ];


tadpoleISPParity[topo_Association, loopIndex_Integer, ispIndex_Integer] := Module[
   {exprs = normalizeISPExprs[topo], expr, reversed},
   If[ispIndex < 1 || ispIndex > Length[exprs], Return[0]];
   expr = exprs[[ispIndex]];
   reversed = reverseLoopScalarProductExpr[expr, loopIndex];
   Which[
    Expand[reversed - expr] === 0, 1,
    Expand[reversed + expr] === 0, -1,
    True, 0
    ]
   ];


tadpoleOddISPIntegralQ[topo_Association, J[_, _, ispList_]] := Module[{loopData, oddPowers},
   loopData = Select[
     tadpoleLoopReversalData[topo],
     MemberQ[{"massiveFull", "masslessFull"}, #["packType"]] && TrueQ[#["exclusiveLoopQ"]] &
     ];
   AnyTrue[
    loopData,
    Function[data,
     oddPowers = Select[
       Table[
        If[tadpoleISPParity[topo, data["loopIndex"], j] === -1, ispList[[j]], Nothing],
        {j, Length[ispList]}
        ],
       IntegerQ
       ];
     oddPowers =!= {} && And @@ (IntegerQ /@ oddPowers) && OddQ[Total[oddPowers]]
     ]
    ]
   ];


tadpoleMassiveSwapNeededQ[topo_Association, lineIndex_Integer, J[_, linePacks_, _]] :=
  lineIndex <= Length[linePacks] &&
    Lookup[topo["lines"][[lineIndex]], "packType", Missing["packType"]] === "massiveFull" &&
    linePacks[[lineIndex, linePackNPositions[topo["lines"][[lineIndex]], "massiveFull"]]] === {1, 0};


tadpoleMasslessZeroQ[topo_Association, lineIndex_Integer, J[_, linePacks_, _]] :=
  lineIndex <= Length[linePacks] &&
    Lookup[topo["lines"][[lineIndex]], "packType", Missing["packType"]] === "masslessFull" &&
    linePacks[[lineIndex, First[linePackNPositions[topo["lines"][[lineIndex]], "masslessFull"]]]] === 1;


tadpoleSwapLinePack[topo_Association, J[aList_, linePacks_, ispList_], lineIndex_Integer] := Module[
   {newPacks = linePacks, nPositions},
   nPositions = linePackNPositions[topo["lines"][[lineIndex]], "massiveFull"];
   newPacks[[lineIndex, nPositions]] = Reverse[newPacks[[lineIndex, nPositions]]];
   J[aList, newPacks, ispList]
   ];


tadpoleSymmetryRules0[topo_Association] := Module[{data, rules = {}},
   data = tadpoleLoopReversalData[topo];
   AppendTo[rules,
    HoldPattern[(int_J /; tadpoleOddISPIntegralQ[topo, int])] :> 0
    ];
   Do[
    Switch[dataItem["packType"],
     "massiveFull",
     With[{lineIndex = dataItem["lineIndex"]},
      AppendTo[rules,
       HoldPattern[(int_J /; tadpoleMassiveSwapNeededQ[topo, lineIndex, int])] :>
        tadpoleSwapLinePack[topo, int, lineIndex]
       ]
      ],
     "masslessFull",
     With[{lineIndex = dataItem["lineIndex"]},
      AppendTo[rules,
       HoldPattern[(int_J /; tadpoleMasslessZeroQ[topo, lineIndex, int])] :> 0
       ]
      ],
     _, Null
     ],
    {dataItem, data}
    ];
   rules
   ];


effectiveSymmetryRules0[topo_Association] := Module[{userRules = repSymmetry0[topo]},
   If[! ListQ[userRules], Return[userRules]];
   DeleteDuplicates@Join[tadpoleSymmetryRules0[topo], userRules]
   ];


tadpoleSymmetryData[topo_Association] := Module[{data = tadpoleLoopReversalData[topo]},
   <|
    "status" -> "generated",
    "loopReversalData" -> data,
    (* 空 tadpole 集合是合法结果；显式返回空列，避免 Lookup[{}] 产生 Missing。 *)
    "massiveFullLineIndices" -> Lookup[Select[data, #["packType"] === "massiveFull" &], "lineIndex", {}],
    "masslessFullLineIndices" -> Lookup[Select[data, #["packType"] === "masslessFull" &], "lineIndex", {}],
    "automaticRuleCount" -> Length[tadpoleSymmetryRules0[topo]],
    "automaticRules" -> tadpoleSymmetryRules0[topo],
    "userRuleCount" -> Length[repSymmetry0[topo]],
    "effectiveRuleCount" -> Length[effectiveSymmetryRules0[topo]]
    |>
   ];


symmetry[expr_, topo_Association] := Module[
   {rules = effectiveSymmetryRules0[topo], internalExpr, result, publicTimeOnlyQ},
   If[
    ! ListQ[rules] || ! And @@ (validDiscreteReplacementRuleQ /@ rules),
    Message[symmetry::badrules];
    Return[$Failed]
    ];
   publicTimeOnlyQ = Lookup[topo, "ibpMode", "full"] === "timeOnly" &&
     ! FreeQ[expr, J[_String, _List, _List]];
   internalExpr = If[
     publicTimeOnlyQ,
     dsTimeOnlyExpressionToInternal020[expr, topo],
     expr
     ];
   If[internalExpr === $Failed, Return[$Failed]];
   result = internalExpr /. rules;
   If[publicTimeOnlyQ, dsTimeOnlyExpressionToPublic020[result, topo], result]
   ];

(* ::Chapter:: *)
(*标量积与 ISP 覆盖性*)

(* 本章将 Q_e^2 展开到 qq/qk/kk，并记录 propagator z_e 与用户 ISP 的结构信息。
   默认不做 MatrixRank/Solve 等解析检查；覆盖性验证应在单独 numeric/specialized 脚本中代数赋值后执行。 *)


scalarProductVariables[topo_Association] := Module[{nL = topo["nL"], nK = topo["nK"]},
   Join[
    Flatten[Table[qq[i, j], {i, nL}, {j, i, nL}]],
    Flatten[Table[qk[i, j], {i, nL}, {j, nK}]]
    ]
   ];


externalInvariantVariables[topo_Association] := Module[{nK = topo["nK"]},
   Flatten[Table[kk[i, j], {i, nK}, {j, i, nK}]]
   ];


externalLegEnergyVariables[topo_Association] := DeleteDuplicates[
   Variables[vertexExternalEnergy[topo, #] & /@ activeAVertexIds[topo]]
   ];


(* 将两个线性动量表达式的点积展开为 qq/qk/kk。 *)
expandDotExpr[p_, r_, topo_Association] := Module[
   {loops = topo["loopMomenta"], exts = topo["effectiveLoopExternalMomenta"],
    nL = topo["nL"], nK = topo["nK"], lcP, lcR, ecP, ecR, result = 0},
   lcP = Coefficient[p, #] & /@ loops;
   lcR = Coefficient[r, #] & /@ loops;
   ecP = Coefficient[p, #] & /@ exts;
   ecR = Coefficient[r, #] & /@ exts;
   Do[
    result += lcP[[i]] lcR[[j]] qqSym[i, j],
    {i, nL}, {j, nL}
    ];
   Do[
    result += lcP[[i]] ecR[[j]] qk[i, j] + ecP[[j]] lcR[[i]] qk[i, j],
    {i, nL}, {j, nK}
    ];
   Do[
    result += ecP[[i]] ecR[[j]] kkSym[i, j],
    {i, nK}, {j, nK}
    ];
   Expand[result]
   ];



externalInvariantSymbolName[i_Integer, j_Integer, count_Integer : 0] := ToExpression[
   "s" <> coordinateIndexString[Min[i, j], count] <> coordinateIndexString[Max[i, j], count]
   ];


defaultLoopKinematicRulesForTopology[topo_Association] := Module[
   {exts = Lookup[topo, "effectiveLoopExternalMomenta", {}], nK = Lookup[topo, "nK", Length[Lookup[topo, "effectiveLoopExternalMomenta", {}]]]},
   Flatten[Table[sp[exts[[i]], exts[[j]]] -> externalInvariantSymbolName[i, j, nK], {i, nK}, {j, i, nK}]]
   ];


normalizeLoopKinematicRulesForTopology[Automatic, topo_Association] := defaultLoopKinematicRulesForTopology[topo];
normalizeLoopKinematicRulesForTopology[rules_Association, topo_Association] := normalizeLoopKinematicRulesForTopology[Normal[rules], topo];
normalizeLoopKinematicRulesForTopology[rules_List, topo_Association] := Module[
   {defaults = defaultLoopKinematicRulesForTopology[topo], validRules},
   validRules = Select[rules, validReplacementRuleQ];
   Normal[Association[Join[defaults, validRules]]]
   ];
normalizeLoopKinematicRulesForTopology[_, topo_Association] := defaultLoopKinematicRulesForTopology[topo];


loopKinematicInternalToUserRules[topo_Association] := Module[
   {rules = Lookup[topo, "resolvedLoopKinematicRules", defaultLoopKinematicRulesForTopology[topo]]},
   rules /. (Rule | RuleDelayed)[lhs_, rhs_] :> Rule[scalarProductSPInputToInternal[lhs, topo], rhs]
   ];


loopKinematicUserToInternalRules[topo_Association] := Module[
   {rules = loopKinematicInternalToUserRules[topo]},
   Cases[rules, Rule[lhs_, rhs_] :> Rule[rhs, lhs]]
   ];


scalarProductInternalToUser[expr_, topo_Association] := Module[
   {loops = topo["loopMomenta"], exts = topo["effectiveLoopExternalMomenta"]},
   Expand[expr /. Join[
      loopKinematicInternalToUserRules[topo],
      {
       HoldPattern[qq[i_Integer, j_Integer]] :> sp[loops[[i]], loops[[j]]],
       HoldPattern[qk[i_Integer, j_Integer]] :> sp[loops[[i]], exts[[j]]]
       }
      ]]
   ];


scalarProductExpressionValidQ[expr_, topo_Association] := Module[
   {internal = scalarProductInputToInternal[expr, topo], allowedVars},
   allowedVars = Join[scalarProductVariables[topo], externalInvariantVariables[topo]];
   FreeQ[internal, sp] && SubsetQ[allowedVars, Variables[internal]]
   ];


normalizeISPExprs[topo_Association] := scalarProductInputToInternal[Lookup[#, "expr"], topo] & /@ topo["ispData"];


normalizeCoefficientRulesForTopology[rules_List, topo_Association] := normalizeNumericRulesForTopology[rules, topo];
normalizeCoefficientRulesForTopology[rules_, topo_Association] := rules;


normalizeCoefficientRulesForLinearData[rules_, linearData_Association] := Module[
   {topo = Lookup[linearData, "topology", Missing["NoTopology"]]},
   If[parsedTopologyQ[topo], normalizeCoefficientRulesForTopology[rules, topo], rules]
   ];


userCoefficientRulesForLinearData[rules_, linearData_Association] := Module[
   {topo = Lookup[linearData, "topology", Missing["NoTopology"]]},
   If[parsedTopologyQ[topo] && ListQ[rules],
    rules /. (Rule | RuleDelayed)[lhs_, rhs_] :> Rule[scalarProductInternalToUser[lhs, topo], rhs],
    rules
    ]
   ];


loopDependentLineIndices[topo_Association] := Flatten@Position[
   Lookup[topo, "loopCoeffMatrix", {}],
   row_List /; AnyTrue[row, ! zeroQ[#] &]
   ];


expandZList[topo_Association] := Module[{indices = loopDependentLineIndices[topo]},
   If[
    indices === {},
    {},
    expandDotExpr[#, #, topo] & /@ Lookup[topo["lines"][[indices]], "momentum", {}]
    ]
   ];


coefficientMatrix[exprs_List, vars_List] := Table[
   Coefficient[Expand[exprs[[r]]], vars[[c]]],
   {r, Length[exprs]}, {c, Length[vars]}
   ];


linearMomentumExpressionData[expr_, basis_List] := Module[
   {atoms, expandedBasis, exprCoefficients, basisMatrix, coefficientVariables,
    equations, solutions, coefficients, residual, atomResiduals},
   atoms = ds016MomentumAtoms[Join[{expr}, basis]];
   expandedBasis = Expand /@ basis;
   exprCoefficients = Coefficient[Expand[expr], #] & /@ atoms;
   atomResiduals = Function[basisExpression,
       Expand[basisExpression - Total[(Coefficient[basisExpression, #] & /@ atoms) atoms]]
       ] /@ expandedBasis;
   If[! And @@ (ds016ZeroQ /@ Join[
        atomResiduals,
        {Expand[expr - Total[exprCoefficients atoms]]}
        ]),
    Return[<|"expr" -> expr, "basis" -> basis, "coefficients" -> ConstantArray[0, Length[basis]],
      "residual" -> expr, "linearQ" -> False|>]
    ];
   If[basis === {},
    residual = Expand[expr];
    Return[<|"expr" -> expr, "basis" -> basis, "coefficients" -> {},
      "residual" -> residual, "linearQ" -> ds016ZeroQ[residual]|>]
    ];
   basisMatrix = Table[Coefficient[expandedBasis[[i]], atoms[[j]]], {i, Length[basis]}, {j, Length[atoms]}];
   coefficientVariables = Array[Unique["momentumCoefficient$"] &, Length[basis]];
   equations = Thread[coefficientVariables . basisMatrix == exprCoefficients];
   solutions = Quiet[Solve[equations, coefficientVariables]];
   If[solutions === {},
    Return[<|"expr" -> expr, "basis" -> basis, "coefficients" -> ConstantArray[0, Length[basis]],
      "residual" -> expr, "linearQ" -> False|>]
   ];
   coefficients = coefficientVariables /. First[solutions];
   coefficients = coefficients /. Thread[coefficientVariables -> 0];
   residual = Expand[expr - Total[coefficients basis]];
   <|"expr" -> expr, "basis" -> basis, "coefficients" -> coefficients,
    "residual" -> residual, "linearQ" -> ds016ZeroQ[residual]|>
   ];


linearMomentumExpressionQ[expr_, basis_List] := TrueQ[linearMomentumExpressionData[expr, basis]["linearQ"]];


lineMomentumLinearityIssues[topo_Association] := Module[
   {basis = Lookup[topo, "momentumDecompositionBasis",
      Join[topo["loopMomenta"], topo["effectiveLoopExternalMomenta"], Lookup[topo, "independentExternalMomenta", {}]]]},
   DeleteCases[
    MapIndexed[
     With[{data = linearMomentumExpressionData[Lookup[#1, "momentum", 0], basis]},
       If[TrueQ[data["linearQ"]],
        Nothing,
        <|"lineIndex" -> First[#2], "lineId" -> Lookup[#1, "id", Missing["id"]], "momentum" -> Lookup[#1, "momentum", Missing["momentum"]], "residual" -> data["residual"], "basis" -> basis|>
        ]
       ] &,
     topo["lines"]
     ],
    Nothing
    ]
   ];


scalarProductArgumentLinearityIssues[topo_Association] := Module[
   {basis = Join[topo["loopMomenta"], topo["effectiveLoopExternalMomenta"]], rawISPExprs, spData},
   rawISPExprs = Lookup[#, "expr"] & /@ topo["ispData"];
   DeleteCases[
    Flatten[
     MapIndexed[
      Function[{expr, pos},
       MapIndexed[
        Function[{pair, pairPos},
         DeleteCases[
          MapIndexed[
           With[{data = linearMomentumExpressionData[#1, basis]},
             If[TrueQ[data["linearQ"]],
              Nothing,
              <|"ispPosition" -> First[pos], "ispName" -> Lookup[topo["ispData"][[First[pos]]], "name", Missing["name"]], "spPosition" -> First[pairPos], "argumentSlot" -> First[#2], "argument" -> #1, "residual" -> data["residual"], "basis" -> basis|>
              ]
             ] &,
           pair
           ],
          Nothing
          ]
         ],
        Cases[expr, HoldPattern[sp[p_, r_]] :> {p, r}, {0, Infinity}]
        ]
       ],
      rawISPExprs
      ],
     2
     ],
    Nothing
    ]
   ];


externalLegEnergySPArgumentIssues[expr_, topo_Association] := Module[
   {basis = Lookup[topo, "momentumDecompositionBasis",
      Join[topo["loopMomenta"], topo["effectiveLoopExternalMomenta"], Lookup[topo, "independentExternalMomenta", {}]]],
    nL = topo["nL"], pairs},
   pairs = Cases[expr, HoldPattern[sp[p_, r_]] :> {p, r}, {0, Infinity}];
   DeleteCases[
    Flatten[
     MapIndexed[
      Function[{pair, pairPos},
       MapIndexed[
        Function[{arg, argPos},
         Module[{data = linearMomentumExpressionData[arg, basis], loopCoeffs},
          loopCoeffs = Take[data["coefficients"], nL];
          If[TrueQ[data["linearQ"]] && TrueQ[loopCoeffs === ConstantArray[0, nL]],
           Nothing,
           <|
            "spPosition" -> First[pairPos],
            "argumentSlot" -> First[argPos],
            "argument" -> arg,
            "residual" -> data["residual"],
            "loopCoefficients" -> loopCoeffs,
            "basis" -> basis
            |>
           ]
          ]
         ],
        pair
        ]
       ],
      pairs
      ],
     1
     ],
    Nothing
    ]
   ];


externalLegEnergyMomentumDependenceIssues[topo_Association] := Module[
   {declared = Join[topo["loopMomenta"], topo["effectiveLoopExternalMomenta"], Lookup[topo, "independentExternalMomenta", {}]],
    loopSPVars = scalarProductVariables[topo], vertices = activeAVertexIds[topo]},
   DeleteCases[
    Table[
     Module[{raw = rawVertexExternalEnergy[topo, vertex], rawNoSP, directMomenta, spArgIssues, internal, loopSPUsed},
      rawNoSP = raw /. HoldPattern[sp[_, _]] -> 0;
      directMomenta = DeleteDuplicates@Cases[rawNoSP, sym_Symbol /; MemberQ[declared, sym], {0, Infinity}];
      spArgIssues = externalLegEnergySPArgumentIssues[raw, topo];
      internal = scalarProductInputToInternal[raw, topo];
      loopSPUsed = Intersection[Variables[internal], loopSPVars];
      If[directMomenta === {} && spArgIssues === {} && loopSPUsed === {},
       Nothing,
       <|
        "vertexId" -> vertex,
        "rawExternalLegEnergy" -> raw,
        "userExternalLegEnergy" -> scalarProductInternalToUser[internal, topo],
        "directMomentumSymbols" -> directMomenta,
        "spArgumentIssues" -> spArgIssues,
        "loopScalarProducts" -> scalarProductInternalToUser[#, topo] & /@ loopSPUsed,
        "comment" -> "vertices.externalLegEnergy values are scalar time-phase energies; use ssij when tied to loopExternalMomenta, sEi when tied to an independentExternalMomenta magnitude, or an independent scalar otherwise"
        |>
       ]
      ],
     {vertex, vertices}
     ],
    Nothing
    ]
   ];



(* 标量积坐标规则只依赖动量路由、ISP 和外部不变量命名。
   缓存键排除 sector、顶点和 seed 配置，使同一 family 的 shrink sectors 共享一次 LinearSolve。 *)
$scalarProductRuleCache = <||>;


scalarProductRuleCacheKey[topo_Association] := With[
   {
    nL = topo["nL"],
    nK = topo["nK"],
    nE = topo["nE"],
    loopMomenta = topo["loopMomenta"],
    effectiveLoopExternalMomenta = topo["effectiveLoopExternalMomenta"],
    lineMomenta = Lookup[topo["lines"], "momentum"],
    ispExprs = Lookup[topo["ispData"], "expr"],
    resolvedLoopKinematicRules = Lookup[topo, "resolvedLoopKinematicRules", {}]
    },
   HoldComplete[nL, nK, nE, loopMomenta, effectiveLoopExternalMomenta, lineMomenta, ispExprs, resolvedLoopKinematicRules]
   ];


clearScalarProductRuleCache[] := ($scalarProductRuleCache = <||>; Null);


scalarProductRuleCacheReport[] := <|
   "entryCount" -> Length[$scalarProductRuleCache],
   "hashes" -> Keys[$scalarProductRuleCache]
   |>;


(* 小规模线性规则核心：直接作为 ISP 给出的标量积会被保留，不强行改写成 rho。
   本 package 假设用户输入的 z_e 与直接 ISP 已经构成闭合坐标；这里不把 dS 图
   默认当成 overcomplete propagator family 处理，避免改变用户定义的函数族。 *)
makeScalarProductRulesUncached[topo_Association] := Module[
   {spVars, zVars, zExprs, rawISPExprs, ispExprs, ispVars, invalidISPPositions, unsupportedISPExprs,
    coordExprs, coordVars, mat, const, rhs, solVec},
   spVars = scalarProductVariables[topo];
   zVars = z /@ loopDependentLineIndices[topo];
   zExprs = expandZList[topo];
   rawISPExprs = Lookup[#, "expr"] & /@ topo["ispData"];
   ispExprs = scalarProductInputToInternal[#, topo] & /@ rawISPExprs;
   ispVars = Table[rho[j], {j, Length[ispExprs]}];
   invalidISPPositions = Flatten@Position[Map[scalarProductExpressionValidQ[#, topo] &, rawISPExprs], False];
   unsupportedISPExprs = rawISPExprs[[invalidISPPositions]];
   If[unsupportedISPExprs =!= {},
    Return[<|
      "status" -> "notComputed",
      "reason" -> "ISP expressions must be scalar products written as sp[p,r] or linear combinations of them",
      "repSP2Z" -> Missing["NotComputedUnsupportedISP"],
      "solveVars" -> Missing["NotComputedUnsupportedISP"],
      "preservedISPVars" -> ispVars,
      "ispCoordinateVars" -> ispVars,
      "unsupportedISPExprs" -> unsupportedISPExprs
      |>]
    ];
   coordExprs = Join[zExprs, ispExprs];
   coordVars = Join[zVars, ispVars];
   If[Length[coordExprs] =!= Length[spVars],
    Return[<|
      "status" -> "notComputed",
      "reason" -> "z plus ISP equation count does not match scalar-product count",
      "repSP2Z" -> Missing["NotComputedCountMismatch"],
      "solveVars" -> spVars,
      "preservedISPVars" -> ispVars,
      "ispCoordinateVars" -> ispVars,
      "zCount" -> Length[zExprs],
      "ispCount" -> Length[ispExprs],
      "scalarProductCount" -> Length[spVars]
      |>]
    ];
   mat = coefficientMatrix[coordExprs, spVars];
   const = coordExprs /. Thread[spVars -> 0];
   rhs = coordVars - const;
   solVec = Check[LinearSolve[mat, rhs], $Failed];
   If[solVec === $Failed,
    Return[<|
      "status" -> "notComputed",
      "reason" -> "LinearSolve failed",
      "repSP2Z" -> Missing["LinearSolveFailed"],
      "solveVars" -> spVars,
      "preservedISPVars" -> ispVars,
      "ispCoordinateVars" -> ispVars,
      "zCount" -> Length[zExprs],
      "ispCount" -> Length[ispExprs],
      "scalarProductCount" -> Length[spVars]
      |>]
    ];
   <|
    "status" -> "computed",
    "repSP2Z" -> Thread[spVars -> (Expand /@ solVec)],
    "userRepSP2Z" -> (Rule[scalarProductInternalToUser[#[[1]], topo], scalarProductInternalToUser[#[[2]], topo]] & /@ Thread[spVars -> (Expand /@ solVec)]),
    "repZ2SP" -> Thread[zVars -> zExprs],
    "userRepZ2SP" -> (Rule[#[[1]], scalarProductInternalToUser[#[[2]], topo]] & /@ Thread[zVars -> zExprs]),
    "repISP2SP" -> Thread[ispVars -> ispExprs],
    "userRepISP2SP" -> (Rule[#[[1]], scalarProductInternalToUser[#[[2]], topo]] & /@ Thread[ispVars -> ispExprs]),
    "solveVars" -> spVars,
    "preservedISPVars" -> ispVars,
    "ispCoordinateVars" -> ispVars,
    "zCount" -> Length[zExprs],
    "ispCount" -> Length[ispExprs],
    "scalarProductCount" -> Length[spVars]
    |>
   ];


(* 公开入口在命中缓存时直接返回；SHA-256 命中后仍用完整键 SameQ 防止哈希冲突。 *)
makeScalarProductRules[topo_Association] := Module[
   {key, hash, entry, result},
   key = scalarProductRuleCacheKey[topo];
   hash = Hash[key, "SHA256"];
   entry = Lookup[$scalarProductRuleCache, hash, Missing["NotCached"]];
   If[AssociationQ[entry] && SameQ[entry["key"], key],
    Return[entry["value"]]
    ];
   result = makeScalarProductRulesUncached[topo];
   AssociateTo[$scalarProductRuleCache, hash -> <|"key" -> key, "value" -> result|>];
   result
   ];


makeScalarProductData[topo_Association] := Module[
   {spVars, zVars, zExprs, rawISPExprs, ispExprs, ispSymbols, invalidISPPositions, unsupportedISPExprs,
    coordExprs, nSP, structuralNeededISPCount, structuralCountQ, coordinateCountQ},
   spVars = scalarProductVariables[topo];
   zVars = z /@ loopDependentLineIndices[topo];
   zExprs = expandZList[topo];
   rawISPExprs = Lookup[#, "expr"] & /@ topo["ispData"];
   ispExprs = scalarProductInputToInternal[#, topo] & /@ rawISPExprs;
   ispSymbols = Table[rho[j], {j, Length[ispExprs]}];
   invalidISPPositions = Flatten@Position[Map[scalarProductExpressionValidQ[#, topo] &, rawISPExprs], False];
   unsupportedISPExprs = rawISPExprs[[invalidISPPositions]];
   coordExprs = Join[zExprs, ispExprs];
   nSP = Length[spVars];
   structuralNeededISPCount = Max[0, nSP - Length[zExprs]];
   structuralCountQ = Length[ispExprs] >= structuralNeededISPCount;
   coordinateCountQ = Length[coordExprs] === nSP;
   <|
    "scalarProducts" -> (scalarProductInternalToUser[#, topo] & /@ spVars),
    "internalScalarProducts" -> spVars,
    "externalInvariants" -> (scalarProductInternalToUser[#, topo] & /@ externalInvariantVariables[topo]),
    "internalExternalInvariants" -> externalInvariantVariables[topo],
    "loopKinematicNamingReport" -> loopKinematicNamingReport[topo],
    "zVars" -> zVars,
    "zExprs" -> (scalarProductInternalToUser[#, topo] & /@ zExprs),
    "internalZExprs" -> zExprs,
    "ispSymbols" -> ispSymbols,
    "ispExprs" -> (scalarProductInternalToUser[#, topo] & /@ ispExprs),
    "internalISPExprs" -> ispExprs,
    "directISPVars" -> ispSymbols,
    "internalDirectISPVars" -> ispSymbols,
    "unsupportedISPExprs" -> unsupportedISPExprs,
    "solveVars" -> spVars,
    "zCount" -> Length[zExprs],
    "spCount" -> nSP,
    "ispCount" -> Length[ispExprs],
    "directISPCount" -> Length[ispExprs],
    "nonISPScalarProductCount" -> nSP - Length[ispExprs],
    "structuralNeededISPCount" -> structuralNeededISPCount,
    "structuralCountQ" -> structuralCountQ,
    "coordinateCountQ" -> coordinateCountQ,
    "coverageQ" -> Missing["NotCheckedSeedOnly"],
    "independentQ" -> Missing["NotCheckedSeedOnly"],
    "repSP2Z" -> Missing["NotComputedSeedOnly"]
    |>
   ];

(* ::Chapter:: *)
(*IBP 生成元枚举*)

(* 本章生成用户圈动量坐标中的 IBP 算子标签，并把它们编译为内部 normalized
   routing 上的 primitive generators。完整集合包括
   momentum IBP: d/dq_l dot q_m, d/dq_l dot k_j。 *)


(* reference routing 满足 qNormalized=R.qInput+r。传播子和 ISP 已写在
   qNormalized 中，因此 d/dqInput_l 必须乘同一个 Jacobian R[[All,l]]。 *)
affineLoopDerivativeMatrix[topo_Association] := Module[{nL, routing, matrix},
   nL = topo["nL"];
   routing = Lookup[topo, "loopMomentumRoutingAudit", <||>];
   matrix = Lookup[routing, "referenceLoopMatrix", {}];
   If[MatrixQ[matrix] && Dimensions[matrix] === {nL, nL}, matrix, IdentityMatrix[nL]]
   ];


(* 用户的 qInput_m 用初始化时保存的仿射反变换写回 normalized variables。 *)
affineInputLoopVector[topo_Association, loopIndex_Integer] := Expand[
   topo["loopMomenta"][[loopIndex]] /.
    Lookup[topo, "loopMomentumInputToNormalizedRules", {}]
   ];


momentumPrimitiveTerm[coefficient_, dLoop_Integer, vectorType_String, vectorIndex_, vector_] := <|
   "coefficient" -> Expand[coefficient],
   "generator" -> <|
     "type" -> "momentum",
     "dLoop" -> dLoop,
     "vectorType" -> vectorType,
     "vectorIndex" -> vectorIndex,
     "vector" -> Expand[vector]
     |>
   |>;


(* 输入：用户生成元的导数圈变量与乘法向量；输出：保留用户标签并附带内部 primitive terms。
   qInput 的平移余项作为一个外向量表达式整体保留，避免人为增加公开 external basis。 *)
makeAffineMomentumGenerator[
   topo_Association, dLoop_Integer, "loop", vectorIndex_Integer
   ] := Module[
   {loops, derivativeMatrix, derivativeCoefficients, inputVector,
    loopCoefficients, externalRemainder, loopTerms, externalTerms},
   loops = topo["loopMomenta"];
   derivativeMatrix = affineLoopDerivativeMatrix[topo];
   derivativeCoefficients = derivativeMatrix[[All, dLoop]];
   inputVector = affineInputLoopVector[topo, vectorIndex];
   loopCoefficients = Coefficient[inputVector, #] & /@ loops;
   externalRemainder = Expand[inputVector - loopCoefficients . loops];
   loopTerms = Flatten@Table[
      If[zeroQ[derivativeCoefficients[[a]] loopCoefficients[[b]]],
       Nothing,
       momentumPrimitiveTerm[
        derivativeCoefficients[[a]] loopCoefficients[[b]], a, "loop", b, loops[[b]]
        ]
       ],
      {a, topo["nL"]}, {b, topo["nL"]}
      ];
   externalTerms = Table[
     If[zeroQ[derivativeCoefficients[[a]]] || zeroQ[externalRemainder],
      Nothing,
      momentumPrimitiveTerm[
       derivativeCoefficients[[a]], a, "externalExpression",
       Missing["AffineShiftVector"], externalRemainder
       ]
      ],
     {a, topo["nL"]}
     ];
   <|
    "type" -> "momentum",
    "dLoop" -> dLoop,
    "vectorType" -> "loop",
    "vectorIndex" -> vectorIndex,
    "vector" -> loops[[vectorIndex]],
    "normalizedPrimitiveTerms" -> Join[loopTerms, externalTerms]
    |>
   ];


makeAffineMomentumGenerator[
   topo_Association, dLoop_Integer, "external", vectorIndex_Integer
   ] := Module[{derivativeCoefficients, vector, terms},
   derivativeCoefficients = affineLoopDerivativeMatrix[topo][[All, dLoop]];
   vector = topo["effectiveLoopExternalMomenta"][[vectorIndex]];
   terms = Table[
     If[zeroQ[derivativeCoefficients[[a]]],
      Nothing,
      momentumPrimitiveTerm[
       derivativeCoefficients[[a]], a, "external", vectorIndex, vector
       ]
      ],
     {a, topo["nL"]}
     ];
   <|
    "type" -> "momentum",
    "dLoop" -> dLoop,
    "vectorType" -> "external",
    "vectorIndex" -> vectorIndex,
    "vector" -> vector,
    "normalizedPrimitiveTerms" -> terms
    |>
   ];


makeIBPGenerators[topo_Association] := Module[
   {timeGenerators, loopGenerators, externalGenerators},
   timeGenerators = <|"type" -> "time", "vertex" -> #|> & /@ Lookup[topo, "activeVertexIds", topo["vertexIds"]];
   loopGenerators = Flatten[
     Table[
       makeAffineMomentumGenerator[topo, l, "loop", m],
      {l, topo["nL"]}, {m, topo["nL"]}
      ],
     1
     ];
   externalGenerators = Flatten[
     Table[
       makeAffineMomentumGenerator[topo, l, "external", j],
      {l, topo["nL"]}, {j, topo["nK"]}
      ],
     1
     ];
   Join[timeGenerators, loopGenerators, externalGenerators]
   ];


expectedMomentumGeneratorCount[topo_Association] := topo["nL"] (topo["nL"] + topo["nK"]);


(* ::Chapter:: *)
(*轻量 momentum IBP seed*)

(* 本章生成 momentum IBP seed 的传播子幂次项、z/ISP 吸收和 massive/massless building-block 导数项。
   输出会经过 EOM 与 per-line massless endpoint canonical 门禁；shrunk pack 的 bS 幂次也使用当前 pack 指标。 *)

linearTerms[expr_] := Module[{expanded = Expand[expr]},
   If[Head[expanded] === Plus, List @@ expanded, {expanded}]
   ];


shiftLineB[J[aList_, linePacks_, ispList_], e_Integer, delta_] := Module[
   {newLinePacks = linePacks},
   newLinePacks[[e, 1]] = newLinePacks[[e, 1]] + delta;
   J[aList, newLinePacks, ispList]
   ];


fixedLineMomentumMagnitude[topo_Association, e_Integer] := Module[
   {momentum = Expand[topo["lines"][[e, "momentum"]]], rules, directSquare,
     externalData, coefficients, effectiveLoopExternalMomenta, square, binding},
   binding = SelectFirst[
     externalLegMagnitudeBindingData[topo],
     SameQ[
       canonicalExternalLegMomentum[Lookup[#, "momentum", 0]],
       canonicalExternalLegMomentum[momentum]
       ] &,
     Missing["NoMagnitudeBinding"]
     ];
   If[AssociationQ[binding], Return[Lookup[binding, "userMagnitudeExpression"]]];
   rules = Join[
     Lookup[topo, "resolvedMagnitudeKinematicRules", {}],
     Lookup[topo, "resolvedLoopKinematicRules", {}]
     ];
   directSquare = sp[momentum, momentum] /. rules;
   If[directSquare =!= sp[momentum, momentum],
     Return[If[MatchQ[directSquare, Power[_, 2]], directSquare[[1]], Sqrt[directSquare]]]
    ];
   effectiveLoopExternalMomenta = Lookup[topo, "effectiveLoopExternalMomenta", {}];
   externalData = linearMomentumExpressionData[momentum, effectiveLoopExternalMomenta];
   If[TrueQ[externalData["linearQ"]],
    coefficients = externalData["coefficients"];
    square = Sum[
       coefficients[[i]] coefficients[[j]] sp[effectiveLoopExternalMomenta[[i]], effectiveLoopExternalMomenta[[j]]],
       {i, Length[effectiveLoopExternalMomenta]}, {j, Length[effectiveLoopExternalMomenta]}
       ] /. rules;
     Return[If[MatchQ[square, Power[_, 2]], square[[1]], Sqrt[Expand[square]]]]
    ];
   Sqrt[sp[momentum, momentum]]
   ];


(* cycle line 的 xi[e] 是积分变量；bridge 的模长由公开 momentum 和运动学规则唯一确定。 *)
lineMomentumMagnitude[topo_Association, e_Integer] := If[
   lineIndexedPowerQ[topo["lines"][[e]]],
   xi[topo["lines"][[e, "id"]]],
   fixedLineMomentumMagnitude[topo, e]
   ];


(* cycle 线把幂移写入 b/bS；bridge 线没有连续幂指标，同一变化必须成为显式动量系数。 *)
shiftLinePower[topo_Association, int_J, e_Integer, delta_] := If[
   lineIndexedPowerQ[topo["lines"][[e]]],
   shiftLineB[int, e, delta],
   fixedLineMomentumMagnitude[topo, e]^(-delta) int
   ];


zeroPointRuleValue[topo_Association, symbol_, default_: 0] := Module[
   {hits},
   hits = Cases[Lookup[topo, "zeroPointRules", {}], (Rule | RuleDelayed)[lhs_, rhs_] /; lhs === symbol :> rhs];
   If[hits === {}, default, Last[hits]]
   ];


vertexZeroPoint[topo_Association, vertexId_] := zeroPointRuleValue[topo, a0[vertexRepresentative[topo, vertexId]]];


lineBZeroPoint[topo_Association, e_Integer] := zeroPointRuleValue[topo, b0[topo["lines"][[e, "id"]]]];


lineBSZeroPoint[topo_Association, e_Integer] := zeroPointRuleValue[topo, bS0[topo["lines"][[e, "id"]]]];


vertexPowerIndex[topo_Association, J[aList_, linePacks_, ispList_], vertexId_] := Module[
   {pos = vertexASlot[topo, vertexId]},
   If[Head[pos] === Missing, 0, aList[[pos]] + vertexZeroPoint[topo, vertexId]]
   ];


linePowerIndex[topo_Association, J[aList_, linePacks_, ispList_], e_Integer] := Module[
   {base, packType = actualLinePackType[topo, e, linePacks[[e]]]},
   base = lineIntegerPowerIndex[topo, J[aList, linePacks, ispList], e];
   base + If[packType === "shrunk", lineBSZeroPoint[topo, e], lineBZeroPoint[topo, e]]
   ];


shiftISPIndex[J[aList_, linePacks_, ispList_], j_Integer, delta_] := Module[
   {newISPList = ispList},
   newISPList[[j]] = newISPList[[j]] + delta;
   J[aList, linePacks, newISPList]
   ];


(* 将一个线性因子吸收到 b 或 ISP 指标中；无法识别的高次项保持为显式系数。 *)
absorbLinearFactorTerm[term_, int_J, topo_Association] := Module[
   {zVars, ispVars, vars, rules, rebuildMonomial},
   zVars = z /@ loopDependentLineIndices[topo];
   ispVars = Table[rho[j], {j, Length[topo["ispData"]]}];
   vars = Join[zVars, ispVars];
   If[Length[vars] == 0, Return[term int]];
   rules = CoefficientRules[term, vars];
   rebuildMonomial[powers_] := Times @@ MapThread[#1^#2 &, {vars, powers}];
   Total[
    rules /. (powers_ -> coeff_) :> Module[
       {degree = Total[powers], pos, var},
       If[degree === 0,
        coeff int,
        If[Count[powers, 1] === 1 && Count[powers, Except[0 | 1]] === 0 && degree === 1,
         pos = First@FirstPosition[powers, 1];
         var = vars[[pos]];
         Which[
          MatchQ[var, z[_Integer]],
          coeff shiftLinePower[topo, int, var[[1]], -2],
          MemberQ[ispVars, var],
          coeff shiftISPIndex[int, First@FirstPosition[ispVars, var], 1],
          True,
          coeff var int
          ],
         coeff rebuildMonomial[powers] int
         ]
        ]
       ]
    ]
   ];


absorbLinearFactor[factor_, int_J, topo_Association] := Total[
   absorbLinearFactorTerm[#, int, topo] & /@ linearTerms[factor]
   ];


(* fixed line 的幂移会先产生显式模长系数；这里逐项拆出唯一积分，
   再复用裸 J 的吸收规则，避免内部 helper 留在公开导数结果中。 *)
absorbLinearFactorExpressionTerm[factor_, term_, topo_Association] := Module[
   {integrals, integral, coefficient},
   integrals = DeleteDuplicates[Cases[term, _J, {0, Infinity}]];
   If[Length[integrals] =!= 1, Return[factor term]];
   integral = First[integrals];
   coefficient = Cancel[term/integral];
   If[
    FreeQ[coefficient, _J],
    coefficient absorbLinearFactor[factor, integral, topo],
    factor term
    ]
   ];


absorbLinearFactor[factor_, expr_, topo_Association] := Total[
   absorbLinearFactorExpressionTerm[factor, #, topo] & /@ linearTerms[Expand[expr]]
   ];


momentumDivergenceTerm[int_J, gen_Association] := If[
   gen["type"] === "momentum" && gen["vectorType"] === "loop" && gen["dLoop"] === gen["vectorIndex"],
   dim int,
   0
   ];


(* 同分支传播子的有序端点方向由 endpoints[[1]] 定义；++ 取 +1，-- 取 -1。 *)
sameBranchSKSign[line_Association] := If[
   StringTake[Lookup[line, "skType", "++"], 1] === "+",
   1,
   -1
   ];


lineEndpointSlotsAtVertex[line_Association, vertexId_] := Flatten @ Position[
   line["endpoints"],
   vertexId
   ];


applyCompiledTimeDerivativeTerm[
   topo_Association, int_J, e_Integer, endpointSlot_Integer, term_Association
   ] := Module[{result, endpointVertex, nPosition, xPower = term["xPower"]},
   endpointVertex = topo["lines"][[e, "endpoints", endpointSlot]];
   nPosition = linePackNPositions[topo["lines"][[e]], actualLinePackType[topo, e, int[[2, e]]]][[endpointSlot]];
   result = setLinePackEntry[int, e, nPosition, term["targetState"]];
   result = shiftVertexA[result, topo, endpointVertex, xPower];
   result = shiftLinePower[topo, result, e, -(xPower + 1)];
   -term["coefficient"] result
   ];


compiledTimeEndpointDerivativeTerms[
   topo_Association, int_J, e_Integer, endpointSlot_Integer
   ] := Module[{state, terms},
   state = int[[2, e, linePackNPositions[topo["lines"][[e]], actualLinePackType[topo, e, int[[2, e]]]][[endpointSlot]]]];
   terms = Lookup[lineCompiledFunctionSystem[topo["lines"][[e]]], "derivativeTerms", {}];
   Total[
    KroneckerDelta[state, Lookup[#, "sourceState", Missing["NoSourceState"]]] *
       applyCompiledTimeDerivativeTerm[topo, int, e, endpointSlot, #] & /@ terms
    ]
   ];


applyCompiledMomentumDerivativeTerm[
   topo_Association, int_J, e_Integer, endpointSlot_Integer, factor_, term_Association
   ] := Module[{result, endpointVertex, nPosition, xPower = term["xPower"]},
   endpointVertex = topo["lines"][[e, "endpoints", endpointSlot]];
   nPosition = linePackNPositions[topo["lines"][[e]], actualLinePackType[topo, e, int[[2, e]]]][[endpointSlot]];
   result = setLinePackEntry[int, e, nPosition, term["targetState"]];
   result = shiftVertexA[result, topo, endpointVertex, xPower + 1];
   result = shiftLinePower[topo, result, e, 1 - xPower];
   term["coefficient"] absorbLinearFactor[factor, result, topo]
   ];


compiledMomentumEndpointDerivativeTerms[
   topo_Association, int_J, e_Integer, endpointSlot_Integer, factor_
   ] := Module[{state, terms},
   state = int[[2, e, linePackNPositions[topo["lines"][[e]], actualLinePackType[topo, e, int[[2, e]]]][[endpointSlot]]]];
   terms = Lookup[lineCompiledFunctionSystem[topo["lines"][[e]]], "derivativeTerms", {}];
   Total[
    KroneckerDelta[state, Lookup[#, "sourceState", Missing["NoSourceState"]]] *
       applyCompiledMomentumDerivativeTerm[topo, int, e, endpointSlot, factor, #] & /@ terms
    ]
   ];


toggleMasslessLineState[topo_Association, J[aList_, linePacks_, ispList_], e_Integer] := Module[
   {newLinePacks = linePacks, nPosition},
   nPosition = First[linePackNPositions[topo["lines"][[e]], "masslessFull"]];
   newLinePacks[[e, nPosition]] = 1 - newLinePacks[[e, nPosition]];
   J[aList, newLinePacks, ispList]
   ];


(* q 导数同时作用 massive Hankel block 与 massless 指数核。
   masslessFull 使用 d_q M_n = i sigma (tau_u-tau_v) M_(1-n)；
   masslessCross 使用两个端点相位符号的和。 *)
momentumBuildingBlockDerivativeTerms[topo_Association, int_J, gen_Association, repSP2ZRules_List] := Module[
   {dLoop, vector, lineMomenta, lines, loopCoeff, vDotQ,
     shiftedInt, packType, sigma},
   dLoop = gen["dLoop"];
   vector = gen["vector"];
   lineMomenta = Lookup[topo["lines"], "momentum"];
   lines = topo["lines"];
   Total[
    Table[
     loopCoeff = Coefficient[lineMomenta[[e]], topo["loopMomenta"][[dLoop]]];
     If[zeroQ[loopCoeff],
      0,
      vDotQ = Expand[expandDotExpr[vector, lineMomenta[[e]], topo] /. repSP2ZRules];
      packType = lines[[e]]["packType"];
      Switch[packType,
       "massiveFull" | "massiveCross",
        Total[
         Table[
          loopCoeff compiledMomentumEndpointDerivativeTerms[topo, int, e, endpointSlot, vDotQ],
          {endpointSlot, 2}
          ]
        ],
       "masslessFull",
       sigma = sameBranchSKSign[lines[[e]]];
        shiftedInt = shiftLinePower[topo, toggleMasslessLineState[topo, int, e], e, 1];
       loopCoeff (
         -I sigma absorbLinearFactor[
           vDotQ,
           shiftVertexA[shiftedInt, topo, lines[[e]]["endpoints"][[1]], 1],
           topo
           ] +
          I sigma absorbLinearFactor[
           vDotQ,
           shiftVertexA[shiftedInt, topo, lines[[e]]["endpoints"][[2]], 1],
           topo
           ]
         ),
       "masslessCross",
       shiftedInt = shiftLinePower[topo, int, e, 1];
       loopCoeff Total[
         Table[
          -I skEndpointPhaseSign[lines[[e]], endpointSlot] absorbLinearFactor[
            vDotQ,
            shiftVertexA[shiftedInt, topo, lines[[e]]["endpoints"][[endpointSlot]], 1],
            topo
            ],
          {endpointSlot, 2}
          ]
         ],
       _,
       0
       ]
      ],
     {e, topo["nE"]}
     ]
    ]
   ];

momentumPropagatorDerivativeTerms[topo_Association, int_J, gen_Association, repSP2ZRules_List] := Module[
   {dLoop, vector, lineMomenta, loopCoeff, vDotQ, shiftedInt},
   dLoop = gen["dLoop"];
   vector = gen["vector"];
   lineMomenta = Lookup[topo["lines"], "momentum"];
   Total[
    Table[
     loopCoeff = Coefficient[lineMomenta[[e]], topo["loopMomenta"][[dLoop]]];
     If[zeroQ[loopCoeff],
      0,
      vDotQ = Expand[expandDotExpr[vector, lineMomenta[[e]], topo] /. repSP2ZRules];
      shiftedInt = shiftLinePower[topo, int, e, 2];
      -loopCoeff linePowerIndex[topo, int, e] absorbLinearFactor[vDotQ, shiftedInt, topo]
      ],
     {e, topo["nE"]}
     ]
    ]
   ];


(* ISP 是 numerator 因子；momentum IBP 必须同时微分 ISP^r。
   内部坐标使用 qq/qk/kk，方向导数后再沿现有 z/rho 坐标吸收到指标。 *)
momentumDirectionalSPDerivative[topo_Association, expr_, gen_Association] := Module[
   {dLoop = gen["dLoop"], vector = gen["vector"], loops = topo["loopMomenta"], exts = topo["effectiveLoopExternalMomenta"]},
   Expand[
    expr /. {
      HoldPattern[qq[i_Integer, j_Integer]] :>
       If[i === dLoop, expandDotExpr[vector, loops[[j]], topo], 0] +
        If[j === dLoop, expandDotExpr[loops[[i]], vector, topo], 0],
      HoldPattern[qk[i_Integer, j_Integer]] :>
       If[i === dLoop, expandDotExpr[vector, exts[[j]], topo], 0],
      HoldPattern[kk[_Integer, _Integer]] :> 0
      }
    ]
   ];


momentumISPDerivativeTerms[topo_Association, int_J, gen_Association, repSP2ZRules_List] := Module[
   {ispExprs, exponent, derivative},
   ispExprs = normalizeISPExprs[topo];
   Total[
    Table[
     exponent = int[[3, j]];
     If[zeroQ[exponent],
      0,
      derivative = Expand[momentumDirectionalSPDerivative[topo, ispExprs[[j]], gen] /. repSP2ZRules];
      exponent absorbLinearFactor[
        derivative,
        shiftISPIndex[int, j, -1],
        topo
        ]
      ],
     {j, Length[ispExprs]}
     ]
    ]
   ];


applyMomentumGeneratorSeed::nosp =
   "拓扑 `1` 的标量积到 z/ISP 规则未生成：`2`。请补充 ISP 配置后再生成 momentum seed。";


primitiveMomentumGeneratorSeed[
   topo_Association, int_J, gen_Association, ruleData_Association
   ] := Expand[
   momentumDivergenceTerm[int, gen] +
    momentumPropagatorDerivativeTerms[topo, int, gen, ruleData["repSP2Z"]] +
    momentumBuildingBlockDerivativeTerms[topo, int, gen, ruleData["repSP2Z"]] +
    momentumISPDerivativeTerms[topo, int, gen, ruleData["repSP2Z"]]
   ];


applyMomentumGeneratorSeed[topo_Association, int_J, gen_Association] := Module[
   {ruleData, primitiveTerms, results},
   ruleData = makeScalarProductRules[topo];
   If[Lookup[ruleData, "status", "notComputed"] =!= "computed",
     Message[applyMomentumGeneratorSeed::nosp, topo["name"], Lookup[ruleData, "reason", Missing["reason"]]];
     Return[$Failed]
     ];
   primitiveTerms = Lookup[
     gen,
     "normalizedPrimitiveTerms",
     {<|"coefficient" -> 1, "generator" -> KeyDrop[gen, {"normalizedPrimitiveTerms"}]|>}
     ];
   results = Lookup[#, "coefficient", 0] primitiveMomentumGeneratorSeed[
        topo, int, Lookup[#, "generator", <||>], ruleData
        ] & /@ primitiveTerms;
   If[! FreeQ[results, $Failed], Return[$Failed]];
   Expand[Total[results]]
   ];


(* ::Chapter:: *)
(*轻量 time IBP core seed*)

(* 本章接入 time-IBP 的通用 core 项：顶点幂次、外部能量、massive building-block 端点导数、massless 端点翻转项和 massive/massless theta 边界缩并项。
   单独 time batch 会把进一步 shrink-sector 生成标为 pending；canonical batch 完整补齐全部 contact-reachable sectors。 *)

rawVertexExternalEnergy[topo_Association, vertexId_] := Module[
   {sectorEnergies = Lookup[topo, "sectorExternalLegEnergyByVertex", Missing["NotSet"]]},
   Which[
    AssociationQ[sectorEnergies] && KeyExistsQ[sectorEnergies, vertexId], sectorEnergies[vertexId],
    ListQ[sectorEnergies], vertexId /. sectorEnergies,
    True,
    ke[vertexId]
    ]
   ];


vertexExternalEnergy[topo_Association, vertexId_] := scalarProductInputToInternal[rawVertexExternalEnergy[topo, vertexId], topo];


externalLegEnergyDependencyData[topo_Association, vertexId_] := Module[
   {internal, vars, externalVars, externalUsed, independentUsed},
   internal = vertexExternalEnergy[topo, vertexId];
   vars = Variables[internal];
   externalVars = externalInvariantVariables[topo];
   externalUsed = Intersection[vars, externalVars];
   independentUsed = Complement[vars, externalVars];
   <|
    "internalExternalInvariantVariables" -> externalUsed,
    "externalInvariantVariables" -> (scalarProductInternalToUser[#, topo] & /@ externalUsed),
    "internalIndependentExternalLegEnergyParameters" -> independentUsed,
    "independentExternalLegEnergyParameters" -> (scalarProductInternalToUser[#, topo] & /@ independentUsed),
    "usesExternalInvariantQ" -> TrueQ[externalUsed =!= {}],
    "usesIndependentExternalLegEnergyQ" -> TrueQ[independentUsed =!= {}],
    "kind" -> Which[
      externalUsed =!= {} && independentUsed === {}, "externalInvariantExpression",
      externalUsed === {} && independentUsed =!= {}, "independentExternalLegEnergyParameter",
      externalUsed =!= {} && independentUsed =!= {}, "mixedExpression",
      True, "constant"
      ]
    |>
   ];


externalLegEnergyNamingReport[topo_Association] := Module[
   {vertices = activeAVertexIds[topo], raw, internal, user, dependencies},
   raw = AssociationThread[vertices -> (rawVertexExternalEnergy[topo, #] & /@ vertices)];
   internal = AssociationThread[vertices -> (vertexExternalEnergy[topo, #] & /@ vertices)];
   user = AssociationThread[vertices -> (scalarProductInternalToUser[vertexExternalEnergy[topo, #], topo] & /@ vertices)];
   dependencies = AssociationThread[vertices -> (externalLegEnergyDependencyData[topo, #] & /@ vertices)];
   <|
    "convention" -> "vertex external energy uses ke[i] for independent absolute-value parameters; expressions built from external invariant names are normalized to the same scalar-product coordinates used by loop momenta",
    "rawExternalLegEnergies" -> raw,
    "internalExternalLegEnergies" -> internal,
    "userExternalLegEnergies" -> user,
    "dependencyData" -> dependencies,
    "message" -> "每个 vertices.externalLegEnergy 表示对应顶点外腿打包后的相位能量。若它绑定 loopExternalMomenta，输入原始 Sqrt[sp[p,p]] 并由缺省 ssij 或 exact 自定义规则输出；若绑定 independentExternalMomenta 的实际无圈模长，则输出为 sEi；否则写独立标量。不要把 |p1+p2| 与 |p1|+|p2| 混同，程序不会自动生成无圈动量之间的点积关系。"
    |>
   ];


timeVertexPowerTerm[topo_Association, J[aList_, linePacks_, ispList_], vertexId_] := Module[
   {pos = vertexASlot[topo, vertexId]},
   If[Head[pos] === Missing, Return[0]];
   -vertexPowerIndex[topo, J[aList, linePacks, ispList], vertexId] shiftVertexA[J[aList, linePacks, ispList], topo, vertexId, -1]
   ];


(* 顶点 + 对应 exp[-i E tau]，顶点 - 对应 exp[+i E tau]。 *)
vertexExternalPhaseDerivativeCoefficient[topo_Association, vertexId_] := If[
   Lookup[topo["vertexSignAssoc"], vertexId, "+"] === "+",
   -I,
   I
   ];


timeExternalEnergyTerm[topo_Association, int_J, vertexId_] :=
   vertexExternalPhaseDerivativeCoefficient[topo, vertexId] vertexExternalEnergy[topo, vertexId] int;


skEndpointPhaseSign[line_Association, endpointSlot_Integer] := Module[
   {skType = Lookup[line, "skType", "++"], chars},
   chars = Characters[skType];
   If[Length[chars] < endpointSlot,
    If[endpointSlot === 1, 1, -1],
    If[chars[[endpointSlot]] === "+", 1, -1]
    ]
   ];


(* massless time regular 原子：masslessFull 在有序端点上执行 {b,n}->{b-1,1-n}，两端系数为 +I sigma/-I sigma；
   若 sector 已使两端点 coincident，endpointSlots 同时包含 1、2，两个 regular 项必须在此处相消。
   masslessCross 没有 n 或 theta，只按各端点 SK 相位移动 b。 *)
timeMasslessEndpointDerivativeTerms[topo_Association, J[aList_, linePacks_, ispList_], vertexId_] := Module[
   {pos, connectedLines, lines = topo["lines"], endpointSlots, endpointSign,
    newLinePacks, shiftedIntegral, nPosition, sigma},
   pos = vertexPosition[topo, vertexId];
   If[Head[pos] === Missing, Return[0]];
   connectedLines = topo["vertexLines"][[pos]][[All, 1]];
   Total[
    Table[
     endpointSlots = lineEndpointSlotsAtVertex[lines[[e]], vertexId];
     Switch[lines[[e]]["packType"],
      "masslessFull",
      sigma = sameBranchSKSign[lines[[e]]];
      Total[
       Table[
        endpointSign = sigma If[endpointSlot === 1, 1, -1];
        newLinePacks = linePacks;
        nPosition = First[linePackNPositions[lines[[e]], "masslessFull"]];
        newLinePacks[[e, nPosition]] = 1 - newLinePacks[[e, nPosition]];
        shiftedIntegral = shiftLinePower[topo, J[aList, newLinePacks, ispList], e, -1];
        I endpointSign shiftedIntegral,
        {endpointSlot, endpointSlots}
        ]
       ],
      "masslessCross",
      Total[
       Table[
        endpointSign = skEndpointPhaseSign[lines[[e]], endpointSlot];
        shiftedIntegral = shiftLinePower[topo, J[aList, linePacks, ispList], e, -1];
        I endpointSign shiftedIntegral,
        {endpointSlot, endpointSlots}
        ]
       ],
      _,
      0
      ],
     {e, connectedLines}
     ]
    ]
   ];


(* 缩并后同一条线的两个原端点可能落到同一 active vertex；
   此时必须同时微分两个 building block，不能只取 FirstPosition。 *)
timeMassiveBuildingBlockDerivativeTerms[topo_Association, J[aList_, linePacks_, ispList_], vertexId_] := Module[
   {pos, connectedLines, lines = topo["lines"], endpointSlots},
   pos = vertexPosition[topo, vertexId];
   If[Head[pos] === Missing, Return[0]];
   connectedLines = topo["vertexLines"][[pos]][[All, 1]];
   Total[
    Table[
     If[! MemberQ[{"massiveFull", "massiveCross"}, actualLinePackType[topo, e, linePacks[[e]]]],
      0,
      endpointSlots = lineEndpointSlotsAtVertex[lines[[e]], vertexId];
       Total[
        Table[
         compiledTimeEndpointDerivativeTerms[topo, J[aList, linePacks, ispList], e, endpointSlot],
         {endpointSlot, endpointSlots}
         ]
       ]
      ],
     {e, connectedLines}
     ]
    ]
   ];


lineCompiledShrinkTerms[line_Association] := Lookup[
   lineCompiledFunctionSystem[line],
   "shrinkTerms",
   {}
   ];


lineShrinkZeroPointShift[line_Association] := Module[
   {compiled = lineCompiledFunctionSystem[line]},
   Lookup[
    line,
    "shrinkZeroPointShift",
    Lookup[compiled, "shrinkZeroPointShift", 0]
    ]
   ];


sectorZeroPointRules[topo_Association, shrunkLines_List, repMap_Association, activeVertices_List] := Module[
   {lineRules, vertexRules, classVertices, classShrunkLines, shiftSum},
   lineRules = Flatten@Table[
      If[MemberQ[shrunkLines, e],
       bS0[topo["lines"][[e, "id"]]] -> lineBZeroPoint[topo, e] + lineShrinkZeroPointShift[topo["lines"][[e]]],
       Switch[actualLinePackType[topo, e, makeLinePack[topo["lines"][[e]]]],
        "shrunk", bS0[topo["lines"][[e, "id"]]] -> lineBSZeroPoint[topo, e],
        _, b0[topo["lines"][[e, "id"]]] -> lineBZeroPoint[topo, e]
        ]
       ],
      {e, topo["nE"]}
      ];
   vertexRules = Table[
     classVertices = Select[topo["vertexIds"], Lookup[repMap, #] === rep &];
     classShrunkLines = Select[shrunkLines, Function[e, And @@ (Function[v, Lookup[repMap, v] === rep] /@ topo["lines"][[e, "endpoints"]])]];
     shiftSum = Total[lineShrinkZeroPointShift[topo["lines"][[#]]] & /@ classShrunkLines];
     a0[rep] -> Total[vertexZeroPoint[topo, #] & /@ classVertices] - shiftSum,
     {rep, activeVertices}
     ];
   DeleteDuplicates@Join[vertexRules, lineRules]
   ];


(* massive Wronskian 缩并带一个额外 1/q 和 1/(-tau)，故 bS=b+1 且 merged a 减 1；
   massless 反对称 theta-delta 不带这些因子，必须保持 bS=b 且 merged a 不移位。 *)
lineShrinkBShift[line_Association] := If[
   Lookup[line, "massType", "massive"] === "massless",
   0,
   Lookup[First[Lookup[lineCompiledFunctionSystem[line], "shrinkTerms", {<|"bShift" -> 1|>}]], "bShift", 1]
   ];


shrinkLineIntegral[topo_Association, J[aList_, linePacks_, ispList_], e_Integer, bShift_: Automatic, aShift_: Automatic] := Module[
   {line = topo["lines"][[e]], uSlot, vSlot, oldActive, newRepMap, newActive, newAList, newLinePacks = linePacks,
     mergedRep, oldSlotsForNewRep, slotValues, effectiveBShift, effectiveAShift, powerCoefficient},
   effectiveBShift = If[bShift === Automatic, lineShrinkBShift[line], bShift];
    effectiveAShift = If[
      aShift === Automatic,
      If[Lookup[line, "massType", "massive"] === "massless", 0, effectiveBShift],
      aShift
      ];
   uSlot = vertexASlot[topo, line["endpoints"][[1]]];
   vSlot = vertexASlot[topo, line["endpoints"][[2]]];
   If[Head[uSlot] === Missing || Head[vSlot] === Missing, Return[$Failed]];
   oldActive = activeAVertexIds[topo];
   newRepMap = vertexRepresentativeMap[topo["vertexIds"], Join[
      ({#, vertexRepresentative[topo, #]} & /@ topo["vertexIds"]),
      {line["endpoints"]}
      ]];
   newActive = DeleteDuplicates[Lookup[newRepMap, topo["vertexIds"]]];
   mergedRep = Lookup[newRepMap, line["endpoints"][[1]]];
   newAList = Table[
     oldSlotsForNewRep = Flatten[Position[Lookup[newRepMap, oldActive], newActive[[i]]]];
     slotValues = aList[[oldSlotsForNewRep]];
     If[newActive[[i]] === mergedRep,
       If[uSlot === vSlot, 2 aList[[uSlot]], Total[slotValues]] - effectiveAShift,
       Total[slotValues]
       ],
     {i, Length[newActive]}
     ];
   powerCoefficient = If[lineIndexedPowerQ[line], 1, fixedLineMomentumMagnitude[topo, e]^(-effectiveBShift)];
   newLinePacks[[e]] = If[
     lineIndexedPowerQ[line],
     {lineIntegerPowerIndex[topo, J[aList, linePacks, ispList], e] + effectiveBShift},
     {}
     ];
   powerCoefficient J[newAList, newLinePacks, ispList]
   ];


shrinkLinesIntegral[
   topo_Association,
   J[aList_, linePacks_, ispList_],
   specs_List
   ] := Module[
   {selectedLines, oldActive, pairs, newRepMap, newActive, newAList,
    newLinePacks = linePacks, oldSlotsForNewRep, selectedShiftForRep, powerCoefficient = 1},
   selectedLines = Lookup[specs, "lineIndex"];
   oldActive = activeAVertexIds[topo];
   pairs = topo["lines"][[#, "endpoints"]] & /@ selectedLines;
   newRepMap = vertexRepresentativeMap[topo["vertexIds"], Join[
      ({#, vertexRepresentative[topo, #]} & /@ topo["vertexIds"]),
      pairs
      ]];
   newActive = DeleteDuplicates[Lookup[newRepMap, topo["vertexIds"]]];
   newAList = Table[
     oldSlotsForNewRep = Flatten[Position[Lookup[newRepMap, oldActive], newActive[[i]]]];
     selectedShiftForRep = Total[MapThread[
        If[SameQ @@ Lookup[newRepMap, #1] && First[Lookup[newRepMap, #1]] === newActive[[i]], #2, 0] &,
        {pairs, Lookup[specs, "aShift"]}
        ]];
     Total[aList[[oldSlotsForNewRep]]] - selectedShiftForRep,
     {i, Length[newActive]}
     ];
   Scan[
     Function[spec,
      If[lineIndexedPowerQ[topo["lines"][[spec["lineIndex"]]]],
       newLinePacks[[spec["lineIndex"]]] = {
         lineIntegerPowerIndex[topo, J[aList, linePacks, ispList], spec["lineIndex"]] + spec["bShift"]
         },
       powerCoefficient *= fixedLineMomentumMagnitude[topo, spec["lineIndex"]]^(-spec["bShift"]);
       newLinePacks[[spec["lineIndex"]]] = {}
       ]
      ],
    specs
    ];
   powerCoefficient J[newAList, newLinePacks, ispList]
   ];


thetaBoundaryAtomicTerms[
   topo_Association,
   J[aList_, linePacks_, ispList_],
   e_Integer,
   vertexId_
   ] := Module[
   {line = topo["lines"][[e]], endpointSlots, endpointSlot, endpointOrientation,
    pack = linePacks[[e]], packType, coeff, shrinkTerms},
   endpointSlots = lineEndpointSlotsAtVertex[line, vertexId];
   If[Length[endpointSlots] =!= 1, Return[{}]];
   endpointSlot = First[endpointSlots];
   endpointOrientation = If[endpointSlot === 1, 1, -1];
   packType = actualLinePackType[topo, e, pack];
   Switch[packType,
    "massiveFull",
    coeff = With[{nPositions = linePackNPositions[line, packType]},
      -sameBranchSKSign[line] KroneckerDelta[Total[pack[[nPositions]]], 1]
        (-1)^pack[[nPositions[[endpointSlot]]]]
      ];
    shrinkTerms = lineCompiledShrinkTerms[line];
    (<|
        "lineIndex" -> e,
       "coefficient" -> coeff Lookup[#, "coefficient", 0],
        "bShift" -> Lookup[#, "bShift", 1],
        "aShift" -> Lookup[#, "bShift", 1]
        |> &) /@ shrinkTerms,
    "masslessFull",
    {<|
      "lineIndex" -> e,
      "coefficient" -> -2 endpointOrientation KroneckerDelta[pack[[First[linePackNPositions[line, packType]]]], 1],
      "bShift" -> 0,
      "aShift" -> 0
      |>},
    _,
    {}
    ]
   ];


thetaBundleKey[topo_Association, e_Integer] := Sort[vertexPosition[topo, #] & /@ topo["lines"][[e, "endpoints"]]];


(* 同一代表顶点对的 full lines 共享一个 theta。共同边界
   Product[A]-Product[B] 在 J=(A+B)/2、D=A-B 基底中只含非空奇数 D 子集，
   k-line contact 的系数为 2^(1-k)。每个子集只合并顶点一次，绝不生成 delta^k。 *)
timeThetaBoundaryShrinkTerms[topo_Association, int : J[aList_, linePacks_, ispList_], vertexId_] := Module[
   {pos, connectedLines, eligibleLines, bundles, oddSubsets, atomicChoices},
   pos = vertexPosition[topo, vertexId];
   If[Head[pos] === Missing, Return[0]];
   connectedLines = topo["vertexLines"][[pos]][[All, 1]];
   eligibleLines = Select[
     connectedLines,
     MemberQ[{"massiveFull", "masslessFull"}, actualLinePackType[topo, #, linePacks[[#]]]] &&
       Length[lineEndpointSlotsAtVertex[topo["lines"][[#]], vertexId]] === 1 &
     ];
   bundles = GatherBy[eligibleLines, thetaBundleKey[topo, #] &];
   Total[Table[
     oddSubsets = Select[Rest[Subsets[bundle]], OddQ[Length[#]] &];
     Total[Table[
       (* 精确为零的端点 contact 不参与组合；否则 0 J 会被错误送入 sector normalization。 *)
       atomicChoices = Select[
           thetaBoundaryAtomicTerms[topo, int, #, vertexId],
           ! zeroQ[Lookup[#, "coefficient", 0]] &
           ] & /@ selected;
       If[AnyTrue[atomicChoices, # === {} &],
        0,
        Total[
         Function[choice,
            normalizeContactTerm018[
             topo,
             int,
             2^(1 - Length[selected]) Times @@ Lookup[choice, "coefficient"] *
              shrinkLinesIntegral[topo, int, KeyDrop[#, "coefficient"] & /@ choice]
             ]
            ] /@ Tuples[atomicChoices]
         ]
        ],
       {selected, oddSubsets}
       ]],
     {bundle, bundles}
     ]]
   ];


seedUnsupportedPendingFeatures[topo_Association] := DeleteDuplicates@Join[
    unsupportedSeedFeaturesForTopology[topo]
    ];


momentumIBPPendingFeatures[topo_Association] := seedUnsupportedPendingFeatures[topo];


timeIBPPendingFeatures[topo_Association] := DeleteDuplicates@Join[
    seedUnsupportedPendingFeatures[topo],
    If[thetaFullLineIndices[topo] =!= {}, {"shrinkSectorSeedGeneration"}, {}]
    ];


timeGeneratorLabel[gen_Association] := {gen["type"], gen["vertex"]};


applyTimeGeneratorSeed::badgen = "time seed 只能使用 time 生成元，收到：`1`。";


applyTimeGeneratorSeed[topo_Association, int_J, gen_Association] := Module[
   {vertexId},
   If[gen["type"] =!= "time",
    Message[applyTimeGeneratorSeed::badgen, gen];
    Return[$Failed]
    ];
   vertexId = gen["vertex"];
   Expand[
    timeVertexPowerTerm[topo, int, vertexId] +
     timeExternalEnergyTerm[topo, int, vertexId] +
     timeMasslessEndpointDerivativeTerms[topo, int, vertexId] +
     timeMassiveBuildingBlockDerivativeTerms[topo, int, vertexId] +
     timeThetaBoundaryShrinkTerms[topo, int, vertexId]
    ]
   ];


(* ::Chapter:: *)
(*符号 seed 的结构索引*)

(* general template 的连续指标只在这里枚举变量；实际目标域由 DSGenerateIBP 接收。 *)

indexVariableQ[x_] := ! TrueQ[IntegerQ[x] || NumericQ[x]];


continuousIndexVariables[J[aList_, linePacks_, ispList_]] := Select[
   Join[aList, Cases[Flatten[linePacks], _b | _bS], ispList],
   indexVariableQ
   ];


Options[selectedDiscreteSeedRules] = {};


selectedDiscreteSeedRules[topo_Association, OptionsPattern[]] := Module[
   {rules = enumerateDiscreteStates[makeBaseIntegral[topo], topo]["rules"]},
   <|"status" -> "generated", "mode" -> "all", "ruleCount" -> Length[rules], "rules" -> rules|>
   ];


momentumGeneratorLabel[gen_Association] := {gen["type"], gen["dLoop"], gen["vectorType"], gen["vectorIndex"]};


(* ::Chapter:: *)
(*独立变量导数 seed*)

(* 本章生成微分方程阶段使用的独立变量导数。顶点能量 ke[i] 只微分顶点相位；
   外部不变量先在约束坐标中解成外动量矢量导数 k_i.d/dk_j 的线性组合，解的 basis 显式返回。 *)

externalVectorDerivativeGenerators[topo_Association] := Flatten[
   Table[
    <|
     "type" -> "externalVector",
     "vectorIndex" -> i,
     "dExternal" -> j,
     "vector" -> topo["effectiveLoopExternalMomenta"][[i]]
     |>,
    {i, topo["nK"]}, {j, topo["nK"]}
    ],
   1
   ];


externalVectorDerivativeGeneratorBasis[topo_Association, Automatic] := externalVectorDerivativeGeneratorBasis[topo, "upperTriangular"];
externalVectorDerivativeGeneratorBasis[topo_Association, "upperTriangular"] := Flatten[
   Table[
    If[i <= j,
     <|
      "type" -> "externalVector",
      "vectorIndex" -> i,
      "dExternal" -> j,
      "vector" -> topo["effectiveLoopExternalMomenta"][[i]]
      |>,
     Nothing
     ],
    {i, topo["nK"]}, {j, topo["nK"]}
    ],
   1
   ];
externalVectorDerivativeGeneratorBasis[topo_Association, "all"] := externalVectorDerivativeGenerators[topo];
externalVectorDerivativeGeneratorBasis[topo_Association, gens_List] := gens;
externalVectorDerivativeGeneratorBasis[topo_Association, _] := externalVectorDerivativeGeneratorBasis[topo, "upperTriangular"];


externalVectorDerivativeLabel[gen_Association] := {gen["type"], gen["vectorIndex"], gen["dExternal"]};


externalVectorSPCoordinateDerivative[topo_Association, coordinate_, gen_Association] := Module[
   {dExternal = gen["dExternal"], vector = gen["vector"], loops = topo["loopMomenta"], exts = topo["effectiveLoopExternalMomenta"]},
   Expand[
    coordinate /. {
      HoldPattern[qq[_Integer, _Integer]] :> 0,
      HoldPattern[qk[i_Integer, j_Integer]] :>
       If[j === dExternal, expandDotExpr[loops[[i]], vector, topo], 0],
      HoldPattern[kk[i_Integer, j_Integer]] :>
       If[i === dExternal, expandDotExpr[vector, exts[[j]], topo], 0] +
        If[j === dExternal, expandDotExpr[exts[[i]], vector, topo], 0]
      }
   ]
   ];


(* 对非线性动力学量表达式必须显式执行链式法则，不能把坐标直接替成其方向导数。 *)
externalVectorDirectionalSPDerivative[topo_Association, expr_, gen_Association] := Module[
   {coordinates},
   coordinates = DeleteDuplicates[Cases[expr, _qq | _qk | _kk, {0, Infinity}]];
   Expand[Total[
     D[expr, #] externalVectorSPCoordinateDerivative[topo, #, gen] & /@ coordinates
     ]]
   ];


externalVectorPropagatorDerivativeTerms[topo_Association, int_J, gen_Association, repSP2ZRules_List] := Module[
   {dExternal, vector, lineMomenta, extCoeff, vDotQ, shiftedInt},
   dExternal = gen["dExternal"];
   vector = gen["vector"];
   lineMomenta = Lookup[topo["lines"], "momentum"];
   Total[
    Table[
     extCoeff = Coefficient[lineMomenta[[e]], topo["effectiveLoopExternalMomenta"][[dExternal]]];
     If[zeroQ[extCoeff] || externalLegCoordinateLineQ[lineMomenta[[e]], topo],
      0,
      vDotQ = Expand[expandDotExpr[vector, lineMomenta[[e]], topo] /. repSP2ZRules];
      shiftedInt = shiftLinePower[topo, int, e, 2];
      -extCoeff linePowerIndex[topo, int, e] absorbLinearFactor[vDotQ, shiftedInt, topo]
      ],
     {e, topo["nE"]}
     ]
    ]
   ];


(* 外部动力学量导数与 qIBP 共用最终 AT 编译结果；WT 只用于 theta shrink。 *)
externalVectorBuildingBlockDerivativeTerms[topo_Association, int_J, gen_Association, repSP2ZRules_List] := Module[
   {dExternal, vector, lineMomenta, lines, extCoeff, vDotQ,
    shiftedInt, packType, sigma},
   dExternal = gen["dExternal"];
   vector = gen["vector"];
   lineMomenta = Lookup[topo["lines"], "momentum"];
   lines = topo["lines"];
   Total[
    Table[
     extCoeff = Coefficient[lineMomenta[[e]], topo["effectiveLoopExternalMomenta"][[dExternal]]];
     If[zeroQ[extCoeff] || externalLegCoordinateLineQ[lineMomenta[[e]], topo],
      0,
      vDotQ = Expand[expandDotExpr[vector, lineMomenta[[e]], topo] /. repSP2ZRules];
      packType = lines[[e]]["packType"];
      Switch[packType,
       "massiveFull" | "massiveCross",
       Total[
        Table[
         extCoeff compiledMomentumEndpointDerivativeTerms[
           topo, int, e, endpointSlot, vDotQ
           ],
         {endpointSlot, 2}
         ]
        ],
       "masslessFull",
       sigma = sameBranchSKSign[lines[[e]]];
       shiftedInt = shiftLinePower[topo, toggleMasslessLineState[topo, int, e], e, 1];
       extCoeff (
         -I sigma absorbLinearFactor[
           vDotQ,
           shiftVertexA[shiftedInt, topo, lines[[e]]["endpoints"][[1]], 1],
           topo
           ] +
          I sigma absorbLinearFactor[
           vDotQ,
           shiftVertexA[shiftedInt, topo, lines[[e]]["endpoints"][[2]], 1],
           topo
           ]
         ),
       "masslessCross",
       shiftedInt = shiftLinePower[topo, int, e, 1];
       extCoeff Total[
         Table[
          -I skEndpointPhaseSign[lines[[e]], endpointSlot] absorbLinearFactor[
            vDotQ,
            shiftVertexA[shiftedInt, topo, lines[[e]]["endpoints"][[endpointSlot]], 1],
            topo
            ],
          {endpointSlot, 2}
          ]
         ],
       _,
       0
       ]
      ],
     {e, topo["nE"]}
     ]
    ]
   ];


externalVectorISPDerivativeTerms[topo_Association, int_J, gen_Association, repSP2ZRules_List] := Module[
   {ispExprs, exponent, derivative},
   ispExprs = normalizeISPExprs[topo];
   Total[
    Table[
     exponent = int[[3, j]];
     If[zeroQ[exponent],
      0,
      derivative = Expand[externalVectorDirectionalSPDerivative[topo, ispExprs[[j]], gen] /. repSP2ZRules];
      exponent absorbLinearFactor[
        derivative,
        shiftISPIndex[int, j, -1],
        topo
        ]
      ],
     {j, Length[ispExprs]}
     ]
    ]
   ];


externalVectorExternalLegEnergyDerivativeTerms[topo_Association, int_J, gen_Association] := Module[
   {vertices = activeAVertexIds[topo], derivative},
   Total[
    Table[
     derivative = externalVectorDirectionalSPDerivative[topo, vertexExternalEnergy[topo, vertexId], gen];
     If[zeroQ[derivative],
      0,
      -vertexExternalPhaseDerivativeCoefficient[topo, vertexId] derivative shiftVertexA[int, topo, vertexId, 1]
      ],
     {vertexId, vertices}
     ]
    ]
   ];


applyExternalVectorDerivativeSeed::badgen = "external-vector seed 只能使用 externalVector 生成元，收到：`1`。";
applyExternalVectorDerivativeSeed::nosp =
   "拓扑 `1` 的标量积到 z/ISP 规则未生成：`2`。请补充 ISP 配置后再生成 external-vector seed。";


applyExternalVectorDerivativeSeed[topo_Association, int_J, gen_Association] := Module[
   {ruleData},
   If[gen["type"] =!= "externalVector",
    Message[applyExternalVectorDerivativeSeed::badgen, gen];
    Return[$Failed]
    ];
   ruleData = makeScalarProductRules[topo];
   If[Lookup[ruleData, "status", "notComputed"] =!= "computed",
    Message[applyExternalVectorDerivativeSeed::nosp, topo["name"], Lookup[ruleData, "reason", Missing["reason"]]];
    Return[$Failed]
    ];
   Expand[
    externalVectorPropagatorDerivativeTerms[topo, int, gen, ruleData["repSP2Z"]] +
     externalVectorBuildingBlockDerivativeTerms[topo, int, gen, ruleData["repSP2Z"]] +
     externalVectorISPDerivativeTerms[topo, int, gen, ruleData["repSP2Z"]] +
     externalVectorExternalLegEnergyDerivativeTerms[topo, int, gen]
    ]
   ];


Options[makeExternalInvariantDerivativeDecomposition] = {
   ExternalInvariantCoordinateVariables -> Automatic,
   ExternalVectorOperatorBasis -> Automatic
   };

makeExternalInvariantDerivativeDecomposition::badvar =
   "变量 `1` 不在当前外部不变量坐标 `2` 中；若它是组合变量，请显式给出 ExternalInvariantCoordinateVariables。";
makeExternalInvariantDerivativeDecomposition::nosol =
   "变量 `1` 无法由所选 external-vector operator basis `2` 解出。";


normalizeExternalInvariantCoordinateList[topo_Association, Automatic] := externalInvariantVariables[topo];
normalizeExternalInvariantCoordinateList[topo_Association, coords_List] := scalarProductInputToInternal[#, topo] & /@ coords;
normalizeExternalInvariantCoordinateList[topo_Association, coord_] := {scalarProductInputToInternal[coord, topo]};


makeExternalInvariantDerivativeDecomposition[topo_Association, var_, OptionsPattern[]] := Module[
   {target, coords, targetPos, gens, matrix, coeffs, equations, solutions, selectedRules,
    selectedCoefficients, residual, fullGens, fullMatrix, rank, nullity, freeCoeffRules},
   target = scalarProductInputToInternal[var, topo];
   coords = normalizeExternalInvariantCoordinateList[topo, OptionValue[ExternalInvariantCoordinateVariables]];
   targetPos = FirstPosition[coords, target, Missing["NotFound"]];
   If[Head[targetPos] === Missing,
    Message[makeExternalInvariantDerivativeDecomposition::badvar, var, scalarProductInternalToUser[#, topo] & /@ coords];
    Return[<|"status" -> "badVariable", "targetVariable" -> scalarProductInternalToUser[target, topo], "coordinateVariables" -> scalarProductInternalToUser[#, topo] & /@ coords|>]
    ];
   targetPos = First[targetPos];
   gens = externalVectorDerivativeGeneratorBasis[topo, OptionValue[ExternalVectorOperatorBasis]];
   matrix = Table[
     externalVectorDirectionalSPDerivative[topo, coords[[r]], gens[[c]]],
     {r, Length[coords]}, {c, Length[gens]}
     ];
   coeffs = Array[cv, Length[gens]];
   equations = Thread[matrix . coeffs == UnitVector[Length[coords], targetPos]];
   solutions = Quiet[Solve[equations, coeffs]];
   If[solutions === {},
    Message[makeExternalInvariantDerivativeDecomposition::nosol, var, externalVectorDerivativeLabel /@ gens];
    Return[<|"status" -> "noSolution", "targetVariable" -> scalarProductInternalToUser[target, topo], "operatorBasis" -> externalVectorDerivativeLabel /@ gens, "matrix" -> matrix|>]
    ];
   selectedRules = First[solutions];
   freeCoeffRules = Thread[Complement[coeffs, selectedRules[[All, 1]]] -> 0];
   selectedCoefficients = Expand[coeffs /. selectedRules /. freeCoeffRules];
   residual = Simplify[Expand[matrix . selectedCoefficients - UnitVector[Length[coords], targetPos]]];
   fullGens = externalVectorDerivativeGenerators[topo];
   fullMatrix = Table[
     externalVectorDirectionalSPDerivative[topo, coords[[r]], fullGens[[c]]],
     {r, Length[coords]}, {c, Length[fullGens]}
     ];
   rank = Quiet[MatrixRank[fullMatrix]];
   nullity = If[IntegerQ[rank], Length[fullGens] - rank, Missing["SymbolicRankNotComputed"]];
   <|
    "status" -> "solved",
    "targetVariable" -> scalarProductInternalToUser[target, topo],
    "internalTargetVariable" -> target,
    "coordinateVariables" -> (scalarProductInternalToUser[#, topo] & /@ coords),
    "internalCoordinateVariables" -> coords,
    "operatorBasis" -> (externalVectorDerivativeLabel /@ gens),
    "operators" -> gens,
    "matrix" -> matrix,
    "coefficientRules" -> Thread[externalVectorDerivativeLabel /@ gens -> selectedCoefficients],
    "coefficients" -> selectedCoefficients,
    "residual" -> residual,
    "solutionFamilyRules" -> solutions,
    "fullOperatorCount" -> Length[fullGens],
    "coordinateCount" -> Length[coords],
    "fullOperatorRank" -> rank,
    "nullity" -> nullity,
    "nonUniqueQ" -> TrueQ[IntegerQ[nullity] && nullity > 0]
    |>
   ];


directExternalLegEnergyVariableDerivativeSeed[topo_Association, int_J, var_] := Module[
   {internalVar, vertices = activeAVertexIds[topo], derivative},
   internalVar = scalarProductInputToInternal[var, topo];
   Total[
    Table[
     derivative = D[vertexExternalEnergy[topo, vertexId], internalVar];
     If[zeroQ[derivative],
      0,
      -vertexExternalPhaseDerivativeCoefficient[topo, vertexId] derivative shiftVertexA[int, topo, vertexId, 1]
      ],
     {vertexId, vertices}
     ]
    ]
   ];


Options[applyExternalInvariantVariableDerivativeSeed] = Options[makeExternalInvariantDerivativeDecomposition];
applyExternalInvariantVariableDerivativeSeed[topo_Association, int_J, var_, opts : OptionsPattern[]] := Module[
   {decomp, terms},
   decomp = makeExternalInvariantDerivativeDecomposition[topo, var, FilterRules[{opts}, Options[makeExternalInvariantDerivativeDecomposition]]];
   If[Lookup[decomp, "status", "failed"] =!= "solved", Return[$Failed]];
   terms = MapThread[
     #1 applyExternalVectorDerivativeSeed[topo, int, #2] &,
     {decomp["coefficients"], decomp["operators"]}
     ];
   Expand[Total[terms]]
   ];


Options[applyIndependentVariableDerivativeSeed] = Options[makeExternalInvariantDerivativeDecomposition];
applyIndependentVariableDerivativeSeed[topo_Association, int_J, var_, opts : OptionsPattern[]] := Module[
   {internalVar},
   internalVar = scalarProductInputToInternal[var, topo];
   If[MemberQ[externalInvariantVariables[topo], internalVar],
    Return[applyExternalInvariantVariableDerivativeSeed[topo, int, internalVar, opts]]
    ];
   Expand[directExternalLegEnergyVariableDerivativeSeed[topo, int, internalVar]]
   ];


(* 独立变量集合按 external invariant 坐标与未被其表达的顶点能量参数组成。*)
independentVariableDerivativeVariables[topo_Association] := DeleteDuplicates@Join[
   externalInvariantVariables[topo],
   Complement[externalLegEnergyVariables[topo], externalInvariantVariables[topo]]
   ];


independentVariableDerivativeKind[topo_Association, var_] := If[
   MemberQ[externalInvariantVariables[topo], var],
   "externalInvariant",
   "externalLegEnergy"
   ];


Options[makeIndependentVariableDerivativeSeedBatch] = Options[makeExternalInvariantDerivativeDecomposition];
makeIndependentVariableDerivativeSeedBatch[topo_Association, int_J, opts : OptionsPattern[]] := Module[
   {generators, results, derivative, decomposition, kind, variable, failureQ},
   generators = makeIndependentVariableDerivativeGenerators[topo];
   results = Table[
     variable = generator["variable"];
     kind = generator["kind"];
     decomposition = If[kind === "externalInvariant",
       makeExternalInvariantDerivativeDecomposition[
        topo,
        variable,
        FilterRules[{opts}, Options[makeExternalInvariantDerivativeDecomposition]]
        ],
       Missing["NotApplicable"]
       ];
     derivative = If[kind === "externalInvariant" && Lookup[decomposition, "status", "failed"] =!= "solved",
       $Failed,
       applyIndependentVariableDerivativeSeed[
        topo,
        int,
        variable,
        Sequence @@ FilterRules[{opts}, Options[makeExternalInvariantDerivativeDecomposition]]
        ]
       ];
     derivative = If[derivative === $Failed, $Failed, applySeedCanonical[Expand[derivative], topo]];
     <|
      "variable" -> variable,
      "userVariable" -> generator["userVariable"],
      "kind" -> kind,
      "decomposition" -> decomposition,
      "derivative" -> derivative,
      "status" -> If[derivative === $Failed, "failed", "generated"],
      "forbiddenNData" -> If[derivative === $Failed, {}, forbiddenNData[topo, derivative]],
      "canonicalQ" -> If[derivative === $Failed, False, ! containsForbiddenNQ[topo, derivative]]
      |>,
     {generator, generators}
     ];
   failureQ = AnyTrue[results, Lookup[#, "status", "failed"] =!= "generated" &];
   <|
    "status" -> If[failureQ, "failed", "generated"],
    "caseName" -> Lookup[topo, "name", Missing["caseName"]],
    "variables" -> generators,
    "variableCount" -> Length[generators],
    "integral" -> int,
    "equationCount" -> Length[results],
    "equations" -> results,
    "allCanonicalQ" -> And @@ Lookup[results, "canonicalQ", True],
    "failedVariables" -> Lookup[Select[results, Lookup[#, "status", "failed"] =!= "generated" &], "userVariable", {}]
    |>
   ];


(* ::Chapter:: *)
(*自动 shrink-sector seed*)

(* 本章从 top-sector 的 massive full lines 派生 shrunk sectors。
   为避免大图爆炸，默认有 sector 数量保护；超过保护时只保留 pending gate，不展开。 *)

vertexRepresentativeMap[vertexIds_List, pairs_List] := Module[
   {parent, pos, find, union},
   parent = AssociationThread[vertexIds -> vertexIds];
   pos = AssociationThread[vertexIds -> Range[Length[vertexIds]]];
   find[x_] := If[parent[x] === x, x, parent[x] = find[parent[x]]];
   union[x_, y_] := Module[{rx = find[x], ry = find[y], rep, other},
     If[rx === ry, Return[]];
     rep = If[pos[rx] <= pos[ry], rx, ry];
     other = If[rep === rx, ry, rx];
     parent[other] = rep;
     ];
   Scan[union[#[[1]], #[[2]]] &, pairs];
   Association[Table[v -> find[v], {v, vertexIds}]]
   ];


remapExtLegsToRepresentatives[extLegs_List, repMap_Association] := Map[
   If[ListQ[#] && Length[#] >= 2,
     Join[#[[1 ;; 1]], {Lookup[repMap, #[[2]], #[[2]]]}, Drop[#, 2]],
     #
     ] &,
   extLegs
   ];


remapExternalLegEnergiesToRepresentatives[sectorEnergies_, repMap_Association] := Module[
   {rules, grouped},
   Which[
    AssociationQ[sectorEnergies],
    rules = Normal[sectorEnergies] /. (v_ -> val_) :> (Lookup[repMap, v, v] -> val);
    grouped = Merge[rules, Total];
    grouped,
    ListQ[sectorEnergies],
    rules = sectorEnergies /. (v_ -> val_) :> (Lookup[repMap, v, v] -> val);
    Normal[Merge[rules, Total]],
    True,
    sectorEnergies
    ]
   ];


shrinkSectorTopology[topo_Association, shrunkLines_List] := Module[
   {pairs, repMap, activeVertices, fixedA, newLines, newExtLegs, newZeroPointRules,
    newBMatrix, newVertexLines, sectorTopo},
   pairs = topo["lines"][[#, "endpoints"]] & /@ shrunkLines;
   repMap = vertexRepresentativeMap[topo["vertexIds"], pairs];
   activeVertices = DeleteDuplicates[Lookup[repMap, topo["vertexIds"]]];
   fixedA = Association[Thread[Complement[topo["vertexIds"], activeVertices] -> 0]];
   newLines = MapIndexed[
     Module[{line = #1, e = First[#2], endpoints},
       endpoints = Lookup[repMap, line["endpoints"]];
       If[MemberQ[shrunkLines, e],
        Join[line, <|"originalEndpoints" -> Lookup[line, "originalEndpoints", line["endpoints"]], "endpoints" -> endpoints, "state" -> "shrunk", "packType" -> "shrunk"|>],
        Join[line, <|"originalEndpoints" -> Lookup[line, "originalEndpoints", line["endpoints"]], "endpoints" -> endpoints|>]
        ]
       ] &,
     topo["lines"]
     ];
   newExtLegs = remapExtLegsToRepresentatives[topo["extLegs"], repMap];
   newZeroPointRules = sectorZeroPointRules[topo, shrunkLines, repMap, activeVertices];
   newBMatrix = Table[
     Which[
      newLines[[e, "endpoints", 1]] === topo["vertexIds"][[v]], 1,
      newLines[[e, "endpoints", 2]] === topo["vertexIds"][[v]], -1,
      True, 0
      ],
     {v, Length[topo["vertexIds"]]}, {e, Length[newLines]}
     ];
   newVertexLines = Table[
     Select[Table[{e, newBMatrix[[v, e]]}, {e, Length[newLines]}], #[[2]] =!= 0 &],
     {v, Length[topo["vertexIds"]]}
     ];
   (* sector 只消费已解析 root topology，并更新 contact 真正改变的内部字段。
      不重新构造用户 case，也不再次调用公开 parser。 *)
   sectorTopo = Join[topo, <|
     "name" -> topo["name"] <> "_sector_" <> StringRiffle["e" <> ToString[#] & /@ shrunkLines, "_"],
     "lines" -> newLines,
     "extLegs" -> newExtLegs,
     "sectorExternalLegEnergyByVertex" -> remapExternalLegEnergiesToRepresentatives[
       topo["sectorExternalLegEnergyByVertex"],
       repMap
       ],
     "bMatrix" -> newBMatrix,
     "vertexLines" -> newVertexLines,
     "zeroPointRules" -> newZeroPointRules,
     "rootZeroPointRules" -> Lookup[topo, "rootZeroPointRules", topo["zeroPointRules"]],
     "sectorVertexRepresentativeMap" -> repMap,
     "activeVertexIds" -> activeVertices,
     "fixedAVertexValues" -> fixedA,
     "sectorShrunkLines" -> shrunkLines
     |>];
   Join[sectorTopo, <|
    "sectorMetadata" -> makeSectorMetadata[sectorTopo],
    "tadpoleSymmetryData" -> tadpoleSymmetryData[sectorTopo],
    "effectiveSymmetryRules" -> effectiveSymmetryRules0[sectorTopo]
    |>]
   ];


massiveFullLineIndices[topo_Association] := Flatten@Position[Lookup[topo["lines"], "packType"], "massiveFull"];


masslessFullLineIndices[topo_Association] := Flatten@Position[Lookup[topo["lines"], "packType"], "masslessFull"];


thetaFullLineIndices[topo_Association] := Sort @ Join[
   massiveFullLineIndices[topo],
   masslessFullLineIndices[topo]
   ];


contactReachableShrinkSubsets[topo_Association] := Module[
   {lines = thetaFullLineIndices[topo], queue = {{}}, seen = <||>, result = {},
    state, pairs, repMap, eligible, bundles, events, newState, key},
   AssociateTo[seen, ToString[{}, InputForm] -> True];
   While[queue =!= {},
    state = First[queue];
    queue = Rest[queue];
    pairs = topo["lines"][[#, "endpoints"]] & /@ state;
    repMap = vertexRepresentativeMap[topo["vertexIds"], pairs];
    eligible = Select[
      Complement[lines, state],
      ! SameQ @@ Lookup[repMap, topo["lines"][[#, "endpoints"]]] &
      ];
    bundles = GatherBy[
      eligible,
      Sort[vertexPosition[topo, #] & /@ Lookup[repMap, topo["lines"][[#, "endpoints"]]]] &
      ];
    events = Flatten[Select[Rest[Subsets[#]], OddQ[Length[#]] &] & /@ bundles, 1];
    Do[
     newState = Sort[Union[state, event]];
     key = ToString[newState, InputForm];
     If[! KeyExistsQ[seen, key],
      AssociateTo[seen, key -> True];
      AppendTo[result, newState];
      AppendTo[queue, newState]
      ],
     {event, events}
     ];
    ];
   SortBy[result, {Length[#], #} &]
   ];


shrinkSectorSubsets[topo_Association] := Module[
   {lines = thetaFullLineIndices[topo], subsets},
   If[lines === {}, Return[<|"status" -> "generated", "subsets" -> {}, "completeCoverageQ" -> True|>]];
   subsets = contactReachableShrinkSubsets[topo];
   <|"status" -> "generated", "subsets" -> subsets, "completeCoverageQ" -> True|>
   ];


Options[makeTopologyData] = {PrecomputeShrinkSectorMetadata -> False};


makeTopologyIndexMaps[topo_Association, metadata_Association] := <|
   "vertexIdToOriginalASlot" -> metadata["vertexIdToOriginalASlot"],
   "vertexIdToCompactASlot" -> metadata["vertexIdToCompactASlot"],
   "lineIdToSlot" -> metadata["lineIdToSlot"],
   "bSymbolToLineSlot" -> metadata["bSymbolToLineSlot"],
   "compactASlots" -> metadata["compactASlots"],
   "lineSlots" -> metadata["lineSlots"],
   "ispSlots" -> metadata["ispSlots"]
   |>;


makeTopologySeedSummary[topo_Association] := <|
   "continuousVariables" -> continuousIndexVariables[makeBaseIntegral[topo]],
   "discreteVariables" -> Flatten[discreteVarsForLine /@ topo["lines"]],
   "discreteStateCount" -> discreteStateCount[topo],
   "momentumGeneratorCount" -> Count[Lookup[makeIBPGenerators[topo], "type"], "momentum"],
   "timeGeneratorCount" -> Count[Lookup[makeIBPGenerators[topo], "type"], "time"],
   "loopKinematicNamingReport" -> loopKinematicNamingReport[topo],
    "externalLegEnergyNamingReport" -> externalLegEnergyNamingReport[topo]
    |>;


vertexPairBundleKey[line_Association] := Sort[Lookup[line, "originalEndpoints", line["endpoints"]]];


masslessBundleCandidates[topo_Association] := Module[
   {indexedLines, masslessFullLines, grouped},
   indexedLines = MapIndexed[Join[#1, <|"lineIndex" -> First[#2]|>] &, topo["lines"]];
   masslessFullLines = Select[indexedLines, #["packType"] === "masslessFull" &];
   grouped = GroupBy[masslessFullLines, vertexPairBundleKey];
   KeyValueMap[
    If[Length[#2] > 1,
      <|
       "vertexPair" -> #1,
       "lineIndices" -> Lookup[#2, "lineIndex"],
       "lineIds" -> Lookup[#2, "id"],
       "packTemplates" -> (makeLinePack /@ #2),
        "bundleMode" -> "commonThetaOddSubsetContacts"
       |>,
      Nothing
      ] &,
    grouped
    ]
   ];


masslessEndpointConventionData[topo_Association] := MapIndexed[
   <|
     "lineIndex" -> First[#2],
     "lineId" -> #1["id"],
     "orderedEndpoints" -> #1["endpoints"],
     "n1ReferenceEndpoint" -> #1["endpoints"][[1]],
     "n1OppositeEndpoint" -> #1["endpoints"][[2]],
     "convention" -> "n=1 is the antisymmetric state defined by the first endpoint; swapping endpoints flips n=1"
     |> &,
   Select[topo["lines"], #["packType"] === "masslessFull" &]
   ];


unsupportedSeedFeaturesForTopology[topo_Association] := DeleteDuplicates@Join[
    {}
    ];


topologyValidationReport[topo_Association] := Module[
   {issues = {}, appendIssue, vertexIds, lineIds, packTypes, allowedPackTypes,
     duplicateLoopMomenta, duplicateExternalMomenta, duplicateExternalLegMomenta,
     loopExternalMomentumOverlap, externalLegMomentumOverlap,
     vertexSigns, activeVertexIds, fixedAVertexIds, badVertexSigns, badActiveVertexIds, badFixedAVertexIds,
    extLegs, badExtLegShapePositions, badExtLegVertexData, sectorExternalLegEnergies, externalLegEnergyKeys, badExternalLegEnergyKeys,
    ispNames, badISPRangeData,
     zeroPointRuleValidationReport,
    badMassTypeLines, badSKTypeLines, badStateLines,
    badEndpointLines, lineMomentumVars, declaredMomentumVars, undeclaredMomentumVars,
    nonLinearLineMomentumData, nonLinearScalarProductArgumentData, externalLegEnergyMomentumDependenceData,
    spData, pendingFeatures, ruleData, kinematicAudit,
      topologyMomentumAudit, momentumIBPRequiredQ},
   appendIssue[severity_, code_, data_: <||>] := AppendTo[issues, Join[<|"severity" -> severity, "code" -> code|>, data]];
   vertexIds = topo["vertexIds"];
   vertexSigns = Lookup[topo["vertices"], "vertexType"];
   activeVertexIds = Lookup[topo, "activeVertexIds", vertexIds];
   fixedAVertexIds = Keys[Lookup[topo, "fixedAVertexValues", <||>]];
   extLegs = Lookup[topo, "extLegs", {}];
   sectorExternalLegEnergies = Lookup[topo, "sectorExternalLegEnergyByVertex", <||>];
   kinematicAudit = Lookup[topo, "kinematicCoordinateAudit", <||>];
   topologyMomentumAudit = <|
     "status" -> Lookup[Lookup[topo, "momentumDeclarationAudit", <||>], "status", "invalid"],
     "issues" -> Lookup[Lookup[topo, "momentumDeclarationAudit", <||>], "issues", {}]
     |>;
   momentumIBPRequiredQ = Lookup[topo, "ibpMode", "full"] === "full";
   Scan[
    appendIssue[Lookup[#, "severity", "error"], Lookup[#, "code", "momentumDeclarationIssue"], KeyDrop[#, {"severity", "code"}]] &,
    Lookup[topologyMomentumAudit, "issues", {}]
    ];
   ispNames = Lookup[topo["ispData"], "name", {}];
   If[Lookup[kinematicAudit, "status", "complete"] === "incomplete",
    appendIssue["error", "incompleteKinematicCoordinates", <|
      "coordinateRank" -> Lookup[kinematicAudit, "coordinateRank", Missing["rank"]],
      "parameterRank" -> Lookup[kinematicAudit, "parameterRank", Missing["parameterRank"]],
      "baseCoordinateCount" -> Lookup[kinematicAudit, "baseCoordinateCount", Missing["count"]],
      "missingDirections" -> Lookup[kinematicAudit, "missingDirections", {}],
      "parameterMissingDirections" -> Lookup[kinematicAudit, "parameterMissingDirections", {}],
      "unsupportedRulePositions" -> Lookup[kinematicAudit, "unsupportedRulePositions", {}]
      |>]
    ];
   If[Lookup[kinematicAudit, "status", "complete"] === "overcomplete",
    appendIssue["warning", "overcompleteKinematicCoordinates", <|
      "constraintResiduals" -> Lookup[kinematicAudit, "constraintResiduals", {}],
      "ruleDependencies" -> Lookup[kinematicAudit, "ruleDependencies", {}],
      "comment" -> "IBP generation may continue, but inverse conversion and ds in the redundant coordinates are disabled until the family definition supplies the constraints."
      |>]
    ];
   lineIds = Lookup[topo["lines"], "id"];
   packTypes = Lookup[topo["lines"], "packType"];
   allowedPackTypes = {"massiveFull", "massiveCross", "masslessFull", "masslessCross", "shrunk"};
   duplicateLoopMomenta = Cases[Tally[topo["loopMomenta"]], {mom_, count_} /; count > 1 :> mom];
   If[duplicateLoopMomenta =!= {},
    appendIssue["error", "duplicateLoopMomenta", <|"loopMomenta" -> topo["loopMomenta"], "duplicates" -> duplicateLoopMomenta|>]
    ];
   duplicateExternalMomenta = Cases[Tally[topo["effectiveLoopExternalMomenta"]], {mom_, count_} /; count > 1 :> mom];
   If[duplicateExternalMomenta =!= {},
    appendIssue["error", "duplicateExternalMomenta", <|"effectiveLoopExternalMomenta" -> topo["effectiveLoopExternalMomenta"], "duplicates" -> duplicateExternalMomenta|>]
    ];
   duplicateExternalLegMomenta = Cases[Tally[Lookup[topo, "independentExternalMomenta", {}]], {mom_, count_} /; count > 1 :> mom];
   If[duplicateExternalLegMomenta =!= {},
    appendIssue["error", "duplicateExternalLegMomenta", <|"independentExternalMomenta" -> Lookup[topo, "independentExternalMomenta", {}], "duplicates" -> duplicateExternalLegMomenta|>]
    ];
   loopExternalMomentumOverlap = Intersection[topo["loopMomenta"], topo["effectiveLoopExternalMomenta"]];
   If[loopExternalMomentumOverlap =!= {},
    appendIssue["error", "loopExternalMomentumOverlap", <|"overlap" -> loopExternalMomentumOverlap, "loopMomenta" -> topo["loopMomenta"], "effectiveLoopExternalMomenta" -> topo["effectiveLoopExternalMomenta"]|>]
    ];
   externalLegMomentumOverlap = Intersection[
     Lookup[topo, "independentExternalMomenta", {}],
     Join[topo["loopMomenta"], topo["effectiveLoopExternalMomenta"]]
     ];
   If[externalLegMomentumOverlap =!= {},
    appendIssue["error", "externalLegMomentumOverlap", <|"overlap" -> externalLegMomentumOverlap|>]
    ];
   If[! DuplicateFreeQ[lineIds],
    appendIssue["error", "duplicateLineIds", <|"lineIds" -> lineIds|>]
    ];
   If[ispNames =!= {} && ! DuplicateFreeQ[ispNames],
    appendIssue["error", "duplicateISPNames", <|"ispNames" -> ispNames|>]
    ];
   badISPRangeData = DeleteCases[
     MapIndexed[
      If[KeyExistsQ[#1, "range"] && ! validIndexRangeSpecQ[#1["range"]],
        <|"ispPosition" -> First[#2], "name" -> Lookup[#1, "name", Missing["name"]], "rangeSpec" -> #1["range"]|>,
        Nothing
        ] &,
      topo["ispData"]
      ],
     Nothing
     ];
   If[badISPRangeData =!= {},
    appendIssue["error", "malformedISPRangeSpecs", <|"isps" -> badISPRangeData, "allowed" -> "integer or nonempty integer list"|>]
    ];
   zeroPointRuleValidationReport = validateCoefficientRules[topo["zeroPointRules"]];
   If[Lookup[zeroPointRuleValidationReport, "status", "ok"] =!= "ok",
    appendIssue["error", "invalidZeroPointRules", KeyDrop[zeroPointRuleValidationReport, {"status"}]]
    ];
   If[! DuplicateFreeQ[vertexIds],
    appendIssue["error", "duplicateVertexIds", <|"vertexIds" -> vertexIds|>]
    ];
   badVertexSigns = DeleteDuplicates[Complement[vertexSigns, {"+", "-"}]];
   If[badVertexSigns =!= {},
    appendIssue["error", "unknownVertexSigns", <|"vertexSigns" -> badVertexSigns, "allowedVertexSigns" -> {"+", "-"}|>]
    ];
   badActiveVertexIds = Complement[activeVertexIds, vertexIds];
   If[badActiveVertexIds =!= {},
    appendIssue["error", "activeVertexIdsNotInVertexData", <|"activeVertexIds" -> badActiveVertexIds, "vertexIds" -> vertexIds|>]
    ];
   badFixedAVertexIds = Complement[fixedAVertexIds, vertexIds];
   If[badFixedAVertexIds =!= {},
    appendIssue["error", "fixedAVertexValuesNotInVertexData", <|"fixedAVertexIds" -> badFixedAVertexIds, "vertexIds" -> vertexIds|>]
    ];
   If[! ListQ[extLegs],
    appendIssue["error", "malformedExtLegs", <|"reason" -> "extLegs must be a list of entries with at least {label, vertexId, energy}"|>],
    badExtLegShapePositions = Flatten @ Position[extLegs, entry_ /; !(ListQ[entry] && Length[entry] >= 3), {1}, Heads -> False];
    If[badExtLegShapePositions =!= {},
     appendIssue["error", "malformedExtLegs", <|"badPositions" -> badExtLegShapePositions|>]
     ];
    badExtLegVertexData = DeleteCases[
      MapIndexed[
       If[ListQ[#1] && Length[#1] >= 2 && ! MemberQ[vertexIds, #1[[2]]],
         <|"extLegPosition" -> First[#2], "vertexId" -> #1[[2]]|>,
         Nothing
         ] &,
       extLegs
       ],
      Nothing
      ];
    If[badExtLegVertexData =!= {},
     appendIssue["error", "extLegVertexNotInVertexData", <|"extLegs" -> badExtLegVertexData, "vertexIds" -> vertexIds|>]
     ]
    ];
   externalLegEnergyKeys = Which[
     AssociationQ[sectorExternalLegEnergies], Keys[sectorExternalLegEnergies],
     ListQ[sectorExternalLegEnergies] && And @@ (MatchQ[#, _Rule | _RuleDelayed] & /@ sectorExternalLegEnergies), Cases[sectorExternalLegEnergies, (Rule | RuleDelayed)[v_, _] :> v],
     sectorExternalLegEnergies === <||>, {},
     True, Missing["MalformedExternalLegEnergies"]
     ];
   If[externalLegEnergyKeys === Missing["MalformedExternalLegEnergies"],
     appendIssue["error", "malformedExternalLegEnergies", <|"reason" -> "internal sectorExternalLegEnergyByVertex must be an Association or list of rules"|>],
    badExternalLegEnergyKeys = Complement[externalLegEnergyKeys, vertexIds];
    If[badExternalLegEnergyKeys =!= {},
      appendIssue["error", "sectorExternalLegEnergyNotInVertices", <|"externalLegEnergyKeys" -> badExternalLegEnergyKeys, "vertexIds" -> vertexIds|>]
     ]
    ];
   externalLegEnergyMomentumDependenceData = externalLegEnergyMomentumDependenceIssues[topo];
   If[externalLegEnergyMomentumDependenceData =!= {},
    appendIssue["error", "invalidExternalLegEnergyMomentumDependence", <|"issues" -> externalLegEnergyMomentumDependenceData, "comment" -> "vertices.externalLegEnergy must use declared loopExternalMomenta/independentExternalMomenta through raw Sqrt[sp[p,p]] expressions, or independent scalar phase parameters"|>]
    ];
   badEndpointLines = Select[
     MapIndexed[Join[#1, <|"lineIndex" -> First[#2]|>] &, topo["lines"]],
     ! SubsetQ[vertexIds, Lookup[#, "endpoints", {}]] &
     ];
   If[badEndpointLines =!= {},
    appendIssue["error", "lineEndpointNotInVertexData", <|"lineIndices" -> Lookup[badEndpointLines, "lineIndex"], "endpoints" -> Lookup[badEndpointLines, "endpoints"]|>]
    ];
   If[Complement[packTypes, allowedPackTypes] =!= {},
    appendIssue["error", "unknownPackTypes", <|"packTypes" -> DeleteDuplicates[Complement[packTypes, allowedPackTypes]]|>]
    ];
   badMassTypeLines = Select[
     MapIndexed[Join[#1, <|"lineIndex" -> First[#2]|>] &, topo["lines"]],
     ! MemberQ[{"massive", "massless"}, Lookup[#, "massType", "massive"]] &
     ];
   If[badMassTypeLines =!= {},
    appendIssue["error", "unknownLineMassTypes", <|"lineIndices" -> Lookup[badMassTypeLines, "lineIndex"], "massTypes" -> Lookup[badMassTypeLines, "massType"], "allowedMassTypes" -> {"massive", "massless"}|>]
    ];
   badSKTypeLines = Select[
     MapIndexed[Join[#1, <|"lineIndex" -> First[#2]|>] &, topo["lines"]],
     ! MemberQ[{"++", "--", "+-", "-+"}, Lookup[#, "skType", "++"]] &
     ];
   If[badSKTypeLines =!= {},
    appendIssue["error", "unknownLineSKTypes", <|"lineIndices" -> Lookup[badSKTypeLines, "lineIndex"], "skTypes" -> Lookup[badSKTypeLines, "skType"], "allowedSKTypes" -> {"++", "--", "+-", "-+"}|>]
    ];
   badStateLines = Select[
     MapIndexed[Join[#1, <|"lineIndex" -> First[#2]|>] &, topo["lines"]],
     ! MemberQ[{"full", "shrunk"}, Lookup[#, "state", "full"]] &
     ];
   If[badStateLines =!= {},
    appendIssue["error", "unknownLineStates", <|"lineIndices" -> Lookup[badStateLines, "lineIndex"], "states" -> Lookup[badStateLines, "state"], "allowedStates" -> {"full", "shrunk"}|>]
    ];
   lineMomentumVars = DeleteDuplicates[Flatten[Variables /@ Lookup[topo["lines"], "momentum"]]];
   declaredMomentumVars = ds016MomentumAtoms@Lookup[
     topo,
     "momentumDecompositionBasis",
     Join[topo["loopMomenta"], topo["effectiveLoopExternalMomenta"], Lookup[topo, "independentExternalMomenta", {}]]
     ];
   undeclaredMomentumVars = Complement[lineMomentumVars, declaredMomentumVars];
   If[undeclaredMomentumVars =!= {},
    appendIssue["error", "undeclaredMomentumVariables", <|"variables" -> undeclaredMomentumVars, "declared" -> declaredMomentumVars|>]
    ];
   nonLinearLineMomentumData = lineMomentumLinearityIssues[topo];
   If[nonLinearLineMomentumData =!= {},
    appendIssue["error", "nonLinearLineMomenta", <|"issues" -> nonLinearLineMomentumData, "comment" -> "line momenta must be linear combinations of loopMomenta, loopExternalMomenta and declared independentExternalMomenta"|>]
    ];
   nonLinearScalarProductArgumentData = scalarProductArgumentLinearityIssues[topo];
   If[nonLinearScalarProductArgumentData =!= {},
    appendIssue["error", "nonLinearScalarProductArguments", <|"issues" -> nonLinearScalarProductArgumentData, "comment" -> "sp[p,r] arguments must be linear momentum combinations before scalar products are expanded"|>]
    ];
   spData = makeScalarProductData[topo];
   If[momentumIBPRequiredQ && spData["unsupportedISPExprs"] =!= {},
    appendIssue["error", "unsupportedISPExpressions", <|"expressions" -> spData["unsupportedISPExprs"], "allowedScalarProducts" -> spData["scalarProducts"]|>]
    ];
   If[momentumIBPRequiredQ && ! TrueQ[spData["structuralCountQ"]],
    appendIssue["error", "insufficientISPData", <|"needed" -> spData["structuralNeededISPCount"], "providedDirect" -> spData["directISPCount"], "provided" -> spData["ispCount"]|>]
    ];
   If[momentumIBPRequiredQ && spData["unsupportedISPExprs"] === {} && TrueQ[spData["structuralCountQ"]] && ! TrueQ[spData["coordinateCountQ"]],
    appendIssue["error", "scalarProductCoordinateCountMismatch", <|
      "zCount" -> spData["zCount"],
      "nonISPScalarProductCount" -> spData["nonISPScalarProductCount"],
      "assumption" -> "topology input must provide a closed z/ISP coordinate system; no automatic redundant propagator subset is selected"
      |>]
    ];
   If[momentumIBPRequiredQ && spData["unsupportedISPExprs"] === {} && TrueQ[spData["structuralCountQ"]] && TrueQ[spData["coordinateCountQ"]],
    ruleData = makeScalarProductRules[topo];
    If[Lookup[ruleData, "status", "notComputed"] =!= "computed",
     appendIssue["error", "scalarProductCoordinateSolveFailed", <|
       "reason" -> Lookup[ruleData, "reason", Missing["reason"]],
       "solveVars" -> Lookup[ruleData, "solveVars", spData["solveVars"]],
       "zCount" -> spData["zCount"],
       "nonISPScalarProductCount" -> spData["nonISPScalarProductCount"]
       |>]
     ]
    ];
   pendingFeatures = unsupportedSeedFeaturesForTopology[topo];
   If[pendingFeatures =!= {},
    appendIssue["pending", "unsupportedSeedFeatures", <|"pendingFeatures" -> pendingFeatures|>]
    ];
   <|
    "status" -> If[Count[Lookup[issues, "severity", {}], "error"] == 0, "ok", "issues"],
    "errorCount" -> Count[Lookup[issues, "severity", {}], "error"],
    "warningCount" -> Count[Lookup[issues, "severity", {}], "warning"],
    "pendingCount" -> Count[Lookup[issues, "severity", {}], "pending"],
    "pendingFeatures" -> pendingFeatures,
    "issues" -> issues
    |>
   ];


topologyValidationErrorQ[report_Association] := TrueQ[Lookup[report, "errorCount", 0] > 0];
topologyValidationErrorQ[_] := False;


makeTopologyData[case_Association, OptionsPattern[]] := Module[
   {topo, topMetadata, subsetData, sectorTopos, sectorMetadataList, inputReport, validationReport},
   If[caseInputPreflightErrorQ[case],
    inputReport = caseInputRequirementReport[case];
    validationReport = caseInputErrorReport[case];
    Return[Join[case, <|
       "status" -> "invalidInput",
       "reason" -> If[inputReport["missingRequiredKeys"] =!= {}, "missingRequiredCaseKeys", "malformedCaseInput"],
       "inputRequirementReport" -> inputReport,
       "validationReport" -> validationReport,
       "sectorMetadataList" -> {},
       "indexMaps" -> <||>,
       "seedSummary" -> <||>,
       "masslessBundleCandidates" -> {},
       "precomputedShrinkSectorSummary" -> <|"status" -> "skipped"|>,
       "precomputedShrinkSectorKeys" -> {}
       |>]]
    ];
   topo = parseTopology[case];
   (* parser 已经负责报告具体输入错误；这里把失败转换为 DSInit 可识别的
      invalidInput 数据，避免 $Failed 继续流入 sector metadata producer。 *)
   If[topo === $Failed,
    validationReport = <|
      "status" -> "issues",
      "errorCount" -> 1,
      "warningCount" -> 0,
      "pendingCount" -> 0,
      "pendingFeatures" -> {},
      "issues" -> {
        <|"severity" -> "error", "code" -> "topologyParseFailed"|>
        }
      |>;
    Return[Join[case, <|
       "status" -> "invalidInput",
       "reason" -> "topologyParseFailed",
       "validationReport" -> validationReport,
       "sectorMetadataList" -> {},
       "indexMaps" -> <||>,
       "seedSummary" -> <||>,
       "masslessBundleCandidates" -> {},
       "precomputedShrinkSectorSummary" -> <|"status" -> "skipped"|>,
       "precomputedShrinkSectorKeys" -> {}
       |>]]
    ];
   topMetadata = makeSectorMetadata[topo];
   subsetData = If[TrueQ[OptionValue[PrecomputeShrinkSectorMetadata]],
     shrinkSectorSubsets[topo],
     <|"status" -> "skipped", "subsets" -> {}, "completeCoverageQ" -> False|>
     ];
   sectorTopos = If[Lookup[subsetData, "status", "skipped"] === "generated",
     shrinkSectorTopology[topo, #] & /@ Lookup[subsetData, "subsets", {}],
     {}
     ];
   sectorMetadataList = Join[{topMetadata}, Lookup[#, "sectorMetadata", makeSectorMetadata[#]] & /@ sectorTopos];
   Join[topo, <|
     "sectorMetadata" -> topMetadata,
     "sectorMetadataList" -> sectorMetadataList,
     "indexMaps" -> makeTopologyIndexMaps[topo, topMetadata],
     "seedSummary" -> makeTopologySeedSummary[topo],
     "validationReport" -> topologyValidationReport[topo],
    "tadpoleSymmetryData" -> tadpoleSymmetryData[topo],
    "effectiveSymmetryRules" -> effectiveSymmetryRules0[topo],
    "masslessBundleCandidates" -> masslessBundleCandidates[topo],
     "masslessEndpointConventions" -> masslessEndpointConventionData[topo],
     "precomputedShrinkSectorSummary" -> KeyDrop[subsetData, "subsets"],
     "precomputedShrinkSectorKeys" -> Lookup[sectorMetadataList, "sectorKey"]
     |>]
   ];

(* ::Chapter:: *)
(*统一 canonical seed 与 Kira 导出门禁*)

(* 本章只合并已经生成的 seed，并给出 Kira 前的 readiness 判断。
   若 time/momentum seed 仍有 pending features 或 forbidden n，Kira exporter 必须返回 notReady，不写文件。 *)

canonicalSeedReadyQ[batch_Association] := TrueQ[
   Lookup[batch, "status", "missing"] === "generated" &&
    Lookup[batch, "completeCanonicalQ", False] &&
    Lookup[batch, "forbiddenNData", {"missing"}] === {}
   ];


(* seed record 必须直接携带定长 sectorKey。只给收缩集合而不给总位宽时无法无歧义恢复 key，
   因此缺字段直接返回 Missing。 *)
seedEntrySourceSectorKey[entry_Association] := Lookup[entry, "sectorKey", Missing["SectorKeyRequired"]];


seedEntryIBPClass[entry_Association] := Module[{source = ToString[Lookup[entry, "source", "unknown"], InputForm]},
   Which[
    StringContainsQ[source, "Momentum", IgnoreCase -> True], "qIBP",
    StringContainsQ[source, "momentum", IgnoreCase -> True], "qIBP",
    StringContainsQ[source, "Time", IgnoreCase -> True], "tIBP",
    StringContainsQ[source, "time", IgnoreCase -> True], "tIBP",
    True, "unknownIBP"
    ]
   ];


classifyCanonicalSeedBatch[batch_Association] := Module[
   {equations = Lookup[batch, "equations", {}], tagged, grouped},
   tagged = Join[#, <|"sectorKey" -> seedEntrySourceSectorKey[#], "ibpClass" -> seedEntryIBPClass[#]|>] & /@ equations;
   grouped = GroupBy[tagged, {#sectorKey &, #ibpClass &}];
   <|
    "status" -> If[Lookup[batch, "status", "missing"] === "generated", "generated", "notGenerated"],
    "caseName" -> Lookup[batch, "caseName", Missing["caseName"]],
    "sectorKeys" -> Sort[DeleteDuplicates[Lookup[tagged, "sectorKey", {}]]],
    "classes" -> Sort[DeleteDuplicates[Lookup[tagged, "ibpClass", {}]]],
    "equationCount" -> Length[tagged],
    "bySectorThenClass" -> grouped,
    "summary" -> Association @ KeyValueMap[#1 -> Map[Length, #2] &, grouped]
    |>
   ];


canonicalGeneratorStrings[list_] := Sort[ToString[#, InputForm] & /@ list];


expectedGeneratorsBySector[batch_Association] := Module[
   {topEntry, topMetadata, topKey, shrinkSummaries, shrinkEntries},
   topMetadata = SelectFirst[
     Lookup[batch, "sectorMetadataList", {}],
     Lookup[#, "sectorShrunkLines", Missing["NoShrunkLines"]] === {} &,
     Missing["TopSectorMetadata"]
     ];
   topKey = If[Head[topMetadata] === Missing, Missing["TopSectorKey"], topMetadata["sectorKey"]];
   topEntry = topKey -> <|
      "qIBP" -> Lookup[Lookup[batch, "momentumSummary", <||>], "generators", {}],
      "tIBP" -> Lookup[Lookup[batch, "timeSummary", <||>], "generators", {}]
      |>;
   shrinkSummaries = Lookup[Lookup[batch, "shrinkSectorSummary", <||>], "sectorSummaries", {}];
   shrinkEntries = Table[
     Lookup[summary, "sectorKey", Missing["SectorKeyRequired"]] -> <|
       "qIBP" -> Lookup[Lookup[summary, "momentumSummary", <||>], "generators", {}],
       "tIBP" -> Lookup[Lookup[summary, "timeSummary", <||>], "generators", {}]
       |>,
     {summary, shrinkSummaries}
     ];
   Association[Join[{topEntry}, shrinkEntries]]
   ];


sectorGeneratorCoverageChecks[equations_List, expectedBySector_Association, sectorKeys_List] := Association @ Table[
    Module[
     {sectorEquations, actualQ, actualT, expected = Lookup[expectedBySector, sectorKey, <|"qIBP" -> {}, "tIBP" -> {}|>]},
     sectorEquations = Select[equations, seedEntrySourceSectorKey[#] === sectorKey &];
     actualQ = DeleteDuplicates @ Lookup[Select[sectorEquations, seedEntryIBPClass[#] === "qIBP" &], "generator", {}];
     actualT = DeleteDuplicates @ Lookup[Select[sectorEquations, seedEntryIBPClass[#] === "tIBP" &], "generator", {}];
     sectorKey -> <|
       "actualQGenerators" -> actualQ,
       "actualTGenerators" -> actualT,
       "expectedQGenerators" -> expected["qIBP"],
       "expectedTGenerators" -> expected["tIBP"],
       "qGeneratorCoverageQ" -> TrueQ[canonicalGeneratorStrings[actualQ] === canonicalGeneratorStrings[expected["qIBP"]]],
       "tGeneratorCoverageQ" -> TrueQ[canonicalGeneratorStrings[actualT] === canonicalGeneratorStrings[expected["tIBP"]]]
       |>
     ],
    {sectorKey, sectorKeys}
    ];


makeCanonicalSeedCoverageReport[batch_Association] := Module[
   {classified, summary, sectorKeys, equations, sectorClassChecks, topEquations, topQGenerators, topTGenerators,
     expectedQGenerators, expectedTGenerators, expectedBySector, sectorGeneratorChecks,
     totalClassifiedCount, forbiddenData, reportQ, momentumRequiredQ, ibpMode, rootLineCount},
   classified = classifyCanonicalSeedBatch[batch];
   summary = Lookup[classified, "summary", <||>];
   sectorKeys = Lookup[Lookup[batch, "sectorMetadataList", {}], "sectorKey", {}];
   equations = Lookup[batch, "equations", {}];
   ibpMode = Lookup[batch, "ibpMode", "full"];
   rootLineCount = Lookup[
     SelectFirst[Lookup[batch, "sectorMetadataList", {}], AssociationQ, <||>],
     "rootLineCount",
     0
     ];
   momentumRequiredQ = ibpMode === "full";
   sectorClassChecks = Association @ Table[
      sectorKey -> <|
        "qIBPCount" -> Lookup[Lookup[summary, sectorKey, <||>], "qIBP", 0],
        "tIBPCount" -> Lookup[Lookup[summary, sectorKey, <||>], "tIBP", 0],
        "hasBothQAndT" -> TrueQ[
          Lookup[Lookup[summary, sectorKey, <||>], "qIBP", 0] > 0 &&
           Lookup[Lookup[summary, sectorKey, <||>], "tIBP", 0] > 0
          ],
        "hasRequiredClasses" -> TrueQ[
          Lookup[Lookup[summary, sectorKey, <||>], "tIBP", 0] > 0 &&
           (! momentumRequiredQ || Lookup[Lookup[summary, sectorKey, <||>], "qIBP", 0] > 0)
          ]
        |>,
      {sectorKey, sectorKeys}
      ];
   topEquations = Select[
     equations,
     sectorKeyTopQ[seedEntrySourceSectorKey[#], ibpMode, rootLineCount] &
     ];
   topQGenerators = DeleteDuplicates @ Lookup[Select[topEquations, seedEntryIBPClass[#] === "qIBP" &], "generator", {}];
   topTGenerators = DeleteDuplicates @ Lookup[Select[topEquations, seedEntryIBPClass[#] === "tIBP" &], "generator", {}];
   expectedQGenerators = Lookup[Lookup[batch, "momentumSummary", <||>], "generators", {}];
   expectedTGenerators = Lookup[Lookup[batch, "timeSummary", <||>], "generators", {}];
   expectedBySector = expectedGeneratorsBySector[batch];
   sectorGeneratorChecks = sectorGeneratorCoverageChecks[equations, expectedBySector, sectorKeys];
   totalClassifiedCount = Total[Cases[summary, _Integer, Infinity]];
   forbiddenData = DeleteCases[Flatten[Lookup[equations, "forbiddenNData", {}]], Null];
   reportQ = TrueQ[
     Lookup[batch, "status", Missing["status"]] === "generated" &&
      Lookup[batch, "pendingFeatures", {"missing"}] === {} &&
      Lookup[batch, "forbiddenNData", {"missing"}] === {} &&
      TrueQ[Lookup[batch, "eomCanonicalQ", False]] &&
      TrueQ[And @@ Lookup[equations, "eomCanonicalQ", {False}]] &&
      forbiddenData === {} &&
      totalClassifiedCount === Lookup[batch, "equationCount", Missing["equationCount"]] &&
      Sort[Lookup[classified, "sectorKeys", {}]] === Sort[sectorKeys] &&
      FreeQ[Lookup[classified, "classes", {}], "unknownIBP"] &&
      And @@ Lookup[Values[sectorClassChecks], "hasRequiredClasses", {False}] &&
      And @@ Flatten[Lookup[Values[sectorGeneratorChecks], {"qGeneratorCoverageQ", "tGeneratorCoverageQ"}, {False}]] &&
      canonicalGeneratorStrings[topQGenerators] === canonicalGeneratorStrings[expectedQGenerators] &&
      canonicalGeneratorStrings[topTGenerators] === canonicalGeneratorStrings[expectedTGenerators]
     ];
   <|
    "status" -> If[reportQ, "ready", "notReady"],
    "caseName" -> Lookup[batch, "caseName", Missing["caseName"]],
    "passQ" -> reportQ,
    "canonicalSeedReadyQ" -> canonicalSeedReadyQ[batch],
    "sectorKeys" -> sectorKeys,
    "classes" -> Lookup[classified, "classes", {}],
    "classSummary" -> summary,
    "sectorClassChecks" -> sectorClassChecks,
    "sectorGeneratorChecks" -> sectorGeneratorChecks,
    "topQGenerators" -> topQGenerators,
    "topTGenerators" -> topTGenerators,
    "expectedQGenerators" -> expectedQGenerators,
    "expectedTGenerators" -> expectedTGenerators,
    "classifiedEquationCount" -> totalClassifiedCount,
    "equationCount" -> Lookup[batch, "equationCount", Missing["equationCount"]],
    "pendingFeatures" -> Lookup[batch, "pendingFeatures", {}],
    "forbiddenNData" -> Lookup[batch, "forbiddenNData", {}],
    "entryForbiddenNData" -> forbiddenData,
    "eomCanonicalQ" -> Lookup[batch, "eomCanonicalQ", False]
    |>
   ];



(* ::Chapter:: *)
(*seed MMA 保存与读取*)

(* seed batch 是前端产物，保存为 Mathematica 表达式；Kira exporter 只读取 linear-system 数据。 *)

Options[writeSeedBatchMMA] = {OutputDirectory -> None, SeedFileBaseName -> Automatic};


validOutputDirectoryOptionQ[value_] := value === None || value === Automatic || (StringQ[value] && StringLength[value] > 0);


validSeedFileBaseNameQ[value_] := value === Automatic || (StringQ[value] && StringLength[value] > 0);


validateSeedMMAOutputOptions[outputDir_, fileBaseName_] := Module[{issues = {}},
   If[! validOutputDirectoryOptionQ[outputDir],
    AppendTo[issues, <|"code" -> "invalidOutputDirectory", "optionKey" -> "OutputDirectory", "value" -> outputDir, "allowed" -> {"None", "Automatic", "non-empty string"}|>]
    ];
   If[! validSeedFileBaseNameQ[fileBaseName],
    AppendTo[issues, <|"code" -> "invalidSeedFileBaseName", "optionKey" -> "SeedFileBaseName", "value" -> fileBaseName, "allowed" -> {"Automatic", "non-empty string"}|>]
    ];
   <|"status" -> If[issues === {}, "ok", "invalidSeedMMAOutputOptions"], "issues" -> issues|>
   ];


validateSeedBatchForSave[batch_Association] := Module[{issues = {}},
   If[Lookup[batch, "status", Missing["status"]] =!= "generated",
    AppendTo[issues, <|"code" -> "seedBatchNotGenerated", "status" -> Lookup[batch, "status", Missing["status"]]|>]
    ];
   If[! KeyExistsQ[batch, "equations"] || ! ListQ[Lookup[batch, "equations", Missing["equations"]]],
    AppendTo[issues, <|"code" -> "seedBatchMissingEquations"|>]
    ];
   If[! KeyExistsQ[batch, "equationCount"],
    AppendTo[issues, <|"code" -> "seedBatchMissingEquationCount"|>]
    ];
   <|"status" -> If[issues === {}, "ok", "invalidSeedBatch"], "issues" -> issues|>
   ];


writeSeedBatchMMA[batch_Association, OptionsPattern[]] := Module[
   {outputDir = OptionValue[OutputDirectory], fileBaseName = OptionValue[SeedFileBaseName], outputOptionReport, batchReport, file},
   outputOptionReport = validateSeedMMAOutputOptions[outputDir, fileBaseName];
   If[outputOptionReport["status"] =!= "ok",
    Return[<|"status" -> "notWritten", "reason" -> "invalidSeedMMAOutputOptions", "outputOptionValidationReport" -> outputOptionReport|>]
    ];
   batchReport = validateSeedBatchForSave[batch];
   If[batchReport["status"] =!= "ok",
    Return[<|"status" -> "notWritten", "reason" -> "invalidSeedBatch", "seedBatchValidationReport" -> batchReport, "outputOptionValidationReport" -> outputOptionReport|>]
    ];
   If[! StringQ[outputDir], Return[<|"status" -> "notWritten", "reason" -> "OutputDirectory was not requested", "outputOptionValidationReport" -> outputOptionReport|>]];
   If[! DirectoryQ[outputDir], CreateDirectory[outputDir, CreateIntermediateDirectories -> True]];
   If[fileBaseName === Automatic, fileBaseName = Lookup[batch, "caseName", "seed_batch"] <> "_canonical_seed"];
   file = FileNameJoin[{outputDir, fileBaseName <> ".m"}];
   Put[batch, file];
   <|"status" -> "written", "file" -> file, "equationCount" -> Lookup[batch, "equationCount", Missing["equationCount"]]|>
   ];


readSeedBatchMMA[file_String] := If[FileExistsQ[file],
   Get[file],
   <|"status" -> "notRead", "reason" -> "fileNotFound", "file" -> file|>
   ];


(* ::Chapter:: *)
(*Kira user-defined system 导出*)

(* 本章按参考 codebubble 的 reduce_user_defined_system 格式导出。
   Kira exporter 的输入是已经撒点/替换参数后的线性系统数据，不直接消费 seed batch。
   userSystem/ibp.kira 每个非零方程一个 block，list 给目标积分编号，jobs.yaml 调用 Kira。
   默认只返回字符串；只有显式给 OutputDirectory 才写文件。 *)

kiraBackendSymbol[name_String] := Symbol["Global`" <> name];


(* Kira/Fermat 只接收小写原子变量；原始 Mathematica 变量可以是 kk[1,1] 这类非原子表达式。 *)
kiraCoefficientVariableMap[variables_List, preferredBackendSymbols_List : {}] := MapIndexed[
   Function[{variable, position},
    <|
     "original" -> variable,
     "backend" -> If[
       MemberQ[preferredBackendSymbols, variable] && Head[variable] === Symbol,
       SymbolName[Unevaluated[variable]],
       "dsc" <> ToString[First[position]]
       ]
     |>
    ],
   variables
   ];


(* Fermat 的 user-system 文本不直接接收 Sqrt 等分数幂。把这些代数生成元先整体视为
   扩大有理函数域中的独立原子；importer 再按 manifest 恢复原 Mathematica 表达式。 *)
kiraCoefficientAlgebraicGenerators[coefficients_List] := SortBy[
   DeleteDuplicates @ Cases[
     coefficients,
     (power : Power[_, exponent_Rational] /; Denominator[exponent] > 1) :> power,
     Infinity
     ],
   ToString[InputForm[#]] &
   ];


kiraBackendCoefficientRules[variableMap_List] := Join[
   (Lookup[#, "original"] -> kiraBackendSymbol[Lookup[#, "backend"]]) & /@ variableMap
   ];


(* Kira 不消费虚数。Complex 原子必须在 kinematic map 与积分 phase gauge 阶段消失；
   这里只映射已经通过实数化门禁的实 coefficient variables。 *)
kiraBackendCoefficientExpression[expr_, variableMap_List] :=
   expr /. kiraBackendCoefficientRules[variableMap];


kiraCoefficientString[expr_, variableMap_List] := StringReplace[
   ToString[kiraBackendCoefficientExpression[expr, variableMap], InputForm],
   WhitespaceCharacter .. -> ""
   ];
kiraCoefficientString[expr_] := kiraCoefficientString[expr, {}];


kiraBackendCoefficientStringQ[text_String] := StringMatchQ[
   text,
   RegularExpression["[a-z0-9_+\\-*/^().]+"]
   ];


kiraBackendCoefficientSyntaxReport[linearEquations_List, variableMap_List] := Module[
   {coefficients, strings, badPositions},
   coefficients = Last /@ Flatten[kiraNonzeroCoefficientRules[#["coefficientRules"]] & /@ linearEquations];
   strings = kiraCoefficientString[#, variableMap] & /@ coefficients;
   badPositions = Flatten@Position[strings, text_ /; ! kiraBackendCoefficientStringQ[text], {1}, Heads -> False];
   <|
    "status" -> If[badPositions === {}, "valid", "invalid"],
    "coefficientCount" -> Length[strings],
    "badPositions" -> badPositions,
    "badStrings" -> Lookup[AssociationThread[Range[Length[strings]] -> strings], badPositions, {}]
    |>
   ];


kiraNonzeroCoefficientRules[coefficientRules_List] := Select[
   (First[#] -> Cancel[Last[#]]) & /@ coefficientRules,
   ! TrueQ[Last[#] === 0] &
   ];


kiraEquationBlock[linearEquation_Association, coefficientRules_List, variableMap_List : {}] := Module[
   {nonzeroRules},
   nonzeroRules = kiraNonzeroCoefficientRules[coefficientRules];
   If[nonzeroRules === {}, Return[""]];
   StringRiffle[ToString[First[#]] <> "*(" <> kiraCoefficientString[Last[#], variableMap] <> ")" & /@ nonzeroRules, "\n"] <> "\n\n"
   ];


(* ::Section::Closed:: *)
(*Kira 内部虚轴运动学坐标*)

(* 顶点相位参数和 fixed massless propagator 模长都沿同一虚轴实数化；角色只从
   topology metadata 判断。cycle line 的 xi 是积分变量，不进入外部坐标列表。 *)
kiraImaginaryAxisKinematicOccurrences[topo_Association] := Module[
   {report, dependencies, vertexOccurrences, fixedMasslessOccurrences},
   report = externalLegEnergyNamingReport[topo];
   dependencies = Lookup[report, "dependencyData", <||>];
   vertexOccurrences = Flatten[
     KeyValueMap[
      Function[{vertexId, data},
       (<|
          "physical" -> #,
          "vertexId" -> vertexId,
          "role" -> "independentExternalLegEnergyParameter",
          "externalLegEnergyKind" -> Lookup[data, "kind", Missing["kind"]]
          |> &) /@ Lookup[data, "internalIndependentExternalLegEnergyParameters", {}]
      ],
      dependencies
      ],
     1
     ];
   fixedMasslessOccurrences = Map[
     Function[e,
      With[{line = topo["lines"][[e]]},
       <|
        "physical" -> fixedLineMomentumMagnitude[topo, e],
        "lineIndex" -> e,
        "lineId" -> Lookup[line, "id", e],
        "role" -> "fixedMasslessPropagatorMomentumMagnitude",
        "massType" -> "massless",
        "linePowerMode" -> Lookup[line, "linePowerMode", "fixedCoefficient"]
        |>
       ]
      ],
     Select[
      Range[Length[topo["lines"]]],
      Lookup[topo["lines"][[#]], "massType", "massive"] === "massless" &&
        ! lineIndexedPowerQ[topo["lines"][[#]]] &
      ]
     ];
   Join[vertexOccurrences, fixedMasslessOccurrences]
   ];


kiraKinematicBackendName[physical_Symbol] := "i" <> ToLowerCase[SymbolName[Unevaluated[physical]]];


kiraBackendKinematicConventionData[linearData_Association, linearEquations_List] := Module[
   {topo, occurrences, grouped, physicalAtoms, nonAtomic, backendNames, duplicateNames,
     coefficientVariables, existingNames, collisions, reservedNames = {"dsii", "ccc"}, records},
   topo = Lookup[linearData, "topology", Missing["topology"]];
   If[! parsedTopologyQ[topo],
    Return[<|"status" -> "notRequired", "scope" -> "KiraBackendOnly",
      "reason" -> "parsedTopologyUnavailable"|>]
    ];
   occurrences = kiraImaginaryAxisKinematicOccurrences[topo];
   If[occurrences === {},
     Return[<|"status" -> "notRequired", "scope" -> "KiraBackendOnly",
       "reason" -> "noImaginaryAxisKinematicAtoms"|>]
     ];
   physicalAtoms = DeleteDuplicates[Lookup[occurrences, "physical", {}]];
   nonAtomic = Select[physicalAtoms, Head[#] =!= Symbol &];
   If[nonAtomic =!= {},
     Return[<|"status" -> "invalid", "scope" -> "KiraBackendOnly",
       "reason" -> "backendKinematicCoordinatesMustBeAtomicSymbols",
       "message" -> "Kira 实数化要求每个顶点相位参数和固定无质量传播子模长都是独立符号。请用 KinematicRules 把列出的表达式绑定到独立坐标后再导出。",
       "nonAtomicPhysicalKinematics" -> nonAtomic,
       "sourceOccurrences" -> occurrences|>]
     ];
   backendNames = kiraKinematicBackendName /@ physicalAtoms;
   duplicateNames = Keys@Select[Counts[backendNames], # > 1 &];
   coefficientVariables = linearCoefficientDiagnostics[linearEquations]["coefficientVariables"];
   existingNames = SymbolName[Unevaluated[#]] & /@ Select[coefficientVariables, Head[#] === Symbol &];
   collisions = DeleteDuplicates@Join[
      duplicateNames,
      Intersection[backendNames, reservedNames],
      Intersection[backendNames, existingNames]
      ];
   If[collisions =!= {} || ! And @@ (StringMatchQ[#, RegularExpression["[a-z][a-z0-9_]*"]] & /@ backendNames),
     Return[<|"status" -> "invalid", "scope" -> "KiraBackendOnly",
       "reason" -> "backendKinematicNameCollisionOrInvalidName",
       "message" -> "Kira 后端坐标名发生冲突或不符合命名规则。请通过 KinematicRules 选用互不重复的字母数字符号。",
       "backendNames" -> backendNames, "conflictingBackendNames" -> collisions,
      "existingCoefficientVariableNames" -> existingNames,
      "sourceOccurrences" -> occurrences|>]
    ];
   grouped = GroupBy[occurrences, Lookup[#, "physical"] &];
   records = MapThread[
     Function[{physical, backendName},
      With[{backend = kiraBackendSymbol[backendName], sources = Lookup[grouped, physical, {}]},
       <|
        "physical" -> physical,
        "backend" -> backend,
         "backendName" -> backendName,
         "sourceVertexIds" -> DeleteDuplicates[Lookup[sources, "vertexId", {}]],
         "sourceLineIds" -> DeleteDuplicates[Lookup[sources, "lineId", {}]],
         "sourceRoles" -> DeleteDuplicates[Lookup[sources, "role", {}]],
        "physicalToBackendRule" -> Rule[physical, -I backend],
        "backendToPhysicalRule" -> Rule[backend, I physical],
        "physicalDerivativeFromBackendFactor" -> I,
        "backendDerivativeFromPhysicalFactor" -> -I,
        "eulerOperatorInvariantQ" -> True
        |>
       ]
      ],
     {physicalAtoms, backendNames}
     ];
   <|
    "status" -> "configured",
    "scope" -> "KiraBackendOnly",
    "physicalToBackendConvention" -> "physicalKinematic == -I backendKinematic",
    "derivativeConvention" -> "D[physicalKinematic] == I D[backendKinematic]",
    "eulerConvention" -> "physicalKinematic D[physicalKinematic] == backendKinematic D[backendKinematic]",
    "records" -> records,
    "physicalToBackendRules" -> Lookup[records, "physicalToBackendRule", {}],
    "backendToPhysicalRules" -> Lookup[records, "backendToPhysicalRule", {}],
    "physicalKinematics" -> physicalAtoms,
    "backendKinematics" -> Lookup[records, "backend", {}],
    "backendNames" -> backendNames,
    "sourceOccurrences" -> occurrences,
    "collisionCheck" -> <|"status" -> "passed", "conflictingBackendNames" -> {}|>
    |>
   ];


kiraKinematicCoefficientRulesData[rules_List, topo_Association, convention_Association] := Module[
   {status, records, physicalToBackend, normalizedRules, effectiveRules, physicalInputRules,
     normalizedPhysicalRules, physicalPointRules, backendRules, invalidNumericRules},
   status = Lookup[convention, "status", "notRequired"];
   normalizedRules = normalizeCoefficientRulesForTopology[rules, topo];
   effectiveRules = DeleteDuplicatesBy[Join[normalizedRules, rules], First];
   If[status =!= "configured",
     Return[<|"status" -> "notRequired", "userRules" -> rules,
       "normalizedPhysicalRules" -> normalizedRules, "backendRules" -> effectiveRules,
       "physicalPointRules" -> effectiveRules, "kinematicNumericRules" -> {},
       "coordinateConsistencyQ" -> True|>]
    ];
   records = Lookup[convention, "records", {}];
   physicalToBackend = Lookup[convention, "physicalToBackendRules", {}];
   invalidNumericRules = Cases[
     rules,
     rule : (Rule | RuleDelayed)[lhs_, rhs_] /;
       AnyTrue[records, SameQ[Lookup[#, "physical"], lhs] &] && ! kiraExactRationalQ[rhs] :> rule
     ];
   If[invalidNumericRules =!= {},
     Return[<|"status" -> "invalid", "reason" -> "backendKinematicNumericValuesMustBeExactRealRationals",
       "message" -> "Kira 实数化坐标的数值必须是精确实有理数。请去掉小数、复数或未定参数后重试。",
        "invalidNumericRules" -> invalidNumericRules|>]
     ];
   physicalInputRules = Map[
     Function[rule,
       With[{lhs = First[rule], rhs = Last[rule],
          record = SelectFirst[records, SameQ[Lookup[#, "physical"], First[rule]] &, Missing["notKinematic"]]},
        If[AssociationQ[record],
         Rule[lhs, -I rhs],
         Rule[lhs, rhs]
         ]
        ]
       ],
      rules
      ];
   normalizedPhysicalRules = normalizeCoefficientRulesForTopology[physicalInputRules, topo];
   physicalPointRules = DeleteDuplicatesBy[
     Join[normalizedPhysicalRules, physicalInputRules],
     First
     ];
   backendRules = Map[
     Function[rule,
       With[{lhs = First[rule], rhs = Last[rule],
          record = SelectFirst[records, SameQ[Lookup[#, "physical"], First[rule]] &, Missing["notKinematic"]]},
        If[AssociationQ[record],
         Rule[Lookup[record, "backend"], Cancel[I rhs]],
         Rule[lhs /. physicalToBackend, rhs /. physicalToBackend]
         ]
        ]
       ],
      physicalPointRules
      ];
   <|
     "status" -> "converted",
     "userRules" -> rules,
     "normalizedPhysicalRules" -> normalizedPhysicalRules,
     "backendRules" -> backendRules,
     "physicalPointRules" -> physicalPointRules,
     "kinematicNumericRules" -> Select[backendRules, MemberQ[Lookup[records, "backend", {}], First[#]] &],
     "coordinateConsistencyQ" -> True
     |>
   ];


kiraRestoreKinematicVariableIdentities[expressions_, convention_Association] := Module[{identityRules},
   identityRules = (Lookup[#, "backend"] -> Lookup[#, "physical"]) & /@ Lookup[convention, "records", {}];
   expressions /. identityRules
   ];


(* ::Section::Closed:: *)
(*Gaussian 相位有理化*)

(* 全数值关系中的系数只允许位于 Gaussian 有理数的实轴或虚轴。逐积分列相位
   与逐方程公共相位共同把它们变成严格有理数；manifest 保存列相位供 importer 反变换。 *)
kiraExactRationalQ[value_] := IntegerQ[value] || Head[value] === Rational;


kiraGaussianPairMultiply[{leftReal_, leftImaginary_}, {rightReal_, rightImaginary_}] := {
   leftReal rightReal - leftImaginary rightImaginary,
   leftReal rightImaginary + leftImaginary rightReal
   };


kiraGaussianPairPower[pair_List, exponent_Integer] := Which[
   exponent === 0, {1, 0},
   exponent > 0, Nest[kiraGaussianPairMultiply[#, pair] &, {1, 0}, exponent],
   exponent < 0,
   With[{positive = kiraGaussianPairPower[pair, -exponent]},
    {
     positive[[1]]/(positive[[1]]^2 + positive[[2]]^2),
     -positive[[2]]/(positive[[1]]^2 + positive[[2]]^2)
     }
    ]
   ];


(* Kira coefficient variables and algebraic generators are real by contract. Only actual
   Complex atoms need pair propagation; this avoids expanding every large rational function. *)
kiraGaussianRealImaginaryPair[expr_] := Which[
   FreeQ[expr, _Complex], {expr, 0},
   Head[expr] === Complex, {Re[expr], Im[expr]},
   Head[expr] === Plus, Total[kiraGaussianRealImaginaryPair /@ List @@ expr],
   Head[expr] === Times, Fold[kiraGaussianPairMultiply, {1, 0}, kiraGaussianRealImaginaryPair /@ List @@ expr],
   MatchQ[expr, Power[_, _Integer]],
   kiraGaussianPairPower[kiraGaussianRealImaginaryPair[First[expr]], Last[expr]],
   True, $Failed
   ];


kiraGaussianAxisData[value_] := Module[
   {normalized = Cancel[value], pair, real, imaginary, realZeroQ, imaginaryZeroQ},
   pair = kiraGaussianRealImaginaryPair[normalized];
   If[pair === $Failed,
    Return[<|"status" -> "invalid", "coefficient" -> value,
      "normalizedCoefficient" -> normalized, "reason" -> "unsupportedComplexCoefficientStructure"|>]
    ];
   real = Cancel[First[pair]];
   imaginary = Cancel[Last[pair]];
   realZeroQ = TrueQ[real === 0];
   imaginaryZeroQ = TrueQ[imaginary === 0];
   Which[
    Xor[realZeroQ, imaginaryZeroQ],
    <|"status" -> "valid", "axisPhase" -> If[realZeroQ, 1, 0],
      "rationalPart" -> If[realZeroQ, imaginary, real]|>,
    True,
    <|"status" -> "invalid", "coefficient" -> value, "normalizedCoefficient" -> normalized,
      "realPart" -> real, "imaginaryPart" -> imaginary|>
    ]
   ];


kiraGaussianPhaseRationalize[linearEquations_List, integralCount_Integer] := Module[
   {equationData, invalidItems = {}, invalidItemCount = 0, buildAxisItem,
    constraints, edges, adjacency, participatingIDs, invalidIDs,
    phaseByID = <||>, componentByID = <||>, conflicts = {}, componentCount = 0,
    queue, current, neighbor, expected, transformedEquations, rowPhases,
    invalidTransformed = {}, invalidTransformedCount = 0, rationalCoefficientCount,
    phaseRules, appliedQ},
   buildAxisItem[rule_] := Module[{data = kiraGaussianAxisData[Last[rule]]},
     If[Lookup[data, "status", "invalid"] =!= "valid",
      invalidItemCount++;
      If[Length[invalidItems] < 20,
       AppendTo[invalidItems, Join[<|"id" -> First[rule], "coefficient" -> Last[rule]|>, data]]
       ]
      ];
     {First[rule], Lookup[data, "axisPhase", Missing["axisPhase"]]}
     ];
   equationData = Map[
     Function[equation,
      buildAxisItem /@ kiraNonzeroCoefficientRules[equation["coefficientRules"]]
      ],
     linearEquations
     ];
   If[invalidItemCount > 0,
    Return[<|"status" -> "notRationalizable", "reason" -> "coefficientsMustBePureAxisRealRationalFunctions",
      "invalidCoefficientCount" -> invalidItemCount, "invalidCoefficients" -> invalidItems|>]
    ];
   participatingIDs = Sort@DeleteDuplicates@Cases[equationData, {id_Integer, _Integer} :> id, Infinity];
   invalidIDs = Select[participatingIDs, ! IntegerQ[#] || # < 1 || # > integralCount &];
   If[invalidIDs =!= {},
    Return[<|"status" -> "notRationalizable", "reason" -> "integralIDOutsideDeclaredRange", "invalidIntegralIDs" -> invalidIDs|>]
    ];
   constraints = Flatten[
     Map[
      Function[items,
       If[Length[items] < 2, {},
        ({First[First[items]], First[#],
            BitXor[Last[First[items]], Last[#]]} &) /@ Rest[items]
        ]
       ],
      equationData
      ],
     1
     ];
   edges = Flatten[
     ({#[[1]] -> {#[[2]], #[[3]]}, #[[2]] -> {#[[1]], #[[3]]}} &) /@ constraints,
     1
     ];
   adjacency = GroupBy[edges, First -> Last];
   Do[
    If[! KeyExistsQ[phaseByID, start],
     componentCount++;
     phaseByID[start] = 0;
     componentByID[start] = componentCount;
     queue = {start};
     While[queue =!= {},
      current = First[queue];
      queue = Rest[queue];
      Do[
       neighbor = edge[[1]];
       expected = BitXor[phaseByID[current], edge[[2]]];
       If[KeyExistsQ[phaseByID, neighbor],
        If[phaseByID[neighbor] =!= expected,
         AppendTo[conflicts, <|"sourceID" -> current, "targetID" -> neighbor,
           "relativePhase" -> edge[[2]], "sourcePhase" -> phaseByID[current],
           "targetPhase" -> phaseByID[neighbor]|>]
         ],
        phaseByID[neighbor] = expected;
        componentByID[neighbor] = componentCount;
        AppendTo[queue, neighbor]
        ],
       {edge, Lookup[adjacency, current, {}]}
       ]
      ]
     ],
    {start, participatingIDs}
    ];
   If[conflicts =!= {},
    Return[<|"status" -> "notRationalizable", "reason" -> "inconsistentIntegralPhaseConstraints",
      "constraintCount" -> Length[constraints], "componentCount" -> componentCount,
      "conflictCount" -> Length[conflicts], "conflicts" -> Take[conflicts, UpTo[20]]|>]
    ];
   Do[
    If[! KeyExistsQ[phaseByID, id], phaseByID[id] = 0],
    {id, Range[integralCount]}
    ];
   rowPhases = Map[
     Function[items,
      If[items === {}, 0,
       Mod[Last[First[items]] + phaseByID[First[First[items]]], 2]
       ]
      ],
     equationData
     ];
   transformedEquations = MapThread[
     Function[{equation, rowPhase},
      Join[equation, <|"coefficientRules" -> Map[
          Function[rule,
           First[rule] -> Cancel[Last[rule] I^phaseByID[First[rule]] I^(-rowPhase)]
           ],
          equation["coefficientRules"]
          ]|>]
      ],
     {linearEquations, rowPhases}
     ];
   transformedEquations = Map[
     Function[equation,
      Join[equation, <|"coefficientRules" -> Map[
          Function[rule,
           If[TrueQ[Last[rule] === 0], rule,
            Module[{data = kiraGaussianAxisData[Last[rule]]},
             If[Lookup[data, "status", "invalid"] === "valid" && Lookup[data, "axisPhase", 1] === 0,
              First[rule] -> Lookup[data, "rationalPart"],
              invalidTransformedCount++;
              If[Length[invalidTransformed] < 20, AppendTo[invalidTransformed, data]];
              rule
              ]
             ]
            ]
           ],
          equation["coefficientRules"]
          ]|>]
     ],
     transformedEquations
     ];
   If[invalidTransformedCount > 0,
    Return[<|"status" -> "notRationalizable", "reason" -> "phaseTransformDidNotProduceRealRationalFunctions",
      "invalidCoefficientCount" -> invalidTransformedCount,
      "invalidCoefficients" -> invalidTransformed|>]
    ];
   rationalCoefficientCount = Total[
     Length[kiraNonzeroCoefficientRules[# ["coefficientRules"]]] & /@ transformedEquations
     ];
   phaseRules = SortBy[Normal[phaseByID], First];
   appliedQ = AnyTrue[Flatten[equationData, 1], MatchQ[#, {_, 1}] &];
   <|
    "status" -> If[appliedQ, "applied", "notRequired"],
    "passQ" -> True,
    "physicalToBackendConvention" -> "J[id] == I^phase[id] Kira[id]",
    "integralCount" -> integralCount,
    "integralPhaseRules" -> phaseRules,
    "equationRowPhases" -> MapIndexed[First[#2] -> #1 &, rowPhases],
    "constraintCount" -> Length[constraints],
    "participatingIntegralCount" -> Length[participatingIDs],
    "componentCount" -> componentCount,
    "conflictCount" -> 0,
    "rationalCoefficientCount" -> rationalCoefficientCount,
    "coefficientDomain" -> "realRationalFunctions",
    "linearEquations" -> transformedEquations
    |>
   ];


kiraPureRationalCoefficientSystemQ[linearEquations_List] := Module[{coefficients},
   coefficients = Last /@ Flatten[kiraNonzeroCoefficientRules[# ["coefficientRules"]] & /@ linearEquations];
   coefficients =!= {} && And @@ (kiraExactRationalQ /@ coefficients)
   ];


kiraNumericCoefficientSystemQ[linearEquations_List] := Module[
   {coefficients},
   coefficients = Last /@ Flatten[kiraNonzeroCoefficientRules[#["coefficientRules"]] & /@ linearEquations];
   coefficients =!= {} && And @@ (TrueQ[NumericQ[#]] & /@ coefficients)
   ];


linearCoefficientDiagnostics[linearEquations_List] := Module[
   {coefficients, algebraicGenerators, ordinaryVariables, coefficientVariables},
   coefficients = Last /@ Flatten[kiraNonzeroCoefficientRules[#["coefficientRules"]] & /@ linearEquations];
   algebraicGenerators = kiraCoefficientAlgebraicGenerators[coefficients];
   ordinaryVariables = DeleteDuplicates[Variables[coefficients]];
   coefficientVariables = DeleteDuplicates[Join[algebraicGenerators, ordinaryVariables]];
   <|
    "numericCoefficientSystemQ" -> TrueQ[coefficients =!= {} && And @@ (TrueQ[NumericQ[#]] & /@ coefficients)],
    "coefficientVariables" -> coefficientVariables,
    "coefficientAlgebraicGenerators" -> algebraicGenerators
    |>
   ];


kiraAppendNumericDummyQ[setting_, numericSystemQ_] := If[setting === Automatic, TrueQ[numericSystemQ], TrueQ[setting]];


kiraDummyCoefficientString[symbol_] := If[StringQ[symbol], symbol, kiraCoefficientString[symbol]];


kiraDummyEquationBlock[id_Integer, symbol_] := ToString[id] <> "*(" <> kiraDummyCoefficientString[symbol] <> ")\n\n";


kiraTargetIntegralItems[targetSpec_] := Which[
   targetSpec === Automatic || targetSpec === All, Automatic,
   ListQ[targetSpec], targetSpec,
   True, {targetSpec}
   ];


kiraResolveTargetIntegralIDs[linearData_Association, targetSpec_, numericDummyIntegralId_] := Module[
   {maxId = linearData["integralCount"], idRules, items, resolvedPairs, missingItems, ids, invalidIds, targetIDs},
   idRules = Association[linearData["integralRules"]];
   items = kiraTargetIntegralItems[targetSpec];
   If[items === Automatic,
    targetIDs = Range[maxId],
    resolvedPairs = ({#, Which[
          IntegerQ[#], #,
          MatchQ[#, _J], Lookup[idRules, #, Missing["targetIntegralNotInSystem", #]],
          True, Missing["unsupportedTargetIntegral", #]
          ]} &) /@ items;
    missingItems = Cases[resolvedPairs, {item_, _Missing} :> item];
    ids = DeleteDuplicates[Cases[resolvedPairs, {_, id_Integer} :> id]];
    invalidIds = Select[ids, # < 1 || # > maxId &];
    If[missingItems =!= {} || invalidIds =!= {},
     Return[<|"status" -> "invalidTargetIntegrals", "targetSpec" -> targetSpec, "missingTargetItems" -> missingItems, "invalidTargetIDs" -> invalidIds, "maxIntegralID" -> maxId|>]
     ];
    targetIDs = ids
    ];
   If[numericDummyIntegralId =!= None && ! MemberQ[targetIDs, numericDummyIntegralId],
    targetIDs = Append[targetIDs, numericDummyIntegralId]
    ];
   If[targetIDs === {},
    Return[<|"status" -> "invalidTargetIntegrals", "targetSpec" -> targetSpec, "missingTargetItems" -> {}, "invalidTargetIDs" -> {}, "reason" -> "empty target list"|>]
    ];
   <|"status" -> "resolved", "targetSpec" -> targetSpec, "targetIDs" -> targetIDs, "targetIntegralCount" -> Length[targetIDs]|>
   ];


applyKiraCoefficientRulesToLinearEquation[linearEquation_Association, coeffRules_List] := Module[
   {rules, const},
   rules = linearEquation["coefficientRules"] /. coeffRules;
   rules = (First[#] -> Cancel[Last[#]]) & /@ rules;
   const = Cancel[linearEquation["constantTerm"] /. coeffRules];
   Join[linearEquation, <|"coefficientRules" -> rules, "constantTerm" -> const|>]
   ];


validReplacementRuleQ[rule_] := MatchQ[Unevaluated[rule], _Rule | _RuleDelayed];


validateCoefficientRules[coeffRules_] := Module[
   {badPositions},
   If[! ListQ[coeffRules],
    Return[<|"status" -> "invalidCoefficientRules", "reason" -> "coefficient rules must be a list of Rule or RuleDelayed entries", "coefficientRules" -> coeffRules|>]
    ];
   badPositions = Flatten @ Position[Unevaluated[coeffRules], rule_ /; ! validReplacementRuleQ[rule], {1}, Heads -> False];
   If[badPositions === {},
    <|"status" -> "ok", "coefficientRules" -> coeffRules|>,
    <|"status" -> "invalidCoefficientRules", "reason" -> "coefficient rules contain non-rule entries", "badPositions" -> badPositions, "coefficientRules" -> coeffRules|>
    ]
   ];


defaultKiraJobOptions[] := <|
   "RunInitiate" -> True,
   "RunTriangular" -> False,
   "RunBackSubstitution" -> False,
   "RunFirefly" -> True,
   "WriteKira2MathJob" -> True,
   "Kira2MathTarget" -> "list",
   "AppendNumericDummyEquation" -> Automatic,
   "NumericDummySymbol" -> "ccc",
   "WriteRunScript" -> False,
   "RunScriptName" -> "run.sh",
   "KiraCommand" -> "kira",
   "KiraParallelJobs" -> 10,
   "RunScriptDos2Unix" -> True,
   "RunScriptCleanup" -> True
   |>;


allowedKiraJobOptionKeys[] := Keys[defaultKiraJobOptions[]];


validKiraJobOptionValueQ[key_, value_] /; MemberQ[{"RunInitiate", "RunTriangular", "RunBackSubstitution", "RunFirefly", "WriteKira2MathJob", "WriteRunScript", "RunScriptDos2Unix", "RunScriptCleanup"}, key] := BooleanQ[value];
validKiraJobOptionValueQ["AppendNumericDummyEquation", value_] := value === Automatic || BooleanQ[value];
validKiraJobOptionValueQ["KiraParallelJobs", value_] := IntegerQ[value] && value > 0;
validKiraJobOptionValueQ[key_, value_] /; MemberQ[{"Kira2MathTarget", "NumericDummySymbol", "RunScriptName", "KiraCommand"}, key] := StringQ[value] && value =!= "";
validKiraJobOptionValueQ[_, _] := False;


validateKiraJobOptions[jobOptions_] := Module[
   {raw, unknownKeys, badValueData},
   If[jobOptions === Automatic,
    Return[<|"status" -> "ok", "issues" -> {}|>]
    ];
   raw = Which[
     AssociationQ[jobOptions], jobOptions,
     ListQ[jobOptions], Check[Association[jobOptions], $Failed],
     True, $Failed
     ];
   If[raw === $Failed,
    Return[<|"status" -> "invalidKiraJobOptions", "reason" -> "KiraJobOptions must be Automatic, Association, or rule list", "issues" -> {<|"code" -> "malformedKiraJobOptions", "jobOptions" -> jobOptions|>}|>]
    ];
   unknownKeys = Complement[Keys[raw], allowedKiraJobOptionKeys[]];
   badValueData = KeyValueMap[
     If[MemberQ[allowedKiraJobOptionKeys[], #1] && ! validKiraJobOptionValueQ[#1, #2],
       <|"optionKey" -> #1, "optionValue" -> #2|>,
       Nothing
       ] &,
     raw
     ];
   If[unknownKeys === {} && badValueData === {},
    <|"status" -> "ok", "issues" -> {}|>,
    <|
     "status" -> "invalidKiraJobOptions",
     "reason" -> "KiraJobOptions contain unknown keys or invalid values",
     "unknownKiraJobOptionKeys" -> unknownKeys,
     "malformedKiraJobOptionValues" -> badValueData,
     "allowedKiraJobOptionKeys" -> allowedKiraJobOptionKeys[]
     |>
    ]
   ];


validateKiraOutputDirectory[outputDir_] := If[validOutputDirectoryOptionQ[outputDir],
   <|"status" -> "ok", "issues" -> {}|>,
   <|"status" -> "invalidOutputDirectory", "issues" -> {<|"code" -> "invalidOutputDirectory", "optionKey" -> "OutputDirectory", "value" -> outputDir, "allowed" -> {"None", "Automatic", "non-empty string"}|>}|>
   ];


normalizeKiraJobOptions[Automatic] := defaultKiraJobOptions[];
normalizeKiraJobOptions[opts_Association] := Join[defaultKiraJobOptions[], opts];
normalizeKiraJobOptions[opts_List] := normalizeKiraJobOptions[Association[opts]];
normalizeKiraJobOptions[_] := defaultKiraJobOptions[];


kiraYAMLBool[value_] := If[TrueQ[value], "true", "false"];


kiraJobsYAML[jobOptions_: Automatic] := Module[
   {opts = normalizeKiraJobOptions[jobOptions], text},
   text = "jobs:\n" <>
     "  - reduce_user_defined_system:\n" <>
     "      input_system: {files: [\"userSystem\"], config: false}\n" <>
     "      select_integrals:\n" <>
      "        select_mandatory_list:\n" <>
      "          - [list]\n" <>
      "      run_initiate: " <> kiraYAMLBool[Lookup[opts, "RunInitiate", True]] <> "\n" <>
      "      run_triangular: " <> kiraYAMLBool[Lookup[opts, "RunTriangular", False]] <> "\n" <>
      "      run_back_substitution: " <> kiraYAMLBool[Lookup[opts, "RunBackSubstitution", False]] <> "\n" <>
      "      run_firefly: " <> kiraYAMLBool[Lookup[opts, "RunFirefly", True]] <> "\n";
   If[TrueQ[Lookup[opts, "WriteKira2MathJob", True]],
    text = text <>
      "  - kira2math:\n" <>
      "      target:\n" <>
      "       - [" <> ToString[Lookup[opts, "Kira2MathTarget", "list"]] <> "]\n"
    ];
   text
   ];


kiraRunCommand[opts_Association] := Module[
   {cmd = ToString[Lookup[opts, "KiraCommand", "kira"]], parallel = Lookup[opts, "KiraParallelJobs", 10]},
   If[IntegerQ[parallel] && parallel > 0,
    cmd <> " --parallel=" <> ToString[parallel] <> " jobs.yaml",
    cmd <> " jobs.yaml"
    ]
   ];


kiraRunScript[jobOptions_: Automatic] := Module[
   {opts = normalizeKiraJobOptions[jobOptions], lines = {}},
   If[TrueQ[Lookup[opts, "RunScriptCleanup", True]],
    lines = Join[lines, {
        "rm -f kira.log",
        "rm -f firefly.log",
        "rm -rf results",
        "rm -rf tmp",
        "rm -rf sectormappings",
        "rm -rf firefly_saves",
        "rm -rf ff_save"
        }]
    ];
   If[TrueQ[Lookup[opts, "RunScriptDos2Unix", True]],
    lines = Join[lines, {
        "dos2unix list",
        "dos2unix userSystem/ibp.kira",
        "dos2unix jobs.yaml"
        }]
    ];
   lines = Append[lines, kiraRunCommand[opts]];
   StringRiffle[lines, "\n"] <> "\n"
   ];


makeKiraInputStrings[linearData_Association, coeffRules_ : {}, jobOptions_: Automatic, targetSpec_: Automatic] := Module[
   {linearEquations, badEquations, exportedEquations, ibpText, listText, jobsText, runScriptText, repKira2JText, repJ2KiraText, metadataText,
    normalizedJobOptions, jobOptionReport, coefficientRuleReport, topologyReport, requiredKeys, missingKeys, coefficientDiagnostics, numericCoefficientSystemQ, appendNumericDummyQ,
    numericDummyIntegralId, targetData, targetIntegralCount, kiraBlockCount, numericDummySymbol,
    rawCoeffRules, normalizedCoeffRules, effectiveCoeffRules, backendEffectiveCoeffRules,
    kinematicConvention, kinematicRuleData, physicalCoefficientVariables, physicalCoefficientAlgebraicGenerators,
    coefficientVariableMap, backendCoefficientVariables, imaginaryUnitUsedQ, backendSyntaxReport,
    gaussianPhaseGauge, gaussianPhaseGaugeManifest, pureRationalBackendQ, backendTextAudit},
   topologyReport = Lookup[linearData, "topologyValidationReport", Missing["NoTopologyValidationReport"]];
   If[Lookup[linearData, "status", "missing"] =!= "generated",
    Return[<|"status" -> "notGenerated", "reason" -> "linear data missing", "topologyValidationReport" -> topologyReport|>]
    ];
   If[topologyValidationErrorQ[topologyReport],
    Return[<|"status" -> "invalidTopology", "reason" -> "topology validation has errors", "topologyValidationReport" -> topologyReport|>]
    ];
   requiredKeys = {"linearEquations", "integralRules", "integralCount", "equationCount"};
   missingKeys = Select[requiredKeys, ! KeyExistsQ[linearData, #] &];
   If[missingKeys =!= {},
    Return[<|"status" -> "notLinearSystem", "reason" -> "Kira exporter expects makeLinearSystemData output", "topologyValidationReport" -> topologyReport, "missingKeys" -> missingKeys|>]
    ];
   jobOptionReport = validateKiraJobOptions[jobOptions];
   If[Lookup[jobOptionReport, "status", "ok"] =!= "ok",
    Return[Join[jobOptionReport, <|"topologyValidationReport" -> topologyReport|>]]
    ];
   rawCoeffRules = coeffRules;
   coefficientRuleReport = validateCoefficientRules[rawCoeffRules];
   If[Lookup[coefficientRuleReport, "status", "ok"] =!= "ok",
    Return[Join[coefficientRuleReport, <|"topologyValidationReport" -> topologyReport|>]]
    ];
   normalizedCoeffRules = normalizeCoefficientRulesForLinearData[rawCoeffRules, linearData];
   (* linearData 可能已投影到公开坐标，也可能仍含内部 Gram 原子；两种规则必须在
      同一点并用，避免 s11 只被规范化成 kk 后反而漏掉公开 Sqrt[s11]。 *)
   effectiveCoeffRules = DeleteDuplicatesBy[
     Join[normalizedCoeffRules, rawCoeffRules],
     First
     ];
   kinematicConvention = kiraBackendKinematicConventionData[linearData, linearData["linearEquations"]];
   If[Lookup[kinematicConvention, "status", "invalid"] === "invalid",
    dsErrorPrint[Lookup[kinematicConvention, "message",
      "Kira 后端运动学坐标没有通过实数化门禁。请检查顶点相位参数和固定无质量传播子模长。"]];
    Return[<|"status" -> "invalidBackendKinematicConvention",
      "reason" -> Lookup[kinematicConvention, "reason", "invalidBackendKinematicConvention"],
      "message" -> Lookup[kinematicConvention, "message", Missing["message"]],
      "backendKinematicConvention" -> kinematicConvention,
      "topologyValidationReport" -> topologyReport|>]
    ];
   kinematicRuleData = kiraKinematicCoefficientRulesData[rawCoeffRules, linearData["topology"], kinematicConvention];
   If[Lookup[kinematicRuleData, "status", "invalid"] === "invalid",
    dsErrorPrint[Lookup[kinematicRuleData, "message",
      "Kira 后端运动学坐标的数值规则无效。请使用精确实有理数。"]];
    Return[<|"status" -> "invalidBackendKinematicNumericRules",
      "reason" -> Lookup[kinematicRuleData, "reason", "invalidBackendKinematicNumericRules"],
      "message" -> Lookup[kinematicRuleData, "message", Missing["message"]],
      "backendKinematicConvention" -> kinematicConvention,
      "backendKinematicRuleData" -> kinematicRuleData,
      "topologyValidationReport" -> topologyReport|>]
    ];
   backendEffectiveCoeffRules = Lookup[kinematicRuleData, "backendRules", effectiveCoeffRules];
   linearEquations = applyKiraCoefficientRulesToLinearEquation[
       #,
       Lookup[kinematicConvention, "physicalToBackendRules", {}]
       ] & /@ linearData["linearEquations"];
   linearEquations = applyKiraCoefficientRulesToLinearEquation[#, backendEffectiveCoeffRules] & /@ linearEquations;
   badEquations = Select[linearEquations, ! TrueQ[#["linearQ"]] || ! TrueQ[#["constantTerm"] === 0] &];
   If[badEquations =!= {},
    Return[<|"status" -> "notLinear", "topologyValidationReport" -> topologyReport, "badEquationCount" -> Length[badEquations], "badEquations" -> badEquations|>]
    ];
   exportedEquations = Select[linearEquations, kiraNonzeroCoefficientRules[#["coefficientRules"]] =!= {} &];
   If[exportedEquations === {},
    Return[<|"status" -> "emptySystem", "reason" -> "all linear equations are zero after coefficient rules", "topologyValidationReport" -> topologyReport|>]
    ];
   normalizedJobOptions = normalizeKiraJobOptions[jobOptions];
   coefficientDiagnostics = linearCoefficientDiagnostics[exportedEquations];
   gaussianPhaseGauge = kiraGaussianPhaseRationalize[
     exportedEquations,
     linearData["integralCount"]
     ];
   If[MemberQ[{"notRationalizable"}, Lookup[gaussianPhaseGauge, "status", "notRationalizable"]],
    Return[Join[
      KeyDrop[gaussianPhaseGauge, {"linearEquations"}],
      <|"status" -> "invalidGaussianPhaseGauge", "topologyValidationReport" -> topologyReport|>
      ]]
    ];
   If[MemberQ[{"applied", "notRequired"}, Lookup[gaussianPhaseGauge, "status", "notApplicable"]],
    exportedEquations = gaussianPhaseGauge["linearEquations"]
    ];
   gaussianPhaseGaugeManifest = KeyDrop[gaussianPhaseGauge, {"linearEquations"}];
   coefficientDiagnostics = linearCoefficientDiagnostics[exportedEquations];
   pureRationalBackendQ = kiraPureRationalCoefficientSystemQ[exportedEquations] &&
     coefficientDiagnostics["coefficientVariables"] === {};
   If[pureRationalBackendQ,
    normalizedJobOptions = Join[normalizedJobOptions, <|
       "RunTriangular" -> True,
       "RunBackSubstitution" -> True,
       "RunFirefly" -> False,
       "AppendNumericDummyEquation" -> False
       |>]
    ];
   physicalCoefficientVariables = kiraRestoreKinematicVariableIdentities[
     coefficientDiagnostics["coefficientVariables"],
     kinematicConvention
     ];
   physicalCoefficientAlgebraicGenerators = kiraRestoreKinematicVariableIdentities[
     coefficientDiagnostics["coefficientAlgebraicGenerators"],
     kinematicConvention
     ];
   coefficientVariableMap = kiraCoefficientVariableMap[
     coefficientDiagnostics["coefficientVariables"],
     Lookup[kinematicConvention, "backendKinematics", {}]
     ];
   imaginaryUnitUsedQ = ! FreeQ[
      Last /@ Flatten[kiraNonzeroCoefficientRules[#["coefficientRules"]] & /@ exportedEquations],
      _Complex
      ];
   If[imaginaryUnitUsedQ,
    Return[<|
      "status" -> "invalidRealBackendCoefficients",
      "reason" -> "Kira coefficients must be real after massless momentum mapping and integral phase gauge",
      "gaussianPhaseGauge" -> gaussianPhaseGaugeManifest,
      "topologyValidationReport" -> topologyReport
      |>]
    ];
   backendCoefficientVariables = Lookup[coefficientVariableMap, "backend", {}];
   backendSyntaxReport = kiraBackendCoefficientSyntaxReport[exportedEquations, coefficientVariableMap];
   If[Lookup[backendSyntaxReport, "status", "invalid"] =!= "valid",
    Return[<|
      "status" -> "invalidBackendCoefficientSyntax",
      "reason" -> "Kira/Fermat coefficients must contain only lowercase atomic variables and rational operators",
      "coefficientVariableMap" -> coefficientVariableMap,
      "backendCoefficientSyntaxReport" -> backendSyntaxReport,
      "topologyValidationReport" -> topologyReport
      |>]
    ];
   numericCoefficientSystemQ = coefficientDiagnostics["numericCoefficientSystemQ"];
   appendNumericDummyQ = If[pureRationalBackendQ,
     False,
     kiraAppendNumericDummyQ[Lookup[normalizedJobOptions, "AppendNumericDummyEquation", Automatic], numericCoefficientSystemQ]
     ];
   numericDummyIntegralId = If[appendNumericDummyQ, linearData["integralCount"] + 1, None];
   targetData = kiraResolveTargetIntegralIDs[linearData, targetSpec, numericDummyIntegralId];
   If[Lookup[targetData, "status", "missing"] =!= "resolved",
    Return[Join[targetData, <|"topologyValidationReport" -> topologyReport|>]]
    ];
   targetIntegralCount = targetData["targetIntegralCount"];
   kiraBlockCount = Length[exportedEquations] + If[appendNumericDummyQ, 1, 0];
   numericDummySymbol = Lookup[normalizedJobOptions, "NumericDummySymbol", "ccc"];
   ibpText = StringJoin[kiraEquationBlock[#, #["coefficientRules"], coefficientVariableMap] & /@ exportedEquations] <>
     If[appendNumericDummyQ, kiraDummyEquationBlock[numericDummyIntegralId, numericDummySymbol], ""];
   backendTextAudit = <|
     "status" -> If[
       ! StringContainsQ[ibpText, "dsii"] && ! StringContainsQ[ibpText, "Complex"] &&
        ! StringContainsQ[ibpText, numericDummySymbol] &&
        If[pureRationalBackendQ, backendCoefficientVariables === {}, True],
       "passed",
       "failed"
       ],
     "pureRationalBackendQ" -> pureRationalBackendQ,
     "containsBackendImaginaryUnitQ" -> StringContainsQ[ibpText, "dsii"],
     "containsComplexTokenQ" -> StringContainsQ[ibpText, "Complex"],
     "containsNumericDummySymbolQ" -> StringContainsQ[ibpText, numericDummySymbol],
     "backendCoefficientVariables" -> backendCoefficientVariables
     |>;
   If[Lookup[backendTextAudit, "status", "failed"] =!= "passed",
    Return[<|"status" -> "invalidRealBackendText",
      "reason" -> "Kira text must be real and must not contain dsii, Complex, or a numeric dummy",
      "gaussianPhaseGauge" -> gaussianPhaseGaugeManifest,
      "backendTextAudit" -> backendTextAudit,
      "topologyValidationReport" -> topologyReport|>]
    ];
   listText = StringRiffle[ToString /@ targetData["targetIDs"], "\n"] <> "\n";
   jobsText = kiraJobsYAML[normalizedJobOptions];
   runScriptText = If[TrueQ[Lookup[normalizedJobOptions, "WriteRunScript", False]], kiraRunScript[normalizedJobOptions], Missing["RunScriptDisabled"]];
   (* backend 只消费整数 ID；反向映射保留 loop J 或带 sector 的 tree token 原对象。 *)
   repKira2JText = ToString[InputForm[
       linearData["integralRules"] /. Rule[integral_, id_Integer] :> Rule[Tuserweight[id], integral]
       ]] <> "\n";
   repJ2KiraText = ToString[InputForm[linearData["integralRules"]]] <> "\n";
   metadataText = ToString[InputForm[
       Join[
        KeyDrop[linearData, {"integralList", "integralRules", "linearEquations", "topology"}],
        <|
         "exportedEquationCount" -> Length[exportedEquations],
         "kiraBlockCount" -> kiraBlockCount,
         "targetIntegralCount" -> targetIntegralCount,
         "targetIntegralIDs" -> targetData["targetIDs"],
         "kiraTargetIntegrals" -> targetSpec,
         "numericCoefficientSystemQ" -> numericCoefficientSystemQ,
          "coefficientVariables" -> physicalCoefficientVariables,
          "coefficientAlgebraicGenerators" -> physicalCoefficientAlgebraicGenerators,
          "backendExpressionVariables" -> coefficientDiagnostics["coefficientVariables"],
          "coefficientVariableMap" -> coefficientVariableMap,
         "backendCoefficientVariables" -> backendCoefficientVariables,
          "backendImaginaryUnit" -> None,
          "backendCoefficientSyntaxReport" -> backendSyntaxReport,
           "gaussianPhaseGauge" -> gaussianPhaseGaugeManifest,
           "backendKinematicConvention" -> kinematicConvention,
           "backendKinematicRuleData" -> kinematicRuleData,
          "pureRationalBackendQ" -> pureRationalBackendQ,
          "backendTextAudit" -> backendTextAudit,
          "numericDummyAppendedQ" -> appendNumericDummyQ,
         "numericDummyIntegralId" -> numericDummyIntegralId,
          "kiraCoefficientRules" -> backendEffectiveCoeffRules,
          "physicalCoefficientRules" -> Lookup[kinematicRuleData, "physicalPointRules", effectiveCoeffRules],
         "normalizedKiraCoefficientRules" -> normalizedCoeffRules,
         "userKiraCoefficientRules" -> rawCoeffRules,
         "kiraJobOptions" -> normalizedJobOptions
         |>
        ]
       ]] <> "\n";
   <|
    "status" -> "generated",
    "topologyValidationReport" -> topologyReport,
    "ibp.kira" -> ibpText,
    "list" -> listText,
    "jobs.yaml" -> jobsText,
    "runScriptName" -> Lookup[normalizedJobOptions, "RunScriptName", "run.sh"],
    "run.sh" -> runScriptText,
    "result/repkira2J.m" -> repKira2JText,
    "result/repJ2kira.m" -> repJ2KiraText,
    "result/kira_export_metadata.m" -> metadataText,
    "equationCount" -> linearData["equationCount"],
    "exportedEquationCount" -> Length[exportedEquations],
    "kiraBlockCount" -> kiraBlockCount,
    "integralCount" -> linearData["integralCount"],
    "targetIntegralCount" -> targetIntegralCount,
    "targetIntegralIDs" -> targetData["targetIDs"],
    "kiraTargetIntegrals" -> targetSpec,
    "numericCoefficientSystemQ" -> numericCoefficientSystemQ,
     "coefficientVariables" -> physicalCoefficientVariables,
     "coefficientAlgebraicGenerators" -> physicalCoefficientAlgebraicGenerators,
     "backendExpressionVariables" -> coefficientDiagnostics["coefficientVariables"],
    "coefficientVariableMap" -> coefficientVariableMap,
    "backendCoefficientVariables" -> backendCoefficientVariables,
     "backendImaginaryUnit" -> None,
     "backendCoefficientSyntaxReport" -> backendSyntaxReport,
      "gaussianPhaseGauge" -> gaussianPhaseGaugeManifest,
      "backendKinematicConvention" -> kinematicConvention,
      "backendKinematicRuleData" -> kinematicRuleData,
     "pureRationalBackendQ" -> pureRationalBackendQ,
     "backendTextAudit" -> backendTextAudit,
     "numericDummyAppendedQ" -> appendNumericDummyQ,
    "numericDummyIntegralId" -> numericDummyIntegralId,
     "kiraCoefficientRulesApplied" -> backendEffectiveCoeffRules,
     "physicalCoefficientRulesApplied" -> Lookup[kinematicRuleData, "physicalPointRules", effectiveCoeffRules],
    "normalizedKiraCoefficientRulesApplied" -> normalizedCoeffRules,
    "userKiraCoefficientRulesApplied" -> rawCoeffRules
    |>
   ];


(* Kira 的 Linux user-system parser 会把 CRLF 空行误读为非空方程；所有后端文本固定写成无 BOM UTF-8/LF。 *)
writeKiraUTF8LFText[path_String, text_String] := Module[{stream, normalized, bytes},
   normalized = StringReplace[text, {"\r\n" -> "\n", "\r" -> "\n"}];
   bytes = ToCharacterCode[normalized, "UTF8"];
   stream = OpenWrite[path, BinaryFormat -> True];
   BinaryWrite[stream, bytes, "Byte"];
   Close[stream];
   path
   ];


writeKiraInputFiles[outputDir_String, strings_Association] := Module[
   {userSystemDir, resultDir, runScriptName, files},
   userSystemDir = FileNameJoin[{outputDir, "userSystem"}];
   resultDir = FileNameJoin[{outputDir, "result"}];
   If[! DirectoryQ[outputDir], CreateDirectory[outputDir, CreateIntermediateDirectories -> True]];
   If[! DirectoryQ[userSystemDir], CreateDirectory[userSystemDir, CreateIntermediateDirectories -> True]];
   If[! DirectoryQ[resultDir], CreateDirectory[resultDir, CreateIntermediateDirectories -> True]];
   runScriptName = Lookup[strings, "runScriptName", "run.sh"];
   files = <|
     FileNameJoin[{userSystemDir, "ibp.kira"}] -> strings["ibp.kira"],
     FileNameJoin[{outputDir, "list"}] -> strings["list"],
     FileNameJoin[{outputDir, "jobs.yaml"}] -> strings["jobs.yaml"],
     FileNameJoin[{resultDir, "repkira2J.m"}] -> strings["result/repkira2J.m"],
     FileNameJoin[{resultDir, "repJ2kira.m"}] -> strings["result/repJ2kira.m"],
     FileNameJoin[{resultDir, "kira_export_metadata.m"}] -> strings["result/kira_export_metadata.m"]
     |>;
   If[StringQ[Lookup[strings, "run.sh", Missing["NoRunScript"]]],
    files = Join[files, <|FileNameJoin[{outputDir, runScriptName}] -> strings["run.sh"]|>]
    ];
   KeyValueMap[writeKiraUTF8LFText, files];
   Keys[files]
   ];


Options[makeKiraExportData] = {
   OutputDirectory -> None,
   KiraCoefficientRules -> {},
   KiraTargetIntegrals -> Automatic,
   KiraJobOptions -> Automatic
   };

makeKiraExportData::notlinearinput = "Kira 导出只接受 linear-system 数据，不直接接受 seed batch：`1`。";
makeKiraExportData::badlinear = "linear-system 不能导出 Kira：`1`。";


makeKiraExportData[linearData_Association, OptionsPattern[]] := Module[
   {linearForExport, strings, outputDir, outputDirReport, filesWritten, topologyReport},
   topologyReport = Lookup[linearData, "topologyValidationReport", Missing["NoTopologyValidationReport"]];
   outputDir = OptionValue[OutputDirectory];
   outputDirReport = validateKiraOutputDirectory[outputDir];
   If[outputDirReport["status"] =!= "ok",
    Return[<|
      "status" -> "notReady",
      "caseName" -> Lookup[linearData, "caseName", Missing["caseName"]],
      "topologyValidationReport" -> topologyReport,
      "reason" -> "invalidOutputDirectory",
      "kiraInput" -> outputDirReport
      |>]
    ];
   If[! KeyExistsQ[linearData, "linearEquations"],
    Message[makeKiraExportData::notlinearinput, Lookup[linearData, "caseName", Missing["caseName"]]];
    Return[<|
      "status" -> "notReady",
      "caseName" -> Lookup[linearData, "caseName", Missing["caseName"]],
      "topologyValidationReport" -> topologyReport,
      "reason" -> "Kira exporter expects linear-system data; save seed batch as MMA first, then call makeLinearSystemData after numeric/sampling choices"
      |>]
    ];
   (* serializer 只读取 linearData 的现有顺序；显式重排必须在此边界之前完成。 *)
   linearForExport = linearData;
   strings = makeKiraInputStrings[linearForExport, OptionValue[KiraCoefficientRules], OptionValue[KiraJobOptions], OptionValue[KiraTargetIntegrals]];
   If[Lookup[strings, "status", "missing"] =!= "generated",
    If[! MemberQ[{"invalidKiraJobOptions", "invalidCoefficientRules"}, Lookup[strings, "status", Missing["status"]]],
     Message[makeKiraExportData::badlinear, Lookup[strings, "status", Missing["status"]]]
     ];
    Return[<|"status" -> "notReady", "caseName" -> linearForExport["caseName"], "topologyValidationReport" -> Lookup[linearForExport, "topologyValidationReport", topologyReport], "reason" -> "linear system is not exportable", "linearSystem" -> linearForExport, "kiraInput" -> strings|>]
    ];
   filesWritten = If[StringQ[outputDir], writeKiraInputFiles[outputDir, strings], {}];
   <|
    "status" -> "ready",
    "caseName" -> linearForExport["caseName"],
    "topologyValidationReport" -> Lookup[linearForExport, "topologyValidationReport", topologyReport],
    "linearSystem" -> linearForExport,
    "kiraInput" -> strings,
    "equationCount" -> Lookup[linearForExport, "equationCount", Missing["equationCount"]],
    "exportedEquationCount" -> Lookup[strings, "exportedEquationCount", Missing["exportedEquationCount"]],
    "kiraBlockCount" -> Lookup[strings, "kiraBlockCount", Missing["kiraBlockCount"]],
    "integralCount" -> Lookup[linearForExport, "integralCount", Missing["integralCount"]],
    "kiraOrderingReport" -> Lookup[linearForExport, "kiraOrderingReport", <||>],
    "manualIntegralOrderReport" -> Lookup[linearForExport, "manualIntegralOrderReport", <||>],
    "targetIntegralCount" -> Lookup[strings, "targetIntegralCount", Missing["targetIntegralCount"]],
    "targetIntegralIDs" -> Lookup[strings, "targetIntegralIDs", Missing["targetIntegralIDs"]],
    "kiraTargetIntegrals" -> Lookup[strings, "kiraTargetIntegrals", Missing["kiraTargetIntegrals"]],
    "numericCoefficientSystemQ" -> Lookup[strings, "numericCoefficientSystemQ", Missing["numericCoefficientSystemQ"]],
    "coefficientVariables" -> Lookup[strings, "coefficientVariables", Missing["coefficientVariables"]],
    "coefficientAlgebraicGenerators" -> Lookup[strings, "coefficientAlgebraicGenerators", {}],
    "backendExpressionVariables" -> Lookup[strings, "backendExpressionVariables", {}],
    "coefficientVariableMap" -> Lookup[strings, "coefficientVariableMap", {}],
     "backendCoefficientVariables" -> Lookup[strings, "backendCoefficientVariables", {}],
     "backendImaginaryUnit" -> Lookup[strings, "backendImaginaryUnit", None],
     "backendCoefficientSyntaxReport" -> Lookup[strings, "backendCoefficientSyntaxReport", <||>],
     "gaussianPhaseGauge" -> Lookup[strings, "gaussianPhaseGauge", <|"status" -> "notApplicable"|>],
     "backendKinematicConvention" -> Lookup[strings, "backendKinematicConvention", <|"status" -> "notRequired"|>],
     "backendKinematicRuleData" -> Lookup[strings, "backendKinematicRuleData", <|"status" -> "notRequired"|>],
     "physicalCoefficientRulesApplied" -> Lookup[strings, "physicalCoefficientRulesApplied", {}],
     "pureRationalBackendQ" -> TrueQ[Lookup[strings, "pureRationalBackendQ", False]],
     "backendTextAudit" -> Lookup[strings, "backendTextAudit", <|"status" -> "notRun"|>],
     "numericDummyAppendedQ" -> Lookup[strings, "numericDummyAppendedQ", Missing["numericDummyAppendedQ"]],
    "numericDummyIntegralId" -> Lookup[strings, "numericDummyIntegralId", Missing["numericDummyIntegralId"]],
    "outputDirectory" -> outputDir,
    "writeFilesQ" -> StringQ[outputDir],
    "filesWritten" -> filesWritten
    |>
   ];


(* ::Chapter:: *)
(*线性系统中间格式*)

(* 本章把 seed 方程转成后端可消费的线性系统数据。
   这里只做积分抽取、编号和逐项系数收集，不做求解、排序优化或 Kira 文件语法假设。 *)

integralObjectsInExpression[expr_] := DeleteDuplicates[Cases[expr, _J, {0, Infinity}]];


integralObjectsInBatch[batch_Association] := DeleteDuplicates[
   Flatten[integralObjectsInExpression /@ Lookup[Lookup[batch, "equations", {}], "equation", {}]]
   ];


makeIntegralIndex[integrals_List] := AssociationThread[integrals, Range[Length[integrals]]];


linearTermData[term_, integralIndex_Association] := Module[
   {integrals, integral, coeff},
   integrals = Cases[term, _J, {0, Infinity}];
   Which[
    Length[integrals] == 0,
    <|"kind" -> "constant", "term" -> term|>,
    Length[DeleteDuplicates[integrals]] == 1,
    integral = First[DeleteDuplicates[integrals]];
    coeff = Cancel[term/integral];
    If[FreeQ[coeff, _J] && KeyExistsQ[integralIndex, integral],
     <|"kind" -> "linear", "id" -> integralIndex[integral], "integral" -> integral, "coefficient" -> coeff|>,
     <|"kind" -> "nonlinear", "term" -> term|>
     ],
    True,
    <|"kind" -> "nonlinear", "term" -> term|>
    ]
   ];


combineLinearCoefficientRules[linearTerms_List] := Normal @ Merge[
    (Rule[#["id"], #["coefficient"]] & /@ linearTerms),
    Total
    ];

reindexLinearEquation[eq_Association, oldIdToNewId_Association] := Module[
   {rules},
   rules = Select[
     (Lookup[oldIdToNewId, First[#], Missing["DroppedIntegral"]] -> Last[#]) & /@ Lookup[eq, "coefficientRules", {}],
     Head[First[#]] =!= Missing &
     ];
   Join[eq, <|"coefficientRules" -> Normal@Merge[rules, Total]|>]
   ];


reorderLinearSystemIntegrals[linearData_Association, integralOrder_List] := Module[
   {oldIntegrals, requested, newIntegrals, oldRules, oldIdToIntegral, newIndex, oldIdToNewId,
    newLinearEquations, metadataList, orderingSpec, manualOrderReport},
   If[Lookup[linearData, "status", "missing"] =!= "generated" || ! KeyExistsQ[linearData, "integralList"], Return[linearData]];
   oldIntegrals = linearData["integralList"];
   requested = normaliseIntegralOrder[integralOrder, oldIntegrals];
   newIntegrals = Join[requested, Select[oldIntegrals, ! MemberQ[requested, #] &]];
   oldRules = Association[linearData["integralRules"]];
   oldIdToIntegral = Association[Reverse /@ Normal[oldRules]];
   newIndex = makeIntegralIndex[newIntegrals];
   oldIdToNewId = Association @ KeyValueMap[#1 -> newIndex[#2] &, oldIdToIntegral];
   newLinearEquations = reindexLinearEquation[#, oldIdToNewId] & /@ linearData["linearEquations"];
   metadataList = Lookup[linearData, "sectorMetadataList", {}];
   orderingSpec = Join[Lookup[linearData, "kiraOrdering", <||>], <|"ManualIntegralOrder" -> requested|>];
   manualOrderReport = <|
     "requestedIntegrals" -> integralOrder,
     "matchedIntegrals" -> requested,
     "missingIntegralOrderItems" -> DeleteDuplicates @ missingIntegralOrderItems[integralOrder, oldIntegrals],
     "allRequestedIntegralsMatchedQ" -> TrueQ[missingIntegralOrderItems[integralOrder, oldIntegrals] === {}]
     |>;
   Join[linearData, <|
     "integralList" -> newIntegrals,
     "integralRules" -> Normal[newIndex],
     "kiraOrdering" -> orderingSpec,
     "kiraOrderingReport" -> kiraOrderingMatchReport[orderingSpec, newIntegrals],
     "manualIntegralOrderReport" -> manualOrderReport,
     "integralSortKeys" -> (integralSortKey[#, orderingSpec, metadataList] & /@ newIntegrals),
     "integralMetadata" -> integralMetadataList[newIntegrals, metadataList, orderingSpec],
     "linearEquations" -> newLinearEquations
     |>]
   ];


linearizeSeedEquation[entry_Association, integralIndex_Association] := Module[
   {terms, termData, linearPieces, constantTerms, nonlinearTerms},
   terms = linearTerms[Lookup[entry, "equation", 0]];
   termData = linearTermData[#, integralIndex] & /@ terms;
   linearPieces = Select[termData, #["kind"] === "linear" &];
   constantTerms = Lookup[Select[termData, #["kind"] === "constant" &], "term", {}];
   nonlinearTerms = Lookup[Select[termData, #["kind"] === "nonlinear" &], "term", {}];
   Join[
    KeyDrop[entry, "equation"],
    <|
     "coefficientRules" -> combineLinearCoefficientRules[linearPieces],
     "constantTerm" -> Total[constantTerms],
     "nonlinearTerms" -> nonlinearTerms,
     "linearQ" -> TrueQ[nonlinearTerms === {}]
     |>
    ]
   ];


(* 来源 metadata 用于 producer coverage；去重键只包含后端真正消费的数学方程。 *)
linearEquationMathematicalKey[entry_Association] := {
   Lookup[entry, "coefficientRules", {}],
   Lookup[entry, "constantTerm", 0],
   Lookup[entry, "nonlinearTerms", {}]
   };


Options[makeLinearSystemData] = {KiraOrdering -> Automatic, AuditLevel -> "standard"};


dsSeedArtifactContractTrustedQ[batch_Association] := Module[
   {contract = Lookup[batch, "artifactContract", <||>]},
   If[! AssociationQ[contract] || ! TrueQ[Lookup[contract, "sealedQ", False]], Return[False]];
   TrueQ[
    StringQ[Lookup[contract, "sourceDigest", None]] &&
     Lookup[batch, "equationCount", Missing["equationCount"]] === Length[Lookup[batch, "equations", {}]]
    ]
   ];


dsSeedArtifactSourceDigestAuditQ[batch_Association] := Module[
   {contract = Lookup[batch, "artifactContract", <||>], context, rangeAudit},
   context = <|"inputHash" -> Lookup[Lookup[batch, "dSIBPContextSummary", <||>], "inputHash", Missing["inputHash"]]|>;
   rangeAudit = <|"rangeRules" -> Lookup[batch, "targetEnvelopeRules", {}]|>;
   TrueQ[
    dsGeneratedIBPSourceDigest[Lookup[batch, "equations", {}], context, rangeAudit] ===
     Lookup[contract, "sourceDigest", Missing["sourceDigest"]]
    ]
   ];


makeLinearSystemData[batch_Association, topoSpec_: Automatic, OptionsPattern[]] := Module[
   {topo, integrals, integralIndex, equations, linearEquations, coefficientDiagnostics,
    metadataList, metadata, orderingSpec, orderingReport, seedCoverageReport, topologyReport,
    sourceContract, trustedProducerQ, completeSystemQ, sourceDigest,
    auditLevel = OptionValue[AuditLevel], sourceDigestAuditQ},
   topo = normalizeTopologySpec[topoSpec];
   topologyReport = Lookup[batch, "topologyValidationReport", Missing["NoTopologyValidationReport"]];
   If[MatchQ[topologyReport, _Missing] && parsedTopologyQ[topo],
    topologyReport = topologyValidationReport[topo]
    ];
   If[Lookup[batch, "status", "missing"] =!= "generated",
    Return[<|"status" -> "notGenerated", "sourceStatus" -> Lookup[batch, "status", Missing["status"]], "topologyValidationReport" -> topologyReport|>]
    ];
   If[topologyValidationErrorQ[topologyReport],
    Return[<|"status" -> "invalidTopology", "caseName" -> Lookup[batch, "caseName", Missing["caseName"]], "topologyValidationReport" -> topologyReport|>]
    ];
   If[Lookup[batch, "pendingFeatures", {}] =!= {},
    Return[<|
      "status" -> "notReady",
      "caseName" -> Lookup[batch, "caseName", Missing["caseName"]],
      "topologyValidationReport" -> topologyReport,
      "reason" -> "pendingFeatures",
      "pendingFeatures" -> Lookup[batch, "pendingFeatures", {}]
      |>]
    ];
   If[Lookup[batch, "forbiddenNData", {}] =!= {},
    Return[<|
      "status" -> "notReady",
      "caseName" -> Lookup[batch, "caseName", Missing["caseName"]],
      "topologyValidationReport" -> topologyReport,
      "reason" -> "forbiddenNData",
      "forbiddenNData" -> Lookup[batch, "forbiddenNData", {}]
      |>]
    ];
   If[! KeyExistsQ[batch, "completeCanonicalQ"],
    Return[<|
      "status" -> "notReady",
      "caseName" -> Lookup[batch, "caseName", Missing["caseName"]],
      "topologyValidationReport" -> topologyReport,
      "reason" -> "notCanonicalSeedBatch",
      "comment" -> "linear/Kira stages require complete DSGenerateIBP output with all-sector qIBP/tIBP coverage"
      |>]
    ];
   sourceContract = Lookup[batch, "artifactContract", <||>];
   trustedProducerQ = dsSeedArtifactContractTrustedQ[batch];
   sourceDigestAuditQ = If[
     trustedProducerQ && auditLevel === "full",
     dsSeedArtifactSourceDigestAuditQ[batch],
     Missing["NotRunAtStandardAuditLevel"]
     ];
   If[trustedProducerQ && auditLevel === "full" && ! TrueQ[sourceDigestAuditQ],
    Return[<|"status" -> "notReady", "reason" -> "sourceDigestMismatch",
      "sourceDigestAuditQ" -> sourceDigestAuditQ|>]
    ];
   seedCoverageReport = If[
     trustedProducerQ && auditLevel === "standard",
     <|
      "status" -> "producerCertified",
      "passQ" -> True,
      "completeSystemQ" -> TrueQ[Lookup[sourceContract, "completeSystemQ", False]],
      "sourceDigest" -> Lookup[sourceContract, "sourceDigest", Missing["sourceDigest"]],
      "auditLevel" -> Lookup[sourceContract, "auditLevel", "standard"]
      |>,
     makeCanonicalSeedCoverageReport[batch]
     ];
   If[trustedProducerQ && auditLevel === "full" &&
     TrueQ[Lookup[sourceContract, "completeSystemQ", False]] &&
     ! TrueQ[Lookup[seedCoverageReport, "passQ", False]],
    Return[<|"status" -> "notReady", "reason" -> "fullCoverageAuditFailed",
      "seedCoverageReport" -> seedCoverageReport|>]
    ];
   completeSystemQ = TrueQ[
     trustedProducerQ && Lookup[sourceContract, "completeSystemQ", False]
     ];
   orderingReport = validateKiraOrderingSpec[OptionValue[KiraOrdering]];
   If[Lookup[orderingReport, "status", "ok"] =!= "ok",
    Return[<|"status" -> "notReady", "caseName" -> Lookup[batch, "caseName", Missing["caseName"]], "topologyValidationReport" -> topologyReport, "reason" -> "invalidKiraOrdering", "kiraOrderingValidationReport" -> orderingReport|>]
    ];
   metadataList = batchSectorMetadataList[batch, topo];
   orderingSpec = resolveKiraOrderingSpec[batch, topo, OptionValue[KiraOrdering]];
   integrals = sortIntegralsForKira[integralObjectsInBatch[batch], orderingSpec, metadataList];
   integralIndex = makeIntegralIndex[integrals];
   equations = If[
     trustedProducerQ,
     Lookup[batch, "equations", {}],
     symmetry[Lookup[batch, "equations", {}], topo]
     ];
   If[equations === $Failed,
    Return[<|"status" -> "notReady", "caseName" -> Lookup[batch, "caseName", Missing["caseName"]], "topologyValidationReport" -> topologyReport, "reason" -> "invalidSymmetryRules"|>]
    ];
   linearEquations = linearizeSeedEquation[#, integralIndex] & /@ equations;
   linearEquations = DeleteDuplicatesBy[linearEquations, linearEquationMathematicalKey];
   coefficientDiagnostics = linearCoefficientDiagnostics[linearEquations];
   sourceDigest = If[
     trustedProducerQ,
     Lookup[sourceContract, "sourceDigest", Missing["sourceDigest"]],
     IntegerString[Hash[Lookup[batch, "equations", {}], "SHA256"], 16, 64]
     ];
   metadata = If[metadataList === {}, Missing["NoSectorMetadata"], First[metadataList]];
   <|
    "status" -> "generated",
    "caseName" -> Lookup[batch, "caseName", Missing["caseName"]],
    "topology" -> topo,
    "integralCount" -> Length[integrals],
    "sourceEquationCount" -> Length[equations],
    "equationCount" -> Length[linearEquations],
    "duplicateEquationCount" -> Length[equations] - Length[linearEquations],
    "integralList" -> integrals,
    "integralRules" -> Normal[integralIndex],
    "kiraOrdering" -> orderingSpec,
    "kiraOrderingReport" -> kiraOrderingMatchReport[orderingSpec, integrals],
    "integralSortKeys" -> (integralSortKey[#, orderingSpec, metadataList] & /@ integrals),
    "integralMetadata" -> integralMetadataList[integrals, metadataList, orderingSpec],
    "sectorMetadata" -> metadata,
    "sectorMetadataList" -> metadataList,
    "topologyValidationReport" -> topologyReport,
    "seedCoverageReport" -> seedCoverageReport,
    "completeSystemQ" -> completeSystemQ,
    "artifactContract" -> <|
      "contractVersion" -> 1,
      "sealedQ" -> trustedProducerQ,
      "sourceDigest" -> sourceDigest,
      "completeSystemQ" -> completeSystemQ,
      "subsetQ" -> ! completeSystemQ,
       "auditLevel" -> If[trustedProducerQ, auditLevel, "consumerFull"],
       "sourceDigestAuditQ" -> sourceDigestAuditQ,
      "capabilities" -> <|
        "serializeQ" -> True,
        "formalReductionQ" -> completeSystemQ,
        "reuseProducerCanonicalQ" -> trustedProducerQ
        |>
      |>,
    "tadpoleSymmetryData" -> tadpoleSymmetryData[topo],
    "effectiveSymmetryRules" -> effectiveSymmetryRules0[topo],
    "linearEquations" -> linearEquations,
    "linearQ" -> And @@ (Lookup[linearEquations, "linearQ"]),
    "nonlinearEquationCount" -> Count[Lookup[linearEquations, "linearQ"], False],
    "numericCoefficientSystemQ" -> coefficientDiagnostics["numericCoefficientSystemQ"],
    "coefficientVariables" -> coefficientDiagnostics["coefficientVariables"]
    |>
   ];


Options[applyCoefficientRulesToLinearSystem] = {CoefficientRules -> {}};


applyCoefficientRulesToLinearSystem[linearData_Association, OptionsPattern[]] := Module[
   {rawRules = OptionValue[CoefficientRules], rules, ruleReport, linearEquations, coefficientDiagnostics},
   If[Lookup[linearData, "status", "missing"] =!= "generated" || ! KeyExistsQ[linearData, "linearEquations"], Return[linearData]];
   ruleReport = validateCoefficientRules[rawRules];
   If[Lookup[ruleReport, "status", "ok"] =!= "ok",
    Return[Join[linearData, <|"status" -> "notReady", "reason" -> "invalidCoefficientRules", "coefficientRuleValidationReport" -> ruleReport|>]]
    ];
   rules = normalizeCoefficientRulesForLinearData[rawRules, linearData];
   linearEquations = applyKiraCoefficientRulesToLinearEquation[#, rules] & /@ linearData["linearEquations"];
   coefficientDiagnostics = linearCoefficientDiagnostics[linearEquations];
   Join[linearData, <|
     "coefficientRulesApplied" -> rules,
     "userCoefficientRulesApplied" -> userCoefficientRulesForLinearData[rules, linearData],
     "linearEquations" -> linearEquations,
     "linearQ" -> And @@ (Lookup[linearEquations, "linearQ"]),
     "nonlinearEquationCount" -> Count[Lookup[linearEquations, "linearQ"], False],
     "numericCoefficientSystemQ" -> coefficientDiagnostics["numericCoefficientSystemQ"],
     "coefficientVariables" -> coefficientDiagnostics["coefficientVariables"]
     |>]
   ];


Options[makeSampledLinearSystemData] = Join[Options[makeLinearSystemData], {CoefficientRules -> Automatic}];


makeSampledLinearSystemData[batch_Association, topoSpec_: Automatic, OptionsPattern[]] := Module[
   {topo, linearData, rules},
   topo = normalizeTopologySpec[topoSpec];
   linearData = makeLinearSystemData[batch, topo,
     KiraOrdering -> OptionValue[KiraOrdering], AuditLevel -> OptionValue[AuditLevel]];
   If[Lookup[linearData, "status", "missing"] =!= "generated", Return[linearData]];
   rules = Replace[OptionValue[CoefficientRules], Automatic -> {}];
   applyCoefficientRulesToLinearSystem[linearData, CoefficientRules -> rules]
   ];

(* ::Chapter:: *)
(*公开原子 API 与表示转换*)

(* J 本身不携带 topology；公开短签名必须使用显式注册的 context，不能从指标形状猜测物理类型。 *)
If[! ValueQ[$dSIBPTopologyContext], $dSIBPTopologyContext = Missing["NotSet"]];


dSIBPPublicAPI::notopo =
   "公开 API 缺少 topology context。请使用显式 topo 参数，或先调用 setIBPTopologyContext[topo]。";
dSIBPPublicAPI::badtopo = "topology context 无效或解析失败：`1`。";
dSIBPPublicAPI::badshape = "表达式中的 J 与 topology context 不兼容：`1`。";
dSIBPPublicAPI::badstate = "IBP 公开算子要求所有 full-line 离散态已显式取 0/1：`1`。";
dSIBPPublicAPI::badgen = "找不到请求的 IBP 生成元：`1`。";
dSIBPPublicAPI::badvar = "变量 `1` 不在当前 topology 初始化的外部独立变量列表 `2` 中。";
dSIBPPublicAPI::ambiguousvar = "变量 `1` 属于过完备动力学坐标；在 family 定义中补充约束并重选独立变量前，ds 已禁用。";
dSIBPPublicAPI::noinverse = "当前动力学规则没有唯一的用户坐标到基础标量积反向映射；rep2innerform 已拒绝。审计：`1`。";
dSIBPPublicAPI::nonlinear = "ds 只接受 J 的线性组合；检测到非线性或非多项式 J 依赖：`1`。";
dSIBPPublicAPI::derivativefailed = "变量 `1` 的积分导数生成失败。";


setIBPTopologyContext[spec_Association] := Module[{topo = normalizeTopologySpec[spec]},
   If[! parsedTopologyQ[topo],
    Message[dSIBPPublicAPI::badtopo, spec];
    Return[$Failed]
    ];
   $dSIBPTopologyContext = topo
   ];


clearIBPTopologyContext[] := ($dSIBPTopologyContext = Missing["NotSet"]);


resolvePublicTopologyContext[spec_: Automatic] := Module[{topo},
   (* 原子接口既接受 parsed topology，也接受 DSInit 返回的完整 context；后者必须先解包，
      否则会被误当作原始 case 再次 parse，并产生缺少 lineInput 等字段的伪错误。 *)
   topo = Which[
     spec === Automatic, $dSIBPTopologyContext,
     parsedTopologyQ[spec], spec,
     AssociationQ[spec] && TrueQ[dsContextQ[spec]], spec["topology"],
     True, normalizeTopologySpec[spec]
     ];
   If[! parsedTopologyQ[topo],
    If[spec === Automatic,
     Message[dSIBPPublicAPI::notopo],
     Message[dSIBPPublicAPI::badtopo, spec]
     ];
    Return[$Failed]
    ];
   topo
   ];


publicExpectedPackLength[line_Association, packType_String] := Switch[packType,
   "massiveFull" | "massiveCross", If[lineIndexedPowerQ[line], 3, 2],
   "masslessFull", If[lineIndexedPowerQ[line], 2, 1],
   "masslessCross" | "shrunk", If[lineIndexedPowerQ[line], 1, 0],
   _, Missing["UnknownPackType", packType]
   ];


publicIntegralShapeIssues[topo_Association, J[aList_, linePacks_, ispList_]] := Module[
   {issues = {}, expected, packType},
   If[Length[aList] =!= Length[activeAVertexIds[topo]],
    AppendTo[issues, <|"slot" -> "aList", "expected" -> Length[activeAVertexIds[topo]], "actual" -> Length[aList]|>]
    ];
   If[Length[linePacks] =!= topo["nE"],
    AppendTo[issues, <|"slot" -> "linePacks", "expected" -> topo["nE"], "actual" -> Length[linePacks]|>],
    Do[
     If[! ListQ[linePacks[[e]]],
      AppendTo[issues, <|"slot" -> "linePack", "lineIndex" -> e, "reason" -> "notList"|>],
      packType = actualLinePackType[topo, e, linePacks[[e]]];
      expected = publicExpectedPackLength[topo["lines"][[e]], packType];
      If[Head[expected] === Missing || Length[linePacks[[e]]] =!= expected,
       AppendTo[issues, <|"slot" -> "linePack", "lineIndex" -> e, "packType" -> packType, "expected" -> expected, "actual" -> Length[linePacks[[e]]]|>]
       ]
      ],
     {e, Length[linePacks]}
     ]
    ];
   If[Length[ispList] =!= Length[topo["ispData"]],
    AppendTo[issues, <|"slot" -> "ispList", "expected" -> Length[topo["ispData"]], "actual" -> Length[ispList]|>]
    ];
   issues
   ];


publicResolvedDiscreteStateIssues[topo_Association, J[aList_, linePacks_, ispList_]] := Module[
   {issues = {}, packType, pack},
   If[Length[linePacks] =!= topo["nE"], Return[issues]];
   Do[
    pack = linePacks[[e]];
    If[ListQ[pack],
     packType = actualLinePackType[topo, e, pack];
     Switch[packType,
      "massiveFull" | "massiveCross",
      Do[If[! MemberQ[{0, 1}, pack[[p]]], AppendTo[issues, <|"lineIndex" -> e, "packType" -> packType, "packPosition" -> p, "value" -> pack[[p]]|>]],
       {p, linePackNPositions[topo["lines"][[e]], packType]}],
      "masslessFull",
      With[{p = First[linePackNPositions[topo["lines"][[e]], packType]]},
       If[! MemberQ[{0, 1}, pack[[p]]], AppendTo[issues, <|"lineIndex" -> e, "packType" -> packType, "packPosition" -> p, "value" -> pack[[p]]|>]]],
      _, Null
      ]
     ],
    {e, Length[linePacks]}
    ];
   issues
   ];


publicExpressionIntegrals[expr_] := DeleteDuplicates[Cases[expr, _J, {0, Infinity}]];


validatePublicExpression[expr_, topo_Association, requireDiscreteQ_: False] := Module[
   {integrals = publicExpressionIntegrals[expr], shapeIssues, stateIssues},
   shapeIssues = Flatten[publicIntegralShapeIssues[topo, #] & /@ integrals];
   If[shapeIssues =!= {}, Message[dSIBPPublicAPI::badshape, shapeIssues]; Return[False]];
   If[TrueQ[requireDiscreteQ],
    stateIssues = Flatten[publicResolvedDiscreteStateIssues[topo, #] & /@ integrals];
    If[stateIssues =!= {}, Message[dSIBPPublicAPI::badstate, stateIssues]; Return[False]]
    ];
   True
   ];


publicLoopIndex[topo_Association, item_] := Which[
   IntegerQ[item] && 1 <= item <= topo["nL"], item,
   MemberQ[topo["loopMomenta"], item], First@FirstPosition[topo["loopMomenta"], item],
   True, Missing["LoopMomentumNotFound", item]
   ];


publicExternalIndex[topo_Association, item_] := Which[
   IntegerQ[item] && 1 <= item <= topo["nK"], item,
   MemberQ[topo["effectiveLoopExternalMomenta"], item], First@FirstPosition[topo["effectiveLoopExternalMomenta"], item],
   True, Missing["ExternalMomentumNotFound", item]
   ];


publicApplyIBPGenerator[expr_, topo_Association, gen_Association] := Module[{internalExpr, result},
   If[gen["type"] === "time" && ! dsTopologyCapabilityQ[topo, "timeIBPUsableQ"],
    dsErrorPrint["当前 topology 未通过 time-IBP capability gate。 The current topology failed the time-IBP capability gate."]; Return[$Failed]
    ];
   If[gen["type"] === "momentum" && ! dsTopologyCapabilityQ[topo, "momentumIBPUsableQ"],
    dsErrorPrint["当前 topology 未通过 momentum-IBP capability gate。 The current topology failed the momentum-IBP capability gate."]; Return[$Failed]
    ];
    internalExpr = If[
      Lookup[topo, "ibpMode", "full"] === "timeOnly",
      dsTimeOnlyExpressionToInternal020[expr, topo],
      expr
      ];
    If[internalExpr === $Failed, Return[$Failed]];
    If[! validatePublicExpression[internalExpr, topo, True], Return[$Failed]];
    result = Expand[internalExpr /. int_J :> If[
         gen["type"] === "time",
         applyTimeGeneratorSeed[topo, int, gen],
         applyMomentumGeneratorSeed[topo, int, gen]
         ]];
    If[! FreeQ[result, $Failed], Return[$Failed]];
    result = applySeedCanonical[result, topo];
    If[result === $Failed, Return[$Failed]];
    dsTimeOnlyExpressionToPublic020[result, topo]
    ];


dtau[vertex_, expr_, topoSpec_Association] := Module[{topo, gen},
   topo = resolvePublicTopologyContext[topoSpec];
   If[topo === $Failed, Return[$Failed]];
   gen = SelectFirst[makeIBPGenerators[topo], #1["type"] === "time" && #1["vertex"] === vertex &, Missing["NotFound"]];
   If[Head[gen] === Missing, Message[dSIBPPublicAPI::badgen, {"dtau", vertex}]; Return[$Failed]];
   publicApplyIBPGenerator[expr, topo, gen]
   ];
dtau[vertex_, expr_] := Module[{topo = resolvePublicTopologyContext[]},
   If[topo === $Failed, $Failed, dtau[vertex, expr, topo]]
   ];


dqq[dLoop_, vectorLoop_, expr_, topoSpec_Association] := Module[{topo, i, j, gen},
   topo = resolvePublicTopologyContext[topoSpec];
   If[topo === $Failed, Return[$Failed]];
   i = publicLoopIndex[topo, dLoop];
   j = publicLoopIndex[topo, vectorLoop];
   If[Head[i] === Missing || Head[j] === Missing, Message[dSIBPPublicAPI::badgen, {"dqq", dLoop, vectorLoop}]; Return[$Failed]];
   gen = makeAffineMomentumGenerator[topo, i, "loop", j];
   publicApplyIBPGenerator[expr, topo, gen]
   ];
dqq[dLoop_, vectorLoop_, expr_] := Module[{topo = resolvePublicTopologyContext[]},
   If[topo === $Failed, $Failed, dqq[dLoop, vectorLoop, expr, topo]]
   ];


dqk[dLoop_, vectorExternal_, expr_, topoSpec_Association] := Module[{topo, i, j, gen},
   topo = resolvePublicTopologyContext[topoSpec];
   If[topo === $Failed, Return[$Failed]];
   i = publicLoopIndex[topo, dLoop];
   j = publicExternalIndex[topo, vectorExternal];
   If[Head[i] === Missing || Head[j] === Missing, Message[dSIBPPublicAPI::badgen, {"dqk", dLoop, vectorExternal}]; Return[$Failed]];
   gen = makeAffineMomentumGenerator[topo, i, "external", j];
   publicApplyIBPGenerator[expr, topo, gen]
   ];
dqk[dLoop_, vectorExternal_, expr_] := Module[{topo = resolvePublicTopologyContext[]},
   If[topo === $Failed, $Failed, dqk[dLoop, vectorExternal, expr, topo]]
   ];


publicProtectJMap[expr_, function_] := Module[{integrals, tokens, forward, backward},
   integrals = publicExpressionIntegrals[expr];
   tokens = Table[Unique["heldJ$"], {Length[integrals]}];
   forward = Thread[integrals -> tokens];
   backward = Thread[tokens -> integrals];
   function[expr /. forward] /. backward
   ];


userISPToInternalRules[topo_Association] := MapIndexed[
   Rule[Lookup[#1, "name", rho[First[#2]]], rho[First[#2]]] &,
   topo["ispData"]
   ];


internalISPToUserRules[topo_Association] := MapIndexed[
   Rule[rho[First[#2]], Lookup[#1, "name", rho[First[#2]]]] &,
   topo["ispData"]
   ];


rep2innerform[expr_, topoSpec_Association] := Module[{topo, audit},
   topo = resolvePublicTopologyContext[topoSpec];
   If[topo === $Failed, Return[$Failed]];
   audit = Lookup[topo, "kinematicCoordinateAudit", <||>];
   (* 反向映射同时受动量声明层和坐标 Jacobian 层约束；局部坐标可逆不能覆盖
      过完备 loop/无圈动量声明已经关闭的 topology capability。 *)
   If[! dsTopologyCapabilityQ[topo, "inverseKinematicsUsableQ"],
    Message[dSIBPPublicAPI::noinverse, <|
      "capabilities" -> Lookup[topo, "capabilities", <||>],
      "coordinateAudit" -> KeyTake[audit, {"status", "constraintResiduals", "parameterDependencies"}]
      |>];
    Return[$Failed]
    ];
   If[AssociationQ[audit] && ! TrueQ[Lookup[audit, "inverseAvailableQ", True]],
    Message[dSIBPPublicAPI::noinverse, KeyTake[audit, {"status", "constraintResiduals", "parameterDependencies"}]];
    Return[$Failed]
    ];
   If[! validatePublicExpression[expr, topo], Return[$Failed]];
   publicProtectJMap[expr, Function[body,
     Expand[scalarProductInputToInternal[body, topo] /. userISPToInternalRules[topo]]
     ]]
   ];
rep2innerform[expr_] := Module[{topo = resolvePublicTopologyContext[]},
   If[topo === $Failed, $Failed, rep2innerform[expr, topo]]
   ];


rep2outform[expr_, topoSpec_Association] := Module[{topo},
   topo = resolvePublicTopologyContext[topoSpec];
   If[topo === $Failed, Return[$Failed]];
   If[! validatePublicExpression[expr, topo], Return[$Failed]];
   publicProtectJMap[expr, Function[body,
     Expand[scalarProductInternalToUser[body /. internalISPToUserRules[topo], topo]]
     ]]
   ];
rep2outform[expr_] := Module[{topo = resolvePublicTopologyContext[]},
   If[topo === $Failed, $Failed, rep2outform[expr, topo]]
   ];


(* ::Section:: *)
(*外部动力学变量总导数*)

(* ds 的变量必须使用 topology 初始化时的外部名称；内部 kk 坐标只在底层导数中出现。 *)
publicIndependentVariableDerivativeData[topo_Association] := makeIndependentVariableDerivativeGenerators[topo];


resolvePublicIndependentVariableDerivativeData[topo_Association, userVariable_] := SelectFirst[
   publicIndependentVariableDerivativeData[topo],
   SameQ[Lookup[#, "userVariable", Missing["userVariable"]], userVariable] &,
   Missing["UnknownIndependentVariable", userVariable]
   ];


(* 先把每个不同的 J 替成惰性 token；这样 D 只作用于显式系数，不会进入指标。 *)
publicLinearIntegralDecomposition[expr_] := Module[
   {expanded, integrals, tokens, forwardRules, backwardRules, heldExpression,
    constantTerm, coefficients, reconstructed},
   expanded = Expand[expr];
   integrals = publicExpressionIntegrals[expanded];
   tokens = Table[Unique["heldJ$"], {Length[integrals]}];
   forwardRules = Thread[integrals -> tokens];
   backwardRules = Thread[tokens -> integrals];
   heldExpression = Expand[expanded /. forwardRules];
   constantTerm = Expand[heldExpression /. Thread[tokens -> 0]];
   coefficients = Coefficient[heldExpression, #] & /@ tokens;
   reconstructed = Expand[constantTerm + Total[MapThread[#1 #2 &, {coefficients, tokens}]]];
   If[! TrueQ[Expand[heldExpression - reconstructed] === 0],
    Return[<|
      "status" -> "nonlinear",
      "heldExpression" -> heldExpression,
      "integrals" -> integrals
      |>]
    ];
   <|
    "status" -> "linear",
    "expression" -> expanded,
    "heldExpression" -> heldExpression,
    "integrals" -> integrals,
    "tokens" -> tokens,
    "coefficients" -> coefficients,
    "constantTerm" -> constantTerm,
    "backwardRules" -> backwardRules
    |>
   ];


(* 乘积法则：显式系数导数与每个 J 的指标导数分别生成，合并后再统一 canonical。 *)
ds[expr_, userVariable_] := Module[{topo = resolvePublicTopologyContext[]},
   If[topo === $Failed, $Failed, ds[expr, userVariable, topo]]
   ];


integrandVertexFactor[topo_Association, int_J] := Times @@ Table[
   With[{vertex = activeAVertexIds[topo][[slot]]},
    (-tau[vertex])^(int[[1, slot]] + vertexZeroPoint[topo, vertex])
     Exp[vertexExternalPhaseDerivativeCoefficient[topo, vertex] vertexExternalEnergy[topo, vertex] tau[vertex]]
    ],
   {slot, Length[activeAVertexIds[topo]]}
   ];


integrandBuildingBlock[line_Association, pack_List, momentumMagnitude_] := Module[
   {lineId = line["id"], endpoints = line["endpoints"], packType = line["packType"], nPositions},
   nPositions = linePackNPositions[line, packType];
   Switch[packType,
    "massiveFull" | "massiveCross",
    Hh[MassiveBlock[lineFunctionPreset018[line], Lookup[line, "nu", nu], Lookup[line, "skType", "++"], lineId, endpoints, momentumMagnitude, pack[[nPositions[[1]]]], pack[[nPositions[[2]]]]]],
    "masslessFull",
    Hh[MasslessBlock[Lookup[line, "skType", "++"], lineId, endpoints, momentumMagnitude, pack[[First[nPositions]]]]],
    "masslessCross",
    Hh[MasslessCrossBlock[Lookup[line, "skType", "+-"], lineId, endpoints, momentumMagnitude]],
    _, 1
    ]
   ];


integrandLineFactor[topo_Association, int_J, e_Integer] := Module[
   {line = topo["lines"][[e]], pack = int[[2, e]], packType, momentumMagnitude, denominator},
   packType = actualLinePackType[topo, e, pack];
   momentumMagnitude = lineMomentumMagnitude[topo, e];
   denominator = momentumMagnitude^(-linePowerIndex[topo, int, e]);
   If[
    packType === "shrunk",
    denominator,
    denominator integrandBuildingBlock[Join[line, <|"packType" -> packType|>], pack, momentumMagnitude]
    ]
   ];


integrandISPFactor[topo_Association, int_J] := Times @@ Table[
   Lookup[topo["ispData"][[j]], "name", rho[j]]^int[[3, j]],
   {j, Length[topo["ispData"]]}
   ];


rep2Integrand[expr_] := Module[{topo = resolvePublicTopologyContext[]},
   If[topo === $Failed, $Failed, rep2Integrand[expr, topo]]
   ];


(* ::Chapter:: *)
(*Tree 积分表示与 family metadata*)

(* 本章只处理单槽 J[vertexPacks]。tree 与 loop 共用 Head，但任何入口都先按 arity 分派；
   tree family 的 massiveLegs 顺序同时固定 binary master 顺序和矩阵 tensor 顺序。 *)

integralKind[J[_List, _List, _List]] := "Loop";
integralKind[J[packs_List]] := If[VectorQ[packs, ListQ], "Tree", $Failed];
integralKind[_] := $Failed;


makeTreeFamilyData::badinput = "tree family 输入无效：`1`。";
treeIntegralShape::badshape = "tree J 的 pack 形状与 family 不一致：`1`。";


treeSignedExternalLegEnergy[vertex_Association] :=
  (* package 的 + 顶点相位为 Exp[-I E tau]，故论文 e^(I k0 tau) 中 k0=-E。 *)
  If[vertex["vertexType"] === "+", -vertex["externalLegEnergy"], vertex["externalLegEnergy"]];


treeVertexIssues[vertex_Association] := Module[{issues = {}, legs, required, missing},
   required = {"id", "nu0", "vertexType", "externalLegEnergy", "massiveLegs"};
   missing = Complement[required, Keys[vertex]];
   If[missing =!= {}, AppendTo[issues, <|"code" -> "missingVertexKeys", "keys" -> missing|>]];
   If[! MemberQ[{"+", "-"}, Lookup[vertex, "vertexType", Missing["vertexType"]]],
    AppendTo[issues, <|"code" -> "badVertexType", "value" -> Lookup[vertex, "vertexType", Missing["vertexType"]]|>]
    ];
   legs = Lookup[vertex, "massiveLegs", {}];
   If[! ListQ[legs] || ! And @@ (AssociationQ /@ legs),
    AppendTo[issues, <|"code" -> "badMassiveLegs"|>],
    Do[
     missing = Complement[{"id", "nu", "momentumMagnitude"}, Keys[legs[[i]]]];
     If[missing =!= {}, AppendTo[issues, <|"code" -> "missingLegKeys", "leg" -> i, "keys" -> missing|>]],
     {i, Length[legs]}
     ]
    ];
   issues
   ];


makeTreeFamilyData[spec_Association] := Module[{vertices, issues, ids, normalized},
   vertices = Lookup[spec, "vertices", Missing["vertices"]];
   If[! ListQ[vertices] || ! And @@ (AssociationQ /@ vertices),
    Message[makeTreeFamilyData::badinput, {<|"code" -> "badVertices"|>}];
    Return[$Failed]
    ];
   issues = Flatten[treeVertexIssues /@ vertices];
   ids = Lookup[vertices, "id", Missing["id"]];
   If[DuplicateFreeQ[ids] =!= True, AppendTo[issues, <|"code" -> "duplicateVertexIds", "ids" -> ids|>]];
   If[issues =!= {}, Message[makeTreeFamilyData::badinput, issues]; Return[$Failed]];
   normalized = Map[
     Join[#, <|
          "signedExternalLegEnergy" -> treeSignedExternalLegEnergy[#],
         "p" -> Length[#"massiveLegs"]
         |>] &,
     vertices
     ];
   <|
    "status" -> "ready",
    "name" -> Lookup[spec, "name", "treeFamily"],
    "vertices" -> normalized,
    "vertexOrder" -> Lookup[normalized, "id"],
    "packLengths" -> (1 + Lookup[normalized, "p"]),
     "sector" -> Lookup[spec, "sector", ""],
    "sourceStructure" -> Lookup[spec, "sourceStructure", {}],
    "topology" -> Lookup[spec, "topology", Missing["NoLoopTopology"]],
    "requiresSourceRules" -> TrueQ[Lookup[spec, "requiresSourceRules", False]],
    "loopBaselinePowers" -> Lookup[spec, "loopBaselinePowers", <||>]
    |>
   ];


treeIntegralShapeIssues[int_J, data_Association] := Module[{packs, expected, issues = {}},
   If[integralKind[int] =!= "Tree", Return[{<|"code" -> "notTreeIntegral", "integral" -> int|>}]];
   packs = First[int];
   expected = data["packLengths"];
   If[Length[packs] =!= Length[expected],
    Return[{<|"code" -> "vertexCount", "expected" -> Length[expected], "actual" -> Length[packs]|>}]
    ];
   Do[
    If[! ListQ[packs[[i]]] || Length[packs[[i]]] =!= expected[[i]],
     AppendTo[issues, <|"code" -> "packLength", "vertex" -> data["vertexOrder"][[i]], "expected" -> expected[[i]], "actual" -> If[ListQ[packs[[i]]], Length[packs[[i]]], Missing["NotList"]]|>]
     ],
    {i, Length[packs]}
    ];
   issues
   ];


treeIntegralQ[int_J, data_Association] := treeIntegralShapeIssues[int, data] === {};


treeBinaryStates[p_Integer?NonNegative] := IntegerDigits[#, 2, p] & /@ Range[0, 2^p - 1];


treeVertexMasterPacks[vertex_Association, aValue_: 0] := Prepend[#, aValue] & /@ treeBinaryStates[vertex["p"]];


treeMasterList[data_Association, aValues_: Automatic] := Module[{values, perVertex},
   values = If[aValues === Automatic, ConstantArray[0, Length[data["vertices"]]], aValues];
   If[! ListQ[values] || Length[values] =!= Length[data["vertices"]], Return[$Failed]];
   perVertex = MapThread[treeVertexMasterPacks, {data["vertices"], values}];
   J[#] & /@ Tuples[perVertex]
   ];


(* ::Chapter:: *)
(*Tree 通用迭代矩阵*)

(* 这里直接实现 2401.00129 Eq. (3.37)、(3.47)、(3.50)。M1 与 M0tilde 的逆
   只由对角元构造；数值零对角元返回 $Failed，符号分母则作为 singularLoci 显式记录。 *)

treePauli[0] = {{1, 0}, {0, 1}};
treePauli[1] = {{0, 1}, {1, 0}};
treePauli[2] = {{0, -I}, {I, 0}};
treePauli[3] = {{1, 0}, {0, -1}};
treeBasisTransform = 1/Sqrt[2] {{1, -I}, {-I, 1}};
treeBasisTransformInverse = 1/Sqrt[2] {{1, I}, {I, 1}};


treeTensorProduct[mats_List] := Which[
   mats === {}, {{1}},
   Length[mats] === 1, First[mats],
   True, KroneckerProduct @@ mats
   ];


treePauliLift[p_Integer?Positive, slot_Integer, component_Integer] := treeTensorProduct[
   Table[If[i === slot, treePauli[component], treePauli[0]], {i, p}]
   ];


treeM1[vertex_Association, effectiveNu0_] := Module[{p, nus, identity},
   p = vertex["p"];
   nus = Lookup[vertex["massiveLegs"], "nu"];
   identity = IdentityMatrix[2^p];
   If[p === 0,
    effectiveNu0 identity,
    Sum[(nus[[i]] + 1/2) treePauliLift[p, i, 3], {i, p}] +
     (effectiveNu0 - p/2 - Total[nus]) identity
    ]
   ];


treeM0[vertex_Association] := Module[{p, momentumMagnitudes, identity},
   p = vertex["p"];
   momentumMagnitudes = Lookup[vertex["massiveLegs"], "momentumMagnitude"];
   identity = IdentityMatrix[2^p];
   If[p === 0,
    I vertex["signedExternalLegEnergy"] identity,
    -I Sum[momentumMagnitudes[[i]] treePauliLift[p, i, 2], {i, p}] + I vertex["signedExternalLegEnergy"] identity
    ]
   ];


treeM0Tilde[vertex_Association] := Module[{p, momentumMagnitudes, identity},
   p = vertex["p"];
   momentumMagnitudes = Lookup[vertex["massiveLegs"], "momentumMagnitude"];
   identity = IdentityMatrix[2^p];
   If[p === 0,
    I vertex["signedExternalLegEnergy"] identity,
    -I Sum[momentumMagnitudes[[i]] treePauliLift[p, i, 3], {i, p}] + I vertex["signedExternalLegEnergy"] identity
    ]
   ];


treeTp[vertex_Association] := treeTensorProduct[ConstantArray[treeBasisTransform, vertex["p"]]];
treeTpInverse[vertex_Association] := treeTensorProduct[ConstantArray[treeBasisTransformInverse, vertex["p"]]];


treeDiagonalInverse::singular = "tree 递推需要除以零，当前参数位于奇异层：`1`。请避开这些参数值，或先加入 regulator 再进行递推。";


treeDiagonalInverse[matrix_List] := Module[{diag, offDiagonal},
   diag = Diagonal[matrix];
   offDiagonal = matrix - DiagonalMatrix[diag];
   If[! TrueQ[offDiagonal === ConstantArray[0, Dimensions[matrix]]], Return[$Failed]];
   If[AnyTrue[diag, TrueQ[# === 0] &], Message[treeDiagonalInverse::singular, Select[diag, TrueQ[# === 0] &]]; Return[$Failed]];
   DiagonalMatrix[1/diag]
   ];


treeAminusMatrix[vertex_Association, shift_: 0] := Module[{m1, inv},
   m1 = treeM1[vertex, vertex["nu0"] + shift];
   inv = treeDiagonalInverse[m1];
   If[inv === $Failed, $Failed, -inv . treeM0[vertex]]
   ];


treeAplusMatrix[vertex_Association, shift_: 0] := Module[{m0t, inv},
   m0t = treeM0Tilde[vertex];
   inv = treeDiagonalInverse[m0t];
   If[inv === $Failed, $Failed,
    -treeTpInverse[vertex] . inv . treeTp[vertex] . treeM1[vertex, vertex["nu0"] + shift + 1]
    ]
   ];


makeTreeRecurrenceData[data_Association] := Module[{vertices, records},
   vertices = data["vertices"];
   records = Map[
     Function[vertex,
      <|
       "vertex" -> vertex["id"],
       "states" -> treeBinaryStates[vertex["p"]],
       "M1" -> treeM1[vertex, vertex["nu0"]],
       "M0" -> treeM0[vertex],
       "M0Tilde" -> treeM0Tilde[vertex],
       "Aminus" -> treeAminusMatrix[vertex, 0],
       "Aplus" -> treeAplusMatrix[vertex, 0],
       "singularLoci" -> Join[
         Thread[Diagonal[treeM1[vertex, vertex["nu0"]]] != 0],
         Thread[Diagonal[treeM0Tilde[vertex]] != 0]
         ]
       |>
      ],
     vertices
     ];
   <|"status" -> If[FreeQ[records, $Failed], "ready", "singular"], "vertices" -> records|>
   ];


treeStateIndex[bits_List] := FromDigits[bits, 2] + 1;


treeLoopIntegralFromTree::unsupported = "tree 到 loop seed 的反投影尚不支持该 line pack：`1`。";


treeLoopIntegralFromTree[int_J, data_Association] := Module[
   {topo = Lookup[data, "topology", Missing["NoLoopTopology"]], packs, aList, linePacks, baseline, line, states, vertexIndex, legIndex},
   If[Head[topo] === Missing || ! treeIntegralQ[int, data], Return[$Failed]];
   packs = First[int];
   aList = packs[[All, 1]];
   baseline = Lookup[data, "loopBaselinePowers", <||>];
   linePacks = Table[
     line = topo["lines"][[e]];
     Switch[line["packType"],
      "massiveFull" | "massiveCross",
      states = Table[
        vertexIndex = FirstPosition[data["vertexOrder"], line["endpoints"][[slot]], Missing["NoVertex"]];
        If[Head[vertexIndex] === Missing, Return[$Failed]];
        vertexIndex = First[vertexIndex];
        legIndex = FirstPosition[Lookup[data["vertices"][[vertexIndex, "massiveLegs"]], "id"], {line["id"], slot}, Missing["NoLeg"]];
        If[Head[legIndex] === Missing, Return[$Failed]];
        packs[[vertexIndex, 1 + First[legIndex]]],
        {slot, 2}
        ];
      If[lineIndexedPowerQ[line], Prepend[states, Lookup[baseline, line["id"], 0]], states],
      "shrunk",
      If[lineIndexedPowerQ[line], {Lookup[baseline, line["id"], 0]}, {}],
      _,
      Message[treeLoopIntegralFromTree::unsupported, <|"line" -> line["id"], "packType" -> line["packType"]|>];
      Return[$Failed]
      ],
     {e, topo["nE"]}
     ];
   J[aList, linePacks, ConstantArray[0, Length[topo["ispData"]]]]
   ];


treeSourceAwareStepFromTopology[int_J, vertexIndex_Integer, endpoint_Integer, data_Association] := Module[
   {packs = First[int], current, seedA, states, vertexId, localIntegrals, loopIntegrals, records, ruleData, rules, result},
   current = packs[[vertexIndex, 1]];
   If[current === endpoint, Return[int]];
   seedA = If[current < endpoint, current + 1, current];
   states = treeBinaryStates[data["vertices"][[vertexIndex, "p"]]];
   vertexId = data["vertexOrder"][[vertexIndex]];
   localIntegrals = J[ReplacePart[packs, vertexIndex -> Prepend[#, seedA]]] & /@ states;
   loopIntegrals = treeLoopIntegralFromTree[#, data] & /@ localIntegrals;
   If[! FreeQ[loopIntegrals, $Failed], Return[$Failed]];
   records = Map[DSTreeSeeds[vertexId, #, data["topology"]] &, loopIntegrals];
   If[! FreeQ[records, $Failed], Return[$Failed]];
   ruleData = makeTreeTimeReductionRules[records, data];
   rules = Lookup[ruleData, If[current < endpoint, "minus", "plus"], {}];
   result = Replace[int, rules];
   If[result === int,
    Message[makeTreeFamilyData::badinput, {<|"code" -> "missingSourceAwareStep", "integral" -> int, "vertex" -> vertexId|>}];
    $Failed,
    result
    ]
   ];


treeSingleStepIntegral[int_J, vertexIndex_Integer, endpoint_Integer, data_Association] := Module[
   {packs, pack, current, bits, vertex, matrix, nextA, states, row, newPacks, directRules, directResult},
   If[! treeIntegralQ[int, data], Return[$Failed]];
   packs = First[int];
   pack = packs[[vertexIndex]];
   current = First[pack];
   bits = Rest[pack];
   vertex = data["vertices"][[vertexIndex]];
   If[TrueQ[Lookup[data, "requiresSourceRules", False]] && AssociationQ[Lookup[data, "topology", Missing["NoLoopTopology"]]],
    Return[treeSourceAwareStepFromTopology[int, vertexIndex, endpoint, data]]
    ];
   directRules = Lookup[Lookup[data, "timeReductionRules", <||>], If[current < endpoint, "minus", "plus"], {}];
   directResult = Replace[int, directRules];
   If[directResult =!= int, Return[directResult]];
   If[TrueQ[Lookup[data, "requiresSourceRules", False]],
    Message[makeTreeFamilyData::badinput, {<|"code" -> "missingSourceAwareStep", "integral" -> int, "vertex" -> vertex["id"]|>}];
    Return[$Failed]
    ];
   Which[
    current < endpoint,
    matrix = treeAminusMatrix[vertex, current + 1];
    nextA = current + 1,
    current > endpoint,
    matrix = treeAplusMatrix[vertex, current - 1];
    nextA = current - 1,
    True,
    Return[int]
    ];
   If[matrix === $Failed, Return[$Failed]];
   states = treeBinaryStates[vertex["p"]];
   row = treeStateIndex[bits];
   Total[MapIndexed[
     Function[{state, index},
      newPacks = ReplacePart[packs, vertexIndex -> Prepend[state, nextA]];
      matrix[[row, First[index]]] J[newPacks]
      ],
     states
     ]]
   ];


treeEndpointData::badend = "tree 迭代终点无效：`1`。";


normalizeTreeEndpoints[end_, data_Association] := Module[{result},
   result = If[end === Automatic, ConstantArray[0, Length[data["vertices"]]], end];
   If[! ListQ[result] || Length[result] =!= Length[data["vertices"]] || ! And @@ (IntegerQ /@ result),
    Message[treeEndpointData::badend, end];
    Return[$Failed]
    ];
   result
   ];


treeReducibleIntegrals[expr_, endpoints_List, data_Association] := Select[
   DeleteDuplicates[Cases[expr, int_J /; treeIntegralQ[int, data], {0, Infinity}]],
   Or @@ MapThread[Unequal, {First[#][[All, 1]], endpoints}] &
   ];


treeEndpointDistance[int_J, endpoints_List] := Total[Abs[First[int][[All, 1]] - endpoints]];


treeProgressKeyLessQ[target : {_, _}, source : {_, _}] :=
  target[[1]] < source[[1]] || (target[[1]] === source[[1]] && target[[2]] < source[[2]]);


treeRecurrenceStateHash[expr_] := IntegerString[Hash[HoldComplete[expr], "SHA256"], 16, 64];


(* 每条递推边都必须严格降低到指定 endpoint 的整数 L1 距离。该局部序保证递推 DAG
   不依赖任意迭代次数；完整表达式哈希另用于诊断实现错误造成的状态循环。 *)
treeSingleFamilyStepProgress[source_J, replacement_, endpoints_List, data_Association] := Module[
   {sourceDistance, targets, invalidTargets, targetDistances},
   sourceDistance = treeEndpointDistance[source, endpoints];
   targets = DeleteDuplicates[Cases[replacement, int_J, {0, Infinity}]];
   invalidTargets = Select[
     targets,
     ! treeIntegralQ[#, data] || ! And @@ (IntegerQ /@ First[#][[All, 1]]) &
     ];
   targetDistances = treeEndpointDistance[#, endpoints] & /@ Complement[targets, invalidTargets];
   <|
    "passQ" -> TrueQ[invalidTargets === {} && And @@ (# < sourceDistance & /@ targetDistances)],
    "source" -> source,
    "sourceDistance" -> sourceDistance,
    "targetIntegrals" -> targets,
    "targetDistances" -> targetDistances,
    "invalidTargets" -> invalidTargets
    |>
   ];


repIterativeData::badindex = "tree a 指标必须是可判定整数：`1`。";
repIterativeData::nosector = "tree 积分无法唯一匹配 sector family：`1`。";
repIterativeData::noprogress = "tree 递推没有严格趋近指定终点：`1`。 Tree recurrence did not strictly approach the requested endpoint: `1`.";
repIterativeData::cycle = "tree 递推检测到重复 canonical 状态：`1`。 Tree recurrence encountered a repeated canonical state: `1`.";


Options[repIterativeData] = {};


repIterativeData[expr_, end_: Automatic, data_Association, OptionsPattern[]] := Module[
   {endpoints, result = Expand[expr], integrals, allA, steps = 0, firstInt, packs, vertexIndex,
    replacement, progressData, seenStates = <||>, stateHash},
   endpoints = normalizeTreeEndpoints[end, data];
   If[endpoints === $Failed, Return[<|"status" -> "error", "result" -> $Failed, "steps" -> 0|>]];
   integrals = Select[DeleteDuplicates[Cases[result, int_J /; treeIntegralQ[int, data], {0, Infinity}]], True &];
   allA = Flatten[First[#][[All, 1]] & /@ integrals];
   If[! And @@ (IntegerQ /@ allA), Message[repIterativeData::badindex, Select[allA, ! IntegerQ[#] &]]; Return[<|"status" -> "error", "result" -> $Failed, "steps" -> 0|>]];
   AssociateTo[seenStates, treeRecurrenceStateHash[result] -> True];
   While[(integrals = treeReducibleIntegrals[result, endpoints, data]) =!= {},
    firstInt = First[integrals];
    packs = First[firstInt];
    vertexIndex = SelectFirst[Range[Length[endpoints]], packs[[#, 1]] =!= endpoints[[#]] &];
    replacement = treeSingleStepIntegral[firstInt, vertexIndex, endpoints[[vertexIndex]], data];
    If[replacement === $Failed, Return[<|"status" -> "singular", "result" -> $Failed, "steps" -> steps|>]];
    progressData = treeSingleFamilyStepProgress[firstInt, replacement, endpoints, data];
    If[! TrueQ[progressData["passQ"]],
     Message[repIterativeData::noprogress, progressData];
     Return[<|"status" -> "error", "reason" -> "nonDecreasingRecurrence", "result" -> $Failed,
       "steps" -> steps, "progressData" -> progressData|>]
     ];
    result = Expand[result /. firstInt -> replacement];
    stateHash = treeRecurrenceStateHash[result];
    If[KeyExistsQ[seenStates, stateHash],
     Message[repIterativeData::cycle, stateHash];
     Return[<|"status" -> "error", "reason" -> "recurrenceCycle", "result" -> $Failed,
       "steps" -> steps, "stateHash" -> stateHash|>]
     ];
    AssociateTo[seenStates, stateHash -> True];
    steps++;
    ];
   <|"status" -> "reduced", "result" -> result, "steps" -> steps, "endpoints" -> endpoints|>
   ];


$dSIBPTreeFamilyContext = None;
repIterative0 = {};


treeFirstStepToZero[int_J, data_Association] := Module[{packs, vertexIndex},
   If[! treeIntegralQ[int, data], Return[int]];
   packs = First[int];
   If[! And @@ (IntegerQ /@ packs[[All, 1]]), Return[$Failed]];
   vertexIndex = SelectFirst[Range[Length[packs]], packs[[#, 1]] =!= 0 &, Missing["AtEndpoint"]];
   If[Head[vertexIndex] === Missing, int, treeSingleStepIntegral[int, vertexIndex, 0, data]]
   ];


makeRepIterative0[data_Association] /; KeyExistsQ[data, "vertices"] && ! KeyExistsQ[data, "families"] := With[{family = data}, {
    HoldPattern[int_J] /; treeIntegralQ[int, family] && And @@ (IntegerQ /@ First[int][[All, 1]]) && AnyTrue[First[int][[All, 1]], # != 0 &] :>
     treeFirstStepToZero[int, family]
    }];


setTreeFamilyContext[data_Association] := Module[{},
   $dSIBPTreeFamilyContext = data;
   repIterative0 = makeRepIterative0[data];
   data
   ];


Options[repIterative] = Options[repIterativeData];


repIterative[expr_, end_: Automatic, OptionsPattern[]] := If[
   AssociationQ[$dSIBPTreeFamilyContext],
    If[KeyExistsQ[$dSIBPTreeFamilyContext, "families"],
     repIterativeSectorData[expr, end, $dSIBPTreeFamilyContext]["result"],
     repIterativeData[expr, end, $dSIBPTreeFamilyContext]["result"]
    ],
   Message[makeTreeFamilyData::badinput, {<|"code" -> "treeContextNotSet"|>}];
   $Failed
   ];


repIterative[expr_, end_, data_Association, OptionsPattern[]] /;
   KeyExistsQ[data, "vertices"] && ! KeyExistsQ[data, "families"] := Module[{activeData},
   (* 显式 family 调用同时刷新公开原始规则，保证随后可直接使用 repIterative0。 *)
   activeData = setTreeFamilyContext[data];
   repIterativeData[expr, end, activeData]["result"]
   ];


(* ::Chapter:: *)
(*Loop time-IBP 到 tree 的显式投影*)

(* 投影器从每个 loop J 的 shrunk pack 重建目标 sector，再按原始 line/endpoint 顺序打包 n。
   所有 theta、WT、zero-point 与 coincident canonical 已在 dtau 中完成；这里不重写传播子边界公式。 *)

loopToTreeProjection::badloop = "loop-to-tree 投影只接受合法三槽 loop J：`1`。";
loopToTreeProjection::mixedcontact = "mixed-sign line 不得产生 theta/contact shrink：`1`。";


loopIntegralShrunkLines[int : J[_, linePacks_List, _], topo_Association] := Select[
   Range[Min[Length[linePacks], topo["nE"]]],
   actualLinePackType[topo, #, linePacks[[#]]] === "shrunk" &&
     ! MemberQ[{"masslessCross", "massiveCross"}, topo["lines"][[#, "packType"]]] &
   ];


loopTreeTargetTopology[int : J[_, linePacks_List, _], topo_Association] := Module[{shrunk, badMixed},
   shrunk = loopIntegralShrunkLines[int, topo];
   badMixed = Select[shrunk, MemberQ[{"+-", "-+"}, Lookup[topo["lines"][[#]], "skType", "++"]] &];
   If[badMixed =!= {}, Message[loopToTreeProjection::mixedcontact, badMixed]; Return[$Failed]];
   If[shrunk === {}, topo, shrinkSectorTopology[topo, shrunk]]
   ];


(* loop 的 a0 留在 tree family 的 nu0 中；b0/bS0 没有 tree 指标槽，必须随完整物理幂次进入显式能量系数。 *)
loopTreeProjectionCoefficient[int : J[_, _, _], targetTopo_Association] := Module[
   {momentumMagnitudePower},
   momentumMagnitudePower = Times @@ Table[
      lineMomentumMagnitude[targetTopo, e]^(-linePowerIndex[targetTopo, int, e]),
      {e, targetTopo["nE"]}
      ];
   momentumMagnitudePower
   ];


loopTreeRelativeProjectionCoefficient[
   int : J[targetA_List, _, _], targetTopo_Association,
   referenceInt : J[referenceA_List, _, _], referenceTopo_Association
   ] := Module[{targetLinePower, referenceLinePower, momentumMagnitude},
   Times @@ Table[
      targetLinePower = linePowerIndex[targetTopo, int, e];
      referenceLinePower = linePowerIndex[referenceTopo, referenceInt, e];
      momentumMagnitude = lineMomentumMagnitude[targetTopo, e];
      momentumMagnitude^(-Expand[targetLinePower - referenceLinePower]),
      {e, targetTopo["nE"]}
      ]
   ];


projectLoopIntegralToTree[int : J[aList_List, linePacks_List, ispList_List], topo_Association, referenceInt : J[_, _, _]] := Module[
   {targetTopo, referenceTopo, active, repMap, packs, legStates, originalEndpoints, line, lineId,
    relativeCoefficient},
   If[Length[linePacks] =!= topo["nE"] || ispList =!= ConstantArray[0, Length[ispList]],
    Message[loopToTreeProjection::badloop, int]; Return[$Failed]
    ];
   targetTopo = loopTreeTargetTopology[int, topo];
   referenceTopo = loopTreeTargetTopology[referenceInt, topo];
   If[MemberQ[{targetTopo, referenceTopo}, $Failed], Return[$Failed]];
   active = activeAVertexIds[targetTopo];
   repMap = Lookup[targetTopo, "sectorVertexRepresentativeMap", AssociationThread[topo["vertexIds"] -> topo["vertexIds"]]];
   If[Length[aList] =!= Length[active], Message[loopToTreeProjection::badloop, int]; Return[$Failed]];
   packs = Table[
     legStates = {};
     Do[
      If[! MemberQ[loopIntegralShrunkLines[int, topo], e] && MemberQ[{"massiveFull", "massiveCross"}, topo["lines"][[e, "packType"]]],
       line = topo["lines"][[e]];
       lineId = line["id"];
       originalEndpoints = Lookup[line, "originalEndpoints", line["endpoints"]];
       Do[
        If[Lookup[repMap, originalEndpoints[[slot]], originalEndpoints[[slot]]] === active[[v]],
          AppendTo[legStates, linePacks[[e, linePackNPositions[line, line["packType"]][[slot]]]]]
         ],
        {slot, 2}
        ]
       ],
      {e, topo["nE"]}
      ];
     Prepend[legStates, aList[[v]]],
     {v, Length[active]}
     ];
   relativeCoefficient = loopTreeRelativeProjectionCoefficient[int, targetTopo, referenceInt, referenceTopo];
   relativeCoefficient J[packs]
   ];


projectLoopIntegralToTree[int : J[_, _, _], topo_Association] := Module[{targetTopo = loopTreeTargetTopology[int, topo]},
   If[targetTopo === $Failed, $Failed, loopTreeProjectionCoefficient[int, targetTopo] projectLoopIntegralToTree[int, topo, int]]
   ];
projectLoopIntegralToTree[int_J, _Association, _J] := (Message[loopToTreeProjection::badloop, int]; $Failed);
projectLoopIntegralToTree[int_J, _Association] := (Message[loopToTreeProjection::badloop, int]; $Failed);


projectLoopTimeEquationToTree[expr_, topoSpec_Association, referenceInt_: Automatic] := Module[{topo, result, ref},
   topo = normalizeTopologySpec[topoSpec];
   If[topo === $Failed, Return[$Failed]];
   ref = If[referenceInt === Automatic, FirstCase[expr, int : J[_, _, _] :> int, Missing["NoReferenceIntegral"], Infinity], referenceInt];
   If[Head[ref] === Missing || integralKind[ref] =!= "Loop", Message[loopToTreeProjection::badloop, ref]; Return[$Failed]];
   result = Expand[expr /. int : J[_, _, _] :> projectLoopIntegralToTree[int, topo, ref]];
   If[FreeQ[result, $Failed], result, $Failed]
   ];


treeContactSourceIntegrals[expr_] := Module[{integrals, maxVertices},
   integrals = DeleteDuplicates[Cases[expr, int : J[_List] :> int, Infinity]];
   If[integrals === {}, Return[{}]];
   maxVertices = Max[Length[First[#]] & /@ integrals];
   Select[integrals, Length[First[#]] < maxVertices &]
   ];


DSTreeSeeds[vertex_, int : J[_, _, _], topoSpec_Association] := Module[
   {topo, loopSeed, treeSeed, contactAudit, shrinkConsumptionTrace, allLoopIntegrals, shrunkByTerm, projectedReference},
   topo = normalizeTopologySpec[topoSpec];
   If[topo === $Failed, Return[$Failed]];
   loopSeed = dtau[vertex, int, topo];
   If[loopSeed === $Failed, Return[$Failed]];
   allLoopIntegrals = DeleteDuplicates[Cases[loopSeed, J[_, _, _], Infinity]];
   shrunkByTerm = Association@Table[
      ToString[item, InputForm] -> loopIntegralShrunkLines[item, topo],
      {item, allLoopIntegrals}
      ];
   contactAudit = Map[
     Function[e,
      <|"lineIndex" -> e, "lineId" -> topo["lines"][[e, "id"]], "skType" -> Lookup[topo["lines"][[e]], "skType", Missing["skType"]],
       "thetaAllowedQ" -> MemberQ[{"++", "--"}, Lookup[topo["lines"][[e]], "skType", "++"]]|>
      ],
     DeleteDuplicates[Flatten[Values[shrunkByTerm]]]
     ];
   shrinkConsumptionTrace = Map[
     Function[e,
      <|"lineIndex" -> e, "lineId" -> topo["lines"][[e, "id"]],
       "skType" -> Lookup[topo["lines"][[e]], "skType", Missing["skType"]],
       "consumer" -> If[Lookup[topo["lines"][[e]], "massType", "massive"] === "massive", "WT/shrinkTerms", "masslessContact"]|>
      ],
     DeleteDuplicates[Flatten[Values[shrunkByTerm]]]
     ];
   If[AnyTrue[contactAudit, ! TrueQ[#"thetaAllowedQ"] &], Message[loopToTreeProjection::mixedcontact, contactAudit]; Return[$Failed]];
   treeSeed = projectLoopTimeEquationToTree[loopSeed, topo, int];
   projectedReference = projectLoopIntegralToTree[int, topo, int];
   <|
    "status" -> If[treeSeed === $Failed, "error", "generated"],
    "generator" -> dtau[vertex],
    "loopSeed" -> loopSeed,
    "treeSeed" -> treeSeed,
    "treeIntegral" -> projectedReference,
    "referenceA" -> First[int],
    "contactAudit" -> contactAudit,
    "shrinkConsumptionTrace" -> shrinkConsumptionTrace,
    "mixedContactQ" -> AnyTrue[contactAudit, MemberQ[{"+-", "-+"}, #"skType"] &],
    "shrunkLinesByIntegral" -> shrunkByTerm
    |>
   ];


DSTreeSeeds[int : J[_, _, _], topoSpec_Association] := Module[{topo, vertices},
   topo = normalizeTopologySpec[topoSpec];
   If[topo === $Failed, Return[$Failed]];
   vertices = Lookup[Select[makeIBPGenerators[topo], #"type" === "time" &], "vertex"];
   DSTreeSeeds[#, int, topo] & /@ vertices
   ];


treeSeedRuleGroupKey[record_Association, data_Association] := Module[
   {treeInt, treeInts, vertexId, vertexIndex, packs},
   treeInts = Cases[record["treeIntegral"], _J, {0, Infinity}];
   treeInt = If[treeInts === {}, Missing["NoTreeIntegral"], First[treeInts]];
   vertexId = record["generator"] /. dtau[v_] :> v;
   vertexIndex = FirstPosition[data["vertexOrder"], vertexId, Missing["NoVertex"]];
   If[Head[treeInt] === Missing || Head[vertexIndex] === Missing, Return[Missing["BadSeedRecord"]]];
   vertexIndex = First[vertexIndex];
   packs = First[treeInt];
   {vertexId, ReplacePart[packs, vertexIndex -> {packs[[vertexIndex, 1]], "state"}]}
   ];


makeTreeTimeReductionRules::incomplete = "tree time seed 状态组不完整：`1`。";


makeTreeTimeReductionRules[records_List, data_Association] := Module[
   {good, grouped, minusRules = {}, plusRules = {}, sourceQ = False, groupData, vertexId, vertexIndex,
    vertex, states, ordered, currentInts, minusInts, equations, m1, m0, invM1, invM0t, tp, tpInv,
    regular, source, lowerRhs, upperRhs},
   good = Select[Flatten[records], AssociationQ[#] && Lookup[#, "status", "error"] === "generated" &];
   grouped = GroupBy[good, treeSeedRuleGroupKey[#, data] &];
   KeyValueMap[
    Function[{key, group},
     If[Head[key] === Missing, Return[]];
     vertexId = key[[1]];
     vertexIndex = First@FirstPosition[data["vertexOrder"], vertexId];
     vertex = data["vertices"][[vertexIndex]];
     states = treeBinaryStates[vertex["p"]];
     ordered = SortBy[group, Function[record,
        treeStateIndex[Rest[First[First[Cases[record["treeIntegral"], _J, {0, Infinity}]]][[vertexIndex]]]]
        ]];
     If[Length[ordered] =!= Length[states], Message[makeTreeTimeReductionRules::incomplete, <|"key" -> key, "expected" -> Length[states], "actual" -> Length[ordered]|>]; Return[]];
     currentInts = First[Cases[#"treeIntegral", _J, {0, Infinity}]] & /@ ordered;
     minusInts = MapThread[
       Function[{int, state}, J[ReplacePart[First[int], vertexIndex -> Prepend[state, First[int][[vertexIndex, 1]] - 1]]]],
       {currentInts, states}
       ];
     equations = Lookup[ordered, "treeSeed"];
     m1 = treeM1[vertex, vertex["nu0"] + First[currentInts[[1]]][[vertexIndex, 1]]];
     m0 = treeM0[vertex];
     invM1 = treeDiagonalInverse[m1];
     invM0t = treeDiagonalInverse[treeM0Tilde[vertex]];
     If[MemberQ[{invM1, invM0t}, $Failed], Return[]];
     tp = treeTp[vertex];
     tpInv = treeTpInverse[vertex];
     regular = m1.minusInts + m0.currentInts;
     source = Expand[equations - regular];
     If[! TrueQ[source === ConstantArray[0, Length[source]]], sourceQ = True];
     lowerRhs = Expand[treeAminusMatrix[vertex, First[currentInts[[1]]][[vertexIndex, 1]]] . currentInts - invM1.source];
     upperRhs = Expand[treeAplusMatrix[vertex, First[currentInts[[1]]][[vertexIndex, 1]] - 1] . minusInts - tpInv.invM0t.tp.source];
     minusRules = Join[minusRules, Thread[minusInts -> lowerRhs]];
     plusRules = Join[plusRules, Thread[currentInts -> upperRhs]];
     ],
    grouped
    ];
   <|
    "status" -> If[Length[grouped] === 0, "empty", "generated"],
    "minus" -> DeleteDuplicates[minusRules],
    "plus" -> DeleteDuplicates[plusRules],
    "sourceQ" -> sourceQ,
    "groupCount" -> Length[grouped]
    |>
   ];


attachTreeTimeReductionRules[data_Association, ruleData_Association] /;
   KeyExistsQ[data, "vertices"] && ! KeyExistsQ[data, "families"] := Join[
   data,
   <|
    "timeReductionRules" -> KeyTake[ruleData, {"minus", "plus"}],
    "requiresSourceRules" -> TrueQ[ruleData["sourceQ"]]
    |>
   ];


makeTreeFamilyDataFromTopology[topoSpec_Association] := Module[
   {topo, active, repMap, vertices, originalClass, signs, legs},
   topo = normalizeTopologySpec[topoSpec];
   If[topo === $Failed, Return[$Failed]];
   active = activeAVertexIds[topo];
   repMap = Lookup[topo, "sectorVertexRepresentativeMap", AssociationThread[topo["vertexIds"] -> topo["vertexIds"]]];
   vertices = Table[
     originalClass = Select[topo["vertexIds"], Lookup[repMap, #, #] === active[[v]] &];
     signs = DeleteDuplicates[Lookup[topo["vertexSignAssoc"], originalClass]];
     If[Length[signs] =!= 1,
      Message[makeTreeFamilyData::badinput, {<|"code" -> "mixedSignMergedVertex", "vertices" -> originalClass|>}];
      Return[$Failed]
      ];
     legs = Flatten@Table[
        If[MemberQ[{"massiveFull", "massiveCross"}, topo["lines"][[e, "packType"]]],
         Table[
          If[topo["lines"][[e, "endpoints", slot]] === active[[v]],
           {<|
             "id" -> {topo["lines"][[e, "id"]], slot},
             "nu" -> topo["lines"][[e, "nu"]],
              "momentumMagnitude" -> lineMomentumMagnitude[topo, e]
             |>},
           {}
           ],
          {slot, 2}
          ],
         {}
         ],
        {e, topo["nE"]}
        ];
     <|
      "id" -> active[[v]],
       "vertexType" -> First[signs],
      "nu0" -> vertexZeroPoint[topo, active[[v]]],
       "externalLegEnergy" -> vertexExternalEnergy[topo, active[[v]]],
      "massiveLegs" -> legs
      |>,
     {v, Length[active]}
     ];
   makeTreeFamilyData[<|
     "name" -> topo["name"],
     "sector" -> sectorKeyFromShrunkLines[topo, Lookup[topo, "sectorShrunkLines", {}]],
     "vertices" -> vertices,
     "topology" -> topo,
     "requiresSourceRules" -> thetaFullLineIndices[topo] =!= {}
     |>]
   ];


makeTreeSectorFamilies[topoSpec_Association] := Module[{topo, subsets, sectorTopos, families},
   topo = normalizeTopologySpec[topoSpec];
   If[topo === $Failed, Return[$Failed]];
   subsets = contactReachableShrinkSubsets[topo];
   sectorTopos = Join[{topo}, shrinkSectorTopology[topo, #] & /@ subsets];
   families = makeTreeFamilyDataFromTopology /@ sectorTopos;
   If[MemberQ[families, $Failed], Return[$Failed]];
   <|
    "status" -> "ready",
    "families" -> families,
    "sectorOrder" -> Lookup[families, "sector"],
    "topFamily" -> First[families]
    |>
   ];


attachTreeTimeReductionRules[context_Association, ruleData_Association] /; KeyExistsQ[context, "families"] := Module[
   {families = context["families"]},
   families[[1]] = attachTreeTimeReductionRules[First[families], ruleData];
   Join[context, <|"families" -> families, "topFamily" -> First[families]|>]
   ];


treeFamilyForIntegral::ambiguous = "tree J 的 pack 形状同时匹配多个 sector：`1`。当前表示无法唯一确定 sector，已拒绝继续约化。";


treeFamilyForIntegral[int_J, context_Association] := Module[{matches},
   matches = Select[context["families"], treeIntegralQ[int, #] &];
   Which[
    matches === {}, Missing["NoTreeSectorFamily", int],
    Length[matches] === 1, First[matches],
    True,
    Message[treeFamilyForIntegral::ambiguous, <|"integral" -> int, "sectors" -> Lookup[matches, "sector"]|>];
    Missing["AmbiguousTreeSectorFamily", int, Lookup[matches, "sector"]]
    ]
   ];


makeRepIterative0[context_Association] /; KeyExistsQ[context, "families"] := With[{sectorContext = context}, {
    HoldPattern[int_J] /; integralKind[int] === "Tree" &&
       AssociationQ[treeFamilyForIntegral[int, sectorContext]] &&
       And @@ (IntegerQ /@ First[int][[All, 1]]) &&
       AnyTrue[First[int][[All, 1]], # != 0 &] :>
     treeFirstStepToZero[int, treeFamilyForIntegral[int, sectorContext]]
    }];


treeSectorEndpoints[end_, family_Association, context_Association] := Which[
   end === Automatic, ConstantArray[0, Length[family["vertices"]]],
   AssociationQ[end], Lookup[end, family["sector"], ConstantArray[0, Length[family["vertices"]]]],
   family["sector"] === context["topFamily"]["sector"], end,
   True, ConstantArray[0, Length[family["vertices"]]]
   ];


treeSectorIntegralDistanceData[int_J, end_, context_Association] := Module[{family, endpoints},
   family = treeFamilyForIntegral[int, context];
   If[Head[family] === Missing, Return[$Failed]];
   endpoints = normalizeTreeEndpoints[treeSectorEndpoints[end, family, context], family];
   If[endpoints === $Failed || ! And @@ (IntegerQ /@ First[int][[All, 1]]), Return[$Failed]];
   <|"integral" -> int, "sector" -> family["sector"], "endpoints" -> endpoints,
    "remainingThetaLines" -> Length[thetaFullLineIndices[family["topology"]]],
    "distance" -> treeEndpointDistance[int, endpoints],
    "progressKey" -> {Length[thetaFullLineIndices[family["topology"]]], treeEndpointDistance[int, endpoints]}|>
   ];


treeSectorStepProgress[source_J, replacement_, end_, context_Association] := Module[
   {sourceData, targets, targetData, invalidTargets},
   sourceData = treeSectorIntegralDistanceData[source, end, context];
   targets = DeleteDuplicates[Cases[replacement, int_J /; integralKind[int] === "Tree", {0, Infinity}]];
   targetData = treeSectorIntegralDistanceData[#, end, context] & /@ targets;
   invalidTargets = Pick[targets, Head[#] === Symbol && # === $Failed & /@ targetData];
   <|
    "passQ" -> TrueQ[sourceData =!= $Failed && invalidTargets === {} &&
       And @@ (treeProgressKeyLessQ[Lookup[#, "progressKey", {Infinity, Infinity}], sourceData["progressKey"]] & /@ targetData)],
    "source" -> sourceData,
    "targets" -> targetData,
    "invalidTargets" -> invalidTargets
    |>
   ];


Options[repIterativeSectorData] = Options[repIterativeData];


repIterativeSectorData[expr_, end_: Automatic, context_Association, OptionsPattern[]] := Module[
   {result = Expand[expr], steps = 0, integrals, item, family, endpoints, packs, vertexIndex, unresolved,
    scanFailure, replacement, progressData, seenStates = <||>, stateHash},
   AssociateTo[seenStates, treeRecurrenceStateHash[result] -> True];
   While[True,
    integrals = DeleteDuplicates[Cases[result, int_J /; integralKind[int] === "Tree", {0, Infinity}]];
    unresolved = {};
    scanFailure = None;
    Do[
     family = treeFamilyForIntegral[item, context];
     If[Head[family] === Missing,
      Message[repIterativeData::nosector, item];
      scanFailure = <|"status" -> "error", "result" -> $Failed, "steps" -> steps|>;
      Break[]
      ];
     If[! And @@ (IntegerQ /@ First[item][[All, 1]]),
      Message[repIterativeData::badindex, First[item][[All, 1]]];
      scanFailure = <|"status" -> "error", "result" -> $Failed, "steps" -> steps|>;
      Break[]
      ];
     endpoints = normalizeTreeEndpoints[treeSectorEndpoints[end, family, context], family];
     If[endpoints === $Failed,
      scanFailure = <|"status" -> "error", "result" -> $Failed, "steps" -> steps|>;
      Break[]
      ];
     If[Or @@ MapThread[Unequal, {First[item][[All, 1]], endpoints}], AppendTo[unresolved, item]],
     {item, integrals}
     ];
    If[AssociationQ[scanFailure], Return[scanFailure]];
    If[unresolved === {}, Break[]];
    item = First[unresolved];
    family = treeFamilyForIntegral[item, context];
    endpoints = normalizeTreeEndpoints[treeSectorEndpoints[end, family, context], family];
    If[endpoints === $Failed, Return[<|"status" -> "error", "result" -> $Failed, "steps" -> steps|>]];
    packs = First[item];
    vertexIndex = SelectFirst[Range[Length[endpoints]], packs[[#, 1]] =!= endpoints[[#]] &];
    replacement = treeSingleStepIntegral[item, vertexIndex, endpoints[[vertexIndex]], family];
    If[replacement === $Failed, Return[<|"status" -> "error", "reason" -> "recurrenceSingular", "result" -> $Failed, "steps" -> steps|>]];
    progressData = treeSectorStepProgress[item, replacement, end, context];
    If[! TrueQ[progressData["passQ"]],
     Message[repIterativeData::noprogress, progressData];
     Return[<|"status" -> "error", "reason" -> "nonDecreasingRecurrence", "result" -> $Failed,
       "steps" -> steps, "progressData" -> progressData|>]
     ];
    result = Expand[result /. item -> replacement];
    stateHash = treeRecurrenceStateHash[result];
    If[KeyExistsQ[seenStates, stateHash],
     Message[repIterativeData::cycle, stateHash];
     Return[<|"status" -> "error", "reason" -> "recurrenceCycle", "result" -> $Failed,
       "steps" -> steps, "stateHash" -> stateHash|>]
     ];
    AssociateTo[seenStates, stateHash -> True];
    steps++;
    ];
   <|"status" -> "reduced", "result" -> result, "steps" -> steps|>
   ];


repIterative[expr_, end_, context_Association, OptionsPattern[]] /; KeyExistsQ[context, "families"] := Module[{activeContext},
   (* 多 sector context 保留既有唯一分派门禁，同时同步本轮可直接替换的单步规则。 *)
   activeContext = setTreeFamilyContext[context];
   repIterativeSectorData[expr, end, activeContext]["result"]
   ];


(* ::Chapter:: *)
(*Tree dlog connection*)

(* dlog 输出始终把 primitive matrix、letter/constant-matrix 对和同序 master list 一起返回。
   多顶点 top block 是按 vertexOrder 的 Kronecker sum；contact source 由 DSTreeSeeds 单独提供。 *)

treeVertexDLogData[vertex_Association] := Module[
   {p, states, momentumMagnitudes, k0, omega0, omegaEx, tp, tpInv, m1, omega, letters, coeffs},
   p = vertex["p"];
   states = treeBinaryStates[p];
   momentumMagnitudes = Lookup[vertex["massiveLegs"], "momentumMagnitude"];
   k0 = vertex["signedExternalLegEnergy"];
   omega0 = -I DiagonalMatrix[Table[
      Log[k0 + Sum[(2 states[[row, i]] - 1) momentumMagnitudes[[i]], {i, p}]],
      {row, Length[states]}
      ]];
   omegaEx = DiagonalMatrix[Table[
     -Sum[states[[row, i]] (2 vertex["massiveLegs"][[i, "nu"]] + 1) Log[momentumMagnitudes[[i]]], {i, p}],
      {row, Length[states]}
      ]];
   tp = treeTp[vertex];
   tpInv = treeTpInverse[vertex];
   m1 = treeM1[vertex, vertex["nu0"] + 1];
   omega = Expand[omegaEx - I tpInv . omega0 . tp . m1];
   (* 每个顶点先列 massive-line momentum magnitudes，再按 binary master order 列 cut letters；
      多顶点输出按 vertexOrder 稳定拼接，保证矩阵序列化可复现。 *)
   letters = DeleteDuplicates@Join[momentumMagnitudes, Cases[omega0, Log[arg_] :> arg, Infinity]];
   coeffs = Association@Table[letter -> Map[Coefficient[#, Log[letter]] &, omega, {2}], {letter, letters}];
   <|"vertex" -> vertex["id"], "states" -> states, "omega" -> omega, "letters" -> letters, "letterMatrices" -> coeffs|>
   ];


treeEmbedVertexMatrix[matrix_List, vertexIndex_Integer, dimensions_List] := treeTensorProduct[
   Table[If[i === vertexIndex, matrix, IdentityMatrix[dimensions[[i]]]], {i, Length[dimensions]}]
   ];


dsTreeDLogBlock018[data_Association] := Module[{vertices, dimensions, omega, letters, letterMatrices, masters},
   vertices = treeVertexDLogData /@ data["vertices"];
   dimensions = 2^Lookup[data["vertices"], "p"];
   omega = Total[MapIndexed[treeEmbedVertexMatrix[#1["omega"], First[#2], dimensions] &, vertices]];
   letters = DeleteDuplicates[Flatten[Lookup[vertices, "letters"]]];
   letterMatrices = Association@Table[
      letter -> Map[Coefficient[#, Log[letter]] &, omega, {2}],
      {letter, letters}
      ];
   masters = treeMasterList[data];
   <|
    "status" -> "generated",
    "sector" -> data["sector"],
    "omega" -> omega,
    "letters" -> letters,
    "letterMatrices" -> letterMatrices,
    "masters" -> masters,
    "masterCount" -> Length[masters],
    "vertexOrder" -> data["vertexOrder"],
    "vertexBlocks" -> vertices,
    "sourceStructure" -> data["sourceStructure"]
    |>
   ];


dsTreeDLogBlock018[data_Association, seedData_] := Join[
   dsTreeDLogBlock018[data],
   <|"sourceEquations" -> seedData|>
   ];
