(* ::Package:: *)
(* mix bubble+tree 示例：bubble 的 cycle momenta 为 l1 与 l1+k1+k2，前者 massive、
   后者 massless；massless bridge 和 bubble 另一外腿均携带 k1+k2，三点顶点另接 k1、k2。
   本例覆盖 kL/kE 角色、参数重定义、massless 三槽 convention、cycle/bridge contact、圈指标全偶
   与坐标完备性；固定 bridge 不属于圈积分，因此不施加 parity。 *)

(* ::Chapter:: *)
(*加载标准 package*)

exampleDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{exampleDir, "..", "load_current_package.wl"}], CharacterEncoding -> "UTF-8"];
kiraInputSummary = Get[FileNameJoin[{exampleDir, "kira_input_summary.wl"}]];
kiraResultSummary = Get[FileNameJoin[{exampleDir, "kira_result_summary.wl"}]];


(* ::Chapter:: *)
(*详细物理输入*)

(* 两条 cycle propagator 分别对完整 {b,n1,n2} 指标取全偶；line 3 是 fixed bridge，
   没有圈积分 b 指标，不进入 parity constraint。 *)
cycleParityConstraints = {
   b[1] + n[1, 1] + n[1, 2] -> 0,
   b[2] + n[2, 1] + n[2, 2] -> 0
   };

(* loopExternalMomenta 只需给 shift-invariant loop 方向 {k1+k2}；其完整 Gram 产生 ss11。
   independentExternalMomenta 只给 loop Gram 尚未覆盖的实际无圈模长 {k1,k2}，产生 sE1,sE2。
   bridge/bubble 外腿的 |k1+k2| 已由 ss11 覆盖，不再生成 sE3 或 sp[k1,k2]。 *)
caseInput = <|
   "name" -> "mixBubbleTreeK1K2",
   "vertices" -> {
     <|"id" -> v1, "vertexType" -> "+", "externalLegEnergy" -> E1|>,
     <|"id" -> v2, "vertexType" -> "+", "externalLegEnergy" -> E2|>,
     <|"id" -> v3, "vertexType" -> "+", "externalLegEnergy" -> E3|>
     },
   "lines" -> {
     <|"id" -> 1, "endpoints" -> {v1, v2}, "momentum" -> l1,
       "nu" -> nu1, "massType" -> "massive"|>,
     <|"id" -> 2, "endpoints" -> {v1, v2}, "momentum" -> l1 + k1 + k2,
       "massType" -> "massless"|>,
     <|"id" -> 3, "endpoints" -> {v2, v3}, "momentum" -> k1 + k2,
       "massType" -> "massless"|>
     },
   "extLegs" -> {
     {bubbleLeg, v1, k1 + k2},
     {treeLeg1, v3, k1},
     {treeLeg2, v3, k2}
     },
   "loopMomenta" -> {l1},
   "loopExternalMomenta" -> {k1 + k2},
   "independentExternalMomenta" -> {k1, k2},
   "ibpMode" -> "full",
   "ispData" -> {},
   "zeroPointRules" -> {
     a0[v1] -> alpha1, a0[v2] -> alpha2, a0[v3] -> alpha3,
     b0[1] -> beta1, b0[2] -> beta2, b0[3] -> beta3
     },
   "symmetryRules" -> {},
   "parityConstraints" -> cycleParityConstraints
   |>;


(* ::Chapter:: *)
(*变量提案与缺省选项*)

kinematicProposal = DSKinematics[caseInput];
kinematicProposal["selectionTemplate"]
kinematicProposal["parameterRedefinitionGuide"]

(* 缺省为 KinematicRules->Automatic，直接采用上一单元格显示的规则。
   其它缺省为不写 init、不预生成微分 metadata、不覆盖不同输入、注册为当前 context。 *)
kinematicRules = Automatic;
initOptions = {
   WriteInitializationFiles -> False,
   GenerateDerivativeMetadata -> True,
   OverwriteInitialization -> False,
   RegisterAsCurrent -> True,
   ProgressReporting -> True,
   KinematicRules -> kinematicRules
   };


(* ::Chapter:: *)
(*初始化、notation 与公开原子操作*)

DSMessagesOn[];
context = DSInit[caseInput, Sequence @@ initOptions];
topology = Lookup[context, "topology", <||>];
contextInfo = DSInfo[context, "Full"];
notation = DSParameterNotation[context];
publicAPI = DSPublicAPI[];

(* cycle/fixed full line 都保留双端点三槽；fixed line 以短字符串 "F" 标记无 b 指标。 *)
integral = J[
   {0, 0, 0},
   {{0, 0, 0}, {0, 0, 0}, {"F", 0, 1}},
   {}
   ];

timeSeed = dtau[v3, integral, topology];
zeroContactIntegral = ReplacePart[integral, {2, 3} -> {"F", 0, 0}];
zeroContactTimeSeed = dtau[v3, zeroContactIntegral, topology];
loopLoopSeed = dqq[1, 1, integral, topology];
loopExternalSeed = dqk[1, 1, integral, topology];
totalDerivative = ds[ss11^2 integral + sE1 integral, sE1, topology];
innerCoordinate = rep2innerform[sp[k1 + k2, k1 + k2], topology];
outerCoordinate = rep2outform[innerCoordinate, topology];
integrandForm = rep2Integrand[integral, topology];
rawSymmetryRules = repSymmetry0[topology];
canonicalIntegral = symmetry[integral, topology];
expectedTotalDerivative = -ReplacePart[integral, {2, 3} -> {"F", 1, 0}];


(* ::Chapter:: *)
(*用户参数重定义*)

(* 左端必须复制 baseCoordinateOrder 中的原始 sp[...]，不能写 ss11->loopScale。
   右端写自定义参数表达式；规则必须覆盖全部基础方向。新 context 的 seed、ds、DSDE
   和 GenerateDerivativeMetadata 都改用 {loopScale,legScale1,legScale2}。 *)
customRules = {
   sp[k1 + k2, k1 + k2] -> loopScale^2,
   sp[k1, k1] -> legScale1^2,
   sp[k2, k2] -> legScale2^2
   };
customContext = DSRedefineParameters[context, customRules, ProgressReporting -> True];
customTopology = Lookup[customContext, "topology", <||>];
customNotation = DSParameterNotation[customContext];
customDerivative = ds[loopScale^2 integral + legScale1 integral, legScale1, customTopology];
customDerivativeVariables = Lookup[Lookup[customContext, "derivatives", <||>], "operators", {}][[All, "userVariable"]];

(* 只构造 general templates，用 sector key 确认 massive cycle、massless cycle 与 massless bridge
   都进入 contact-reachable 数据；本例不连续撒点，也不运行 reduction。 *)
templateSeedData = DSSeeds[customContext, ProgressReporting -> False];
templateRecords = DSAllSeeds[templateSeedData];
templateSectorKeys = DeleteDuplicates[Lookup[templateRecords, "sectorKey", {}]];


(* ::Chapter:: *)
(*两个外腿模长和的用户坐标绑定*)

(* 保留 k1、k2 两条外腿，并令 E0=|k1|+|k2|。规则左端仍写原始平方不变量；
   当前参数化选择 |k1|=E0-sE2 的物理支，因此应用区域应满足 E0>sE2>0。 *)
boundEnergyInput = Join[
   caseInput,
   <|"name" -> "bubbleTreeBoundEnergy", "vertices" -> ReplacePart[caseInput["vertices"], {3, "externalLegEnergy"} -> E0]|>
   ];
boundEnergyContext0 = DSInit[
   boundEnergyInput,
   WriteInitializationFiles -> False,
   GenerateDerivativeMetadata -> True,
   RegisterAsCurrent -> False,
   ProgressReporting -> False
   ];
boundEnergyContext = DSRedefineParameters[
   boundEnergyContext0,
   {
    sp[k1 + k2, k1 + k2] -> ss11^2,
    sp[k1, k1] -> (E0 - sE2)^2,
    sp[k2, k2] -> sE2^2
    },
   ProgressReporting -> False
   ];
boundEnergyVariables = Lookup[
   Lookup[boundEnergyContext, "derivatives", <||>],
   "operators",
   {}
   ][[All, "userVariable"]];


(* ::Chapter:: *)
(*单一有效外腿的独立拓扑输入*)

(* v3 不再保留 k1、k2 两条外腿，而改为一条 p0 外腿；E0 是 |p0| 及 v3 相位能量。
   这是新的 topology metadata，不是原 family 的纯坐标重命名。 *)
singleLegInput = Join[
   caseInput,
   <|
    "name" -> "bubbleTreeSingleEffectiveLeg",
    "extLegs" -> {
      {bubbleLeg, v1, k1 + k2},
      {treeLeg0, v3, p0}
      },
    "independentExternalMomenta" -> {p0},
    "vertices" -> ReplacePart[caseInput["vertices"], {3, "externalLegEnergy"} -> E0]
    |>
   ];
singleLegContext0 = DSInit[
   singleLegInput,
   WriteInitializationFiles -> False,
   GenerateDerivativeMetadata -> True,
   RegisterAsCurrent -> False,
   ProgressReporting -> False
   ];
singleLegContext = DSRedefineParameters[
   singleLegContext0,
   {
    sp[k1 + k2, k1 + k2] -> ss11^2,
    sp[p0, p0] -> E0^2
    },
   ProgressReporting -> False
   ];
singleLegVariables = Lookup[
   Lookup[singleLegContext, "derivatives", <||>],
   "operators",
   {}
   ][[All, "userVariable"]];


(* ::Chapter:: *)
(*动量列表与自定义坐标的过完备/欠完备输入*)

(* loop 列表过完备：多给独立方向 k1-k2；红色 warning 后可生成 symbolic IBP，
   但 ds、DSDE 与唯一 rep2innerform 关闭。 *)
overLoopInput = Join[caseInput, <|"loopExternalMomenta" -> {k1 + k2, k1 - k2}|>];
overLoopAudit = DSKinematics[overLoopInput];
overLoopContext = DSInit[overLoopInput, RegisterAsCurrent -> False, ProgressReporting -> False];

(* loop 列表欠完备：{k1} 不覆盖 k1+k2；红色 error 给出 missingDirections 并拒绝初始化。 *)
underLoopInput = Join[caseInput, <|"loopExternalMomenta" -> {k1}|>];
underLoopAudit = DSKinematics[underLoopInput];
underLoopContext = DSInit[underLoopInput, RegisterAsCurrent -> False, ProgressReporting -> False];

(* independent 列表过完备：k1-k2 的模长已由 {|k1+k2|,|k1|,|k2|} 决定。 *)
overIndependentInput = Join[caseInput, <|"independentExternalMomenta" -> {k1, k2, k1 - k2}|>];
overIndependentAudit = DSKinematics[overIndependentInput];
overIndependentContext = DSInit[overIndependentInput, RegisterAsCurrent -> False, ProgressReporting -> False];

(* independent 列表欠完备：只给 k1 时缺少 |k2|。 *)
underIndependentInput = Join[caseInput, <|"independentExternalMomenta" -> {k1}|>];
underIndependentAudit = DSKinematics[underIndependentInput];
underIndependentContext = DSInit[underIndependentInput, RegisterAsCurrent -> False, ProgressReporting -> False];

(* 自定义坐标过完备：同一个基础平方再绑定 alternateLegScale1；欠完备：删除 |k2| 规则。 *)
overCustomRules = Append[customRules, sp[k1, k1] -> alternateLegScale1^2];
overCustomAudit = DSKinematics[caseInput, overCustomRules];
overCustomContext = DSRedefineParameters[context, overCustomRules, ProgressReporting -> False];
underCustomRules = Most[customRules];
underCustomAudit = DSKinematics[caseInput, underCustomRules];
underCustomContext = DSRedefineParameters[context, underCustomRules, ProgressReporting -> False];


(* ::Chapter:: *)
(*结果检查*)

bridgePacks = Cases[
   {timeSeed, loopLoopSeed, loopExternalSeed, totalDerivative, customDerivative},
   J[_, packs_, _] :> packs[[3]],
   Infinity
   ];

coverage = Lookup[notation, "requiredMagnitudeCoverage", {}];
coverageRules = Lookup[coverage, "squaredExpression", {}] -> Lookup[coverage, "userSquaredExpression", {}];
defaultDerivativeVariables = Lookup[
   Lookup[context, "derivatives", <||>],
   "operators",
   {}
   ][[All, "userVariable"]];
linePackTypes = Lookup[Lookup[topology, "lines", {}], "packType", {}];
storedParityConstraints = Lookup[topology, "parityConstraints", {}];

DSMessagesOff[];
messagesDisabledQ = ! DSMessagesQ[];
DSMessagesOn[];

summary = <|
   "version" -> dSIBP`$dSIBPVersion,
   "initStatus" -> Lookup[context, "status", "missing"],
   "graphLoopCount" -> Lookup[topology, "graphLoopCount", Missing["loopCount"]],
   "loopExternalMomenta" -> Lookup[notation, "loopExternalMomenta", {}],
   "independentExternalMomenta" -> Lookup[notation, "independentExternalMomenta", {}],
   "notationVariables" -> Lookup[notation, "selectedUserVariables", {}],
   "linePackTypes" -> linePackTypes,
   "parityConstraints" -> storedParityConstraints,
   "cycleParityIsEven" -> (storedParityConstraints === cycleParityConstraints),
   "fixedBridgeExcludedFromParity" -> FreeQ[storedParityConstraints, b[3] | n[3, _]],
   "requiredMagnitudeCoverage" -> coverage,
   "defaultDerivativeVariables" -> defaultDerivativeVariables,
   "customVariables" -> Lookup[customNotation, "selectedUserVariables", {}],
   "customDerivativeVariables" -> customDerivativeVariables,
   "templateStatus" -> Lookup[templateSeedData, "dSIBPStatus", Lookup[templateSeedData, "status", "missing"]],
   "templateSectorKeys" -> templateSectorKeys,
   "boundEnergyStatus" -> Lookup[boundEnergyContext, "status", "failed"],
   "boundEnergyVariables" -> boundEnergyVariables,
   "singleLegStatus" -> Lookup[singleLegContext, "status", "failed"],
   "singleLegMomenta" -> Lookup[Lookup[singleLegContext, "topology", <||>], "independentExternalMomenta", {}],
   "singleLegVariables" -> singleLegVariables,
   "redefinitionUsesOriginalSP" -> StringContainsQ[Lookup[Lookup[notation, "parameterRedefinitionGuide", <||>], "commandExample", ""], "sp[k1 + k2, k1 + k2]"],
   "bridgePackHasNoB" -> And @@ (FreeQ[#, _b | _bS] & /@ bridgePacks),
   "bridgePacks" -> bridgePacks,
   "timeSeedUsesBridgeMagnitude" -> ! FreeQ[timeSeed, ss11],
   "timeSeedContractsBridge" -> ! FreeQ[timeSeed, J[_, {_, _, {"F"}}, _]],
   "zeroContactHasNoFailure" -> FreeQ[zeroContactTimeSeed, $Failed],
   "momentumSeedsUseCanonicalBridgePack" -> And @@ (SameQ[#, {"F", 1, 0}] & /@ Cases[{loopLoopSeed, loopExternalSeed}, J[_, packs_, _] :> packs[[3]], Infinity]),
   "totalDerivativeCoefficientTermExact" -> SameQ[totalDerivative, expectedTotalDerivative],
   "coordinateRoundTrip" -> Together[outerCoordinate - ss11^2] === 0,
   "integrandUsesFixedBridgeMagnitude" -> ! FreeQ[integrandForm, ss11] && FreeQ[integrandForm, xi[3]],
   "rawSymmetryUnionIsList" -> ListQ[rawSymmetryRules],
   "canonicalIntegral" -> canonicalIntegral,
   "publicFunctionCount" -> Length[Lookup[publicAPI, "functions", {}]],
   "overLoopStatus" -> Lookup[overLoopAudit, "status", "missing"],
   "underLoopStatus" -> Lookup[underLoopAudit, "status", "missing"],
   "overIndependentStatus" -> Lookup[overIndependentAudit, "status", "missing"],
   "underIndependentStatus" -> Lookup[underIndependentAudit, "status", "missing"],
   "overCustomStatus" -> Lookup[overCustomAudit, "status", "missing"],
   "underCustomStatus" -> Lookup[underCustomAudit, "status", "missing"],
   "overContextsContinue" -> And @@ (Lookup[#, "status", "failed"] === "initialized" & /@ {overLoopContext, overIndependentContext, overCustomContext}),
   "overDerivativesDisabled" -> And @@ (! TrueQ[Lookup[Lookup[#, "capabilities", <||>], "derivativeUsableQ", True]] & /@ {overLoopContext, overIndependentContext, overCustomContext}),
   "underContextsRejected" -> And @@ (Lookup[#, "status", "initialized"] === "failed" & /@ {underLoopContext, underIndependentContext, underCustomContext}),
   "messagesDisabledQ" -> messagesDisabledQ,
   "contextInfoAvailable" -> AssociationQ[contextInfo],
   "kiraInputSummaryCase" -> Lookup[kiraInputSummary, "case", "missing"],
   "kiraResultSummaryStatus" -> Lookup[kiraResultSummary, "status", "missing"],
   "kiraResultMasterCount" -> Lookup[kiraResultSummary, "masterCount", Missing["NotAvailable"]]
   |>;

Print[InputForm[summary]];
If[! And[
    ToString[summary["version"]] === currentVersion,
    summary["initStatus"] === "initialized",
    summary["graphLoopCount"] === 1,
    summary["loopExternalMomenta"] === {k1 + k2},
    summary["independentExternalMomenta"] === {k1, k2},
    summary["notationVariables"] === {ss11, sE1, sE2},
    summary["linePackTypes"] === {"massiveFull", "masslessFull", "masslessFull"},
    summary["cycleParityIsEven"],
    summary["fixedBridgeExcludedFromParity"],
    summary["defaultDerivativeVariables"] === {ss11, sE1, sE2, E1, E2, E3},
    summary["customVariables"] === {loopScale, legScale1, legScale2},
    summary["customDerivativeVariables"] === {
      loopScale, legScale1, legScale2, E1, E2, E3
      },
    summary["templateStatus"] === "generated",
    ContainsAll[summary["templateSectorKeys"], {"top", "e1", "e2", "e3"}],
    summary["boundEnergyStatus"] === "initialized",
    summary["boundEnergyVariables"] === {ss11, E0, sE2, E1, E2},
    summary["singleLegStatus"] === "initialized",
    summary["singleLegMomenta"] === {p0},
    summary["singleLegVariables"] === {ss11, E0, E1, E2},
    summary["redefinitionUsesOriginalSP"],
    summary["bridgePackHasNoB"],
    summary["timeSeedUsesBridgeMagnitude"],
    summary["timeSeedContractsBridge"],
    summary["zeroContactHasNoFailure"],
    summary["momentumSeedsUseCanonicalBridgePack"],
    summary["totalDerivativeCoefficientTermExact"],
    summary["coordinateRoundTrip"],
    summary["integrandUsesFixedBridgeMagnitude"],
    summary["overLoopStatus"] === "overcomplete",
    summary["underLoopStatus"] === "undercomplete",
    summary["overIndependentStatus"] === "overcomplete",
    summary["underIndependentStatus"] === "undercomplete",
    summary["overCustomStatus"] === "overcomplete",
    summary["underCustomStatus"] === "incomplete",
    summary["overContextsContinue"],
    summary["overDerivativesDisabled"],
    summary["underContextsRejected"],
    summary["messagesDisabledQ"],
    summary["contextInfoAvailable"],
    summary["kiraInputSummaryCase"] === "mixBubbleTree",
    summary["kiraResultSummaryStatus"] === "passed",
    summary["kiraResultMasterCount"] === 81
    ], Exit[1]];
