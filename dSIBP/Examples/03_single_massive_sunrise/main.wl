(* ::Package:: *)
(* Single-massive sunrise 示例：两顶点、三条平行边形成两圈，第一条线为 massive h，
   其余两条线为 massless exponential。本例只生成 general IBP seeds 与 general 参数
   微分算符，不撒连续指标点，不生成 linearData，也不进入 Kira、DE 或 scaling。 *)


(* ::Chapter:: *)
(*加载标准 package 与本例 convention*)

exampleDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{exampleDir, "..", "load_current_package.wl"}], CharacterEncoding -> "UTF-8"];
Get[FileNameJoin[{exampleDir, "family_conventions.wl"}]];


(* ::Chapter:: *)
(*详细物理输入*)

(* 两个 ISP 分别以两条 massless line 的动量与 massive line 动量作标量积；
   在 k321 <-> l3-k321-kL 下它们严格互换，因此可与 line 交换共同 canonicalize。 *)
case = <|
   "name" -> "singleMassiveSunriseExample",
   "vertices" -> {
     <|"id" -> v1, "vertexType" -> "+", "externalLegEnergy" -> kE|>,
     <|"id" -> v2, "vertexType" -> "+", "externalLegEnergy" -> kE|>
     },
   "lines" -> {
     <|"id" -> e1, "endpoints" -> {v1, v2}, "momentum" -> l3,
       "massType" -> "massive", "nu" -> nuM|>,
     <|"id" -> e2, "endpoints" -> {v1, v2}, "momentum" -> k321,
       "massType" -> "massless"|>,
     <|"id" -> e3, "endpoints" -> {v1, v2}, "momentum" -> l3 - k321 - kL,
       "massType" -> "massless"|>
     },
   "loopMomenta" -> {l3, k321},
   "loopExternalMomenta" -> {kL},
   "independentExternalMomenta" -> {},
   "ibpMode" -> "full",
   "ispData" -> {
     <|"name" -> rhoMassless2, "expr" -> sp[k321, l3], "range" -> {0, 1}|>,
     <|"name" -> rhoMassless3, "expr" -> sp[l3 - k321 - kL, l3], "range" -> {0, 1}|>
     },
   "zeroPointRules" -> {
     a0[v1] -> alpha, a0[v2] -> alpha,
     b0[e1] -> betaM, b0[e2] -> beta0, b0[e3] -> beta0
     },
   "parityConstraints" -> sunriseParityConstraints0,
   "symmetryRules" -> sunriseSymmetryRules0
   |>;


(* ::Chapter:: *)
(*General seeds 与参数微分算符*)

kinematicProposal = DSKinematics[case];
context = DSInit[
   case,
   WriteInitializationFiles -> False,
   GenerateDerivativeMetadata -> True
   ];
seedData = DSSeeds[context, ProgressReporting -> False];
allSeeds = DSAllSeeds[seedData];

topSeeds = Select[allSeeds, Lookup[#, "sectorKey", None] === "top" &];
topGenerators = DeleteDuplicates[Lookup[topSeeds, "generator", {}]];
topTimeGenerators = DeleteDuplicates @ Lookup[
    Select[topSeeds, Lookup[#, "ibpClass", ""] === "tIBP" &],
    "generator",
    {}
    ];
topMomentumGenerators = DeleteDuplicates @ Lookup[
    Select[topSeeds, Lookup[#, "ibpClass", ""] === "qIBP" &],
    "generator",
    {}
    ];

derivativeData = Lookup[context, "derivatives", <||>];
derivativeOperators = Lookup[derivativeData, "operators", {}];
derivativeVariables = Lookup[derivativeOperators, "userVariable", {}];
generalTopIntegral = SelectFirst[
   Lookup[topSeeds, "sourceIntegral", {}],
   Head[#] === J &,
   Missing["NoTopIntegral"]
   ];
parameterDerivativeWitnesses = AssociationThread[
   derivativeVariables,
   ds[generalTopIntegral, #, Lookup[context, "topology", <||>]] & /@
    derivativeVariables
   ];


(* ::Chapter:: *)
(*结果检查*)

topology = Lookup[context, "topology", <||>];
routingAudit = Lookup[topology, "loopMomentumRoutingAudit", <||>];
canonicalWitness = symmetry[
   J[
    {1, 0},
    {{0, 0, 1}, {2, 0, 1}, {1, 1, 0}},
    {3, 1}
    ],
   topology
   ];
expectedCanonicalWitness = J[
   {0, 1},
   {{0, 1, 0}, {1, 0, 1}, {2, 1, 0}},
   {1, 3}
   ];

summary = <|
   "initStatus" -> Lookup[context, "status", "missing"],
   "kinematicStatus" -> Lookup[kinematicProposal, "status", "missing"],
   "graphLoopCount" -> Lookup[topology, "graphLoopCount", Missing["loopCount"]],
   "cycleLineIndices" -> Lookup[topology, "cycleLineIndices", {}],
   "bridgeLineIndices" -> Lookup[topology, "bridgeLineIndices", {}],
   "loopRoutingRank" -> Lookup[routingAudit, "loopCoefficientRank", Missing["routingRank"]],
   "externalLegEnergyByVertex" -> Lookup[topology, "sectorExternalLegEnergyByVertex", <||>],
   "massiveLineCount" -> Count[Lookup[Lookup[topology, "lines", {}], "massType", {}], "massive"],
   "ispCount" -> Length[Lookup[topology, "ispData", {}]],
   "topGeneratorCount" -> Length[topGenerators],
   "symmetryRuleCount" -> Length[repSymmetry0[topology]],
   "parityConstraints" -> Lookup[topology, "parityConstraints", {}],
   "canonicalWitness" -> canonicalWitness,
   "seedStatus" -> Lookup[seedData, "dSIBPStatus", "missing"],
   "seedTemplateCount" -> Length[allSeeds],
   "seedSectorKeys" -> DeleteDuplicates[Lookup[allSeeds, "sectorKey", {}]],
   "topTimeGenerators" -> topTimeGenerators,
   "topMomentumGenerators" -> topMomentumGenerators,
   "derivativeStatus" -> Lookup[derivativeData, "status", "missing"],
   "derivativeVariables" -> derivativeVariables,
   "parameterDerivativeWitnessCount" -> Length[parameterDerivativeWitnesses],
   "parameterDerivativeFailureQ" ->
    ! FreeQ[Values[parameterDerivativeWitnesses], $Failed],
   "numericSeedResidualCheck" -> "notApplicableForGeneralSymbolicSeeds",
   "containsSampledOrLinearData" ->
    Or @@ (ValueQ /@ {generatedIBP, linearData, kiraPlan, deData, scalingData})
   |>;

Print[summary];
If[! And[
    ToString[dSIBP`$dSIBPVersion] === currentVersion,
    summary["initStatus"] === "initialized",
    summary["kinematicStatus"] === "complete",
    summary["graphLoopCount"] === 2,
    summary["cycleLineIndices"] === {1, 2, 3},
    summary["bridgeLineIndices"] === {},
    summary["loopRoutingRank"] === 2,
    summary["externalLegEnergyByVertex"] === <|v1 -> kE, v2 -> kE|>,
    summary["massiveLineCount"] === 1,
    summary["ispCount"] === 2,
    summary["topGeneratorCount"] === 8,
    summary["symmetryRuleCount"] === 1,
    summary["parityConstraints"] === sunriseParityConstraints0,
    summary["canonicalWitness"] === expectedCanonicalWitness,
    summary["seedStatus"] === "generated",
    summary["seedTemplateCount"] > 0,
    Length[summary["seedSectorKeys"]] > 1,
    Length[summary["topTimeGenerators"]] === 2,
    Length[summary["topMomentumGenerators"]] === 6,
    summary["derivativeStatus"] === "generated",
    summary["derivativeVariables"] === {ss11, kE},
    summary["parameterDerivativeWitnessCount"] === 2,
    ! summary["parameterDerivativeFailureQ"],
    ! summary["containsSampledOrLinearData"]
    ], Exit[1]];
