(* ::Package:: *)
(* 本模块从闭合 DE 构造 Euler residual，并按显式次数、通用 loop topology
   或 pure-massive-bubble reference convention 生成 master 的齐次次数。 *)

(* ::Chapter:: *)
(*标度关系检查*)

(* 标度证书的范围由实际固定的 Euler 变量决定；postDerivative 也可能只固定非 DE 参数，
   此时 DE 仍可提供全符号证书。 *)

Options[DSScaleCheck] = {
   ScalingRelation -> "Custom",
   ScalingVariables -> Automatic,
   ScalingWeights -> Automatic,
   ScalingDegrees -> Automatic,
   ProgressReporting -> Automatic
   };

DSScaleCheck::badde = "DSScaleCheck 需要 DSDE 返回的 generated DE 数据。";
DSScaleCheck::badspec = "标度 relation/variables/weights/degrees 不完整或长度不一致：`1`。";

dsIntegralPhysicalPowers[int : J[aList_List, linePacks_List, _List], context_Association] := Module[
   {sectorTopo, activeVertices, aPowers, bPowers},
   sectorTopo = sectorTopologyForIntegral018[context["topology"], int];
   If[Head[sectorTopo] === Missing, Return[$Failed]];
   activeVertices = activeAVertexIds[sectorTopo];
   aPowers = MapThread[#1 + vertexZeroPoint[sectorTopo, #2] &, {aList, activeVertices}];
   bPowers = Table[linePowerIndex[sectorTopo, int, e], {e, sectorTopo["nE"]}];
   <|"aPowers" -> aPowers, "bPowers" -> bPowers,
    "sectorKey" -> sectorKeyFromShrunkLines[sectorTopo, Lookup[sectorTopo, "sectorShrunkLines", {}]]|>
   ];

dsPureMassiveBubbleDegree[int_J, context_Association] := Module[{powers, vertexCount, offset},
   powers = dsIntegralPhysicalPowers[int, context];
   If[powers === $Failed, Return[$Failed]];
   vertexCount = Length[powers["aPowers"]];
   offset = Switch[vertexCount, 2, 2, 1, 1, _, Return[$Failed]];
   dim - Total[powers["bPowers"]] - Total[powers["aPowers"]] - offset
   ];


(* ::Section::Closed:: *)
(*通用 loop topology 的 normalized J 次数*)

(* ISP 是 loop scalar product，因而每个幂次贡献两个动量次数。time-only 公开 J 必须先经
   唯一边界转换恢复内部 sector 表示；sector prefactor 再从完整 N_s 结构读取。 *)
dsLoopTopologyDegree[
   int_J,
   variables_List,
   weights_List,
   context_Association
   ] := Module[
   {powers, rootTopo, internalInt, ispPowers, loopCount, vertexCount, prefactorData,
     prefactor, prefactorDegree},
   rootTopo = context["topology"];
   internalInt = If[
     Lookup[rootTopo, "ibpMode", "full"] === "timeOnly",
     dsTimeOnlyPublicIntegralToInternal020[int, rootTopo],
     int
     ];
   If[! MatchQ[internalInt, J[_List, _List, _List]], Return[$Failed]];
   ispPowers = internalInt[[3]];
   powers = dsIntegralPhysicalPowers[internalInt, context];
   If[powers === $Failed, Return[$Failed]];
   loopCount = Lookup[rootTopo, "graphLoopCount", Length[Lookup[rootTopo, "loopMomenta", {}]]];
   vertexCount = Length[powers["aPowers"]];
   prefactorData = sectorPrefactorDataForIntegral018[rootTopo, internalInt];
   prefactor = materializeSectorPrefactor018[prefactorData];
   If[prefactor === $Failed || TrueQ[prefactor === 0], Return[$Failed]];
   prefactorDegree = Together[
     Total[MapThread[#1 #2 D[prefactor, #2] &, {weights, variables}]]/prefactor
     ];
   If[! And @@ (dsScaleZeroQ[D[prefactorDegree, #]] & /@ variables), Return[$Failed]];
   Together[
    loopCount dim - Total[powers["bPowers"]] - Total[powers["aPowers"]] -
     vertexCount + 2 Total[ispPowers] + prefactorDegree
    ]
   ];


dsLoopTopologyExpressionDegree[
   expr_,
   variables_List,
   weights_List,
   context_Association
   ] := Module[
   {linearData, termDegrees, coefficientDegree, integralDegree, referenceDegree},
   linearData = publicLinearIntegralDecomposition[expr];
   If[
    Lookup[linearData, "status", "failed"] =!= "linear" ||
     ! TrueQ[linearData["constantTerm"] === 0],
    Return[$Failed]
    ];
   termDegrees = MapThread[
     Function[{coefficient, int},
      If[
       TrueQ[coefficient === 0],
       Nothing,
       coefficientDegree = Together[
         Total[MapThread[#1 #2 D[coefficient, #2] &, {weights, variables}]]/coefficient
         ];
       integralDegree = dsLoopTopologyDegree[int, variables, weights, context];
       If[
        integralDegree === $Failed ||
         ! And @@ (dsScaleZeroQ[D[coefficientDegree, #]] & /@ variables),
        $Failed,
        Together[coefficientDegree + integralDegree]
        ]
       ]
      ],
     {linearData["coefficients"], linearData["integrals"]}
     ];
   If[termDegrees === {} || MemberQ[termDegrees, $Failed], Return[$Failed]];
   referenceDegree = First[termDegrees];
   If[
    And @@ (dsScaleZeroQ[# - referenceDegree] & /@ Rest[termDegrees]),
    referenceDegree,
    $Failed
    ]
   ];

(* 线性组合的次数同时包含显式动力学系数；所有非零项必须具有同一 Euler 次数。 *)
dsPureMassiveBubbleExpressionDegree[expr_, variables_List, weights_List, context_Association] := Module[
   {linearData, termDegrees, coefficientDegree, integralDegree, referenceDegree},
   linearData = publicLinearIntegralDecomposition[expr];
   If[Lookup[linearData, "status", "failed"] =!= "linear" || ! TrueQ[linearData["constantTerm"] === 0], Return[$Failed]];
   termDegrees = MapThread[
     Function[{coefficient, int},
      If[TrueQ[coefficient === 0],
       Nothing,
       coefficientDegree = Together[Total[MapThread[#1 #2 D[coefficient, #2] &, {weights, variables}]]/coefficient];
       integralDegree = dsPureMassiveBubbleDegree[int, context];
       If[integralDegree === $Failed || ! And @@ (dsScaleZeroQ[D[coefficientDegree, #]] & /@ variables),
        $Failed,
        Together[coefficientDegree + integralDegree]
        ]
       ]
      ],
     {linearData["coefficients"], linearData["integrals"]}
     ];
   If[termDegrees === {} || MemberQ[termDegrees, $Failed], Return[$Failed]];
   referenceDegree = First[termDegrees];
   If[And @@ (dsScaleZeroQ[# - referenceDegree] & /@ Rest[termDegrees]), referenceDegree, $Failed]
   ];

dsScaleZeroQ[expr_] := TrueQ[Together[Expand[expr]] === 0];

DSScaleCheck[deData_Association, spec_: <||>, OptionsPattern[]] := Module[
   {relation, variables, weights, degrees, masters, context, matrices, sources, missingVariables,
    declaredDegrees, sourceManifest, kiraPlan, numericStage, certificateScope, fixedScalingVariables,
    postDerivativeRules, physicalPostDerivativeRules, degreeRules, evaluationRules,
    evaluatedVariables, eulerMatrix, eulerSource,
     matrixResidual, sourceResidual, checks, status, matrixResidualCount, sourceResidualCount},
   If[Lookup[deData, "status", "missing"] =!= "generated",
    Message[DSScaleCheck::badde]; dsErrorPrint["DE 尚未闭合，不能宣称标度检查通过。 The DE is not closed, so the scaling check cannot be reported as passed."]; Return[<|"status" -> "failed", "reason" -> "deNotGenerated"|>]
    ];
   relation = Lookup[spec, "relation", OptionValue[ScalingRelation]];
   variables = Replace[Lookup[spec, "variables", OptionValue[ScalingVariables]], Automatic -> deData["variables"]];
   weights = Replace[Lookup[spec, "weights", OptionValue[ScalingWeights]], Automatic -> ConstantArray[1, Length[variables]]];
   masters = deData["masters"];
   context = deData["context"];
   declaredDegrees = Lookup[Lookup[deData, "activeBasis", <||>], "scalingDegrees", Automatic];
   sourceManifest = Lookup[deData, "sourceManifest", <||>];
   kiraPlan = Lookup[sourceManifest, "kiraPlan", <||>];
   numericStage = Lookup[kiraPlan, "numericStage", "symbolic"];
   postDerivativeRules = If[
     numericStage === "postDerivative",
     Lookup[kiraPlan, "coefficientRules", {}],
     {}
     ];
   physicalPostDerivativeRules = If[postDerivativeRules === {},
     {},
     Lookup[sourceManifest, "physicalCoefficientRulesApplied", postDerivativeRules]
     ];
   fixedScalingVariables = If[
     ListQ[physicalPostDerivativeRules],
     Select[variables, MemberQ[First /@ physicalPostDerivativeRules, #] &],
     {}
     ];
   certificateScope = Which[
     fixedScalingVariables === {}, "symbolic",
     Length[fixedScalingVariables] === Length[variables], "exactPoint",
     True, "fixedSection"
     ];
   degrees = Replace[
     Lookup[spec, "degrees", OptionValue[ScalingDegrees]],
     Automatic :> Which[
       ListQ[declaredDegrees], declaredDegrees,
       relation === "PureMassiveBubble", dsPureMassiveBubbleExpressionDegree[#, variables, weights, context] & /@ masters,
       relation === "LoopTopology", dsLoopTopologyExpressionDegree[#, variables, weights, context] & /@ masters,
       True, $Failed
       ]
     ];
   (* 定点规则中的非微分参数进入齐次次数；DE 坐标本身只用于构造同一点的 Euler 算符。 *)
   degreeRules = If[
     ListQ[physicalPostDerivativeRules],
     Select[physicalPostDerivativeRules, ! MemberQ[variables, First[#]] &],
     {}
     ];
   If[ListQ[degrees], degrees = degrees /. degreeRules];
   If[! ListQ[variables] || ! ListQ[weights] || Length[variables] =!= Length[weights] ||
     ! ListQ[degrees] || Length[degrees] =!= Length[masters] || MemberQ[degrees, $Failed],
    Message[DSScaleCheck::badspec, dsScaleSpecSentence[<|"relation" -> relation, "variables" -> variables, "weights" -> weights, "degrees" -> degrees|>]];
    dsErrorPrint["标度检查规格无效。 The scaling-check specification is invalid."]; Return[<|"status" -> "failed", "reason" -> "invalidScalingSpecification"|>]
    ];
   matrices = deData["matrices"];
   sources = deData["sources"];
   missingVariables = Select[variables, ! KeyExistsQ[matrices, #] &];
   If[missingVariables =!= {},
    Message[DSScaleCheck::badspec, missingVariables]; dsErrorPrint["DE 缺少 Euler 算符所需变量。 The DE lacks variables required by the Euler operator."]; Return[<|"status" -> "failed", "reason" -> "missingDEVariables", "missingVariables" -> missingVariables|>]
    ];
   evaluationRules = If[ListQ[physicalPostDerivativeRules], physicalPostDerivativeRules, {}];
   evaluatedVariables = variables /. evaluationRules;
   eulerMatrix = Total@MapThread[
     #1 #2 (#3 /. evaluationRules) &,
     {weights, evaluatedVariables, Lookup[matrices, #] & /@ variables}
     ];
   eulerSource = Total@MapThread[
     #1 #2 (#3 /. evaluationRules) &,
     {weights, evaluatedVariables, Lookup[sources, #] & /@ variables}
     ];
   matrixResidual = Map[Together[Expand[#]] &, eulerMatrix - DiagonalMatrix[degrees], {2}];
   sourceResidual = Together[Expand[#]] & /@ eulerSource;
   checks = Switch[
     certificateScope,
     "symbolic",
       <|
        "symbolicMatrixRelation" -> And @@ (dsScaleZeroQ /@ Flatten[matrixResidual]),
        "symbolicSourceRelation" -> And @@ (dsScaleZeroQ /@ sourceResidual)
        |>,
     "exactPoint",
       <|
        "exactPointMatrixRelation" -> And @@ (dsScaleZeroQ /@ Flatten[matrixResidual]),
        "exactPointSourceRelation" -> And @@ (dsScaleZeroQ /@ sourceResidual)
        |>,
     _,
       <|
        "fixedSectionMatrixRelation" -> And @@ (dsScaleZeroQ /@ Flatten[matrixResidual]),
        "fixedSectionSourceRelation" -> And @@ (dsScaleZeroQ /@ sourceResidual)
        |>
     ];
    status = If[And @@ Values[checks], "passed", "failed"];
    If[status === "failed",
     matrixResidualCount = Count[Flatten[matrixResidual], item_ /; ! dsScaleZeroQ[item]];
     sourceResidualCount = Count[sourceResidual, item_ /; ! dsScaleZeroQ[item]];
     Message[DSScaleCheck::failed, certificateScope, matrixResidualCount, sourceResidualCount];
     dsErrorPrint[
      "当前微分方程不满足所声明的标度关系，结果已保留为失败诊断而不是通过证书。请检查返回的 matrixResidual 和 sourceResidual。 " <>
       "The current differential equations do not satisfy the declared scaling relation. The result is retained as a failed diagnostic, not as a passed certificate. Inspect matrixResidual and sourceResidual."
      ]
     ];
    <|
    "status" -> status,
    "numericStage" -> numericStage,
    "certificateScope" -> certificateScope,
    "fixedScalingVariables" -> fixedScalingVariables,
    "relation" -> relation,
    "variables" -> variables,
    "weights" -> weights,
    "evaluatedVariables" -> evaluatedVariables,
    "evaluationPointRules" -> evaluationRules,
    "degrees" -> degrees,
    "eulerMatrix" -> eulerMatrix,
    "eulerSource" -> eulerSource,
    "matrixResidual" -> matrixResidual,
    "sourceResidual" -> sourceResidual,
    "checks" -> checks,
    "symbolicQ" -> (certificateScope === "symbolic")
    |>
   ];

DSScaleCheck[deData_, spec_: <||>, OptionsPattern[]] := (Message[DSScaleCheck::badde]; dsErrorPrint["DSScaleCheck 输入必须是 DE Association。 DSScaleCheck input must be a DE Association."]; <|"status" -> "failed", "reason" -> "inputNotAssociation"|>);
