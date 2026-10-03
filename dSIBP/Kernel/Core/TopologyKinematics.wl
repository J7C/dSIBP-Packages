(* ::Package:: *)
(* 本模块提供图论与动量声明审计。它不生成 IBP，只把 topology 的结构圈数、
   bridge/cycle line、圈动量 routing 以及两类用户外动量列表归一为可供 DSInit 门禁读取的 metadata。 *)

(* ::Chapter:: *)
(*基础线性代数工具*)

(* 只在精确代数意义下判零；本模块不使用数值容差决定 topology。 *)
ds016ZeroQ[expr_] := TrueQ[Together[Expand[expr]] === 0];


ds016MatrixRank[rows_List, width_Integer] := Which[
   rows === {}, 0,
   width === 0, 0,
   True, MatrixRank[rows]
   ];


ds016NonzeroRows[rows_List] := Select[rows, ! And @@ (ds016ZeroQ /@ #) &];


ds016RowBasis[rows_List, width_Integer] := If[
   rows === {} || width === 0,
   {},
   ds016NonzeroRows[RowReduce[rows]]
   ];


(* 按输入顺序选择使行空间增秩的对象；返回位置，便于同时保留来源表达式。 *)
ds016IndependentRowPositions[rows_List, initialRows_List : {}] := Module[
   {selected = {}, current = initialRows, rank, nextRank},
   rank = If[current === {}, 0, MatrixRank[current]];
   Do[
    nextRank = MatrixRank[Append[current, rows[[i]]]];
    If[nextRank > rank,
     AppendTo[selected, i];
     AppendTo[current, rows[[i]]];
     rank = nextRank
     ],
    {i, Length[rows]}
    ];
   selected
   ];


(* momentum 字段只允许由已识别向量原子及精确有理系数组成，避免把符号系数误猜成新向量。 *)
ds016MomentumAtoms[expressions_List, excluded_List : {}] := DeleteDuplicates@Cases[
    Unevaluated[expressions],
    symbol_Symbol /; Context[symbol] =!= "System`" && ! MemberQ[excluded, symbol],
    Infinity,
    Heads -> False
    ];


ds016LinearVectorData[expr_, atoms_List] := Module[{expanded, coefficients, residual, rationalQ},
   expanded = Expand[expr];
   coefficients = Coefficient[expanded, #] & /@ atoms;
   residual = Expand[expanded - coefficients . atoms];
   rationalQ = And @@ (MatchQ[#, _Integer | _Rational] & /@ coefficients);
   <|
    "expression" -> expanded,
    "coefficients" -> coefficients,
    "residual" -> residual,
    "linearQ" -> TrueQ[ds016ZeroQ[residual] && rationalQ]
    |>
   ];


ds016RowsForExpressions[expressions_List, atoms_List] := Lookup[
   ds016LinearVectorData[#, atoms] & /@ expressions,
   "coefficients",
   {}
   ];


ds016DirectionExpressions[rows_List, atoms_List] := Expand[# . atoms] & /@ rows;


(* Lookup[{},key,default] 会把空表误解为规则集合并返回 Missing；线性审计中的空数据必须是零行集。 *)
ds016DataColumn[data_List, key_String] := If[data === {}, {}, Lookup[data, key, {}]];


(* ::Chapter:: *)
(*多重图圈数与 bridge 分类*)

(* 连通分量只使用不同端点间的边；自环不改变连通性，但在 E-V+C 中自然贡献一圈。
   两顶点间的一条普通边给 L=1-2+1=0；真正的 tadpole 必须显式写成 endpoints->{v,v}，给 L=1-1+1=1。 *)
ds016ComponentCount[vertexIds_List, endpoints_List] := Module[{edges, graph},
   edges = (UndirectedEdge @@ # &) /@ Select[endpoints, Length[#] === 2 && ! SameQ @@ # &];
   graph = Graph[vertexIds, edges];
   Length[ConnectedComponents[graph]]
   ];


ds016TopologyGraphAudit[vertexIds_List, lines_List] := Module[
   {endpoints, activeLineIndices, activeEndpoints, malformed, unknown, componentCount, loopCount,
    bridgeLines, cycleLines, incidence},
   endpoints = Lookup[lines, "endpoints", Missing["endpoints"]];
   activeLineIndices = Select[
     Range[Length[lines]],
     Lookup[lines[[#]], "state", "full"] =!= "shrunk" && Lookup[lines[[#]], "packType", Automatic] =!= "shrunk" &
     ];
   activeEndpoints = If[activeLineIndices === {}, {}, endpoints[[activeLineIndices]]];
   malformed = Flatten@Position[endpoints, item_ /; ! ListQ[item] || Length[item] =!= 2, {1}, Heads -> False];
   unknown = If[malformed === {}, Complement[DeleteDuplicates@Flatten[endpoints], vertexIds], {}];
   If[malformed =!= {} || unknown =!= {},
    Return[<|
      "status" -> "invalid",
      "issues" -> DeleteCases[{
         If[malformed === {}, Nothing, <|"code" -> "malformedEndpoints", "lineIndices" -> malformed|>],
         If[unknown === {}, Nothing, <|"code" -> "unknownEndpointVertices", "vertices" -> unknown|>]
         }, Nothing]
      |>]
    ];
   componentCount = ds016ComponentCount[vertexIds, activeEndpoints];
   loopCount = Length[activeLineIndices] - Length[vertexIds] + componentCount;
   bridgeLines = Select[
     activeLineIndices,
     ! SameQ @@ endpoints[[#]] &&
       ds016ComponentCount[vertexIds, endpoints[[DeleteCases[activeLineIndices, #]]]] > componentCount &
     ];
   cycleLines = Complement[activeLineIndices, bridgeLines];
   incidence = Table[
     Which[
      ! MemberQ[activeLineIndices, e], 0,
      SameQ @@ endpoints[[e]], 0,
      vertexIds[[v]] === endpoints[[e, 1]], 1,
      vertexIds[[v]] === endpoints[[e, 2]], -1,
      True, 0
      ],
     {v, Length[vertexIds]}, {e, Length[lines]}
     ];
   <|
    "status" -> "valid",
    "vertexCount" -> Length[vertexIds],
    "inputLineCount" -> Length[lines],
    "internalLineCount" -> Length[activeLineIndices],
    "activeLineIndices" -> activeLineIndices,
    "shrunkLineIndices" -> Complement[Range[Length[lines]], activeLineIndices],
    "connectedComponentCount" -> componentCount,
    "graphLoopCount" -> loopCount,
    "bridgeLineIndices" -> bridgeLines,
    "cycleLineIndices" -> cycleLines,
    "selfLoopLineIndices" -> Select[activeLineIndices, SameQ @@ endpoints[[#]] &],
    "incidenceMatrix" -> incidence,
    "cycleSpaceDimension" -> Length[activeLineIndices] - MatrixRank[incidence],
    "issues" -> {}
    |>
   ];


ds016ResolveIBPMode[case_Association, graphLoopCount_Integer] := Module[{requested},
   requested = Lookup[case, "ibpMode", Automatic];
   Which[
    requested === Automatic && graphLoopCount === 0, "timeOnly",
    requested === Automatic, "full",
    MemberQ[{"full", "timeOnly"}, requested], requested,
    True, "invalid"
    ]
   ];


(* ::Chapter:: *)
(*圈动量 routing 与 affine shift 商空间*)

ds016RawISPExpressions[case_Association] := DeleteCases[
   Map[
    Which[
      AssociationQ[#], Lookup[#, "expr", Nothing],
      ListQ[#] && Length[#] >= 2, #[[2]],
      True, Nothing
      ] &,
    Lookup[case, "ispData", {}]
    ],
   Nothing
   ];


ds016ArgumentRoutingData[argument_, loopMomenta_List, referenceMatrix_, referenceResiduals_] := Module[
   {coefficients, residual, transformedResidual, inverse},
   coefficients = Coefficient[Expand[argument], #] & /@ loopMomenta;
   residual = Expand[argument - coefficients . loopMomenta];
   transformedResidual = residual;
   If[loopMomenta =!= {} && MatrixQ[referenceMatrix] && Length[referenceMatrix] === Length[loopMomenta],
    inverse = Inverse[referenceMatrix];
    transformedResidual = Expand[residual - coefficients . inverse . referenceResiduals]
    ];
   <|
    "argument" -> argument,
    "loopCoefficients" -> coefficients,
    "externalResidual" -> residual,
    "shiftInvariantResidual" -> transformedResidual,
    "containsLoopMomentumQ" -> AnyTrue[coefficients, ! ds016ZeroQ[#] &],
    "linearInLoopMomentaQ" -> FreeQ[residual, Alternatives @@ loopMomenta]
    |>
   ];


ds016ISPShiftInvariantDirections[case_Association, loopMomenta_List, referenceMatrix_, referenceResiduals_] := Module[
   {pairs, pairData},
   pairs = Cases[
     ds016RawISPExpressions[case],
     HoldPattern[sp[left_, right_]] :> {left, right},
     {0, Infinity}
     ];
   DeleteDuplicates@Flatten[
     Map[
      Function[pair,
       pairData = ds016ArgumentRoutingData[#, loopMomenta, referenceMatrix, referenceResiduals] & /@ pair;
       If[AnyTrue[pairData, TrueQ[Lookup[#, "containsLoopMomentumQ", False]] &],
        Select[Lookup[pairData, "shiftInvariantResidual", {}], ! ds016ZeroQ[#] &],
        {}
        ]
       ],
      pairs
      ],
     1
     ]
   ];


ds016LoopRoutingAudit[case_Association, lines_List, graphAudit_Association] := Module[
   {mode, loopMomenta, lineMomenta, matrix, activeLineIndices, activeMatrix, residuals, rank, incidenceResidual, bridgeResidual,
    referencePositions = {}, referenceMatrix = {}, referenceResiduals = {}, shiftResiduals,
    lineLinearQ, routingCoefficientQ, fullChecksQ, issues = {}, graphLoopCount},
   graphLoopCount = Lookup[graphAudit, "graphLoopCount", 0];
   mode = ds016ResolveIBPMode[case, graphLoopCount];
   loopMomenta = Lookup[case, "loopMomenta", {}];
   lineMomenta = Lookup[lines, "momentum", 0];
   If[mode === "invalid", AppendTo[issues, <|"severity" -> "error", "code" -> "invalidIBPMode", "value" -> Lookup[case, "ibpMode", Automatic]|>]];
   matrix = Table[Coefficient[Expand[lineMomenta[[e]]], loopMomenta[[l]]], {e, Length[lines]}, {l, Length[loopMomenta]}];
   activeLineIndices = Lookup[graphAudit, "activeLineIndices", Range[Length[lines]]];
   activeMatrix = If[activeLineIndices === {}, {}, matrix[[activeLineIndices]]];
   residuals = Table[Expand[lineMomenta[[e]] - matrix[[e]] . loopMomenta], {e, Length[lines]}];
   lineLinearQ = If[loopMomenta === {}, ConstantArray[True, Length[lines]], FreeQ[#, Alternatives @@ loopMomenta] & /@ residuals];
   routingCoefficientQ = And @@ (MemberQ[{-1, 0, 1}, #] & /@ Flatten[activeMatrix]);
   rank = ds016MatrixRank[activeMatrix, Length[loopMomenta]];
   (* 每条传播子动量可整体反号，endpoints 顺序不固定其代数方向。
      因而用 GF(2) cycle support 检查流守恒，再用有理秩检查独立圈数。 *)
   incidenceResidual = If[
     Length[loopMomenta] === 0 || ! routingCoefficientQ,
     {},
     Mod[Lookup[graphAudit, "incidenceMatrix", {}] . Mod[matrix, 2], 2]
     ];
   bridgeResidual = If[
     Length[loopMomenta] === 0 || Lookup[graphAudit, "bridgeLineIndices", {}] === {},
     {},
     matrix[[Lookup[graphAudit, "bridgeLineIndices", {}]]]
     ];
   fullChecksQ = mode === "full";
   If[fullChecksQ && Length[loopMomenta] =!= graphLoopCount,
    AppendTo[issues, <|"severity" -> "error", "code" -> "loopMomentumCountMismatch", "expected" -> graphLoopCount, "actual" -> Length[loopMomenta]|>]
    ];
   If[fullChecksQ && rank =!= graphLoopCount,
    AppendTo[issues, <|"severity" -> "error", "code" -> "loopRoutingRankMismatch", "expected" -> graphLoopCount, "actual" -> rank|>]
    ];
   If[fullChecksQ && ! routingCoefficientQ,
    AppendTo[issues, <|"severity" -> "error", "code" -> "unsupportedLoopRoutingCoefficients", "allowed" -> {-1, 0, 1}, "matrix" -> matrix|>]
    ];
   If[fullChecksQ && incidenceResidual =!= {} && ! And @@ (ds016ZeroQ /@ Flatten[incidenceResidual]),
    AppendTo[issues, <|"severity" -> "error", "code" -> "loopRoutingOutsideCycleSpace", "residual" -> incidenceResidual|>]
    ];
   If[fullChecksQ && bridgeResidual =!= {} && ! And @@ (ds016ZeroQ /@ Flatten[bridgeResidual]),
    AppendTo[issues, <|"severity" -> "error", "code" -> "bridgeCarriesLoopMomentum", "lineIndices" -> Lookup[graphAudit, "bridgeLineIndices", {}], "coefficients" -> bridgeResidual|>]
    ];
   If[fullChecksQ && ! And @@ lineLinearQ[[activeLineIndices]],
    AppendTo[issues, <|"severity" -> "error", "code" -> "nonlinearLoopMomentumRouting", "lineIndices" -> Select[activeLineIndices, ! TrueQ[lineLinearQ[[#]]] &]|>]
    ];
   If[fullChecksQ && rank === Length[loopMomenta] && Length[loopMomenta] > 0,
    referencePositions = activeLineIndices[[Take[ds016IndependentRowPositions[activeMatrix], UpTo[Length[loopMomenta]]]]];
    referenceMatrix = matrix[[referencePositions]];
    referenceResiduals = residuals[[referencePositions]];
    shiftResiduals = Expand /@ (residuals - matrix . Inverse[referenceMatrix] . referenceResiduals),
    shiftResiduals = residuals
    ];
   <|
    "status" -> If[AnyTrue[issues, Lookup[#, "severity", ""] === "error" &], "invalid", "valid"],
    "ibpMode" -> mode,
    "loopMomenta" -> loopMomenta,
    "loopCoefficientMatrix" -> matrix,
    "loopCoefficientRank" -> rank,
    "lineExternalResiduals" -> residuals,
    "referenceLineIndices" -> referencePositions,
    "referenceLoopMatrix" -> referenceMatrix,
    "referenceExternalResiduals" -> referenceResiduals,
    "shiftInvariantLineResiduals" -> shiftResiduals,
    "incidenceCycleResidual" -> incidenceResidual,
    "issues" -> issues
    |>
   ];


(* ::Chapter:: *)
(*用户 loop 外动量列表完备性*)

ds016SpanAudit[requiredExpressions_List, userExpressions_List, atoms_List] := Module[
   {requiredData, userData, requiredRows, userRows, requiredBasis, userBasis,
    missingPositions, extraPositions, missingRows, extraRows, dependencies,
    invalidRequired, invalidUser, requiredRank, userRank, unionRank, status},
   requiredData = ds016LinearVectorData[#, atoms] & /@ requiredExpressions;
   userData = ds016LinearVectorData[#, atoms] & /@ userExpressions;
   invalidRequired = Flatten@Position[ds016DataColumn[requiredData, "linearQ"], False];
   invalidUser = Flatten@Position[ds016DataColumn[userData, "linearQ"], False];
   requiredRows = ds016DataColumn[requiredData, "coefficients"];
   userRows = ds016DataColumn[userData, "coefficients"];
   requiredBasis = ds016RowBasis[requiredRows, Length[atoms]];
   userBasis = ds016RowBasis[userRows, Length[atoms]];
   requiredRank = Length[requiredBasis];
   userRank = Length[userBasis];
   unionRank = ds016MatrixRank[Join[requiredBasis, userBasis], Length[atoms]];
   missingPositions = ds016IndependentRowPositions[requiredBasis, userBasis];
   extraPositions = ds016IndependentRowPositions[userBasis, requiredBasis];
   missingRows = If[missingPositions === {}, {}, requiredBasis[[missingPositions]]];
   extraRows = If[extraPositions === {}, {}, userBasis[[extraPositions]]];
   dependencies = If[userRows === {} || Length[userRows] <= userRank, {}, NullSpace[Transpose[userRows]]];
   status = Which[
     invalidRequired =!= {} || invalidUser =!= {}, "invalid",
     missingRows =!= {}, "undercomplete",
     extraRows =!= {} || Length[userExpressions] > userRank, "overcomplete",
     True, "exact"
     ];
   <|
    "status" -> status,
    "atoms" -> atoms,
    "requiredExpressions" -> requiredExpressions,
    "userExpressions" -> userExpressions,
    "requiredBasisDirections" -> ds016DirectionExpressions[requiredBasis, atoms],
    "userBasisDirections" -> ds016DirectionExpressions[userBasis, atoms],
    "missingDirections" -> ds016DirectionExpressions[missingRows, atoms],
    "extraDirections" -> ds016DirectionExpressions[extraRows, atoms],
    "userDependencyVectors" -> dependencies,
    "requiredRank" -> requiredRank,
    "userRank" -> userRank,
    "unionRank" -> unionRank,
    "invalidRequiredPositions" -> invalidRequired,
    "invalidUserPositions" -> invalidUser
    |>
   ];


ds016RequiredLoopExternalDirections[case_Association, graphAudit_Association, routingAudit_Association] := Module[
   {cycleLines, lineDirections, ispDirections, mode},
   mode = Lookup[routingAudit, "ibpMode", "invalid"];
   If[mode =!= "full", Return[{}]];
   cycleLines = Lookup[graphAudit, "cycleLineIndices", {}];
   lineDirections = Lookup[routingAudit, "shiftInvariantLineResiduals", {}];
   lineDirections = If[cycleLines === {}, {}, lineDirections[[cycleLines]]];
   ispDirections = ds016ISPShiftInvariantDirections[
     case,
     Lookup[routingAudit, "loopMomenta", {}],
     Lookup[routingAudit, "referenceLoopMatrix", {}],
     Lookup[routingAudit, "referenceExternalResiduals", {}]
     ];
   DeleteDuplicates@Select[Join[lineDirections, ispDirections], ! ds016ZeroQ[#] &]
   ];


(* ::Chapter:: *)
(*独立无圈动量模长完备性*)

ds016CanonicalMomentumSign[expr_] := Module[{expanded = Expand[expr], opposite},
   opposite = Expand[-expanded];
   First@SortBy[{expanded, opposite}, ToString[InputForm[#]] &]
   ];


ds016MagnitudeMomentaInExpression[expr_] := DeleteDuplicates@Cases[
    expr,
    HoldPattern[Power[sp[left_, right_], Rational[1, 2]]] /; Expand[left - right] === 0 :> ds016CanonicalMomentumSign[left],
    {0, Infinity}
    ];


ds016RequiredIndependentMomentumMagnitudes[case_Association, lines_List, graphAudit_Association, routingAudit_Association] := Module[
   {candidateLines, lineCandidates, extLegCandidates, phaseCandidates, loopMomenta, mode, zeroLoopQ},
   loopMomenta = Lookup[routingAudit, "loopMomenta", {}];
   mode = Lookup[routingAudit, "ibpMode", "invalid"];
   zeroLoopQ[expr_] := And @@ (ds016ZeroQ[Coefficient[Expand[expr], #]] & /@ loopMomenta);
   (* full 模式只把 graph bridge 当作独立模长；timeOnly 中所有 active line 都退出
      loop-IBP 表示，因此它们的模长必须由用户的独立外动量列表覆盖。 *)
   candidateLines = If[
     mode === "timeOnly",
     Lookup[graphAudit, "activeLineIndices", {}],
     Lookup[graphAudit, "bridgeLineIndices", {}]
     ];
   lineCandidates = If[candidateLines === {}, {}, Lookup[lines[[candidateLines]], "momentum", {}]];
   extLegCandidates = Cases[Lookup[case, "extLegs", {}], entry_List /; Length[entry] >= 3 :> entry[[3]]];
   phaseCandidates = Flatten[
     ds016MagnitudeMomentaInExpression /@ Lookup[Lookup[case, "vertices", {}], "externalLegEnergy", {}]
     ];
   DeleteDuplicates@Select[
     ds016CanonicalMomentumSign /@ Join[lineCandidates, extLegCandidates, phaseCandidates],
     ! ds016ZeroQ[#] && (mode === "timeOnly" || zeroLoopQ[#]) &
     ]
   ];


ds016GramPairs[count_Integer] := Flatten[Table[{i, j}, {i, count}, {j, i, count}], 1];


ds016BilinearGramRow[left_List, right_List] := Map[
   Function[pair,
    If[pair[[1]] === pair[[2]],
     left[[pair[[1]]]] right[[pair[[1]]]],
     left[[pair[[1]]]] right[[pair[[2]]]] + left[[pair[[2]]]] right[[pair[[1]]]]
     ]
    ],
   ds016GramPairs[Length[left]]
   ];


ds016SquaredGramRow[row_List] := ds016BilinearGramRow[row, row];


ds016LoopGramRows[loopRows_List] := Flatten[
   Table[ds016BilinearGramRow[loopRows[[i]], loopRows[[j]]], {i, Length[loopRows]}, {j, i, Length[loopRows]}],
   1
   ];


ds016QuadraticSpanAudit[requiredMomenta_List, userMomenta_List, loopMomenta_List, atoms_List] := Module[
   {requiredData, userData, loopData, requiredRows, userRows, loopRows, baseRows,
    missingPositions, extraPositions, missingRows, extraRows, userNewPositions, redundantUserPositions,
    invalidRequired, invalidUser, invalidLoop, requiredNewRank, userNewRank, quadraticDependencies, status},
   requiredData = ds016LinearVectorData[#, atoms] & /@ requiredMomenta;
   userData = ds016LinearVectorData[#, atoms] & /@ userMomenta;
   loopData = ds016LinearVectorData[#, atoms] & /@ loopMomenta;
   invalidRequired = Flatten@Position[ds016DataColumn[requiredData, "linearQ"], False];
   invalidUser = Flatten@Position[ds016DataColumn[userData, "linearQ"], False];
   invalidLoop = Flatten@Position[ds016DataColumn[loopData, "linearQ"], False];
   requiredRows = ds016SquaredGramRow /@ ds016DataColumn[requiredData, "coefficients"];
   userRows = ds016SquaredGramRow /@ ds016DataColumn[userData, "coefficients"];
   loopRows = ds016DataColumn[loopData, "coefficients"];
   baseRows = ds016LoopGramRows[loopRows];
   missingPositions = ds016IndependentRowPositions[requiredRows, Join[baseRows, userRows]];
   extraPositions = ds016IndependentRowPositions[userRows, Join[baseRows, requiredRows]];
   userNewPositions = ds016IndependentRowPositions[userRows, baseRows];
   missingRows = If[missingPositions === {}, {}, requiredRows[[missingPositions]]];
   extraRows = If[extraPositions === {}, {}, userRows[[extraPositions]]];
   requiredNewRank = Length[ds016IndependentRowPositions[requiredRows, baseRows]];
   userNewRank = Length[userNewPositions];
   redundantUserPositions = Complement[Range[Length[userRows]], userNewPositions];
   quadraticDependencies = If[
     Join[baseRows, userRows] === {},
     {},
     NullSpace[Transpose[Join[baseRows, userRows]]]
     ];
   status = Which[
     invalidRequired =!= {} || invalidUser =!= {} || invalidLoop =!= {}, "invalid",
     missingRows =!= {}, "undercomplete",
     extraRows =!= {} || Length[userMomenta] > userNewRank, "overcomplete",
     True, "exact"
     ];
   <|
    "status" -> status,
    "atoms" -> atoms,
    "loopGramRank" -> ds016MatrixRank[baseRows, Length[ds016GramPairs[Length[atoms]]]],
    "requiredMomenta" -> requiredMomenta,
    "userMomenta" -> userMomenta,
    "missingMagnitudeSquares" -> If[missingPositions === {}, {}, sp[#, #] & /@ requiredMomenta[[missingPositions]]],
    "extraMagnitudeSquares" -> If[extraPositions === {}, {}, sp[#, #] & /@ userMomenta[[extraPositions]]],
    "redundantUserPositions" -> redundantUserPositions,
    "redundantUserMomenta" -> If[redundantUserPositions === {}, {}, userMomenta[[redundantUserPositions]]],
    "quadraticDependencyOrder" -> Join[
      Table["loopGram" <> ToString[i], {i, Length[baseRows]}],
      Table["userMagnitude" <> ToString[i], {i, Length[userRows]}]
      ],
    "quadraticDependencies" -> quadraticDependencies,
    "requiredIndependentMagnitudeCount" -> requiredNewRank,
    "userIndependentMagnitudeCount" -> userNewRank,
    "invalidRequiredPositions" -> invalidRequired,
    "invalidUserPositions" -> invalidUser,
    "invalidLoopPositions" -> invalidLoop,
    "missingQuadraticRows" -> missingRows,
    "extraQuadraticRows" -> extraRows
    |>
   ];


(* ::Chapter:: *)
(*统一声明审计与 capability gate*)

ds016MomentumDeclarationAudit[case_Association, lines_List, graphAudit_Association, routingAudit_Association] := Module[
   {loopExternal, independentExternal, requiredLoop, requiredIndependent, atoms,
    loopAudit, independentAudit, mode, status, capabilities, issues = {}},
   loopExternal = Lookup[case, "loopExternalMomenta", Lookup[case, "effectiveLoopExternalMomenta", {}]];
   independentExternal = Lookup[case, "independentExternalMomenta", Lookup[case, "independentExternalMomenta", {}]];
   requiredLoop = ds016RequiredLoopExternalDirections[case, graphAudit, routingAudit];
   requiredIndependent = ds016RequiredIndependentMomentumMagnitudes[case, lines, graphAudit, routingAudit];
   atoms = ds016MomentumAtoms[
     Join[requiredLoop, loopExternal, requiredIndependent, independentExternal],
     Lookup[routingAudit, "loopMomenta", {}]
     ];
   mode = Lookup[routingAudit, "ibpMode", "invalid"];
   loopAudit = If[mode === "full",
     ds016SpanAudit[requiredLoop, loopExternal, atoms],
     <|"status" -> "notRequired", "requiredExpressions" -> {}, "userExpressions" -> loopExternal,
       "missingDirections" -> {}, "extraDirections" -> {}, "requiredRank" -> 0, "userRank" -> 0|>
     ];
   independentAudit = ds016QuadraticSpanAudit[requiredIndependent, independentExternal, loopExternal, atoms];
   status = Which[
     Lookup[graphAudit, "status", "invalid"] =!= "valid" || Lookup[routingAudit, "status", "invalid"] =!= "valid", "invalid",
     MemberQ[Lookup[{loopAudit, independentAudit}, "status"], "invalid"], "invalid",
     MemberQ[Lookup[{loopAudit, independentAudit}, "status"], "undercomplete"], "undercomplete",
     MemberQ[Lookup[{loopAudit, independentAudit}, "status"], "overcomplete"], "overcomplete",
     True, "exact"
     ];
   If[Lookup[loopAudit, "status", ""] === "undercomplete",
    AppendTo[issues, <|"severity" -> "error", "code" -> "undercompleteLoopExternalMomenta", "missingDirections" -> Lookup[loopAudit, "missingDirections", {}]|>]
    ];
   If[Lookup[independentAudit, "status", ""] === "undercomplete",
    AppendTo[issues, <|"severity" -> "error", "code" -> "undercompleteIndependentExternalMomenta", "missingMagnitudeSquares" -> Lookup[independentAudit, "missingMagnitudeSquares", {}]|>]
    ];
   If[Lookup[loopAudit, "status", ""] === "overcomplete",
    AppendTo[issues, <|"severity" -> "warning", "code" -> "overcompleteLoopExternalMomenta", "extraDirections" -> Lookup[loopAudit, "extraDirections", {}], "dependencies" -> Lookup[loopAudit, "userDependencyVectors", {}]|>]
    ];
   If[Lookup[independentAudit, "status", ""] === "overcomplete",
    AppendTo[issues, <|
      "severity" -> "warning",
      "code" -> "overcompleteIndependentExternalMomenta",
      "extraMagnitudeSquares" -> Lookup[independentAudit, "extraMagnitudeSquares", {}],
      "redundantUserMomenta" -> Lookup[independentAudit, "redundantUserMomenta", {}],
      "quadraticDependencyOrder" -> Lookup[independentAudit, "quadraticDependencyOrder", {}],
      "quadraticDependencies" -> Lookup[independentAudit, "quadraticDependencies", {}]
      |>]
    ];
   capabilities = <|
     "initializationUsableQ" -> MemberQ[{"exact", "overcomplete"}, status],
     "timeIBPUsableQ" -> MemberQ[{"exact", "overcomplete"}, status],
     "momentumIBPUsableQ" -> TrueQ[mode === "full" && MemberQ[{"exact", "overcomplete"}, status]],
     "derivativeUsableQ" -> TrueQ[status === "exact"],
     "inverseKinematicsUsableQ" -> TrueQ[status === "exact"],
     "backendExportUsableQ" -> TrueQ[status === "exact"]
     |>;
   <|
    "status" -> status,
    "ibpMode" -> mode,
    "loopExternalMomenta" -> loopExternal,
    "independentExternalMomenta" -> independentExternal,
    "requiredLoopExternalDirections" -> requiredLoop,
    "requiredIndependentMomentumMagnitudes" -> requiredIndependent,
    "momentumAtoms" -> atoms,
    "loopExternalAudit" -> loopAudit,
    "independentExternalAudit" -> independentAudit,
    "capabilities" -> capabilities,
    "issues" -> Join[Lookup[routingAudit, "issues", {}], issues]
    |>
   ];


ds016TopologyAndMomentumAudit[case_Association, lines_List, vertexIds_List] := Module[
   {graphAudit, routingAudit, declarationAudit, activeVertexIds},
   activeVertexIds = Lookup[case, "activeVertexIds", vertexIds];
   graphAudit = ds016TopologyGraphAudit[activeVertexIds, lines];
   If[Lookup[graphAudit, "status", "invalid"] =!= "valid",
    Return[<|"status" -> "invalid", "graph" -> graphAudit, "routing" -> <||>, "declarations" -> <||>|>]
    ];
   routingAudit = ds016LoopRoutingAudit[case, lines, graphAudit];
   declarationAudit = ds016MomentumDeclarationAudit[case, lines, graphAudit, routingAudit];
   <|
    "status" -> Lookup[declarationAudit, "status", "invalid"],
    "graph" -> graphAudit,
    "routing" -> routingAudit,
    "declarations" -> declarationAudit,
    "capabilities" -> Lookup[declarationAudit, "capabilities", <||>],
    "issues" -> Join[Lookup[graphAudit, "issues", {}], Lookup[declarationAudit, "issues", {}]]
    |>
   ];
