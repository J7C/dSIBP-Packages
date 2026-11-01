(* ::Package:: *)
(* 用户主积分只定义 J 线性空间中的有序坐标，不建立与 J 并行的物理积分表示。 *)

(* ::Chapter:: *)
(*userMI basis 构造与查询*)

DSUserMI::badlinear = "DSUserMI 需要 DSLinear 返回且尚未附加 userMI 的 linearData。";
DSUserMI::badbasis = "userMI basis 无效：`1`。";


(* ::Section::Closed:: *)
(*精确线性坐标*)

dsUserMIFirstNonzeroPosition[row_List] := FirstCase[
   Range[Length[row]],
   index_ /; ! TrueQ[Together[row[[index]]] === 0],
   Missing["NoPivot"]
   ];


dsUserMICoordinateData[expressions_List, activeIndices_List, integralList_List] := Module[
   {support, outsideSupport, matrix, residuals, reduced, pivotColumns, rank,
    activeRank, spectatorColumns, pivotMatrix, spectatorMatrix, tokens,
    pivotIntegrals, spectatorIntegrals, reverseExpressions, forwardRules,
    reverseRules, forwardResiduals, reverseResiduals, payload},
   support = DeleteDuplicates@Cases[expressions, _J, Infinity];
   outsideSupport = Complement[support, integralList];
   If[support === {} || outsideSupport =!= {},
    Return[<|"status" -> "failed", "reason" -> "basisSupportOutsideLinearData",
      "supportIntegrals" -> support, "outsideSupportIntegrals" -> outsideSupport|>]
    ];
   matrix = Table[
     Coefficient[Expand[expressions[[i]]], support[[j]]],
     {i, Length[expressions]}, {j, Length[support]}
     ];
   residuals = MapThread[Together[#1 - #2.#3] &, {
      expressions,
      matrix,
      ConstantArray[support, Length[expressions]]
      }];
   If[! And @@ (TrueQ[# === 0] & /@ residuals),
    Return[<|"status" -> "failed", "reason" -> "basisMustBeHomogeneousLinearInJ",
      "reconstructionResiduals" -> residuals|>]
    ];
   reduced = Quiet@Check[RowReduce[matrix], $Failed];
   If[reduced === $Failed,
    Return[<|"status" -> "failed", "reason" -> "basisRankComputationFailed"|>]
    ];
   pivotColumns = DeleteMissing[dsUserMIFirstNonzeroPosition /@ reduced];
   rank = Length[pivotColumns];
   activeRank = MatrixRank[matrix[[activeIndices]]];
   If[rank =!= Length[expressions] || activeRank =!= Length[activeIndices],
    Return[<|"status" -> "failed", "reason" -> "basisRowsMustBeIndependent",
      "rank" -> rank, "basisCount" -> Length[expressions],
      "activeRank" -> activeRank, "activeCount" -> Length[activeIndices]|>]
    ];
   spectatorColumns = Complement[Range[Length[support]], pivotColumns];
   pivotMatrix = matrix[[All, pivotColumns]];
   spectatorMatrix = matrix[[All, spectatorColumns]];
   tokens = userMI /@ Range[Length[expressions]];
   pivotIntegrals = support[[pivotColumns]];
   spectatorIntegrals = support[[spectatorColumns]];
   reverseExpressions = Together /@ (Inverse[pivotMatrix].(
        tokens - spectatorMatrix.spectatorIntegrals
        ));
   forwardRules = Thread[tokens -> expressions];
   reverseRules = Thread[pivotIntegrals -> reverseExpressions];
   forwardResiduals = Together /@ (tokens - (expressions /. reverseRules));
   reverseResiduals = Together /@ (pivotIntegrals - (reverseExpressions /. forwardRules));
   payload = <|
     "status" -> "configured",
     "count" -> Length[expressions],
     "tokens" -> tokens,
     "expressions" -> expressions,
     "activeIndices" -> activeIndices,
     "activeTokens" -> tokens[[activeIndices]],
     "activeExpressions" -> expressions[[activeIndices]],
     "supportIntegrals" -> support,
     "supportCount" -> Length[support],
     "coefficientMatrix" -> matrix,
     "rank" -> rank,
     "activeRank" -> activeRank,
     "pivotColumns" -> pivotColumns,
     "pivotIntegrals" -> pivotIntegrals,
     "spectatorColumns" -> spectatorColumns,
     "spectatorIntegrals" -> spectatorIntegrals,
     "forwardRules" -> forwardRules,
     "reverseRules" -> reverseRules,
     "forwardRoundTripResiduals" -> forwardResiduals,
     "reverseRoundTripResiduals" -> reverseResiduals,
     "reversibleQ" -> TrueQ[And @@ (TrueQ[# === 0] & /@ Join[forwardResiduals, reverseResiduals])],
     "sourceIntegralOrderDigest" -> dsKiraExpressionDigest[integralList]
     |>;
   Join[payload, <|"mappingDigest" -> dsKiraExpressionDigest[payload]|>]
   ];


(* ::Section:: *)
(*公开构造与查询*)

DSUserMI[linearData_Association, expressions_List, spec_Association : <||>] := Module[
   {names, activeIndices, coordinateData, setting, prepared},
   If[Lookup[linearData, "dSIBPStatus", "failed"] =!= "generated" ||
     ! ListQ[Lookup[linearData, "integralList", Missing["integralList"]]] ||
     Lookup[Lookup[linearData, "activeBasis", <||>], "status", "disabled"] === "configured",
    Message[DSUserMI::badlinear];
    Return[<|"status" -> "failed", "reason" -> "notUnpreparedLinearData"|>]
    ];
   names = dsKiraActiveBasisNames[Lookup[spec, "names", Automatic], Length[expressions]];
   activeIndices = Replace[Lookup[spec, "activeIndices", Automatic], (Automatic | All) -> Range[Length[expressions]]];
   If[names === $Failed || Length[names] =!= Length[expressions] ||
     ! DuplicateFreeQ[names] || ! And @@ (StringQ[#] && # =!= "" & /@ names) ||
     ! ListQ[activeIndices] || activeIndices === {} || ! DuplicateFreeQ[activeIndices] ||
     ! And @@ (IntegerQ[#] && 1 <= # <= Length[expressions] & /@ activeIndices),
    Message[DSUserMI::badbasis, "spec 的 names 或 activeIndices 无效：names 必须是与表达式等长且互不重复的非空字符串，activeIndices 必须是 1 到 basis 元素个数之间互不重复的整数位置；请修正后重新调用 DSUserMI。 The spec names or activeIndices are invalid: names must be unique nonempty strings matching the expression count, and activeIndices must be unique integer positions between 1 and the number of basis elements; correct them and call DSUserMI again."];
    Return[<|"status" -> "failed", "reason" -> "invalidNamesOrActiveIndices"|>]
    ];
   coordinateData = dsUserMICoordinateData[expressions, activeIndices, linearData["integralList"]];
   If[Lookup[coordinateData, "status", "failed"] =!= "configured" ||
     ! TrueQ[Lookup[coordinateData, "reversibleQ", False]],
    Message[DSUserMI::badbasis, dsReasonSentence[Lookup[coordinateData, "reason", "coordinate map is not reversible"],
      "userMI 坐标映射不可逆，无法把 userMI token 换回 J 线性组合；请检查 basis 表达式的独立性与 linearData 积分列表后重试。 The userMI coordinate map is not reversible, so userMI tokens cannot be mapped back to J linear combinations; check the independence of the basis expressions and the linearData integral list, then retry."]];
    Return[coordinateData]
    ];
   setting = <|
     "names" -> names,
     "expressions" -> expressions,
     "activeIndices" -> activeIndices,
     "derivativeVariables" -> Lookup[spec, "derivativeVariables", Automatic],
     "scalingDegrees" -> Lookup[spec, "scalingDegrees", Automatic],
     "userMIData" -> coordinateData
     |>;
   prepared = dsKiraAttachActiveBasis[linearData, setting];
   If[Lookup[prepared, "status", "failed"] =!= "generated",
    Message[DSUserMI::badbasis, dsReasonSentence[Lookup[prepared, "reason", "active-basis preparation failed"],
      "userMI active basis 未能附加到 linearData；请检查 basis 表达式与 linearData 的同源性后重试。 The userMI active basis could not be attached to linearData; check the basis expressions and that they come from the same linearData source, then retry."]]
    ];
   prepared
   ];


DSUserMI[data_Association] := Lookup[
   Lookup[data, "activeBasis", <||>],
   "userMI",
   Missing["UserMINotConfigured"]
   ];


DSUserMI[data_Association, key_String] := Lookup[DSUserMI[data], key, Missing["UnknownUserMIKey", key]];


DSUserMI[_, ___] := (Message[DSUserMI::badbasis, "DSUserMI 调用格式无效：需要 DSLinear 返回的 linearData、按顺序排列的 basis 表达式列表和可选的 spec Association；请按该格式重新调用。 The DSUserMI call is malformed: it requires linearData from DSLinear, an ordered basis-expression list, and an optional spec Association; call it again in that form."]; <|"status" -> "failed", "reason" -> "invalidArguments"|>);
