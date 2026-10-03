(* ::Package:: *)
(* Single-massive sunrise 的离散对称性与 odd/odd parity convention。顶点交换
   反转每条 full line 的端点指标；两条 massless 平行线交换时，同时交换对应 line
   pack 与成对 ISP 指标。本例只构造 general seeds 和 general 参数微分算符。 *)


(* ::Chapter:: *)
(*Sunrise odd/odd parity*)

sunriseParityConstraints0 = {
   b[e1] + n[e1, 1] + n[e1, 2] -> 1,
   b[e2] + n[e2, 1] + n[e2, 2] +
     b[e3] + n[e3, 1] + n[e3, 2] -> 1
   };


(* ::Chapter:: *)
(*Sunrise 顶点与 massless-line 交换*)

sunriseIntegralQ[J[_, linePacks_, ispList_]] :=
  Length[linePacks] === 3 && Length[ispList] === 2;


sunriseReverseEndpointPack[pack_List] := If[
   Length[pack] === 3,
   {pack[[1]], pack[[3]], pack[[2]]},
   pack
   ];


sunriseSwapVertices[J[aList_, linePacks_, ispList_]] := J[
   Reverse[aList],
   sunriseReverseEndpointPack /@ linePacks,
   ispList
   ];


sunriseSwapMasslessLines[J[aList_, linePacks_, ispList_]] := J[
   aList,
   ReplacePart[linePacks, {2 -> linePacks[[3]], 3 -> linePacks[[2]]}],
   Reverse[ispList]
   ];


sunriseNumericLexGreaterQ[left_List, right_List] := TrueQ[
   VectorQ[Flatten[{left, right}], NumericQ] &&
    left =!= right && OrderedQ[{right, left}]
   ];


sunriseEndpointState[pack_List, slot_Integer] := If[
   Length[pack] === 3,
   pack[[slot + 1]],
   0
   ];


(* 顶点键先无序化两条同类 massless line 的端点态，所以不受 line 交换影响。 *)
sunriseVertexSwapNeededQ[
   J[{a1_, a2_}, {massivePack_, masslessPack2_, masslessPack3_}, _]
   ] := sunriseNumericLexGreaterQ[
   {a1, sunriseEndpointState[massivePack, 1],
    Sort[{sunriseEndpointState[masslessPack2, 1], sunriseEndpointState[masslessPack3, 1]}]},
   {a2, sunriseEndpointState[massivePack, 2],
    Sort[{sunriseEndpointState[masslessPack2, 2], sunriseEndpointState[masslessPack3, 2]}]}
   ];
sunriseVertexSwapNeededQ[_] := False;


(* line 键无序化一条线的两个端点态，所以不受顶点交换影响；ISP 指标随对应线一起比较。 *)
sunriseMasslessLineSwapNeededQ[
   J[_, {_, masslessPack2_, masslessPack3_}, {isp2_, isp3_}]
   ] := sunriseNumericLexGreaterQ[
   Join[{masslessPack2[[1]]}, Sort[Rest[masslessPack2]], {isp2}],
   Join[{masslessPack3[[1]]}, Sort[Rest[masslessPack3]], {isp3}]
   ];
sunriseMasslessLineSwapNeededQ[_] := False;


(* general 连续指标保持符号时不强行排序；本例只用一个显式整数积分检查规则方向，
   不展开连续 seed 点。 *)
sunriseCanonical[int_J] := Module[{result = int},
   If[sunriseVertexSwapNeededQ[result], result = sunriseSwapVertices[result]];
   If[sunriseMasslessLineSwapNeededQ[result], result = sunriseSwapMasslessLines[result]];
   result
   ];


sunriseSymmetryRules0 = {
   HoldPattern[(int_J /; sunriseIntegralQ[int])] :> sunriseCanonical[int]
   };
