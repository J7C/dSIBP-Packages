(* ::Package:: *)
(* pure massive bubble 的 reference parity、canonical symmetry 与逐生成元 seed 范围。
   主例和 fresh-session 交叉检查共同加载本文件，避免初始化 metadata 引用未定义 helper。 *)

(* ::Chapter:: *)
(*Reference parity 与 symmetry*)

exampleTopQ[J[_, {pack1_, pack2_}, {}]] := Length[pack1] === 3 && Length[pack2] === 3;
exampleR1Q[J[_, {pack1_, pack2_}, {}]] := Length[pack1] === 1 && Length[pack2] === 3;
exampleR2Q[J[_, {pack1_, pack2_}, {}]] := Length[pack1] === 3 && Length[pack2] === 1;

exampleSwapVertices[J[{a1_, a2_}, {{b1_, n1_, n2_}, {b2_, n3_, n4_}}, {}]] :=
  J[{a2, a1}, {{b1, n2, n1}, {b2, n4, n3}}, {}];
exampleSwapLines[J[{a1_, a2_}, {{b1_, n1_, n2_}, {b2_, n3_, n4_}}, {}]] :=
  J[{a1, a2}, {{b2, n3, n4}, {b1, n1, n2}}, {}];

exampleTopCanonical[int_J] := Module[{result = int},
   If[TrueQ[result[[1, 1]] > result[[1, 2]]], result = exampleSwapVertices[result]];
   If[TrueQ[result[[2, 1, 1]] > result[[2, 2, 1]]], result = exampleSwapLines[result]];
   If[TrueQ[result[[1, 1]] === result[[1, 2]] && result[[2, 1, 2]] > result[[2, 1, 3]]], result = exampleSwapVertices[result]];
   If[TrueQ[result[[2, 1, 1]] === result[[2, 2, 1]] && result[[2, 1, 2]] > result[[2, 2, 2]]], result = exampleSwapLines[result]];
   result
   ];

exampleR1Canonical[J[aList_, {{bS1_}, {b2_, n3_, n4_}}, {}]] := If[
   TrueQ[n3 > n4],
   J[aList, {{bS1}, {b2, n4, n3}}, {}],
   J[aList, {{bS1}, {b2, n3, n4}}, {}]
   ];
exampleR2ToR1[J[aList_, {{b1_, n1_, n2_}, {bS2_}}, {}]] :=
  exampleR1Canonical[J[aList, {{bS2}, {b1, n1, n2}}, {}]];

(* top 要求 n1+n2+b1、n3+n4+b2 为偶数；R1/R2 对应保留相同 reference parity。 *)
exampleParityZeroQ[J[_, {{b1_, n1_, n2_}, {b2_, n3_, n4_}}, {}]] :=
  TrueQ[OddQ[n1 + n2 + b1] || OddQ[n3 + n4 + b2]];
exampleParityZeroQ[J[_, {{bS1_}, {b2_, n3_, n4_}}, {}]] :=
  TrueQ[OddQ[bS1] || OddQ[n3 + n4 + b2]];
exampleParityZeroQ[J[_, {{b1_, n1_, n2_}, {bS2_}}, {}]] :=
  TrueQ[OddQ[bS2] || OddQ[n1 + n2 + b1]];
exampleParityZeroQ[_] := False;

exampleSymmetryRules0 = {
   HoldPattern[(int_J /; exampleParityZeroQ[int])] :> 0,
   HoldPattern[(int_J /; exampleR2Q[int])] :> exampleR2ToR1[int],
   HoldPattern[(int_J /; exampleTopQ[int])] :> exampleTopCanonical[int],
   HoldPattern[(int_J /; exampleR1Q[int])] :> exampleR1Canonical[int]
   };


(* ::Chapter:: *)
(*Reference top 目标包络*)

(* Reference 四组 top seed 的共同外包范围是 a in [-1,4]、b in [-2,5]。
   这里把它作为最终关系的共同目标包络，同时用于 lower sectors；DSGenerateIBP 会按
   每组 shift 反推出更窄的 seed 点域。只有 DSKiraPlan/DE target closure 明确报告
   缺口时才扩张对应边界；旧脚本给 R1 随手放大的 [-4,8] 不作为缺省输入。 *)
referenceTopTargetEnvelope = {
   {a[v1], -1, 4},
   {a[v2], -1, 4},
   {b[1], -2, 5},
   {b[2], -2, 5}
   };

(* DSMetaSeedRange 的声明集合包含 top 与 shrink sector 实际出现的全部连续指标；
   bS[e] 的目标包络在 DSGenerateIBP 中自动继承同一 root line 的 b[e]。 *)
referenceSeedIndices = {a[v1], a[v2], b[1], b[2], bS[1], bS[2]};
