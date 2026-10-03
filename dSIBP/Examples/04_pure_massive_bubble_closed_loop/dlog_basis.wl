(* ::Package:: *)
(* pure massive bubble 的 21 个 reference-readable dlog 候选关系。公式来自外部 bubble dlog 参考推导的实际输出，
   先保留 reference merged energy P1+P2 的可读形式；闭环中按 P_ref=-P_pkg 取 P1=P2=-P0。
   前 19 个是外部 bubble DE 参考固定使用的 active basis。formal plan 应加载同目录
   reference_user_mi_basis.wl，再由 package DSUserMI 构造；不能直接把本文件的旧平方坐标表达式当作 physical ks basis。 *)

(* ::Chapter:: *)
(*Reference 临时记号到统一 J*)

referenceGToJ[refG[nList_List, aList_List, bList_List]] :=
   J[aList, {{bList[[1]], nList[[1]], nList[[2]]}, {bList[[2]], nList[[3]], nList[[4]]}}, {}];

referenceR1ToJ[refR1[nList_List, {aPower_}, bList_List]] :=
   J[{aPower}, {{bList[[1]]}, {bList[[2]], nList[[1]], nList[[2]]}}, {}];


(* ::Chapter:: *)
(*21 个候选关系与 19 维 active 子序列*)

referenceDlogCandidates0 = {
   -4 (2 nu refG[{0, 0, 0, 0}, {0, 1}, {0, 2}] + refG[{0, 0, 0, 1}, {0, 2}, {0, 1}] + refG[{0, 0, 0, 1}, {1, 1}, {0, 1}]),
   2 (refG[{0, 0, 0, 0}, {0, 1}, {0, 2}] + 2 refG[{0, 0, 0, 1}, {0, 0}, {0, 3}] - refG[{0, 0, 1, 1}, {0, 1}, {0, 2}] - 2 nu refG[{0, 1, 0, 0}, {0, 0}, {1, 2}] - refG[{0, 1, 0, 1}, {0, 1}, {1, 1}] - refG[{0, 1, 1, 0}, {0, 1}, {1, 1}]),
   2 (refG[{0, 0, 0, 1}, {1, 1}, {0, 1}] + refG[{0, 0, 1, 0}, {0, 2}, {0, 1}] + 2 (1 + nu) refG[{0, 0, 1, 1}, {0, 1}, {0, 2}] - 2 nu refG[{1, 1, 0, 0}, {0, 1}, {0, 2}] - refG[{1, 1, 0, 1}, {0, 2}, {0, 1}] - refG[{1, 1, 0, 1}, {1, 1}, {0, 1}]),
   2 (refG[{0, 1, 1, 0}, {0, 1}, {1, 1}] + 2 (1 + nu) refG[{0, 1, 1, 1}, {0, 0}, {1, 2}] + refG[{1, 0, 1, 0}, {0, 1}, {1, 1}] + refG[{1, 1, 0, 0}, {0, 1}, {0, 2}] + 2 refG[{1, 1, 0, 1}, {0, 0}, {0, 3}] - refG[{1, 1, 1, 1}, {0, 1}, {0, 2}]),
   4 (refG[{1, 1, 0, 1}, {1, 1}, {0, 1}] + refG[{1, 1, 1, 0}, {0, 2}, {0, 1}] + 2 (1 + nu) refG[{1, 1, 1, 1}, {0, 1}, {0, 2}]),
   refG[{0, 0, 0, 0}, {0, 0}, {2, 2}],
   refG[{0, 0, 1, 1}, {0, 0}, {2, 2}],
   refG[{1, 1, 1, 1}, {0, 0}, {2, 2}],
   refG[{0, 1, 0, 0}, {0, 1}, {1, 2}],
   refG[{0, 1, 1, 1}, {0, 1}, {1, 2}],
   refG[{1, 0, 0, 0}, {0, 1}, {1, 2}],
   refG[{1, 0, 1, 1}, {0, 1}, {1, 2}],
   2 refG[{0, 1, 0, 1}, {0, 0}, {1, 3}],
   2 refG[{0, 1, 1, 0}, {0, 0}, {1, 3}],
   Sqrt[s11] refR1[{0, 0}, {0}, {2, 2}],
   Sqrt[s11] refR1[{1, 1}, {0}, {2, 2}],
   4 I (1 + 2 nu) Sqrt[s11] refR1[{0, 1}, {0}, {2, 1}]/(P1 + P2),
   Sqrt[s11] refR1[{0, 1}, {1}, {2, 1}],
   2 (refR1[{0, 0}, {1}, {0, 2}] + 2 refR1[{0, 1}, {0}, {0, 3}] - refR1[{1, 1}, {1}, {0, 2}]),
   2 I ((1 + 2 nu) (nu refR1[{0, 0}, {0}, {0, 2}] + refR1[{0, 1}, {1}, {0, 1}]) + 2 I (P1 + P2) (nu refR1[{0, 0}, {1}, {0, 2}] + refR1[{0, 1}, {2}, {0, 1}]) + (1 + 2 nu) (refR1[{0, 1}, {1}, {0, 1}] + (1 + nu) refR1[{1, 1}, {0}, {0, 2}]))/(P1 + P2),
   -(I/2) (-4 nu (1 + 2 nu) refR1[{0, 0}, {0}, {0, 2}] - 8 (1 + 2 nu) refR1[{0, 1}, {1}, {0, 1}] + 4 I (2 (P1 + P2) refR1[{0, 1}, {2}, {0, 1}] + (1 + nu) (I (1 + 2 nu) refR1[{1, 1}, {0}, {0, 2}] + 2 (P1 + P2) refR1[{1, 1}, {1}, {0, 2}])))/(P1 + P2)
   };

referenceDlogCandidates = Expand[
   referenceDlogCandidates0 /. {
      integral_refG :> referenceGToJ[integral],
      integral_refR1 :> referenceR1ToJ[integral]
      }
   ];
referenceDlogActiveIndices = Range[19];
referenceDlogNames = "dlog" <> ToString[#] & /@ Range[Length[referenceDlogCandidates]];

If[Length[referenceDlogCandidates] =!= 21 || Length[referenceDlogActiveIndices] =!= 19,
 Print[Style["Reference dlog basis length mismatch.", Red, Bold]];
 Abort[]
 ];
