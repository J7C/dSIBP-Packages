(* ::Package:: *)
(* pure massive bubble 的 19 个 reference-readable dlog 候选关系。公式来自外部 bubble dlog
   参考推导的实际输出：19 项整体构成外部 bubble DE 参考的完整 active basis，没有辅助关系。
   保留 reference 的可读记号（ks=Sqrt[s11]、merged energy P1+P2）；闭环中按 P_ref=-P_pkg 取
   P1=P2=-P0。formal plan 应加载同目录 reference_user_mi_basis.wl，再由 package DSUserMI
   构造；不能直接把本文件的平方坐标表达式当作 physical ks basis。 *)

(* ::Chapter:: *)
(*Reference 临时记号到统一 J*)

referenceGToJ[refG[nList_List, aList_List, bList_List]] :=
   J[aList, {{bList[[1]], nList[[1]], nList[[2]]}, {bList[[2]], nList[[3]], nList[[4]]}}, {}];

referenceR1ToJ[refR1[nList_List, {aPower_}, bList_List]] :=
   J[{aPower}, {{bList[[1]]}, {bList[[2]], nList[[1]], nList[[2]]}}, {}];


(* ::Chapter:: *)
(*19 个候选关系（全部 active）*)

referenceDlogCandidates0 = {
   2 (refG[{0, 0, 0, 0}, {0, 1}, {0, 2}] + 2 refG[{0, 0, 0, 1}, {0, 0}, {0, 3}] - refG[{0, 0, 1, 1}, {0, 1}, {0, 2}] - 2 nu refG[{0, 1, 0, 0}, {0, 0}, {1, 2}] - refG[{0, 1, 0, 1}, {0, 1}, {1, 1}] - refG[{0, 1, 1, 0}, {0, 1}, {1, 1}]),
   3 refG[{0, 0, 0, 1}, {1, 1}, {0, 1}] + 3 refG[{0, 0, 1, 0}, {0, 2}, {0, 1}] + 4 refG[{0, 0, 1, 1}, {0, 1}, {0, 2}] + 4 nu refG[{0, 0, 1, 1}, {0, 1}, {0, 2}] + refG[{0, 1, 0, 0}, {1, 1}, {-1, 2}] - s11 refG[{0, 1, 0, 0}, {1, 1}, {1, 2}] + 2 refG[{0, 1, 1, 0}, {0, 1}, {-1, 3}] - 2 s11 refG[{0, 1, 1, 0}, {0, 1}, {1, 3}] - refG[{0, 1, 1, 1}, {0, 2}, {-1, 2}] + s11 refG[{0, 1, 1, 1}, {0, 2}, {1, 2}] - refG[{0, 1, 1, 1}, {1, 1}, {-1, 2}] + s11 refG[{0, 1, 1, 1}, {1, 1}, {1, 2}] + refG[{1, 0, 0, 0}, {0, 2}, {-1, 2}] - s11 refG[{1, 0, 0, 0}, {0, 2}, {1, 2}] + 2 refG[{1, 0, 0, 1}, {0, 1}, {-1, 3}] - 2 s11 refG[{1, 0, 0, 1}, {0, 1}, {1, 3}] - 4 nu refG[{1, 1, 0, 0}, {0, 1}, {0, 2}] - 3 refG[{1, 1, 0, 1}, {0, 2}, {0, 1}] - 3 refG[{1, 1, 0, 1}, {1, 1}, {0, 1}],
   -8 nu refG[{0, 0, 0, 0}, {0, 1}, {0, 2}] - 6 refG[{0, 0, 0, 1}, {0, 2}, {0, 1}] - 4 refG[{0, 0, 0, 1}, {1, 1}, {0, 1}] - 2 refG[{0, 1, 0, 0}, {0, 2}, {-1, 2}] + 2 s11 refG[{0, 1, 0, 0}, {0, 2}, {1, 2}] - 4 refG[{0, 1, 0, 1}, {0, 1}, {-1, 3}] + 4 s11 refG[{0, 1, 0, 1}, {0, 1}, {1, 3}] + 2 refG[{0, 1, 1, 1}, {1, 1}, {-1, 2}] - 2 s11 refG[{0, 1, 1, 1}, {1, 1}, {1, 2}] + 2 refG[{1, 1, 0, 1}, {1, 1}, {0, 1}],
   2 (refG[{0, 1, 1, 0}, {0, 1}, {1, 1}] + (2 + 2 nu) refG[{0, 1, 1, 1}, {0, 0}, {1, 2}] + refG[{1, 0, 1, 0}, {0, 1}, {1, 1}] + refG[{1, 1, 0, 0}, {0, 1}, {0, 2}] + 2 refG[{1, 1, 0, 1}, {0, 0}, {0, 3}] - refG[{1, 1, 1, 1}, {0, 1}, {0, 2}]),
   -2 refG[{0, 0, 0, 1}, {1, 1}, {0, 1}] - 2 refG[{0, 1, 0, 0}, {1, 1}, {-1, 2}] + 2 s11 refG[{0, 1, 0, 0}, {1, 1}, {1, 2}] - 4 refG[{1, 0, 1, 0}, {0, 1}, {-1, 3}] + 4 s11 refG[{1, 0, 1, 0}, {0, 1}, {1, 3}] + 2 refG[{1, 0, 1, 1}, {0, 2}, {-1, 2}] - 2 s11 refG[{1, 0, 1, 1}, {0, 2}, {1, 2}] + 4 refG[{1, 1, 0, 1}, {1, 1}, {0, 1}] + 6 refG[{1, 1, 1, 0}, {0, 2}, {0, 1}] + 8 refG[{1, 1, 1, 1}, {0, 1}, {0, 2}] + 8 nu refG[{1, 1, 1, 1}, {0, 1}, {0, 2}],
   Sqrt[s11] refG[{0, 0, 0, 0}, {0, 0}, {2, 2}],
   Sqrt[s11] refG[{0, 0, 1, 1}, {0, 0}, {2, 2}],
   Sqrt[s11] refG[{1, 1, 1, 1}, {0, 0}, {2, 2}],
   Sqrt[s11] refG[{0, 1, 0, 0}, {0, 1}, {1, 2}],
   Sqrt[s11] refG[{0, 1, 1, 1}, {0, 1}, {1, 2}],
   Sqrt[s11] refG[{1, 0, 0, 0}, {0, 1}, {1, 2}],
   Sqrt[s11] refG[{1, 0, 1, 1}, {0, 1}, {1, 2}],
   2 Sqrt[s11] refG[{0, 1, 0, 1}, {0, 0}, {1, 3}],
   2 Sqrt[s11] refG[{0, 1, 1, 0}, {0, 0}, {1, 3}],
   refR1[{0, 0}, {0}, {2, 2}],
   refR1[{1, 1}, {0}, {2, 2}],
   4 I (1 + 2 nu) refR1[{0, 1}, {0}, {2, 1}]/(P1 + P2),
   refR1[{0, 1}, {1}, {2, 1}],
   2 (refR1[{0, 0}, {1}, {0, 2}] + 2 refR1[{0, 1}, {0}, {0, 3}] - refR1[{1, 1}, {1}, {0, 2}])/Sqrt[s11]
   };

referenceDlogCandidates = Expand[
   referenceDlogCandidates0 /. {
      integral_refG :> referenceGToJ[integral],
      integral_refR1 :> referenceR1ToJ[integral]
      }
   ];
referenceDlogActiveIndices = Range[19];
referenceDlogNames = "dlog" <> ToString[#] & /@ Range[Length[referenceDlogCandidates]];

If[Length[referenceDlogCandidates] =!= 19 || Length[referenceDlogActiveIndices] =!= 19,
 Print[Style["Reference dlog basis length mismatch.", Red, Bold]];
 Abort[]
 ];
