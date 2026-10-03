(* ::Package:: *)
(* pure massive bubble 的 reference 候选数据；basis 构造、秩审计和映射统一由 DSUserMI 完成。 *)

(* ::Chapter:: *)
(*Reference 候选与 physical convention*)

userMIBasisDirectory = DirectoryName[$InputFileName];
If[! ValueQ[referenceDlogCandidates],
 Get[FileNameJoin[{userMIBasisDirectory, "dlog_basis.wl"}]]
 ];

pureMassiveBubbleUserMIActiveIndices = Range[19];
pureMassiveBubbleUserMINames = "dlog" <> ToString[#] & /@ Range[21];
pureMassiveBubbleUserMIDerivativeVariables = {ss11, P0};

(* DEscaleCheck 的 stored degrees 在第 15--18 项恢复显式 ks 后各增加 1。
   用 ep=(3-dim)/2 消去 reference epsilon，得到 package 的 dim 表示。 *)
pureMassiveBubbleUserMIScalingDegrees = Join[
   ConstantArray[dim - 5, 5],
   ConstantArray[dim - 6, 9],
   ConstantArray[dim - 4, 5]
   ];

pureMassiveBubbleUserMIExpressions = Expand[
   referenceDlogCandidates /. {
      Sqrt[s11] -> ss11,
      P1 -> -P0,
      P2 -> -P0
      }
   ];

If[
 Length[pureMassiveBubbleUserMIExpressions] =!= 21 ||
  Length[pureMassiveBubbleUserMIActiveIndices] =!= 19 ||
  Length[pureMassiveBubbleUserMIScalingDegrees] =!= 19,
 Print[Style["Pure massive bubble userMI input length mismatch.", Red, Bold]];
 Abort[]
 ];

pureMassiveBubbleUserMIExpressions
