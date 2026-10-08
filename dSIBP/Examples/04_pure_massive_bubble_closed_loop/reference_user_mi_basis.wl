(* ::Package:: *)
(* pure massive bubble 的 reference 候选数据；basis 构造、秩审计和映射统一由 DSUserMI 完成。
   本文件对应当前外部 bubble reference 的 19 项 basis：十九项全部 active，没有辅助关系。 *)

(* ::Chapter:: *)
(*Reference 候选与 physical convention*)

userMIBasisDirectory = DirectoryName[$InputFileName];
If[! ValueQ[referenceDlogCandidates],
 Get[FileNameJoin[{userMIBasisDirectory, "dlog_basis.wl"}]]
 ];

pureMassiveBubbleUserMIActiveIndices = Range[19];
pureMassiveBubbleUserMINames = "dlog" <> ToString[#] & /@ Range[19];
pureMassiveBubbleUserMIDerivativeVariables = {ss11, P0};

(* 新 basis 的 DEscaleCheck 给出统一标度权重 -2-2 ep；用 ep=(3-dim)/2 换到 package 的
   dim 表示即 dim-5，在 parameterProbeRules 点等于 -18/11。旧版按混合次数分组已不再适用。 *)
pureMassiveBubbleUserMIScalingDegrees = ConstantArray[dim - 5, 19];

(* 读入的 reference 记号含 ks=Sqrt[s11] 与 merged energy P1+P2；先替换 1/Sqrt[s11]
   （它解析为 Power[s11, Rational[-1, 2]]，Sqrt 规则匹配不到），再替换 Sqrt[s11] 与 s11=ks^2。 *)
pureMassiveBubbleUserMIExpressions = Expand[
   referenceDlogCandidates /. {
      Power[s11, Rational[-1, 2]] -> 1/ss11,
      Sqrt[s11] -> ss11,
      s11 -> ss11^2,
      P1 -> -P0,
      P2 -> -P0
      }
   ];

If[
 Length[pureMassiveBubbleUserMIExpressions] =!= 19 ||
  Length[pureMassiveBubbleUserMIActiveIndices] =!= 19 ||
  Length[pureMassiveBubbleUserMIScalingDegrees] =!= 19,
 Print[Style["Pure massive bubble userMI input length mismatch.", Red, Bold]];
 Abort[]
 ];

pureMassiveBubbleUserMIExpressions
