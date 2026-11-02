(* ::Package:: *)
(* 本文件为模块源码与成品 examples 提供统一加载入口，两种入口统一核对 VERSION_INDEX.md
   登记的版本基线：运行时版本串必须等于登记的代码内版本串。正式交付要求恰有一个当前程序
   文件，且其文件名版本 token 与当前发布号一致；手册不带版本号，不参与加载校验。
   维护者可用环境变量 DSIBP_PACKAGE_FILE 检查候选字节。 *)

(* ::Chapter:: *)
(*解析唯一交付件*)

exampleLoaderDirectory = DirectoryName[$InputFileName];
packageDeliveryDirectory = DirectoryName[exampleLoaderDirectory];
modulePackageQ = FileExistsQ[FileNameJoin[{packageDeliveryDirectory, "dSIBP.m"}]] &&
   DirectoryQ[FileNameJoin[{packageDeliveryDirectory, "Kernel"}]];
packageOverride = Quiet[Environment["DSIBP_PACKAGE_FILE"]];
packageOverrideQ = StringQ[packageOverride] && StringLength[StringTrim[packageOverride]] > 0;


(* ::Chapter:: *)
(*版本登记基线*)

versionIndexCandidates = Select[
   NestList[DirectoryName, packageDeliveryDirectory, 8],
   FileExistsQ[FileNameJoin[{#, "VERSION_INDEX.md"}]] &
   ];
If[versionIndexCandidates === {},
  Print[Style["ERROR: 未在加载目录上层找到 VERSION_INDEX.md；模块与成品 examples 都必须在带版本索引的目录树内运行。", Red]];
  Abort[]
  ];
versionIndexText = Import[
   FileNameJoin[{First[versionIndexCandidates], "VERSION_INDEX.md"}],
   "Text", CharacterEncoding -> "UTF-8"
   ];
versionIndexValue[label_String] := FirstCase[
   StringCases[versionIndexText, label <> "：`" ~~ value : Except["`"] .. ~~ "`" -> value],
   value_String,
   Missing["VersionIndexEntry"]
   ];
currentReleaseToken = versionIndexValue["当前正式发布"];
runtimeVersionToken = versionIndexValue["代码内版本串"];
If[Head[currentReleaseToken] === Missing || Head[runtimeVersionToken] === Missing,
  Print[Style["ERROR: VERSION_INDEX.md 缺少 当前正式发布 或 代码内版本串 登记项。", Red]];
  Abort[]
  ];


(* ::Chapter:: *)
(*解析唯一交付件*)

If[modulePackageQ && ! packageOverrideQ,
  PrependTo[$Path, packageDeliveryDirectory];
  Needs["dSIBP`"];
  currentVersion = ToString[dSIBP`$dSIBPVersion];
  If[currentVersion =!= runtimeVersionToken,
    Print[Style["ERROR: 模块运行时版本串与 VERSION_INDEX.md 代码内版本串不一致。", Red]];
    Abort[]
    ];
  currentPackagePath = FileNameJoin[{packageDeliveryDirectory, "dSIBP.m"}],
 currentPackageCandidates = If[packageOverrideQ,
  {ExpandFileName[packageOverride]},
  Select[
   FileNames["package_*.wl", packageDeliveryDirectory],
   StringMatchQ[FileBaseName[#], RegularExpression["package_[0-9]{3}(\\.[0-9]+)?"]] &
   ]
  ];
 If[Length[currentPackageCandidates] =!= 1,
  Print[Style["ERROR: package/ 必须恰有一个当前程序文件。", Red]];
  Abort[]
  ];
 currentPackagePath = First[currentPackageCandidates];
 currentVersion = FirstCase[
  StringCases[FileBaseName[currentPackagePath], RegularExpression["package_([0-9]{3}(?:\\.[0-9]+)?)"] -> "$1"],
  value_String,
  Missing["VersionToken"]
  ];
 If[Head[currentVersion] === Missing,
  Print[Style["ERROR: 当前程序文件名缺少版本 token。", Red]];
  Abort[]
  ];
 If[currentVersion =!= currentReleaseToken,
  Print[Style["ERROR: 交付文件名版本 token 与 VERSION_INDEX.md 当前发布号不一致。", Red]];
  Abort[]
  ]
 ];


(* ::Chapter:: *)
(*加载并核对运行时版本*)

Get[currentPackagePath, CharacterEncoding -> "UTF-8"];

If[ToString[dSIBP`$dSIBPVersion] =!= runtimeVersionToken,
  Print[Style["ERROR: 运行时 package 版本与 VERSION_INDEX.md 代码内版本串不一致。", Red]];
  Abort[]
  ];
