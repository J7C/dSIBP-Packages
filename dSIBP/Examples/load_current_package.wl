(* ::Package:: *)

(***
File: load_current_package.wl
Purpose: Single load entry for every example in this directory. It resolves the dSIBP
   source tree that sits next to Examples/, loads it through the standard Needs entry, and
   exposes the identity values that each example self-checks against.
Interface: Get this file with an explicit UTF-8 encoding, then use the public dSIBP
   context. The caller receives currentVersion, currentPackagePath and currentManualPath.
***)


(* ::Chapter:: *)
(* Resolve the source tree next to Examples *)

exampleLoaderDirectory = DirectoryName[$InputFileName];
packageDeliveryDirectory = DirectoryName[exampleLoaderDirectory];
modulePackageQ = FileExistsQ[FileNameJoin[{packageDeliveryDirectory, "dSIBP.m"}]] &&
   DirectoryQ[FileNameJoin[{packageDeliveryDirectory, "Kernel"}]];
If[! TrueQ[modulePackageQ],
  Print["ERROR: no dSIBP source tree (dSIBP.m plus Kernel/) was found next to this Examples directory."];
  Abort[]
  ];


(* ::Chapter:: *)
(* Load the package and publish its identity *)

PrependTo[$Path, packageDeliveryDirectory];
Needs["dSIBP`"];

currentVersion = ToString[dSIBP`$dSIBPVersion];
currentPackagePath = FileNameJoin[{packageDeliveryDirectory, "dSIBP.m"}];
currentManualPath = FileNameJoin[{packageDeliveryDirectory, "Documentation", "dSIBP_user_manual.pdf"}];
Print["dSIBP ", currentVersion, " loaded from ", packageDeliveryDirectory];
