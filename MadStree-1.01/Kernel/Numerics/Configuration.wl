(* ::Package:: *)

(***
File: Configuration.wl
Purpose: Centralizes resolution of the MadStree version directory and the relative location of the built-in FlintNDE Python package.
Conventions: $MSFlintNDERelativePath is the only modifiable backend-directory variable; its value is always relative to the current MadStree version directory, and public functions and the adapter must not hardcode external repository paths in their bodies.
***)

(* ::Chapter:: *)
(* Relative path configuration *)

$MadStreePackageDirectory = DirectoryName[$MadStreeKernelDirectory];
$MadStreeCollectionDirectory = DirectoryName[DirectoryName[$MadStreePackageDirectory]];
$MadStreeProjectDirectory = DirectoryName[$MadStreeCollectionDirectory];
$MSFlintNDERelativePath = FileNameJoin[{
  "Vendor", "FlintNDE"
}];


(* ::Section:: *)
(* Path validation and public queries *)

(* Only ordinary relative paths are accepted, so that configuration cannot accidentally escape the user-specified project root. *)
msRelativePathQ[path_String] := ! AnyTrue[{"/", "\\"}, StringStartsQ[path, #] &] &&
  FreeQ[FileNameSplit[path], ".."] && ! StringMatchQ[path, LetterCharacter ~~ ":" ~~ ___];

msResolvedFlintNDEPath[] := ExpandFileName[
  FileNameJoin[{$MadStreePackageDirectory, $MSFlintNDERelativePath}]
];

MSFlintNDEConfiguration[] := Module[{resolved = msResolvedFlintNDEPath[]},
  <|
    "version" -> $MadStreeVersion,
    "versionDirectory" -> $MadStreePackageDirectory,
    "collectionDirectory" -> $MadStreeCollectionDirectory,
    "projectDirectory" -> $MadStreeProjectDirectory,
    "pathBase" -> $MadStreePackageDirectory,
    "relativePath" -> $MSFlintNDERelativePath,
    "resolvedPath" -> resolved,
    "packageFile" -> FileNameJoin[{resolved, "flintnde", "__init__.py"}],
    "availableQ" -> FileExistsQ[FileNameJoin[{resolved, "flintnde", "__init__.py"}]]
  |>
];

(* Validation runs immediately after modification; on failure the original configuration is kept and no half-updated state is produced. *)
MSSetFlintNDERelativePath[path_String] := Module[{old = $MSFlintNDERelativePath, candidate},
  If[! msRelativePathQ[path],
    Return[msFailure["RelativeFlintNDEPathRequired", <|"path" -> path|>]]
  ];
  candidate = ExpandFileName[FileNameJoin[{$MadStreePackageDirectory, path, "flintnde", "__init__.py"}]];
  If[! FileExistsQ[candidate],
    Return[msFailure["FlintNDEPackageNotFound", <|"path" -> path, "packageFile" -> candidate|>]]
  ];
  $MSFlintNDERelativePath = path;
  If[TrueQ[MSFlintNDEConfiguration[]["availableQ"]],
    MSFlintNDEConfiguration[],
    $MSFlintNDERelativePath = old;
    msFailure["FlintNDEConfigurationFailed", <|"path" -> path|>]
  ]
];

MSSetFlintNDERelativePath[other_] := msFailure[
  "RelativeFlintNDEPathRequired",
  <|"path" -> HoldForm[other]|>
];
