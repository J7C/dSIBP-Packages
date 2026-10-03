(* ::Package:: *)

(***
File: Paths.wl
Purpose: Centralizes calling-directory resolution, absolute-path detection and controlled directory creation, shared by Artifacts.wl and FlintNDE.wl; each consumer keeps its own resolver for the Automatic semantics.
Conventions: This module depends only on System` and references no other MadStree internal symbols.
***)

(* ::Chapter:: *)
(* Calling-directory resolution *)

msRuntimeDirectory[] := Module[{inputFile = $InputFileName, notebookDirectory},
  If[StringQ[inputFile] && inputFile =!= "",
    Return[DirectoryName[ExpandFileName[inputFile]]]
  ];
  notebookDirectory = Quiet@Check[NotebookDirectory[], $Failed];
  If[StringQ[notebookDirectory], ExpandFileName[notebookDirectory], Directory[]]
];


(* ::Section:: *)
(* Absolute-path detection *)

(* Absolute paths are expanded directly; relative paths are always resolved against the calling script directory to avoid depending on the current process directory. *)
msAbsolutePathQ[path_String] := StringMatchQ[
  StringReplace[path, "\\" -> "/"],
  Alternatives[LetterCharacter ~~ ":/" ~~ ___, "/" ~~ ___]
];


(* ::Section:: *)
(* Controlled directory creation *)

(* Returns the directory path on success and $Failed on failure; the caller decides the Failure tag and diagnostic fields. *)
msEnsureDirectory[path_String] := If[
  DirectoryQ[path],
  path,
  Quiet@Check[
    CreateDirectory[path, CreateIntermediateDirectories -> True],
    $Failed
  ]
];


(* ::Section:: *)
(*Windows 安全路径门禁*)

(* 传统 Win32/Wolfram 文件接口在完整路径达到 260 字符时可能失败。该门禁只分类路径过长，
   不尝试创建目录或启动后端，确保它不会被误报为 Python 或文件写入错误。 *)
msRuntimePathLengthFailure[paths_List] := Module[{maximum = 259, strings, overlong, longest},
  If[$OperatingSystem =!= "Windows", Return[None]];
  strings = Select[paths, StringQ];
  overlong = Select[strings, StringLength[#] > maximum &];
  If[overlong === {}, Return[None]];
  longest = First@MaximalBy[overlong, StringLength];
  msFailure["RuntimePathTooLong", <|
    "path" -> longest,
    "pathLength" -> StringLength[longest],
    "safeMaximum" -> maximum,
    "suggestion" ->
      "请通过 MSRuntimeDirectory 指定更短的临时运行目录；该错误发生在 Python 启动之前。"
  |>]
];
