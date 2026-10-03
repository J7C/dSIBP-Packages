(* ::Package:: *)

(***
File: MadStree.m
Purpose: Standard version-root entry point for MadStree.
Interface: Add this version directory to $Path, then call Needs["MadStree`"].
***)

Get[
  FileNameJoin[{DirectoryName[$InputFileName], "Kernel", "MadStree.wl"}],
  CharacterEncoding -> "UTF-8"
];
