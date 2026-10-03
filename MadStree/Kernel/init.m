(* ::Package:: *)

(***
File: init.m
Purpose: Standard MadStree loading entry point.
Interface: Add the current version directory to $Path, then call Needs["MadStree`"].
***)

Get[
  FileNameJoin[{DirectoryName[$InputFileName], "MadStree.wl"}],
  CharacterEncoding -> "UTF-8"
];
