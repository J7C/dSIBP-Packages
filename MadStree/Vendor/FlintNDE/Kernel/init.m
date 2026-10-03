(* ::Package:: *)
(* 文件用途：支持 Wolfram 标准 Kernel/init.m 发现协议。
   该入口只加载版本根 Mathematica/FlintNDE.wl，不复制任何公开定义。 *)

(* ::Chapter:: *)
(*Kernel 标准入口*)

Get[FileNameJoin[{
  DirectoryName[DirectoryName[$InputFileName]],
  "Mathematica",
  "FlintNDE.wl"
}], CharacterEncoding -> "UTF-8"];
