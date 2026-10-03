(* ::Package:: *)
(* 文件用途：FlintNDE 的标准顶层 Wolfram Language 加载入口。
   公开接口统一由 Mathematica/FlintNDE.wl 定义，本文件不保存第二套实现。 *)

(* ::Chapter:: *)
(*标准程序包入口*)

Get[FileNameJoin[{
  DirectoryName[$InputFileName],
  "Mathematica",
  "FlintNDE.wl"
}], CharacterEncoding -> "UTF-8"];
