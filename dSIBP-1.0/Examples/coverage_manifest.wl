(* ::Package:: *)
(* example 覆盖清单。键是相对于 examples/ 的文件路径，值是该文件实际演示的
   公开函数。正式检查会同时核对文件中的调用文本和 package 的 DSPublicAPI[]。 *)

(* ::Chapter:: *)
(*公开函数与成品例子的覆盖关系*)

<|
 "version" -> "1.0",
 "coverage" -> <|
   "01_mixed_bubble_workflow.wl" -> {
     "DSInit", "DSSeeds", "DSAllSeeds", "DSGenerateIBP", "DSLinear", "DSReorderIntegrals", "DSKiraPlan"
     },
   "02_function_system_hankel.wl" -> {
     "DSInit", "DSInfo", "DSSeeds"
     },
   "03_single_massive_sunrise/main.wl" -> {
     "DSKinematics", "DSInit", "DSSeeds", "DSAllSeeds", "ds",
     "symmetry", "repSymmetry0"
     },
   "04_pure_massive_bubble_closed_loop/main.wl" -> {
     "DSInit", "DSInfo", "DSSeeds", "DSAllSeeds", "DSSeedGroups",
      "DSSeedGroupMetadata", "DSMetaSeedRange",
     "DSGenerateIBP", "DSLinear", "DSUserMI", "DSKiraPlan", "DSKiraExport",
     "DSKiraImport", "DSDE", "DSScaleCheck"
     },
   (* 同一成品展示原生 time-only J、massive-only 公式适配与 massless PendingRederivation。 *)
   "05_tree_two_vertex_time_ibp/main.wl" -> {
     "DSMessagesOn", "DSInit", "DSSeeds", "DSLinear", "DSTreeSeeds",
     "repIterative", "DSTreeNaiveIBP", "DSTreeNaiveDE", "DSTreeDLogDE"
     },
   "06_mix_bubble_tree/main.wl" -> {
     "DSKinematics", "DSMessagesOn", "DSMessagesOff", "DSMessagesQ", "DSInit",
     "DSInfo", "DSParameterNotation", "DSRedefineParameters", "DSPublicAPI",
     "DSSeeds", "DSAllSeeds",
     "dtau", "dqq", "dqk", "ds", "rep2innerform", "rep2outform",
     "rep2Integrand", "symmetry", "repSymmetry0"
     }
   |>,
 "runtimeSmokeExamples" -> {
   "01_mixed_bubble_workflow.wl",
   "02_function_system_hankel.wl",
   "03_single_massive_sunrise/main.wl",
   "05_tree_two_vertex_time_ibp/main.wl",
   "06_mix_bubble_tree/main.wl"
   },
 "externalWorkflowExample" -> "04_pure_massive_bubble_closed_loop/main.wl",
 "longTermRepresentativeExamples" -> <|
   "singleMassiveSunrise" -> "03_single_massive_sunrise/main.wl",
   "pureMassiveBubble" -> "04_pure_massive_bubble_closed_loop/main.wl",
   "mixBubbleTree" -> "06_mix_bubble_tree/main.wl"
   |>
 |>
