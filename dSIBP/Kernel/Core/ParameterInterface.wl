(* ::Package:: *)
(* 本模块提供用户参数 notation 与重定义入口。所有新规则都重新经过 DSKinematics/DSInit，
   因而 seed、ds、DSDE 与序列化 metadata 不会持有彼此不一致的坐标状态。 *)

(* ::Chapter:: *)
(*参数 notation*)

dsParameterNotation[topo_Association] := Module[
   {audit = Lookup[topo, "kinematicCoordinateAudit", <||>], declarationAudit},
   declarationAudit = Lookup[topo, "momentumDeclarationAudit", <||>];
   <|
     "loopExternalMomenta" -> Lookup[topo, "loopExternalMomenta", {}],
     "effectiveLoopExternalMomenta" -> Lookup[topo, "effectiveLoopExternalMomenta", Lookup[topo, "effectiveLoopExternalMomenta", {}]],
     "independentExternalMomenta" -> Lookup[topo, "independentExternalMomenta", {}],
     "kEIndices" -> Range[Length[Lookup[topo, "independentExternalMomenta", {}]]],
     "kEMomenta" -> Lookup[topo, "independentExternalMomenta", {}],
     "kEParameterExpressions" -> independentExternalMagnitudeExpressions018[topo],
    "defaultRules" -> Lookup[audit, "defaultRules", {}],
    "selectedRules" -> Lookup[audit, "selectedRules", {}],
    "selectedUserVariables" -> Lookup[audit, "selectedUserVariables", {}],
    "dependentMagnitudeBindings" -> Lookup[audit, "dependentMagnitudeBindings", {}],
    "requiredLoopExternalDirections" -> Lookup[declarationAudit, "requiredLoopExternalDirections", {}],
    "requiredNoLoopMagnitudeMomenta" -> Lookup[declarationAudit, "requiredIndependentMomentumMagnitudes", {}],
    "requiredMagnitudeCoverage" -> kinematicRequiredMagnitudeCoverage[topo],
    "parameterRedefinitionGuide" -> kinematicParameterRedefinitionGuide[audit],
    "coordinateStatus" -> Lookup[audit, "status", "unknown"],
    "capabilities" -> Lookup[topo, "capabilities", <||>]
    |>
   ];


DSParameterNotation[context_Association] := Module[{resolved = dsResolveContext[context], result, guide, guideText},
   If[Head[resolved] === Missing,
    dsErrorPrint["DSParameterNotation 需要有效的 DSInit context。 DSParameterNotation requires a valid DSInit context."]; Return[$Failed]
    ];
   result = dsParameterNotation[resolved["topology"]];
   guide = Lookup[result, "parameterRedefinitionGuide", <||>];
   guideText = If[
     StringQ[Lookup[guide, "commandExample", None]],
     Lookup[guide, "commandExample", ""],
     Lookup[guide, "defaultBehavior", ""]
     ];
   dsInfoPrint[
    "当前参数 " <> ToString[Lookup[result, "selectedUserVariables", {}], InputForm] <>
     "。可选重定义示例：" <> guideText <>
     ". Current parameters: " <> ToString[Lookup[result, "selectedUserVariables", {}], InputForm] <>
     ". Optional redefinition example: " <> guideText
    ];
   result
   ];


DSParameterNotation[] := Module[{context = dsResolveContext[Automatic]},
   If[Head[context] === Missing,
    dsErrorPrint["请先成功调用 DSInit。 Run DSInit successfully first."]; Return[$Failed]
    ];
   DSParameterNotation[context]
   ];


(* ::Chapter:: *)
(*参数重定义*)

Options[DSRedefineParameters] = {ProgressReporting -> Automatic};


DSRedefineParameters[context_Association, rules_, OptionsPattern[]] := Module[
   {resolved = dsResolveContext[context], input, result, generateDerivativeMetadataQ},
   If[Head[resolved] === Missing,
    dsErrorPrint["DSRedefineParameters 需要有效的 DSInit context。 DSRedefineParameters requires a valid DSInit context."]; Return[$Failed]
    ];
   If[! ListQ[rules] && ! AssociationQ[rules],
    dsErrorPrint["参数重定义规则必须是 Rule 列表或 Association。 Parameter redefinition rules must be a Rule list or an Association."]; Return[$Failed]
    ];
   input = KeyDrop[resolved["input"], {"kinematicRules"}];
   generateDerivativeMetadataQ = AssociationQ[Lookup[resolved, "derivatives", Missing["NotGenerated"]]];
   result = DSInit[
     input,
     KinematicRules -> If[AssociationQ[rules], Normal[rules], rules],
     RegisterAsCurrent -> False,
     WriteInitializationFiles -> False,
     GenerateDerivativeMetadata -> generateDerivativeMetadataQ,
     ProgressReporting -> OptionValue[ProgressReporting]
     ];
   If[Lookup[result, "status", "failed"] =!= "initialized",
    dsErrorPrint["参数重定义未通过完备性或 topology 门禁。 Parameter redefinition failed the completeness or topology gate."]; Return[result]
    ];
   Join[result, <|"parameterRedefinition" -> <|
      "sourceInputHash" -> Lookup[resolved, "inputHash", Missing["inputHash"]],
      "rules" -> If[AssociationQ[rules], Normal[rules], rules]
      |>|>]
   ];


DSRedefineParameters[rules_, OptionsPattern[]] := Module[{current, result},
   current = dsResolveContext[Automatic];
   If[Head[current] === Missing,
    dsErrorPrint["请先成功调用 DSInit。 Run DSInit successfully first."]; Return[$Failed]
    ];
   result = DSRedefineParameters[current, rules, ProgressReporting -> OptionValue[ProgressReporting]];
   If[AssociationQ[result] && Lookup[result, "status", "failed"] === "initialized",
    $dSIBPCurrentContext = result;
    setIBPTopologyContext[result["topology"]]
    ];
   result
   ];
