(* ::Package:: *)

(* ::Chapter:: *)
(*Kira 导出边界*)

Options[DSKiraExport] = Join[Options[makeKiraExportData], {
   KiraActiveBasis -> None,
   KiraRequireCompleteSystem -> True,
   KiraNumericStage -> "symbolic",
   ProgressReporting -> Automatic
   }];

DSKiraExport::badlinear = "DSKiraExport 需要 DSLinear 返回的 backend-neutral linearData。";
DSKiraExport::failed = "Kira 输入未生成：`1`。";
DSKiraExport::badbasis = "KiraActiveBasis 未通过验证：`1`。";
DSKiraExport::capability = "linearData 未携带通过 DSLinear 的同源能力门禁。";
DSKiraExport::devarrules = "数值/系数规则与微分阶段合同冲突，Kira 导出已拒绝：`1`。";
DSKiraExport::badstage = "KiraNumericStage 只允许 \"symbolic\" 或 \"postDerivative\"，收到 `1`。";


(* ::Section::Closed:: *)
(*Active basis 与导数 target closure*)

(* active basis 只占用 backend ID；用户侧和物理层仍只使用 J，避免建立平行积分 Head。 *)
dsKiraActiveBasisVariables[Automatic, topo_Association] :=
   scalarProductInternalToUser[#, topo] & /@ independentVariableDerivativeVariables[topo];
dsKiraActiveBasisVariables[variables_List, _Association] := variables;
dsKiraActiveBasisVariables[variable_, _Association] := {variable};

dsKiraActiveBasisNames[Automatic, count_Integer] := "active" <> ToString[#] & /@ Range[count];
dsKiraActiveBasisNames[names_List, _Integer] := names;
dsKiraActiveBasisNames[_, _Integer] := $Failed;

dsKiraLinearizeActiveBasisExpression[expr_, activeID_Integer, integralIndex_Association, name_String] := Module[
   {terms, termData, linearPieces, constantTerms, nonlinearTerms, rules},
   terms = linearTerms[Expand[expr]];
   termData = linearTermData[#, integralIndex] & /@ terms;
   linearPieces = Select[termData, Lookup[#, "kind", "missing"] === "linear" &];
   constantTerms = Lookup[Select[termData, Lookup[#, "kind", "missing"] === "constant" &], "term", {}];
   nonlinearTerms = Lookup[Select[termData, Lookup[#, "kind", "missing"] === "nonlinear" &], "term", {}];
   rules = Join[{activeID -> 1}, (First[#] -> -Last[#]) & /@ combineLinearCoefficientRules[linearPieces]];
   <|
    "activeBasisName" -> name,
    "activeBasisID" -> activeID,
    "coefficientRules" -> rules,
    "constantTerm" -> Total[constantTerms],
    "nonlinearTerms" -> nonlinearTerms,
    "linearQ" -> TrueQ[nonlinearTerms === {} && Total[constantTerms] === 0 && linearPieces =!= {}]
    |>
   ];

dsKiraAttachActiveBasis[linearData_Association, Automatic] /;
   Lookup[Lookup[linearData, "activeBasis", <||>], "status", "disabled"] === "configured" := linearData;


(* active-basis 导数只能序列化显式系数与 J。若执行 helper 仍留在结果中，Put/Get 可能触发
   二次求值并改变 artifact 身份，因此必须在 manifest 写出前拒绝。 *)
dsKiraDerivativeInternalHelperResidues[expr_] := DeleteDuplicates[
   Cases[
    HoldComplete[expr],
    HoldPattern[
     (shiftVertexA | shiftLinePower | shiftLineB | shiftLinePackEntry | setLinePackEntry |
       shiftISPIndex | absorbLinearFactor | absorbLinearFactorExpressionTerm |
       applyEOM | applyEOMToIntegral | eomReduceIntegralAt |
       applyMasslessEndpointCanonical | applyMassiveCoincidentCanonical |
       canonicalizeMassiveCoincidentIntegral | canonicalizeMasslessIntegral018 |
       sectorAwareCanonical018 | sectorAwareCanonicalTerm018 |
       sectorAwareSymmetry018 | sectorAwareSymmetryTerm018 |
       applySeedCanonicalWithRules018)[___]
     ],
    Infinity
    ],
   SameTest -> SameQ
   ];


dsKiraAttachActiveBasis[linearData_Association, setting_] /; setting === None || setting === Automatic :=
   Join[linearData, <|"activeBasis" -> <|"status" -> "disabled", "count" -> 0|>|>];

dsKiraAttachActiveBasis[linearData_Association, setting_Association] := Module[
   {expressions, count, names, activeIndices, activeCount, activeExpressions, activeNames, topo, variables, allowedVariables, badVariables, degrees,
     oldCount, oldRules, idShift, shiftedRules, shiftedEquations, integralIndex,
     basisEquations, badEquations, rawDerivatives, derivativeHelperResidues,
     derivativeIntegrals, missingDerivativeIntegrals,
     relationIDs, activeIDs, auxiliaryIDs, derivativeTargetIDs, targetIDs, userMIData, activeData},
   If[Lookup[linearData, "representation", None] === "sectorTaggedJ[vertexPacks]",
    Return[<|"status" -> "invalidActiveBasis", "reason" -> "treeActiveBasisNotSupported",
      "comment" -> "tree Kira IDs already include sector identity; active-basis derivatives require a separate tagged closure."|>]
    ];
   expressions = Lookup[setting, "expressions", Missing["expressions"]];
   If[! ListQ[expressions] || expressions === {},
    Return[<|"status" -> "invalidActiveBasis", "reason" -> "expressionsMustBeNonemptyList", "activeBasisInput" -> setting|>]
    ];
   count = Length[expressions];
   names = dsKiraActiveBasisNames[Lookup[setting, "names", Automatic], count];
   If[names === $Failed || Length[names] =!= count || ! And @@ (StringQ[#] && # =!= "" & /@ names) || ! DuplicateFreeQ[names],
    Return[<|"status" -> "invalidActiveBasis", "reason" -> "namesMustBeUniqueNonemptyStrings", "activeBasisInput" -> setting|>]
    ];
   activeIndices = Replace[Lookup[setting, "activeIndices", Automatic], (Automatic | All) -> Range[count]];
   If[! ListQ[activeIndices] || activeIndices === {} || ! DuplicateFreeQ[activeIndices] ||
     ! And @@ (IntegerQ[#] && 1 <= # <= count & /@ activeIndices),
    Return[<|"status" -> "invalidActiveBasis", "reason" -> "activeIndicesMustBeUniqueValidPositions", "activeIndices" -> activeIndices|>]
    ];
   activeCount = Length[activeIndices];
   activeExpressions = expressions[[activeIndices]];
   activeNames = names[[activeIndices]];
   topo = Lookup[linearData, "topology", Missing["topology"]];
   If[! parsedTopologyQ[topo],
    Return[<|"status" -> "invalidActiveBasis", "reason" -> "missingParsedTopology"|>]
    ];
   variables = dsKiraActiveBasisVariables[Lookup[setting, "derivativeVariables", Automatic], topo];
   allowedVariables = dsKiraActiveBasisVariables[Automatic, topo];
   badVariables = Complement[variables, allowedVariables];
   If[variables === {} || badVariables =!= {},
    Return[<|"status" -> "invalidActiveBasis", "reason" -> "invalidDerivativeVariables", "derivativeVariables" -> variables, "badVariables" -> badVariables, "allowedVariables" -> allowedVariables|>]
    ];
   degrees = Lookup[setting, "scalingDegrees", Automatic];
   If[degrees =!= Automatic && (! ListQ[degrees] || Length[degrees] =!= activeCount),
    Return[<|"status" -> "invalidActiveBasis", "reason" -> "scalingDegreesLengthMismatch", "scalingDegrees" -> degrees|>]
    ];
   oldCount = Lookup[linearData, "integralCount", 0];
   oldRules = Lookup[linearData, "integralRules", {}];
   idShift = AssociationThread[Range[oldCount], Range[oldCount] + count];
   shiftedRules = oldRules /. (integral_J -> id_Integer) :> (integral -> Lookup[idShift, id, Missing["unknownIntegralID", id]]);
   If[Cases[shiftedRules, _Missing, Infinity] =!= {},
    Return[<|"status" -> "invalidActiveBasis", "reason" -> "integralMapShiftFailed"|>]
    ];
   shiftedEquations = reindexLinearEquation[#, idShift] & /@ Lookup[linearData, "linearEquations", {}];
   integralIndex = Association[shiftedRules];
   relationIDs = Range[count];
   activeIDs = relationIDs[[activeIndices]];
   auxiliaryIDs = Complement[relationIDs, activeIDs];
   basisEquations = MapThread[
     dsKiraLinearizeActiveBasisExpression,
     {expressions, relationIDs, ConstantArray[integralIndex, count], names}
     ];
   badEquations = Select[basisEquations, ! TrueQ[Lookup[#, "linearQ", False]] &];
   If[badEquations =!= {},
    Return[<|"status" -> "invalidActiveBasis", "reason" -> "basisExpressionsMustBeHomogeneousLinearCombinationsOfMappedJ", "badBasisEquations" -> badEquations|>]
    ];
   rawDerivatives = Table[
     dsSectorAwareDerivative[
      activeExpressions[[i]],
      variables[[j]],
      <|"topology" -> topo, "sectors" -> Lookup[linearData, "sectorMetadataList", {}]|>
      ],
     {i, activeCount}, {j, Length[variables]}
     ];
   If[! FreeQ[rawDerivatives, $Failed],
    Return[<|"status" -> "invalidActiveBasis", "reason" -> "basisDerivativeFailed", "derivativeVariables" -> variables|>]
    ];
   derivativeHelperResidues = dsKiraDerivativeInternalHelperResidues[rawDerivatives];
   If[derivativeHelperResidues =!= {},
    Return[<|
      "status" -> "invalidActiveBasis",
      "reason" -> "basisDerivativeContainsInternalHelpers",
      "derivativeVariables" -> variables,
      "internalHelperResidueCount" -> Length[derivativeHelperResidues],
      "internalHelperResidues" -> Take[derivativeHelperResidues, UpTo[20]]
      |>]
    ];
   derivativeIntegrals = DeleteDuplicates[Cases[rawDerivatives, _J, Infinity]];
   missingDerivativeIntegrals = Select[derivativeIntegrals, ! KeyExistsQ[integralIndex, #] &];
   If[missingDerivativeIntegrals =!= {},
    Return[<|"status" -> "invalidActiveBasis", "reason" -> "derivativeTargetOutsideLinearSystem", "missingDerivativeIntegrals" -> missingDerivativeIntegrals, "derivativeVariables" -> variables|>]
    ];
   derivativeTargetIDs = Lookup[integralIndex, derivativeIntegrals, {}];
   targetIDs = DeleteDuplicates[Join[activeIDs, derivativeTargetIDs]];
   userMIData = Lookup[setting, "userMIData", None];
   If[AssociationQ[userMIData],
    userMIData = Join[userMIData, <|
       "backendIDs" -> relationIDs,
       "backendTokens" -> (Tuserweight /@ relationIDs),
       "activeBackendIDs" -> activeIDs,
       "activeBackendTokens" -> (Tuserweight /@ activeIDs),
       "userMIToBackendRules" -> Thread[(userMI /@ relationIDs) -> (Tuserweight /@ relationIDs)],
       "backendToUserMIRules" -> Thread[(Tuserweight /@ relationIDs) -> (userMI /@ relationIDs)]
       |>]
    ];
   activeData = <|
     "status" -> "configured",
     "count" -> count,
     "relationCount" -> count,
     "activeCount" -> activeCount,
     "names" -> names,
     "expressions" -> expressions,
     "ids" -> relationIDs,
     "tokens" -> (Tuserweight /@ relationIDs),
     "activeIndices" -> activeIndices,
     "activeNames" -> activeNames,
     "activeExpressions" -> activeExpressions,
     "activeIDs" -> activeIDs,
     "activeTokens" -> (Tuserweight /@ activeIDs),
     "auxiliaryIDs" -> auxiliaryIDs,
     "equationConvention" -> "Tuserweight[id] == expressions[[i]]",
     "derivativeVariables" -> variables,
     "rawDerivatives" -> rawDerivatives,
     "derivativeTargetIntegrals" -> derivativeIntegrals,
     "derivativeTargetIDs" -> derivativeTargetIDs,
     "targetIntegralIDs" -> targetIDs,
     "scalingDegrees" -> degrees,
     "sourceIntegralCount" -> oldCount,
     "userMI" -> userMIData
     |>;
   Join[linearData, <|
     "integralRules" -> shiftedRules,
     "integralCount" -> oldCount + count,
     "equationCount" -> Lookup[linearData, "equationCount", Length[shiftedEquations]] + count,
     "linearEquations" -> Join[basisEquations, shiftedEquations],
     "activeBasis" -> activeData
     |>]
   ];

dsKiraAttachActiveBasis[_Association, setting_] := <|"status" -> "invalidActiveBasis", "reason" -> "KiraActiveBasisMustBeNoneAutomaticOrAssociation", "activeBasisInput" -> setting|>;

dsKiraShiftExplicitTargetItem[item_Integer, activeData_Association] := Module[{oldCount, offset},
   oldCount = Lookup[activeData, "sourceIntegralCount", 0];
   offset = Lookup[activeData, "count", 0];
   If[1 <= item <= oldCount, item + offset, item]
   ];
dsKiraShiftExplicitTargetItem[item_, _Association] := item;

dsKiraEffectiveTargets[linearData_Association, targetSpec_] := Module[{activeData, activeIDs, shifted},
   activeData = Lookup[linearData, "activeBasis", <|"status" -> "disabled"|>];
   If[Lookup[activeData, "status", "disabled"] =!= "configured", Return[targetSpec]];
   activeIDs = activeData["activeIDs"];
   Which[
    targetSpec === Automatic, activeData["targetIntegralIDs"],
    targetSpec === All, All,
    ListQ[targetSpec], shifted = dsKiraShiftExplicitTargetItem[#, activeData] & /@ targetSpec; DeleteDuplicates[Join[activeIDs, shifted]],
    True, DeleteDuplicates[Append[activeIDs, dsKiraShiftExplicitTargetItem[targetSpec, activeData]]]
    ]
   ];


(* ::Section::Closed:: *)
(*DE 变量符号保留门禁*)

(* active-basis 导数是后续 DSDE 的坐标合同。seed、linearData 或 serializer 任一层的
   替换规则若消去这些变量，外部 reduction 已不足以重建微分方程，必须在写文件前拒绝。 *)
dsKiraDEVariableRuleAudit[linearData_Association, kiraRules_, numericStage_] := Module[
   {topo, activeData, variables, audit, baseData, squaredExpressions, protectedInternal,
    imaginaryAxisOccurrences, imaginaryAxisPhysical, imaginaryAxisProtected,
    preAppliedRawRules, preAppliedNormalizedRules, preAppliedLHS, preAppliedRHS,
    preAppliedImaginaryLHS, preAppliedImaginaryRHS,
    exportRawRules, exportNormalizedRules, exportLHS, exportRHS, exportBadLHS, exportBadRHS,
    rawRules, normalizedRules, lhsRules, rhsRules, touchesProtectedQ, touchesImaginaryAxisQ,
    badLHS, badRHS, passQ},
   topo = Lookup[linearData, "topology", <||>];
   activeData = Lookup[linearData, "activeBasis", <||>];
   variables = Lookup[activeData, "derivativeVariables", {}];
   If[! MemberQ[{"symbolic", "postDerivative"}, numericStage],
    Return[<|"status" -> "failed", "passQ" -> False, "reason" -> "invalidNumericStage", "numericStage" -> numericStage|>]
    ];
   If[Lookup[activeData, "status", "disabled"] =!= "configured" || variables === {},
    Return[<|
      "status" -> "notApplicable",
      "passQ" -> True,
      "deVariables" -> {},
      "reason" -> "activeBasisDerivativesNotConfigured",
      "numericStage" -> numericStage
      |>]
    ];
   audit = Lookup[topo, "kinematicCoordinateAudit", <||>];
   baseData = Lookup[audit, "baseCoordinateData", {}];
   squaredExpressions = Lookup[audit, "baseSquaredUserExpressions", {}];
   protectedInternal = DeleteDuplicates@Join[
      variables,
      scalarProductInputToInternal[#, topo] & /@ variables,
      If[Length[baseData] === Length[squaredExpressions],
       MapThread[
        Function[{data, expression},
         If[
          AnyTrue[variables, Function[variable, ! FreeQ[expression, variable]]],
          Lookup[data, "internalVariable", Nothing],
          Nothing
          ]
         ],
        {baseData, squaredExpressions}
        ],
       {}
       ]
      ];
   imaginaryAxisOccurrences = kiraImaginaryAxisKinematicOccurrences[topo];
   imaginaryAxisPhysical = DeleteDuplicates@Lookup[imaginaryAxisOccurrences, "physical", {}];
   imaginaryAxisProtected = DeleteDuplicates@Join[
      imaginaryAxisPhysical,
      scalarProductInputToInternal[#, topo] & /@ imaginaryAxisPhysical,
      If[Length[baseData] === Length[squaredExpressions],
       MapThread[
        Function[{data, expression},
         If[
          AnyTrue[imaginaryAxisPhysical, Function[variable, ! FreeQ[expression, variable]]],
          Lookup[data, "internalVariable", Nothing],
          Nothing
          ]
         ],
        {baseData, squaredExpressions}
        ],
       {}
       ]
      ];
   preAppliedRawRules = Cases[Lookup[linearData, "coefficientRulesApplied", {}], _Rule | _RuleDelayed];
   preAppliedNormalizedRules = normalizeCoefficientRulesForTopology[preAppliedRawRules, topo];
   preAppliedLHS = Cases[preAppliedNormalizedRules, (Rule | RuleDelayed)[lhs_, _] :> lhs];
   preAppliedRHS = Cases[preAppliedNormalizedRules, (Rule | RuleDelayed)[_, rhs_] :> rhs];
   touchesImaginaryAxisQ[expr_] := AnyTrue[imaginaryAxisProtected, ! FreeQ[Unevaluated[expr], #] &];
   preAppliedImaginaryLHS = Pick[preAppliedRawRules, touchesImaginaryAxisQ /@ preAppliedLHS];
   preAppliedImaginaryRHS = Pick[preAppliedRawRules, touchesImaginaryAxisQ /@ preAppliedRHS];
   exportRawRules = Cases[Replace[kiraRules, Automatic -> {}], _Rule | _RuleDelayed];
   exportNormalizedRules = normalizeCoefficientRulesForTopology[exportRawRules, topo];
   exportLHS = Cases[exportNormalizedRules, (Rule | RuleDelayed)[lhs_, _] :> lhs];
   exportRHS = Cases[exportNormalizedRules, (Rule | RuleDelayed)[_, rhs_] :> rhs];
   rawRules = Join[preAppliedRawRules, exportRawRules];
   normalizedRules = normalizeCoefficientRulesForTopology[rawRules, topo];
   lhsRules = Cases[normalizedRules, (Rule | RuleDelayed)[lhs_, _] :> lhs];
   rhsRules = Cases[normalizedRules, (Rule | RuleDelayed)[_, rhs_] :> rhs];
   touchesProtectedQ[expr_] := AnyTrue[protectedInternal, ! FreeQ[Unevaluated[expr], #] &];
   badLHS = Pick[rawRules, touchesProtectedQ /@ lhsRules];
   badRHS = Pick[rawRules, touchesProtectedQ /@ rhsRules];
   exportBadLHS = Pick[exportRawRules, touchesProtectedQ /@ exportLHS];
   exportBadRHS = Pick[exportRawRules, touchesProtectedQ /@ exportRHS];
   If[numericStage === "postDerivative" &&
     (Lookup[activeData, "rawDerivatives", {}] === {} || Lookup[activeData, "derivativeTargetIntegrals", Missing["closure"]] === Missing["closure"]),
    Return[<|"status" -> "failed", "passQ" -> False, "reason" -> "analyticDerivativeClosureMissing",
      "numericStage" -> numericStage, "deVariables" -> variables|>]
    ];
   If[numericStage === "postDerivative" &&
     (preAppliedImaginaryLHS =!= {} || preAppliedImaginaryRHS =!= {}),
    Return[<|
      "status" -> "failed",
      "passQ" -> False,
      "reason" -> "imaginaryAxisKinematicsAlreadyNumerical",
      "numericStage" -> numericStage,
      "deVariables" -> variables,
      "imaginaryAxisPhysicalKinematics" -> imaginaryAxisPhysical,
      "imaginaryAxisProtectedAtoms" -> imaginaryAxisProtected,
      "preAppliedRuleLHSIntersection" -> preAppliedImaginaryLHS,
      "preAppliedRuleRHSDependencies" -> preAppliedImaginaryRHS
      |>]
    ];
   passQ = If[numericStage === "postDerivative", True, badLHS === {} && badRHS === {}];
   <|
    "status" -> If[passQ, "passed", "failed"],
    "passQ" -> TrueQ[passQ],
    "numericStage" -> numericStage,
    "analyticDerivativeConstructedBeforeRulesQ" -> TrueQ[numericStage === "postDerivative"],
    "deVariablesNumericalizedAfterDerivativeQ" -> TrueQ[
      numericStage === "postDerivative" && (exportBadLHS =!= {} || exportBadRHS =!= {})
      ],
    "deVariables" -> variables,
    "protectedInternalAtoms" -> protectedInternal,
    "imaginaryAxisPhysicalKinematics" -> imaginaryAxisPhysical,
    "imaginaryAxisProtectedAtoms" -> imaginaryAxisProtected,
    "preAppliedRuleLHSIntersection" -> preAppliedImaginaryLHS,
    "preAppliedRuleRHSDependencies" -> preAppliedImaginaryRHS,
    "numericRuleLHSIntersection" -> badLHS,
    "numericRuleRHSDependencies" -> badRHS,
    "rulesAudited" -> rawRules,
    "comment" -> If[numericStage === "postDerivative",
      "rules are applied only after raw active-basis derivatives and derivative target closure were constructed",
      "differential variables remain symbolic"]
    |>
   ];

dsStableTadpoleSymmetryData[data_Association] := KeyDrop[data, {"automaticRules"}];
dsStableTadpoleSymmetryData[_] := <||>;


dsKiraExpressionDigest[expr_] := IntegerString[Hash[expr, "SHA256"], 16, 64];


dsKiraRelativePath[path_String, root_String] := FileNameDrop[
   ExpandFileName[path],
   FileNameDepth[ExpandFileName[root]]
   ];


dsKiraExportFileDigests[exportData_Association] := Module[{root, files},
   root = Lookup[exportData, "outputDirectory", None];
   files = Select[Lookup[exportData, "filesWritten", {}], StringQ[#] && FileExistsQ[#] &];
   If[! StringQ[root], Return[{}]];
   Map[
    <|"path" -> dsKiraRelativePath[#, root],
      "sha256" -> IntegerString[FileHash[#, "SHA256"], 16, 64]|> &,
    files
    ]
   ];


(* 从已写出的 artifact 摘要中读取指定文件，避免再次遍历完整 MMA 线性系统。 *)
dsKiraExportFileDigest[fileDigests_List, relativePath_String] := Lookup[
   SelectFirst[
    fileDigests,
    StringReplace[ToString[Lookup[#, "path", ""]], "\\" -> "/"] === relativePath &,
    <||>
    ],
   "sha256",
   Missing["ArtifactFileDigest", relativePath]
   ];


dsKiraArtifactIdentity[exportData_Association, linearData_Association] := Module[
   {contract = Lookup[linearData, "artifactContract", <||>], active, fileDigests, identity},
   active = Lookup[linearData, "activeBasis", <|"status" -> "disabled", "count" -> 0|>];
   fileDigests = dsKiraExportFileDigests[exportData];
   identity = <|
     "identityVersion" -> 2,
     "linearSourceDigest" -> Lookup[contract, "sourceDigest", Missing["sourceDigest"]],
     "serializedKiraInputDigest" -> dsKiraExportFileDigest[fileDigests, "userSystem/ibp.kira"],
     "integralMapDigest" -> dsKiraExpressionDigest[Lookup[linearData, "integralRules", {}]],
     "targetDigest" -> dsKiraExpressionDigest[Lookup[exportData, "targetIntegralIDs", {}]],
     "coefficientRulesDigest" -> dsKiraExpressionDigest[Lookup[linearData, "coefficientRulesApplied", {}]],
     "activeBasisDigest" -> dsKiraExpressionDigest[active],
     "exportFiles" -> fileDigests
     |>;
   Join[identity, <|"exportContentDigest" -> dsKiraExpressionDigest[identity]|>]
   ];

dsKiraExportManifest[exportData_Association, linearData_Association] := <|
   "status" -> "exported",
   "packageVersion" -> $dSIBPVersion,
   "caseName" -> Lookup[linearData, "caseName", Missing["caseName"]],
   "context" -> Lookup[linearData, "dSIBPContextSummary", <||>],
   "linearArtifactContract" -> Lookup[linearData, "artifactContract", <||>],
   "artifactIdentity" -> dsKiraArtifactIdentity[exportData, linearData],
   "equationCount" -> Lookup[exportData, "exportedEquationCount", Missing["equationCount"]],
   "integralCount" -> Lookup[exportData, "integralCount", Missing["integralCount"]],
   "targetIntegralIDs" -> Lookup[exportData, "targetIntegralIDs", {}],
   "numericDummyIntegralId" -> Lookup[exportData, "numericDummyIntegralId", None],
   "numericDummySymbol" -> Lookup[Lookup[exportData, "kiraInput", <||>], "numericDummySymbol", Missing["numericDummySymbol"]],
   "coefficientVariables" -> Lookup[exportData, "coefficientVariables", {}],
   "coefficientAlgebraicGenerators" -> Lookup[exportData, "coefficientAlgebraicGenerators", {}],
   "backendExpressionVariables" -> Lookup[exportData, "backendExpressionVariables", {}],
   "coefficientVariableMap" -> Lookup[exportData, "coefficientVariableMap", {}],
   "backendCoefficientVariables" -> Lookup[exportData, "backendCoefficientVariables", {}],
   "backendImaginaryUnit" -> Lookup[exportData, "backendImaginaryUnit", None],
   "backendCoefficientSyntaxReport" -> Lookup[exportData, "backendCoefficientSyntaxReport", <||>],
   "gaussianPhaseGauge" -> Lookup[exportData, "gaussianPhaseGauge", <|"status" -> "notApplicable"|>],
   "backendKinematicConvention" -> Lookup[exportData, "backendKinematicConvention", <|"status" -> "notRequired"|>],
   "backendKinematicRuleData" -> Lookup[exportData, "backendKinematicRuleData", <|"status" -> "notRequired"|>],
   "physicalCoefficientRulesApplied" -> Lookup[exportData, "physicalCoefficientRulesApplied", {}],
   "pureRationalBackendQ" -> TrueQ[Lookup[exportData, "pureRationalBackendQ", False]],
   "backendTextAudit" -> Lookup[exportData, "backendTextAudit", <|"status" -> "notRun"|>],
   "numericDummyAppendedQ" -> TrueQ[Lookup[exportData, "numericDummyAppendedQ", False]],
   "integralList" -> Lookup[linearData, "integralList", {}],
   "integralRules" -> Lookup[linearData, "integralRules", {}],
   "kiraOrdering" -> Lookup[linearData, "kiraOrdering", <||>],
    "activeBasis" -> Lookup[linearData, "activeBasis", <|"status" -> "disabled", "count" -> 0|>],
     "deVariableNumericRuleAudit" -> Lookup[linearData, "deVariableNumericRuleAudit", <|"status" -> "notRun"|>],
   "coefficientRulesApplied" -> Lookup[linearData, "coefficientRulesApplied", {}],
   "userCoefficientRulesApplied" -> Lookup[linearData, "userCoefficientRulesApplied", {}],
   "zeroPointRules" -> Lookup[Lookup[linearData, "topology", <||>], "zeroPointRules", {}],
   "symmetryRules" -> Lookup[Lookup[linearData, "topology", <||>], "symmetryRules", {}],
   "tadpoleSymmetryData" -> dsStableTadpoleSymmetryData[Lookup[linearData, "tadpoleSymmetryData", <||>]],
   "loopTreeProjectionConvention" -> Lookup[Lookup[linearData, "dSIBPContextSummary", <||>], "loopTreeProjectionConvention", <||>]
   |>;

DSKiraExport[linearData_Association, opts : OptionsPattern[]] := Module[
   {preparedLinearData, activeSetting = OptionValue[KiraActiveBasis], effectiveTargets,
      makeOptions, exportData, manifest, deVariableRuleAudit, numericStage = OptionValue[KiraNumericStage], outputDirectory = OptionValue[OutputDirectory], manifestPath,
     progress = OptionValue[ProgressReporting]},
    If[! KeyExistsQ[linearData, "linearEquations"],
     Message[DSKiraExport::badlinear]; dsErrorPrint["输入缺少 linearEquations。 The input does not contain linearEquations."]; Return[<|"status" -> "failed", "reason" -> "notLinearData"|>]
     ];
    If[Lookup[linearData, "dSIBPStatus", "failed"] =!= "generated" ||
      ! TrueQ[Lookup[Lookup[linearData, "contextCapabilities", <||>], "timeIBPUsableQ", False]],
     Message[DSKiraExport::capability]; dsErrorPrint["请传入 DSLinear 返回且同源门禁通过的 linearData。 Pass linearData returned by DSLinear with a valid same-source gate."]; Return[<|
       "status" -> "failed", "reason" -> "capabilityGate"
       |>]
     ];
    If[! TrueQ[Lookup[Lookup[linearData, "contextCapabilities", <||>], "backendExportUsableQ", False]],
     Message[DSKiraExport::capability];
     dsErrorPrint[
      "当前动力学参数是欠完备或过完备的，只能保留 symbolic linearData 供诊断，不能导出到 Kira。请先提供一组完备且独立的动力学参数。 " <>
       "The current kinematic parameters are undercomplete or overcomplete. Symbolic linearData may be retained for diagnostics, but it cannot be exported to Kira. Provide a complete independent parameter set first."
      ];
     Return[<|
       "status" -> "failed", "reason" -> "backendExportCapabilityGate",
       "capabilities" -> Lookup[linearData, "contextCapabilities", <||>]
       |>]
     ];
    If[TrueQ[OptionValue[KiraRequireCompleteSystem]] &&
      ! TrueQ[Lookup[linearData, "completeSystemQ", False]],
     Message[DSKiraExport::incomplete];
     dsErrorPrint[
      "当前方程集尚未覆盖正式约化所需的全部关系。请扩大 seed 范围或补齐关系，再调用 DSLinear 和 DSKiraExport。 " <>
       "The current equation set does not cover every relation required for a formal reduction. Enlarge the seed range or add the missing relations, then call DSLinear and DSKiraExport again."
      ];
     Return[<|"status" -> "failed", "reason" -> "incompleteSystemForFormalReduction",
       "completeSystemQ" -> Lookup[linearData, "completeSystemQ", False]|>]
     ];
   preparedLinearData = dsKiraAttachActiveBasis[linearData, activeSetting];
    If[Lookup[preparedLinearData, "status", "missing"] =!= "generated",
    Message[DSKiraExport::badbasis, dsReasonSentence[Lookup[preparedLinearData, "reason", "unknown"]]];
    dsErrorPrint["active basis 或其导数 target closure 未通过导出门禁。 The active basis or its derivative target closure failed the export gate."]; Return[preparedLinearData]
     ];
    If[! MemberQ[{"symbolic", "postDerivative"}, numericStage],
     Message[DSKiraExport::badstage, numericStage]; Return[<|"status" -> "failed", "reason" -> "invalidNumericStage"|>]
     ];
    deVariableRuleAudit = dsKiraDEVariableRuleAudit[preparedLinearData, OptionValue[KiraCoefficientRules], numericStage];
    If[! TrueQ[Lookup[deVariableRuleAudit, "passQ", False]],
     Message[DSKiraExport::devarrules, dsDeVarRulesSentence[KeyTake[deVariableRuleAudit, {"deVariables", "numericRuleLHSIntersection", "numericRuleRHSDependencies"}]]];
     If[
      Lookup[deVariableRuleAudit, "reason", None] === "imaginaryAxisKinematicsAlreadyNumerical",
      dsErrorPrint[
       "需要转换到 Kira 虚轴的运动学坐标已在 DSLinear 阶段被数值替换，无法再保证方程与 manifest 使用同一个点。请用 symbolic 模式生成 linearData，并只在 formal Kira plan 的 postDerivative 阶段提供这些数值。 " <>
        "Kinematic coordinates that require the Kira imaginary-axis map were already replaced numerically by DSLinear, so the equations and manifest cannot be guaranteed to use the same point. Build symbolic linearData and provide these values only at the postDerivative stage of the formal Kira plan."
       ],
      dsErrorPrint["symbolic 阶段必须保留 DE 变量；postDerivative 只允许在解析一阶导数与 closure 已构造后使用。 The symbolic stage must preserve every DE variable; postDerivative is allowed only after analytic first derivatives and their closure have been constructed."]
      ];
     Return[<|
       "status" -> "failed",
       "reason" -> Lookup[deVariableRuleAudit, "reason", "differentialVariablesWouldBeNumerical"],
       "deVariableNumericRuleAudit" -> deVariableRuleAudit
       |>]
     ];
    preparedLinearData = Join[preparedLinearData, <|"deVariableNumericRuleAudit" -> deVariableRuleAudit|>];
    effectiveTargets = dsKiraEffectiveTargets[preparedLinearData, OptionValue[KiraTargetIntegrals]];
    makeOptions = DeleteCases[
      FilterRules[{opts}, Options[makeKiraExportData]],
       HoldPattern[(KiraTargetIntegrals | KiraNumericStage) -> _]
     ];
   exportData = dsStageRun[
     "序列化 Kira 基础输入 / Serializing basic Kira input",
     makeKiraExportData[
      preparedLinearData,
      Sequence @@ makeOptions,
       KiraTargetIntegrals -> effectiveTargets
      ],
     progress
     ];
   If[Lookup[exportData, "status", "missing"] =!= "ready",
    Message[DSKiraExport::failed, dsReasonSentence[Lookup[exportData, "reason", Lookup[exportData, "status", Missing["status"]]],
      "Kira 输入序列化未通过导出门禁；请检查返回结果中的 kiraInput 字段。 Kira-input serialization failed the export gate; inspect the kiraInput field of the returned result."]];
    dsErrorPrint["package 未运行 Kira；当前只报告导出门禁失败。 The package did not run Kira; only the failed export gate is reported."]; Return[exportData]
    ];
   manifest = dsKiraExportManifest[exportData, Lookup[exportData, "linearSystem", preparedLinearData]];
    If[StringQ[outputDirectory],
     manifestPath = FileNameJoin[{outputDirectory, "dsibp-export-manifest.wl"}];
     Quiet[Check[Put[manifest, manifestPath], manifestPath = $Failed]],
     manifestPath = Missing["NotWritten"]
     ];
    If[StringQ[outputDirectory] && manifestPath === $Failed,
     Message[DSKiraExport::writefailed, FileNameJoin[{outputDirectory, "dsibp-export-manifest.wl"}]];
     dsErrorPrint[
      "Kira 输入文件可能已经部分写出，但缺少 manifest，不能用于导入。请修正输出目录后从 DSKiraExport 重新生成整套文件。 " <>
       "Some Kira input files may have been written, but the manifest is missing, so they cannot be imported. Fix the output directory and regenerate the complete file set with DSKiraExport."
      ];
     Return[Join[exportData, <|
       "status" -> "failed",
       "reason" -> "manifestWriteFailed",
       "dSIBPExportManifest" -> manifest,
       "dSIBPExportManifestPath" -> $Failed
       |>]]
     ];
     Join[exportData, <|
       "deVariableNumericRuleAudit" -> deVariableRuleAudit,
       "dSIBPExportManifest" -> manifest,
      "dSIBPExportManifestPath" -> manifestPath
      |>]
    ];
