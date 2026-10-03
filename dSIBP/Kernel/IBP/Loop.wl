(* ::Package:: *)

(* ::Chapter:: *)
(*统一三槽 seed 与 linearData 高层入口*)

(* DSSeeds 只生产符号 general templates；连续指标域由 DSGenerateIBP 持有，
   数值系数由 DSLinear 的 CoefficientRules 持有。 *)
Options[DSSeeds] = {ProgressReporting -> Automatic};
Options[DSLinear] = {
   LinearSystemMode -> "symbolic",
   CoefficientRules -> Automatic,
   KiraOrdering -> Automatic,
   AuditLevel -> "standard",
   ProgressReporting -> Automatic
   };

DSSeeds::noinit = "DSSeeds 需要有效的 DSInit context。";
DSSeeds::failed = "canonical seed 生成未通过门禁：`1`。";
DSSeeds::capability = "当前 context 不具备 seed 生成所需能力：`1`。";
DSLinear::noinit = "DSLinear 需要有效的 DSInit context。";
DSLinear::badseed = "DSLinear 需要 DSSeeds 返回的 canonical seed Association。";
DSLinear::badmode = "LinearSystemMode 只允许 \"symbolic\" 或 \"numeric\"，收到 `1`。";
DSLinear::failed = "linearData 生成未通过门禁：`1`。";
DSLinear::capability = "当前 context 不具备 linearData 生成所需能力：`1`。";
DSLinear::context = "seedData 与 context 不是同一次初始化的产物。";


(* loop seed 的底层生成器使用 qq/qk/kk；公开高层入口必须在序列化前统一投影到
   当前 context 的用户坐标，使 DSLinear 与后续 backend 不再接触内部 Gram 原子。 *)
dsLoopSeedExpressionToPublicCoordinates[expr_, topo_Association] := publicProtectJMap[
   expr,
   Function[body,
    Expand[scalarProductInternalToUser[body /. internalISPToUserRules[topo], topo]]
    ]
   ];


dsPublicLoopSeedEntry[entry_Association, topo_Association] := If[
   KeyExistsQ[entry, "equation"],
   Join[entry, <|
     "equation" -> dsLoopSeedExpressionToPublicCoordinates[entry["equation"], topo]
     |>],
   entry
   ];


dsPublicLoopSeedBatch[seedData_Association, topo_Association] := If[
   Lookup[seedData, "status", "missing"] === "generated",
   Join[seedData, <|
     "equations" -> (dsPublicLoopSeedEntry[#, topo] & /@ Lookup[seedData, "equations", {}]),
     "allSeeds" -> (dsPublicLoopSeedEntry[#, topo] & /@ Lookup[seedData, "allSeeds", {}]),
     "coordinateRepresentation" -> "user"
     |>],
   seedData
   ];


(* 公开 seed 只有 J[aList,linePacks,ispList]。论文 vertex basis 仍可由 tree
   公式模块内部构造，但不得再决定 DSSeeds/DSLinear 的公开积分形状。 *)


DSLinear[seedData_Association, context_: Automatic, opts : OptionsPattern[]] := Module[
   {resolved, workingSeedData, workingContract, internalDigest, linearData, mode = OptionValue[LinearSystemMode],
    auditLevel = OptionValue[AuditLevel], progress = OptionValue[ProgressReporting]},
   resolved = dsResolveContext[context];
   If[Head[resolved] === Missing,
    Message[DSLinear::noinit]; dsErrorPrint["请传入与 seed 同源的 DSInit context。 Pass the DSInit context from which the seeds originated."]; Return[<|"status" -> "failed", "reason" -> "missingContext"|>]
    ];
    If[! KeyExistsQ[seedData, "completeCanonicalQ"],
     Message[DSLinear::badseed]; dsErrorPrint["输入不是 canonical seed batch。 The input is not a canonical seed batch."]; Return[<|"status" -> "failed", "reason" -> "notCanonicalSeedBatch"|>]
     ];
    If[Lookup[Lookup[seedData, "dSIBPContextSummary", <||>], "inputHash", Missing["seedHash"]] =!=
      Lookup[resolved, "inputHash", Missing["contextHash"]],
     Message[DSLinear::context]; dsErrorPrint["seed 与 context 的 inputHash 不一致。 The seed and context inputHash values differ."]; Return[<|
       "status" -> "failed", "reason" -> "contextMismatch"
       |>]
     ];
    If[Lookup[seedData, "dSIBPStatus", "failed"] =!= "generated" ||
      ! dsContextCapabilityQ[resolved, "timeIBPUsableQ"] ||
      (Lookup[resolved["topology"], "ibpMode", "full"] === "full" &&
        ! dsContextCapabilityQ[resolved, "momentumIBPUsableQ"]),
     Message[DSLinear::capability, dsContextCapabilities[resolved]];
     dsErrorPrint["seed 或 context 未通过 linearData 能力门禁。 The seed or context failed the linearData capability gate."]; Return[<|
       "status" -> "failed", "reason" -> "capabilityGate",
       "capabilities" -> dsContextCapabilities[resolved]
       |>]
     ];
   If[! MemberQ[{"symbolic", "numeric"}, mode],
    Message[DSLinear::badmode, mode]; dsErrorPrint["linearData 模式无效。 The linearData mode is invalid."]; Return[<|"status" -> "failed", "reason" -> "invalidLinearSystemMode", "mode" -> mode|>]
    ];
    If[! MemberQ[{"standard", "full"}, auditLevel],
    Return[<|"status" -> "failed", "reason" -> "invalidAuditLevel",
      "auditLevel" -> auditLevel, "allowedAuditLevels" -> {"standard", "full"}|>]
     ];
    workingSeedData = If[
      Lookup[resolved["topology"], "ibpMode", "full"] === "timeOnly",
      dsTimeOnlyDataToInternal020[seedData, resolved],
      seedData
      ];
    If[workingSeedData === $Failed,
     Return[<|"status" -> "failed", "reason" -> "timeOnlyPublicInputRequired"|>]
     ];
    If[Lookup[resolved["topology"], "ibpMode", "full"] === "timeOnly",
     workingContract = Lookup[workingSeedData, "artifactContract", <||>];
     internalDigest = dsGeneratedIBPSourceDigest[
       Lookup[workingSeedData, "equations", {}],
       resolved,
       <|"rangeRules" -> Lookup[workingSeedData, "targetEnvelopeRules", {}]|>
       ];
     workingSeedData = Join[workingSeedData, <|
        "representation" -> "J[aList,linePacks,ispList]",
        "artifactContract" -> Join[workingContract, <|"sourceDigest" -> internalDigest|>]
        |>]
     ];
    linearData = dsStageRun[
     "转换 backend-neutral linearData / Converting to backend-neutral linearData",
     If[mode === "numeric",
       makeSampledLinearSystemData[
         workingSeedData,
        resolved["topology"],
        KiraOrdering -> OptionValue[KiraOrdering],
        AuditLevel -> auditLevel,
        CoefficientRules -> OptionValue[CoefficientRules]
        ],
       makeLinearSystemData[
         workingSeedData,
        resolved["topology"],
        KiraOrdering -> OptionValue[KiraOrdering],
        AuditLevel -> auditLevel
        ]
     ],
     progress
     ];
   If[Lookup[linearData, "status", "missing"] =!= "generated",
    Message[DSLinear::failed, Lookup[linearData, "status", Missing["status"]]];
    dsErrorPrint["linearData 未通过 canonical/linearity 门禁。 linearData failed the canonical or linearity gate."];
    Return[Join[linearData, <|"dSIBPStatus" -> "failed", "dSIBPContextSummary" -> dsContextSummary[resolved]|>]]
    ];
    If[Lookup[resolved["topology"], "ibpMode", "full"] === "timeOnly",
     linearData = dsTimeOnlyDataToPublic020[linearData, resolved];
     If[linearData === $Failed,
      Return[<|"status" -> "failed", "reason" -> "timeOnlyPublicConversionFailed"|>]
      ];
     linearData = Join[linearData, <|"representation" -> "J[sectorKey,timeShifts,stateBits]"|>]
     ];
    Join[linearData, <|
      "dSIBPStatus" -> "generated",
      "dSIBPContextSummary" -> dsContextSummary[resolved],
       "contextCapabilities" -> dsContextCapabilities[resolved],
        "contextInputHash" -> Lookup[resolved, "inputHash", Missing["inputHash"]]
      |>]
   ];
