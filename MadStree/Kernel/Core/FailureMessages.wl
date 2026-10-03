(* ::Package:: *)

(***
文件：FailureMessages.wl
用途：为 MadStree 的结构化 Failure 补充简洁、完整且可直接行动的中英文说明。
接口：msFailure[tag, details] 保留原 tag 与全部诊断字段，并增加 message、messageCN 和 MessageTemplate。
边界：本模块只负责用户反馈，不改变任何数学条件、控制流、状态代码或返回数据字段。
***)

(* ::Chapter:: *)
(*自然语言句子*)

msNaturalSentence[text_String, terminator_String : "."] := Module[{sentence, first},
  sentence = StringTrim[text];
  If[sentence === "", Return[sentence]];
  first = StringTake[sentence, 1];
  If[StringMatchQ[first, LetterCharacter] && first === ToLowerCase[first],
    sentence = ToUpperCase[first] <> StringDrop[sentence, 1]
  ];
  If[! StringMatchQ[sentence, ___ ~~ ("." | "?" | "!" | "。" | "？" | "！")],
    sentence = sentence <> terminator
  ];
  sentence
];


(* ::Chapter:: *)
(*按门禁类别给出的缺省说明*)

msFailureDefaultText[tag_String] := Switch[tag,
  "BoundaryGammaPole",
    {
      "This boundary integral is undefined because its Gamma argument is zero or a negative integer. Add an analytic regulator to the listed component and try again",
      "该边界积分的 Gamma 参数是零或负整数，因此积分未定义。请对列出的分量加入解析正规化后重试"
    },
  "TimeIBPMasterRequiresRegularization",
    {
      "MadStree cannot initialize this time-integral family because a zero-shift master has a non-positive-integer Gamma argument. Add an analytic regulator to the listed vertex time power and initialize the family again",
      "MadStree 无法初始化这个时间积分函数族，因为某个零平移主积分的 Gamma 参数是非正整数。请对列出的顶点时间幂加入解析正规化后重新初始化"
    },
  "TimeIBPShiftHitsDivergentLayer",
    {
      "The requested time-index shift makes a small-time branch Gamma argument a non-positive integer, so this reduction would use a divergent integral. Stay within the reported safe shift range or add an analytic regulator",
      "请求的时间指标平移使某个小时间分支的 Gamma 参数成为非正整数，因此该约化会使用发散积分。请把指标限制在所列安全范围内，或加入解析正规化"
    },
  "UnsupportedHSmallTLeadingPowers" | "UnsupportedHState",
    {
      "MadStree cannot certify the small-time leading powers for one h state in this convention. Use the documented h conventions or add an analytic regulator before initializing the family",
      "MadStree 无法认证当前约定下某个 h 状态的小时间领头幂次。请使用手册支持的 h 约定，或在初始化前加入解析正规化"
    },
  "UnsupportedTensorSlotKind" | "TensorM0DiagonalizationNotCertified",
    {
      "The local tensor basis does not exactly diagonalize M0 for the listed sector. MadStree stopped instead of using an unverified dense inverse; review the slot definition or derive a supported local basis",
      "局部张量基不能在所列区段精确对角化 M0。MadStree 已停止，未使用未经验证的稠密求逆；请检查 slot 定义或推导受支持的局部基"
    },
  "InitializedContextRequired" | "NotMSIntegral",
    {
      "This operation needs a valid MadStree context and an integral created from that context. Initialize the family again and pass the returned objects directly",
      "当前操作需要有效的 MadStree 上下文及由该上下文生成的积分。请重新初始化函数族，并直接传入返回对象"
    },
  "InvalidTreeInputFields" | "InvalidTreeTopology" | "MalformedTreeSpec" |
  "InvalidVertexFamilyInput" | "InvalidVertexFamilyInputFields" | "MalformedVertexFamilySpec" |
  "VertexFamilyKiNuiLengths" | "MalformedContactBundles",
    {
      "The topology or function-family input is incomplete or inconsistent. Correct the listed fields and initialize the family again",
      "拓扑或函数族输入不完整，或者各字段彼此不一致。请修正列出的字段后重新初始化"
    },
  "UnknownSector" | "UnknownComponent" | "SectorIdentityCollision",
    {
      "The requested sector or component does not belong uniquely to this context. Use a sector and component returned by the current initialization",
      "请求的区段或分量不能唯一归属于当前上下文。请使用本次初始化实际返回的区段和分量"
    },
  "InvalidRecurrenceInput" | "BadTimeShifts" | "ZeroShiftAtComponent" |
  "RecurrenceCycle" | "MasterBasisMustBeCompletePermutation",
    {
      "The requested recurrence or master basis is not valid for this context. Check the listed shifts and use each current master exactly once",
      "请求的递推或主积分基底不适用于当前上下文。请检查列出的指标移动，并确保当前每个主积分恰好出现一次"
    },
  "UnsupportedBasisConversion" | "UnsupportedBasisDirection" | "UnsupportedDSIBPJ" |
  "UnknownNuConvention",
    {
      "The requested convention or basis conversion is not supported for these objects. Choose one of the documented conventions and conversion directions",
      "当前对象不支持所请求的约定或基底转换。请选择手册中列出的约定和转换方向"
    },
  "CertifiedDLogRequired" | "InvalidDLogData" | "FormulaMatrixGenerationFailed" |
  "FormulaArtifactWriteFailed" | "EmptyDifferentialEquationSupport",
    {
      "A certified analytic differential equation is not available. Resolve the reported formula or write error before starting boundary generation or numerical transport",
      "当前没有可用的已认证解析微分方程。请先解决所列公式或写出错误，再生成边界或启动数值输运"
    },
  "NumericalRulesRequired" | "IncompleteNumericalPoint" | "ParameterRulesRequired" |
  "IncompleteEvaluationParameters" | "IncompleteEvaluationPoint",
    {
      "The numerical substitutions are missing or incomplete. Supply numerical values for every listed symbol before starting the numerical solve",
      "数值替换缺失或不完整。请为列出的每个符号提供数值后再启动数值求解"
    },
  "EvaluationPointSequenceRequiresPoint" | "EvaluationCoordinateHeader" |
  "EvaluationPointTag" | "EvaluationCoordinateRowWidth" |
  "DuplicateEvaluationCoordinate" | "DuplicateParameterRuleSymbol" |
  "NonDifferentialEvaluationCoordinate" | "FixedAndRunningCoordinateOverlap" |
  "DifferentComplexAffinePlane" | "ZeroAffineDirection" | "IncompleteKinematicPath",
    {
      "The pointSequence table is not a valid sequence of complex-affine coordinate points. Correct the header, row widths, duplicate fields, or fixed/running overlap shown below",
      "pointSequence 表不是有效的复仿射坐标点列。请按下方信息修正表头、行宽、重复字段或固定参数与跑动坐标的重叠"
    },
  "UnknownOption" | "MessageLanguage" | "SingularityMode" |
  "FlintNDEPathPlanningBooleanRequired" | "DirectTransportAvoidModeRequired" |
  "InvalidRankOrder" | "BoundaryScaleMustExceedOne" |
  "BoundaryScaleGreaterThanOneRequired" | "BoundarySeriesOrderMustBeNonNegative" |
  "BoundaryVectorDimension" | "InvalidExportFormats" | "InvalidSignificantDigits" |
  "ParallelTaskCountPositiveIntegerRequired",
    {
      "One option value is invalid for this operation. Replace the listed value with one of the allowed values and try again",
      "当前操作的某个选项值无效。请把列出的值改为允许值后重试"
    },
  "BoundaryChartNotCertified" | "VanishingComponentBoundaryEnergy" |
  "RegularSingularPullbackRequired",
    {
      "The automatic infinity boundary could not be certified for this point. Review the listed boundary component or chart and choose a nonsingular, supported parameter point",
      "当前点的自动无穷远边界未能通过认证。请检查列出的边界分量或坐标图，并选择非奇异且受支持的参数点"
    },
  "SingularBoundaryAnchor" | "DirectPathCrossesSingularity" |
  "SingularPathOnUserPolyline",
    {
      "The requested boundary anchor or path lies on a differential-equation singularity. Change the point sequence or enable the documented singularity route before retrying",
      "请求的边界锚点或路径落在微分方程奇点上。请修改点列，或启用手册说明的奇点路线后重试"
    },
  "RelativeFlintNDEPathRequired" | "FlintNDEConfigurationFailed" |
  "FlintNDEPackageNotFound" | "FlintNDENotAvailable" | "PythonFlintUnavailable",
    {
      "MadStree cannot find a usable FlintNDE installation. Check the configured relative package path and the Python environment, then try again",
      "MadStree 找不到可用的 FlintNDE。请检查配置的相对程序包路径和 Python 环境后重试"
    },
  "RuntimePathTooLong",
    {
      "The runtime file path is too long for the current Windows environment. Move the project or runtime directory to a shorter path and try again",
      "运行文件路径超过当前 Windows 环境可处理的长度。请把项目或运行目录移到更短的路径后重试"
    },
  "RuntimeDirectoryCreationFailed" | "RuntimeInputWriteFailed" |
  "OutputDirectoryCreationFailed" | "EvaluationExportFailed",
    {
      "MadStree could not create or write the requested output. Check the listed path, permissions, and available disk space before retrying",
      "MadStree 无法创建或写入请求的输出。请检查列出的路径、权限和可用磁盘空间后重试"
    },
  "FlintNDELaunchFailed" | "FlintNDEProcessFailed" | "FlintNDEOutputMissing" |
  "FlintNDEOutputInvalid" | "FlintNDEPathEvaluationFailed" | "FlintNDEEpBatchFailed",
    {
      "The FlintNDE backend did not complete successfully. Read the attached backend message, correct the reported Python or numerical input problem, and retry",
      "FlintNDE 后端未能成功完成。请阅读附带的后端说明，修正所列 Python 环境或数值输入问题后重试"
    },
  "SavedPointResultsRequired",
    {
      "There are no saved ordinary-point results to export. Compute at least one saved point and call the export function again",
      "当前没有可导出的普通保存点结果。请先计算至少一个保存点，再调用导出函数"
    },
  "EpValueListEmpty" | "MaximumEpPowerIntegerRequired" |
  "MaximumEpPowerBelowLeadingPower" | "EpGoalDigitsPositiveIntegerRequired" |
  "EpFitExtraOrderNonNegativeIntegerRequired" | "EpFitOrderIncrementPositiveIntegerRequired" |
  "EpFitMaximumRoundsPositiveIntegerRequired" | "EpInitialInternalMaximumPowerRequired" |
  "EpSampleAngleRangePairRequired" | "EpSampleAngleRangeOpenRealIntervalRequired" |
  "EpSamplePointsAndAngleRangeMutuallyExclusive" |
  "EpSamplePointsExactDistinctNonzeroListRequired" |
  "EpValidationPointsExactDistinctNonzeroListRequired" |
  "EpProductionAndValidationPointsMustBeDisjoint",
    {
      "The regulator sampling or fitting options are inconsistent. Correct the listed order, range, or point pool and keep all production and validation points distinct and nonzero",
      "正规化采样或拟合选项彼此不一致。请修正列出的阶数、范围或点池，并确保生产点与验证点互不重复且都不为零"
    },
  "EpDependentPathCoordinatesNotCertified" | "EpNonLaurentExpression" |
  "EpLaurentSeriesUnresolved" | "EpNonIntegerPower" |
  "EpSingularDifferentialEquation" | "EpLeadingPowerAboveRequestedRange",
    {
      "The regulator dependence cannot be certified as the requested Laurent problem. Review the listed expression or path coordinates and choose a supported regulator setup",
      "当前正规化因子依赖无法认证为所请求的 Laurent 问题。请检查列出的表达式或路径坐标，并改用受支持的正规化设置"
    },
  "EpSeriesPointEvaluationFailed" | "EpSeriesFitResultMissing" |
  "EpBatchExecutionConfigurationMismatch",
    {
      "The regulator reconstruction did not produce a usable numerical result. Check the failed sample or worker configuration, then retry with valid nonzero points",
      "正规化重建未得到可用的数值结果。请检查失败的采样点或工作进程配置，并使用有效的非零点重试"
    },
  "ModuleLoadFailed" | "ModuleContractMissing" | "module" | "symbol",
    {
      "MadStree did not load all required modules. Reinstall the complete current version and load the package again",
      "MadStree 未能加载全部必需模块。请重新安装完整的当前版本，然后再次加载程序包"
    },
  _,
    {
      "MadStree stopped because an input or capability check failed. Review the attached fields, correct the reported input, and try again",
      "MadStree 因输入或能力检查未通过而停止。请查看附带字段，修正所列输入后重试"
    }
];


(* ::Chapter:: *)
(*结构化 Failure 构造器*)

msFailure[tag_String, details_Association : <||>] := Module[
  {defaults, english, chinese},
  defaults = msFailureDefaultText[tag];
  english = msNaturalSentence[Lookup[details, "message", First[defaults]], "."];
  chinese = msNaturalSentence[Lookup[details, "messageCN", Last[defaults]], "。"];
  System`Failure[
    tag,
    Join[
      details,
      <|
        "message" -> english,
        "messageCN" -> chinese,
        "MessageTemplate" -> (chinese <> " " <> english)
      |>
    ]
  ]
];


msFailure[tag_String, details_] := msFailure[
  tag,
  <|"details" -> HoldForm[details]|>
];
