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
  "BadStateBits",
    {
      "The stateBits of this MSIntegral do not occur in the state order of the requested sector. Use a state-bit vector listed in that sector's stateOrder and build the integral again",
      "该 MSIntegral 的 stateBits 不在所请求区段的 stateOrder 中。请使用该区段 stateOrder 列出的状态位向量重新构造积分"
    },
  "ComponentTimeCount",
    {
      "MSConvertBasis received a componentTimes list whose length differs from the number of vertex components of the requested sector. Pass exactly one component time per vertex component of that sector and try again",
      "MSConvertBasis 收到的 componentTimes 长度与所请求区段的顶点分量个数不一致。请为该区段的每个顶点分量各传入一个分量时间后重试"
    },
  "StateVectorDimension",
    {
      "MSConvertBasis received a state vector whose length differs from the master count of the requested sector. Pass a vector with exactly one entry per master of that sector and try again",
      "MSConvertBasis 收到的状态向量长度与所请求区段的主积分个数不一致。请传入长度恰等于该区段主积分个数的向量后重试"
    },
  "DSIBPExpressionConversionFailed",
    {
      "MSFromDSIBPExpression found dSIBP J integrals in the expression that are not valid in this MadStree context. Inspect the attached per-integral failures, correct the sector keys, time shifts, or state bits, and convert the expression again",
      "MSFromDSIBPExpression 在表达式中发现无法归属于当前 MadStree 上下文的 dSIBP J 积分。请检查附带的逐积分失败记录，修正 sectorKey、时间平移或状态位后重新转换表达式"
    },
  "ComputedEvaluationRequired",
    {
      "MSExportEvaluationData needs a completed MSEvaluatePath result whose status is computed. Run MSEvaluatePath successfully first and then pass its full result association to the export function again",
      "MSExportEvaluationData 需要状态为 computed 的完整 MSEvaluatePath 结果。请先成功运行 MSEvaluatePath，再把其完整结果 Association 重新传入导出函数"
    },
  "OutputDirectoryRequired",
    {
      "The output directory option must be Automatic or a path string. Replace the attached value with Automatic or a valid directory path and call the function again",
      "输出目录选项必须是 Automatic 或路径字符串。请把附带的值改为 Automatic 或有效目录路径后重新调用该函数"
    },
  "TimePowerRulesRequired",
    {
      "The TimePowerRules option must be Automatic or a list of substitution rules for the vertex time powers. Replace the attached value and call MSFormulaData again",
      "TimePowerRules 选项必须是 Automatic 或顶点时间幂的替换规则列表。请替换附带的值后重新调用 MSFormulaData"
    },
  "ExactSingularConnectionRequired",
    {
      "The pullback of the dlog connection along the automatic boundary curve must stay an exact Q(i)(t) matrix, but an approximate real number entered through the target rules. Provide exact rational or Gaussian-rational parameter values and generate the boundary data again",
      "dlog 联络沿自动边界曲线的回拉必须保持为精确的 Q(i)(t) 矩阵，但目标替换规则引入了近似实数。请提供精确的有理或高斯有理参数值后重新生成边界数据"
    },
  "SectorLeadingSystemFailed",
    {
      "The linear system that fixes the Frobenius leading vector of the attached sector could not be solved, so the infinity boundary is not certified. Check this sector's leading exponent and root index, or choose another boundary target point",
      "无法求解固定所附区段 Frobenius 领头向量的线性方程组，因此无穷远边界未通过认证。请检查该区段的领头指数与根指标，或改选其它边界目标点"
    },
  "SectorLeadingSystemResidual",
    {
      "The solved Frobenius leading vector of the attached sector left a nonzero residual, so this boundary solution is not certified. Review the sector residue matrix and the requested exponent, or choose another boundary target point",
      "所附区段求解出的 Frobenius 领头向量残差不为零，因此该边界解未通过认证。请检查该区段的残差矩阵与所请求的指数，或改选其它边界目标点"
    },
  "FlintNDEExactPathRequired",
    {
      "FlintNDE transport requires the dlog letters and path parameters along each complex-affine segment to be exact affine Q(i) data, but the attached letter or point parameters are not. Use exact rational or Gaussian-rational coordinate values along the path and evaluate again",
      "FlintNDE 输运要求每个复仿射段上的 dlog 字母与路径参数是精确的仿射 Q(i) 数据，但所附字母或点参数不满足。请沿路径使用精确的有理或高斯有理坐标值后重新求值"
    },
  "FlintNDEExactFrobeniusBoundaryRequired",
    {
      "The regular-singular Frobenius boundary sent to FlintNDE must consist of exact Q(i)(t) connection and branch data, but approximate numbers entered the attached boundary. Provide exact target parameter rules and generate the boundary again",
      "送往 FlintNDE 的正则奇点 Frobenius 边界必须由精确的 Q(i)(t) 联络与分支数据构成，但所附边界中出现了近似数。请提供精确的目标参数替换规则后重新生成边界"
    },
  "RuntimeDirectoryRequired",
    {
      "The MSRuntimeDirectory setting must be Automatic or a path string. Replace the attached value with Automatic or a valid directory path and try again",
      "MSRuntimeDirectory 设置必须是 Automatic 或路径字符串。请把附带的值改为 Automatic 或有效目录路径后重试"
    },
  "InvalidEpSeriesArguments",
    {
      "MSReconstructEpSeries was called with arguments that do not match its usage contract. Call it as MSReconstructEpSeries[context, ep, pointSequence, MaximumEpPower->n] with an initialized context, the regulator symbol, and a valid pointSequence table",
      "MSReconstructEpSeries 的调用参数不符合其用法约定。请按 MSReconstructEpSeries[context, ep, pointSequence, MaximumEpPower->n] 的形式调用，并传入已初始化的上下文、正规化符号和有效的 pointSequence 表"
    },
  "VertexFamilyHankelBranches",
    {
      "The hankelBranches list of this compact vertex-family input must contain one entry per h block, that is Length[ki]-1 entries, each equal to 1 or 2. Correct the attached list and initialize the family again",
      "该紧凑顶点函数族输入的 hankelBranches 列表必须为每个 h 块提供一个条目，即共 Length[ki]-1 个、每个取值 1 或 2。请修正所附列表后重新初始化函数族"
    },
  "InvalidContactBundleLines",
    {
      "The explicit thetaBundles must list each full-line id at most once and may only use ids of actual full lines. Correct the attached explicit line-id list against the reported full-line ids and initialize the tree again",
      "显式 thetaBundles 中每个 full 线 id 至多出现一次，且只能使用实际 full 线的 id。请对照所报告的 full 线 id 修正附带的显式线 id 列表后重新初始化树图"
    },
  "ContactBundleMustShareThetaArgument",
    {
      "Every explicit theta bundle must be nonempty and its full lines must share one common pair of endpoints, which fixes their shared theta argument. Correct the attached bundle groups and initialize the tree again",
      "每个显式 theta 束必须非空，且束内 full 线必须共享同一对端点，以固定其共享的 theta 参数。请修正所附的束分组后重新初始化树图"
    },
  _,
    {
      "MadStree stopped because an input or capability check failed. Review the attached fields, correct the reported input, and try again",
      "MadStree 因输入或能力检查未通过而停止。请查看附带字段，修正所列输入后重试"
    }
];


(* ::Chapter:: *)
(*issue Association 的双语句子渲染*)

(* 只渲染 Message 文本；返回数据中的 issue Association 与 code 键保持不变。 *)
msIssueValueText[value_] := ToString[value, InputForm];

msIssueValueListText[values_List] := StringRiffle[
  Map[msIssueValueText, values],
  ", "
];

msIssueSentence[text_String] := text;

msIssueSentence[issue_Association] := Module[
  {code = Lookup[issue, "code", None], defaults, allowedBlockKeys},
  Switch[code,
    "missingTreeFields",
      "树图输入缺少必需字段 " <> msIssueValueListText[Lookup[issue, "fields", {}]] <>
        "；请补齐这些字段后重新初始化。 The tree input is missing the required fields " <>
        msIssueValueListText[Lookup[issue, "fields", {}]] <>
        "; add them and initialize again.",
    "missingVertexFields",
      "第 " <> msIssueValueText[Lookup[issue, "position", Missing["Absent"]]] <>
        " 个顶点缺少必需字段 " <> msIssueValueListText[Lookup[issue, "fields", {}]] <>
        "；请补齐这些字段后重新初始化。 Vertex " <>
        msIssueValueText[Lookup[issue, "position", Missing["Absent"]]] <>
        " is missing the required fields " <>
        msIssueValueListText[Lookup[issue, "fields", {}]] <>
        "; add them and initialize again.",
    "missingLineFields",
      "第 " <> msIssueValueText[Lookup[issue, "position", Missing["Absent"]]] <>
        " 条线缺少必需字段 " <> msIssueValueListText[Lookup[issue, "fields", {}]] <>
        "；请补齐这些字段后重新初始化。 Line " <>
        msIssueValueText[Lookup[issue, "position", Missing["Absent"]]] <>
        " is missing the required fields " <>
        msIssueValueListText[Lookup[issue, "fields", {}]] <>
        "; add them and initialize again.",
    "vertexIdsMustBePresentAndUnique",
      "顶点 id 必须全部存在且互不重复，当前 id 列表为 " <>
        msIssueValueListText[Lookup[issue, "ids", {}]] <>
        "；请修正后重新初始化。 Vertex ids must all be present and unique; the current id list is " <>
        msIssueValueListText[Lookup[issue, "ids", {}]] <>
        ". Correct it and initialize again.",
    "vertexTypeMustBePlusOrMinus",
      "顶点 " <> msIssueValueText[Lookup[issue, "vertex", Missing["Absent"]]] <>
        " 的 vertexType 值 " <> msIssueValueText[Lookup[issue, "value", Missing["Absent"]]] <>
        " 无效，只允许 \"+\" 或 \"-\"；请修正后重新初始化。 Vertex " <>
        msIssueValueText[Lookup[issue, "vertex", Missing["Absent"]]] <>
        " has the invalid vertexType " <>
        msIssueValueText[Lookup[issue, "value", Missing["Absent"]]] <>
        "; only \"+\" or \"-\" are allowed. Correct it and initialize again.",
    "unsupportedLineType",
      "第 " <> msIssueValueText[Lookup[issue, "line", Missing["Absent"]]] <>
        " 条线的 type 值 " <> msIssueValueText[Lookup[issue, "type", Missing["Absent"]]] <>
        " 不受支持，只允许 " <> msIssueValueListText[Lookup[issue, "allowed", {}]] <>
        "；请修正后重新初始化。 Line " <>
        msIssueValueText[Lookup[issue, "line", Missing["Absent"]]] <>
        " has the unsupported type " <>
        msIssueValueText[Lookup[issue, "type", Missing["Absent"]]] <>
        "; only " <> msIssueValueListText[Lookup[issue, "allowed", {}]] <>
        " are allowed. Correct it and initialize again.",
    "lineMustHaveOneOrTwoEndpoints",
      "第 " <> msIssueValueText[Lookup[issue, "line", Missing["Absent"]]] <>
        " 条线的 endpoints 含有 " <>
        msIssueValueText[Lookup[issue, "actual", Missing["Absent"]]] <>
        " 个端点，只允许 1 个（外腿）或 2 个（内线）；请修正后重新初始化。 Line " <>
        msIssueValueText[Lookup[issue, "line", Missing["Absent"]]] <> " has " <>
        msIssueValueText[Lookup[issue, "actual", Missing["Absent"]]] <>
        " endpoints; only one (external leg) or two (internal line) are allowed. Correct it and initialize again.",
    "unknownEndpoint",
      "第 " <> msIssueValueText[Lookup[issue, "line", Missing["Absent"]]] <>
        " 条线的端点 " <> msIssueValueListText[Lookup[issue, "endpoints", {}]] <>
        " 中含有未定义的顶点 id；请定义这些顶点或修正端点后重新初始化。 Line " <>
        msIssueValueText[Lookup[issue, "line", Missing["Absent"]]] <>
        " references undefined vertex ids in its endpoints " <>
        msIssueValueListText[Lookup[issue, "endpoints", {}]] <>
        "; define these vertices or correct the endpoints and initialize again.",
    "missingNuMagnitude",
      "第 " <> msIssueValueText[Lookup[issue, "line", Missing["Absent"]]] <>
        " 条 massive 线缺少 nu（Hankel 阶）；请补充该字段后重新初始化。 Massive line " <>
        msIssueValueText[Lookup[issue, "line", Missing["Absent"]]] <>
        " is missing its nu (Hankel order); add this field and initialize again.",
    "missingMomentum",
      "第 " <> msIssueValueText[Lookup[issue, "line", Missing["Absent"]]] <>
        " 条线缺少 momentum；请补充该字段后重新初始化。 Line " <>
        msIssueValueText[Lookup[issue, "line", Missing["Absent"]]] <>
        " is missing its momentum; add this field and initialize again.",
    "invalidHankelBranches",
      "第 " <> msIssueValueText[Lookup[issue, "line", Missing["Absent"]]] <>
        " 条 massive 线的 hankelBranches 值 " <>
        msIssueValueText[Lookup[issue, "value", Missing["Absent"]]] <>
        " 无效；massive 外腿线需要 1 个、内线需要 2 个取值 1 或 2 的分支指标。请修正后重新初始化。 Line " <>
        msIssueValueText[Lookup[issue, "line", Missing["Absent"]]] <>
        " has the invalid hankelBranches " <>
        msIssueValueText[Lookup[issue, "value", Missing["Absent"]]] <>
        "; a massive external line needs one and a massive internal line needs two branch labels, each equal to 1 or 2. Correct it and initialize again.",
    "unknownGraphMode",
      "graphMode 的值 " <> msIssueValueText[Lookup[issue, "value", Missing["Absent"]]] <>
        " 不受支持，只允许 \"Tree\" 或 \"TimeOnly\"；请修正后重新初始化。 The graphMode value " <>
        msIssueValueText[Lookup[issue, "value", Missing["Absent"]]] <>
        " is not supported; only \"Tree\" or \"TimeOnly\" are allowed. Correct it and initialize again.",
    "internalGraphMustBeConnected",
      "内线构成的图必须连通，当前内线数为 " <>
        msIssueValueText[Lookup[issue, "edgeCount", Missing["Absent"]]] <>
        "；请检查各内线的端点后重新初始化。 The graph formed by the internal lines must be connected; the current internal-line count is " <>
        msIssueValueText[Lookup[issue, "edgeCount", Missing["Absent"]]] <>
        ". Check the internal-line endpoints and initialize again.",
    "internalGraphMustBeATree",
      "在 \"Tree\" 模式下内线图必须是树：当前有 " <>
        msIssueValueText[Lookup[issue, "vertexCount", Missing["Absent"]]] <>
        " 个顶点和 " <>
        msIssueValueText[Lookup[issue, "edgeCount", Missing["Absent"]]] <>
        " 条内线，而树要求内线数恰为顶点数减一。请修正拓扑后重新初始化。 In \"Tree\" mode the internal-line graph must be a tree: there are " <>
        msIssueValueText[Lookup[issue, "vertexCount", Missing["Absent"]]] <>
        " vertices and " <>
        msIssueValueText[Lookup[issue, "edgeCount", Missing["Absent"]]] <>
        " internal lines, while a tree requires exactly vertexCount-1 internal lines. Correct the topology and initialize again.",
    "MalformedContactBundles" | "InvalidContactBundleLines" |
    "ContactBundleMustShareThetaArgument",
      defaults = msFailureDefaultText[code];
      Last[defaults] <> " " <> First[defaults],
    "unknownVertexFamilyFields",
      "输入含有当前模式（" <> msIssueValueText[Lookup[issue, "inputMode", Missing["Absent"]]] <>
        "）不支持的字段 " <> msIssueValueListText[Lookup[issue, "fields", {}]] <>
        "；请删除这些字段或改用对应模式后重新初始化。 Input contains fields " <>
        msIssueValueListText[Lookup[issue, "fields", {}]] <> " not allowed in the " <>
        msIssueValueText[Lookup[issue, "inputMode", Missing["Absent"]]] <>
        " mode; remove them or switch mode and initialize again.",
    "unknownVertexFamilyBlockFields",
      allowedBlockKeys = Switch[
        Lookup[issue, "blockKind", None],
        "h", $msVertexFamilyHBlockKeys,
        "exponential", $msVertexFamilyExponentialBlockKeys,
        _, {}
      ];
      "第 " <> msIssueValueText[Lookup[issue, "position", Missing["Absent"]]] <> " 个 " <>
        msIssueValueText[Lookup[issue, "blockKind", Missing["Absent"]]] <>
        " 块含有未知字段 " <> msIssueValueListText[Lookup[issue, "fields", {}]] <>
        "；允许的字段为 " <> msIssueValueListText[allowedBlockKeys] <>
        "。请修正后重新初始化。 Block " <>
        msIssueValueText[Lookup[issue, "position", Missing["Absent"]]] <> " of kind " <>
        msIssueValueText[Lookup[issue, "blockKind", Missing["Absent"]]] <>
        " contains unknown fields " <> msIssueValueListText[Lookup[issue, "fields", {}]] <>
        "; the allowed fields are " <> msIssueValueListText[allowedBlockKeys] <>
        ". Correct it and initialize again.",
    _,
      "输入未通过 MadStree 的模式检查；请查看返回 Failure 中 issues 字段列出的详情，修正输入后重新初始化。 The input failed a MadStree schema check; see the issues field of the returned Failure for details, correct the input, and initialize again."
  ]
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
