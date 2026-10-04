(* ::Package:: *)
(* 文件用途：FlintNDE 的标准 Wolfram Language 程序包接口。
   功能范围：构造通用单变量有理矩阵微分方程、从原始点生成一次路径计划，并直接执行
   已有计划。Python bridge 负责 exact Q(i) 奇点发现、多项式加简单极点结构认证和
   通用有理矩阵输运；本文件不复制数值算法，也不提供一体化规划兼执行入口。 *)

BeginPackage["FlintNDE`"];


$FlintNDEVersion::usage = "$FlintNDEVersion is the loaded FlintNDE package version.";
FlintNDERationalSystem::usage =
  "FlintNDERationalSystem[matrix,x] constructs an exact rational matrix system. FlintNDE discovers and classifies its singularities internally.";
FlintNDEPartialFractionSystem::usage =
  "FlintNDEPartialFractionSystem[{P0,P1,...},residues,poles] constructs A(x)=Sum Pk x^k+Sum Rj/(x-pj); a constant polynomial is written as {P0}.";
FlintNDEPlanPath::usage =
  "FlintNDEPlanPath[system,start,points] plans a one-variable multipoint jump path from raw user points. WorkingPrecisionDigits defaults to 200. The default singularity mode is Avoid; SingularityJump must be selected explicitly when the path may cross a singularity.";
FlintNDEExecutePath::usage =
  "FlintNDEExecutePath[system,initialVector,plan] executes a plan returned by FlintNDEPlanPath without replanning. WorkingPrecisionDigits defaults to 200.";
FlintNDEEvaluateEpBatch::usage =
  "FlintNDEEvaluateEpBatch[jobs] plans and executes independent fixed-ep jobs in a bounded Python process pool. Each job contains ep, system, start, points and initialVector. WorkingPrecisionDigits->200 and ParallelTaskCount->12 are the defaults; the effective count is Min[ParallelTaskCount,Length[jobs]], and queued jobs start automatically as workers finish.";
ParallelTaskCount::usage =
  "ParallelTaskCount specifies the maximum number of independent ep tasks run concurrently. The default is 12; it is distinct from python-flint ctx.threads.";
MessageLanguage::usage =
  "MessageLanguage selects the language of runtime notices and diagnostics; the default is \"EN\" and \"CN\" selects Chinese.";
SingularityMode::usage =
  "SingularityMode selects path treatment: \"Avoid\" (default) refuses a segment crossing a singularity, while \"SingularityJump\" explicitly permits a singularity jump whose multivalued branch must be confirmed by the user.";

FlintNDEBridgeError::usage = "FlintNDEBridgeError represents a Python bridge failure.";

FlintNDEBridgeError::error = "FlintNDE bridge failure: `1`";
FlintNDEBridgeError::pythonFlintMissing = "Python started but cannot import python-flint: `1`";
FlintNDEBridgeError::launchFailed = "The Python bridge could not be launched: `1`";
FlintNDEBridgeError::outputMissing = "The Python bridge exited without creating its output file: `1`";


Begin["`Private`"];


(* ::Chapter:: *)
(*程序包位置与通用参数*)

$FlintNDEVersion = "1.0";
$FlintNDEMathematicaDirectory = DirectoryName[$InputFileName];
$FlintNDEVersionDirectory = DirectoryName[$FlintNDEMathematicaDirectory];
$FlintNDEPythonModule = "flintnde.mathematica_bridge";
$FlintNDERequestSchema = "flintnde_mathematica_request_v1";


flintNDEResolvePython[Automatic] := "python";
flintNDEResolvePython[command_String] := command;
flintNDEResolvePython[other_] := Failure[
  "InvalidPythonExecutable",
  <|"value" -> other|>
];


flintNDEResolveWorkDirectory[Automatic] := FileNameJoin[{Directory[], "results_temp"}];
flintNDEResolveWorkDirectory[path_String] := ExpandFileName[path];
flintNDEResolveWorkDirectory[other_] := Failure[
  "InvalidWorkDirectory",
  <|"value" -> other|>
];


(* 短 token 只负责运行文件去重；独立 helper 允许测试在不启动 Python 的情况下固定写入目标。 *)
flintNDECreateToken[] := StringTake[StringDelete[CreateUUID[], "-"], 16];


(* ::Section:: *)
(*Windows 安全路径门禁*)

(* 路径过长在任何目录创建或 Python 启动前独立返回，不与 bridge 错误混合。 *)
flintNDEPathLengthFailure[paths_List] := Module[{maximum = 259, strings, overlong, longest},
  If[$OperatingSystem =!= "Windows", Return[None]];
  strings = Select[paths, StringQ];
  overlong = Select[strings, StringLength[#] > maximum &];
  If[overlong === {}, Return[None]];
  longest = First@MaximalBy[overlong, StringLength];
  Failure["RuntimePathTooLong", <|
    "path" -> longest,
    "pathLength" -> StringLength[longest],
    "safeMaximum" -> maximum,
    "suggestion" ->
      "请通过 WorkDirectory 指定更短的临时运行目录；该错误发生在 Python 启动之前。"
  |>]
];


flintNDENormalizeLanguage[value_] := If[
  MemberQ[{"EN", "CN"}, value],
  value,
  Failure[
    "InvalidMessageLanguage",
    <|"value" -> value, "acceptedValues" -> {"EN", "CN"}|>
  ]
];


(* 公开入口先比较调用参数与当前 Options，任何额外规则都结构化拒绝。 *)
flintNDEUnknownOptionNames[rawOptions_List, allowedOptions_List] := Complement[
  First /@ rawOptions,
  First /@ allowedOptions
];


flintNDENormalizeMode["Avoid"] := "avoid";
flintNDENormalizeMode["SingularityJump"] := "singularity_jump";
flintNDENormalizeMode[value_String] := Failure[
  "InvalidSingularityMode",
  <|"value" -> value, "acceptedValues" -> {"Avoid", "SingularityJump"}|>
];
flintNDENormalizeMode[other_] := Failure[
  "InvalidSingularityMode",
  <|"value" -> other|>
];


(* ::Chapter:: *)
(*精确系统构造*)

(* 实现思路：Wolfram 侧只把每个有理函数化为升幂 numerator/denominator 系数；
   是否属于多项式加简单极点快速情形由 Python 内部重构恒等式认证。 *)
flintNDERationalFunctionRecord[expression_, variable_Symbol] := Module[
  {combined, numerator, denominator, coefficients},
  combined = Together[expression];
  numerator = Numerator[combined];
  denominator = Denominator[combined];
  If[! PolynomialQ[numerator, variable] || ! PolynomialQ[denominator, variable],
    Return[Failure[
      "NonRationalMatrixEntry",
      <|"entry" -> expression, "variable" -> HoldForm[variable]|>
    ]]
  ];
  coefficients = Join[
    CoefficientList[numerator, variable],
    CoefficientList[denominator, variable]
  ];
  If[! FreeQ[coefficients, _Real],
    Return[Failure[
      "InexactRationalMatrixEntry",
      <|"entry" -> expression|>
    ]]
  ];
  <|
    "numerator" -> If[TrueQ[numerator === 0], {0}, CoefficientList[numerator, variable]],
    "denominator" -> CoefficientList[denominator, variable]
  |>
];


Options[FlintNDERationalSystem] = {"Name" -> "Mathematica-rational-matrix-system"};

FlintNDERationalSystem[
  matrix_?MatrixQ,
  variable_Symbol,
  OptionsPattern[]
] := Module[{dimension, records, failure},
  dimension = Length[matrix];
  If[dimension == 0 || ! AllTrue[matrix, Length[#] == dimension &],
    Return[Failure["SquareMatrixRequired", <|"dimensions" -> Dimensions[matrix]|>]]
  ];
  records = Map[flintNDERationalFunctionRecord[#, variable] &, matrix, {2}];
  failure = FirstCase[records, _Failure, None, Infinity];
  If[failure =!= None, Return[failure]];
  <|
    "type" -> "rationalMatrix",
    "variable" -> SymbolName[Unevaluated[variable]],
    "name" -> OptionValue["Name"],
    "matrix" -> records
  |>
];

FlintNDERationalSystem[___] := Failure[
  "InvalidRationalSystemArguments",
  <|"usage" -> "FlintNDERationalSystem[matrix, variable]"|>
];


flintNDEPolynomialMatricesQ[matrices_List] :=
  matrices =!= {} && AllTrue[matrices, MatrixQ];

flintNDEPartialFractionRecord[
  polynomialCoefficients_List,
  residues_List,
  poles_List
] := Module[{dimension, dimensions},
  If[! flintNDEPolynomialMatricesQ[polynomialCoefficients],
    Return[Failure["PolynomialMatrixListRequired", <||>]]
  ];
  dimension = Length[First[polynomialCoefficients]];
  dimensions = Dimensions /@ Join[polynomialCoefficients, residues];
  If[
    dimension == 0 ||
    Dimensions[First[polynomialCoefficients]] =!= {dimension, dimension} ||
    ! AllTrue[dimensions, # === {dimension, dimension} &],
    Return[Failure[
      "ConsistentSquareMatricesRequired",
      <|"dimensions" -> dimensions|>
    ]]
  ];
  If[Length[residues] =!= Length[poles],
    Return[Failure[
      "PoleResidueLengthMismatch",
      <|"residueCount" -> Length[residues], "poleCount" -> Length[poles]|>
    ]]
  ];
  <|
    "type" -> "partialFraction",
    "polynomialCoefficients" -> polynomialCoefficients,
    "residues" -> residues,
    "poles" -> poles
  |>
];


FlintNDEPartialFractionSystem[
  polynomialCoefficients_List /; flintNDEPolynomialMatricesQ[polynomialCoefficients],
  residues_List,
  poles_List
] := flintNDEPartialFractionRecord[polynomialCoefficients, residues, poles];


FlintNDEPartialFractionSystem[___] := Failure[
  "InvalidPartialFractionSystemArguments",
  <|"usage" -> "FlintNDEPartialFractionSystem[{P0,P1,...}, residues, poles]"|>
];


(* ::Chapter:: *)
(*JSON 编码与 Python bridge*)

(* 实现思路：所有 exact Rational/Complex 和任意精度 Real 都先写成字符串，
   避免 JSON 导出降为 machine precision；Association 递归保持 schema 结构。 *)
flintNDEEncode[expression_Association] :=
  Association[flintNDEEncode /@ Normal[expression]];
flintNDEEncode[rule_Rule] := rule[[1]] -> flintNDEEncode[rule[[2]]];
flintNDEEncode[expression_List] := flintNDEEncode /@ expression;
flintNDEEncode[value_Complex] := <|
  "re" -> ToString[Re[value], InputForm],
  "im" -> ToString[Im[value], InputForm]
|>;
flintNDEEncode[value_Rational] :=
  ToString[Numerator[value]] <> "/" <> ToString[Denominator[value]];
flintNDEEncode[value_Real] := ToString[value, InputForm];
flintNDEEncode[value_] := value;


flintNDEInvoke[
  request_Association,
  pythonOption_,
  workDirectoryOption_
] := Module[
  {python, workDirectory, bridgeDirectory, requestFile, outputFile, logFile, token, command,
    processResult, exitCode, standardOutput, standardError, logText, process, result,
    launchFailure, pathFailure, requestWrite, outputLoad},
  python = flintNDEResolvePython[pythonOption];
  If[Head[python] === Failure, Return[python]];
  workDirectory = flintNDEResolveWorkDirectory[workDirectoryOption];
  If[Head[workDirectory] === Failure, Return[workDirectory]];
  bridgeDirectory = FileNameJoin[{workDirectory, "bridge"}];
  token = flintNDECreateToken[];
  requestFile = FileNameJoin[{bridgeDirectory, "i-" <> token <> ".json"}];
  outputFile = FileNameJoin[{bridgeDirectory, "o-" <> token <> ".m"}];
  logFile = FileNameJoin[{bridgeDirectory, "l-" <> token <> ".txt"}];
  pathFailure = flintNDEPathLengthFailure[{
    workDirectory, bridgeDirectory, requestFile, outputFile, logFile
  }];
  If[Head[pathFailure] === Failure, Return[pathFailure]];
  If[! DirectoryQ[bridgeDirectory],
    Quiet@Check[
      CreateDirectory[bridgeDirectory, CreateIntermediateDirectories -> True],
      Return[Failure["WorkDirectoryCreationFailed", <|"path" -> bridgeDirectory|>]]
    ]
  ];
  requestWrite = Quiet@Check[Export[requestFile, flintNDEEncode[request], "JSON"], $Failed];
  If[requestWrite === $Failed || ! FileExistsQ[requestFile],
    Return[Failure["RuntimeInputWriteFailed", <|
      "path" -> requestFile, "pathLength" -> StringLength[requestFile]
    |>]]
  ];
  command = {
    python,
    "-m",
    $FlintNDEPythonModule,
    requestFile,
    outputFile
  };
  (* 参数列表启动不经过 shell；连续加载 FLINT DLL 时不保留旧 Run/重定向入口。 *)
  launchFailure = False;
  processResult = Quiet@Check[
    RunProcess[command, All, ProcessDirectory -> $FlintNDEVersionDirectory],
    launchFailure = True;
    $Failed
  ];
  If[AssociationQ[processResult],
    exitCode = Lookup[processResult, "ExitCode", $Failed];
    standardOutput = Lookup[processResult, "StandardOutput", ""];
    standardError = Lookup[processResult, "StandardError", ""],
    launchFailure = True;
    exitCode = $Failed;
    standardOutput = "";
    standardError = ""
  ];
  logText = StringRiffle[Select[{standardOutput, standardError}, StringLength[#] > 0 &], "\n"];
  Quiet@Check[Export[logFile, logText, "Text", CharacterEncoding -> "UTF-8"], Null];
  process = <|
    "ExitCode" -> exitCode,
    "Command" -> command,
    "WorkingDirectory" -> $FlintNDEVersionDirectory,
    "StandardOutput" -> standardOutput,
    "StandardError" -> standardError
  |>;
  If[launchFailure,
    Return[Failure["BridgeLaunchFailed", <|"process" -> process|>]]
  ];
  If[! FileExistsQ[outputFile],
    If[StringContainsQ[logText, Alternatives[
        "No module named 'flint'", "No module named \"flint\"",
        "ModuleNotFoundError: No module named 'flint'"
      ]],
      Message[FlintNDEBridgeError::pythonFlintMissing, logText];
      Return[Failure["PythonFlintUnavailable", <|
        "process" -> process, "requestFile" -> requestFile
      |>]]
    ];
    If[exitCode =!= 0,
      Message[FlintNDEBridgeError::launchFailed, logText];
      Return[Failure["BridgeLaunchFailed", <|
        "process" -> process, "requestFile" -> requestFile
      |>]]
    ];
    Message[FlintNDEBridgeError::outputMissing, logText];
    Return[Failure["BridgeOutputMissing", <|
      "process" -> process, "requestFile" -> requestFile
    |>]]
  ];
  Clear[Global`FlintNDEBridgeResult];
  outputLoad = Quiet@Check[Get[outputFile, CharacterEncoding -> "UTF-8"], $Failed];
  If[outputLoad === $Failed,
    Return[Failure["BridgeOutputInvalid", <|"path" -> outputFile, "process" -> process|>]]
  ];
  result = Global`FlintNDEBridgeResult;
  Quiet[DeleteFile /@ Select[{requestFile, outputFile, logFile}, FileExistsQ]];
  Clear[Global`FlintNDEBridgeResult];
  If[! AssociationQ[result],
    Return[Failure["InvalidBridgeResult", <|"result" -> result|>]]
  ];
  If[Lookup[result, "status", None] === "error",
    Message[FlintNDEBridgeError::error, Lookup[result, "message", "unknown bridge failure"]];
    Return[Failure[
      "BridgeFailure",
      <|"message" -> Lookup[result, "message", "unknown bridge failure"],
        "process" -> process|>
    ]]
  ];
  result
];


(* ::Chapter:: *)
(*两阶段路径规划*)

Options[FlintNDEPlanPath] = {
  "Python" -> Automatic,
  "WorkDirectory" -> Automatic,
  "WorkingPrecisionDigits" -> 200,
  "OutputDigits" -> 40,
  MessageLanguage -> "EN",
  SingularityMode -> "Avoid",
  "RadiusFraction" -> 0.60,
  "MaxStepOverRadius" -> 0.45,
  "SingularityJumpThreshold" -> 0.5,
  "MatchFraction" -> 0.6,
  "MaxSingularityJumps" -> 16
};


flintNDEPlanningNotice[result_Association, mode_String, language_String] := Module[
  {status = Lookup[result, "status", "error"], text},
  text = Which[
    status === "singularPathRefused",
      Lookup[result, "message", "singular path refused"],
    language === "CN" && mode === "singularity_jump",
      "FlintNDEPlanPath：已按输入的原始点完成路径规划；当前显式使用奇点折跃。多值分支等价于某一绕行路径，必须由用户确认。",
    language === "CN",
      "FlintNDEPlanPath：已按输入的原始点完成路径规划；当前使用避开奇点模式（缺省）。请把返回计划交给 FlintNDEExecutePath；执行时不会再次规划。",
    mode === "singularity_jump",
      "FlintNDEPlanPath: the supplied raw points were planned in explicit singularity-jump mode. The selected multivalued branch is equivalent to a detour path and must be confirmed by the user.",
    True,
      "FlintNDEPlanPath: the supplied raw points were planned in avoid-singularity mode (default). Pass the returned plan to FlintNDEExecutePath; execution does not replan."
  ];
  Print[text]
];


FlintNDEPlanPath[
  system_Association,
  start_,
  points_List,
  opts : OptionsPattern[]
] := Module[{language, mode, request, result, unknownOptions},
  unknownOptions = flintNDEUnknownOptionNames[{opts}, Options[FlintNDEPlanPath]];
  If[unknownOptions =!= {},
    Return[Failure[
      "UnknownOption",
      <|"function" -> "FlintNDEPlanPath", "options" -> unknownOptions|>
    ]]
  ];
  language = flintNDENormalizeLanguage[OptionValue[MessageLanguage]];
  If[Head[language] === Failure, Return[language]];
  mode = flintNDENormalizeMode[OptionValue[SingularityMode]];
  If[Head[mode] === Failure, Return[mode]];
  request = <|
    "schema" -> $FlintNDERequestSchema,
    "action" -> "plan",
    "system" -> system,
    "start" -> start,
    "points" -> points,
    "workingPrecisionDigits" -> OptionValue["WorkingPrecisionDigits"],
    "outputDigits" -> OptionValue["OutputDigits"],
    "messageLanguage" -> language,
    "singularityMode" -> mode,
    "radiusFraction" -> OptionValue["RadiusFraction"],
    "maxStepOverRadius" -> OptionValue["MaxStepOverRadius"],
    "singularityJumpThreshold" -> OptionValue["SingularityJumpThreshold"],
    "matchFraction" -> OptionValue["MatchFraction"],
    "maxSingularityJumps" -> OptionValue["MaxSingularityJumps"]
  |>;
  result = flintNDEInvoke[
    request,
    OptionValue["Python"],
    OptionValue["WorkDirectory"]
  ];
  If[AssociationQ[result], flintNDEPlanningNotice[result, mode, language]];
  result
];

FlintNDEPlanPath[___] := Failure[
  "InvalidPlanArguments",
  <|"usage" -> "FlintNDEPlanPath[system, start, {point1,...}]"|>
];


(* ::Chapter:: *)
(*已有计划直接执行*)

Options[FlintNDEExecutePath] = {
  "Python" -> Automatic,
  "WorkDirectory" -> Automatic,
  "WorkingPrecisionDigits" -> 200,
  "OutputDigits" -> 40,
  "PrimaryOrder" -> 40,
  "ReferenceOrder" -> 48,
  "TargetRelativeError" -> "1e-30",
  "CertificationMode" -> "embedded",
  "RadiusFraction" -> 0.60,
  MessageLanguage -> "EN"
};


(* 后端把零中点写成字符串 "0"，而 Wolfram 的 "`digits" 精度标记对零不成立：
   ToExpression["0.0`40"] 得到机器精度零，I 乘该零又会把同一复数的高精度实部一起拖到
   MachinePrecision。因此这里只按字面识别全零尾数并保留精确零，不设阈值容差，
   真正非零的小分量仍按 digits 位任意精度解码。 *)
flintNDEDecodeDecimal[text_String, digits_Integer] := Module[
  {parts, mantissa, exponent},
  parts = StringSplit[ToLowerCase[StringTrim[text]], "e", 2];
  mantissa = First[parts];
  If[mantissa =!= "" &&
    StringFreeQ[StringDelete[mantissa, {"-", "+", "."}], Except["0"]],
    Return[0]
    ];
  exponent = If[Length[parts] === 2, "*^" <> Last[parts], ""];
  If[! StringContainsQ[mantissa, "."], mantissa = mantissa <> ".0"];
  ToExpression[mantissa <> "`" <> ToString[digits] <> exponent]
];

flintNDEDecodeComplex[record_Association, digits_Integer] :=
  flintNDEDecodeDecimal[record["real"], digits] +
  I flintNDEDecodeDecimal[record["imag"], digits];

flintNDEDecodeVector[records_List, digits_Integer] :=
  flintNDEDecodeComplex[#, digits] & /@ records;


flintNDEDecodeSingularValue[record_Association, digits_Integer] := If[
  KeyExistsQ[record, "text"],
  record["text"],
  flintNDEDecodeComplex[record, digits]
];


flintNDEDecodeExecutionResult[result_Association, digits_Integer] := Module[
  {decoded = result, samples, singularTargets},
  If[KeyExistsQ[result, "primaryFinalVector"],
    decoded = Join[decoded, <|
      "primaryFinalVectorRecords" -> result["primaryFinalVector"],
      "primaryFinalVector" -> flintNDEDecodeVector[
        result["primaryFinalVector"], digits
      ]
    |>]
  ];
  If[KeyExistsQ[result, "referenceFinalVector"],
    decoded = Join[decoded, <|
      "referenceFinalVectorRecords" -> result["referenceFinalVector"],
      "referenceFinalVector" -> flintNDEDecodeVector[
        result["referenceFinalVector"], digits
      ]
    |>]
  ];
  If[KeyExistsQ[result, "samplePoints"],
    samples = Map[
      Join[#, <|
        "valueRecords" -> #["value"],
        "value" -> flintNDEDecodeVector[#["value"], digits]
      |>] &,
      result["samplePoints"]
    ];
    decoded = Join[decoded, <|"samplePoints" -> samples|>]
  ];
  If[KeyExistsQ[result, "singularTargets"],
    singularTargets = Map[
      Join[#, <|
        "valueRecords" -> #["values"],
        "values" -> (flintNDEDecodeSingularValue[#, digits] & /@ #["values"])
      |>] &,
      result["singularTargets"]
    ];
    decoded = Join[decoded, <|"singularTargets" -> singularTargets|>]
  ];
  decoded
];


FlintNDEExecutePath[
  system_Association,
  initialVector_List,
  plan_Association,
  opts : OptionsPattern[]
] := Module[
  {planRecord, request, result, digits, language, planningDigits,
   workingDigits, notice, unknownOptions},
  unknownOptions = flintNDEUnknownOptionNames[{opts}, Options[FlintNDEExecutePath]];
  If[unknownOptions =!= {},
    Return[Failure[
      "UnknownOption",
      <|"function" -> "FlintNDEExecutePath", "options" -> unknownOptions|>
    ]]
  ];
  language = flintNDENormalizeLanguage[OptionValue[MessageLanguage]];
  If[Head[language] === Failure, Return[language]];
  If[
    Lookup[plan, "schema", None] =!= "flintnde_mathematica_bridge_v1" ||
    Lookup[plan, "status", None] =!= "complete" ||
    Lookup[plan, "operation", None] =!= "plan" ||
    ! KeyExistsQ[plan, "plan"],
    Return[Failure[
      "InvalidExecutionPlan",
      <|"reason" -> "expected the complete result returned by FlintNDEPlanPath"|>
    ]]
  ];
  planRecord = plan["plan"];
  planningDigits = Lookup[
    planRecord, "planningPrecisionDigits", Missing["Absent"]
  ];
  workingDigits = OptionValue["WorkingPrecisionDigits"];
  If[! IntegerQ[planningDigits],
    Return[Failure[
      "PlannedPathPrecisionMissing",
      <|"requestedPrecisionDigits" -> workingDigits,
        "reason" -> "the plan does not record its planning precision; replan"|>
    ]]
  ];
  If[TrueQ[workingDigits > planningDigits],
    notice = If[
      language === "CN",
      "执行请求 " <> ToString[workingDigits] <> " 位十进制精度，但该路径只按 " <>
        ToString[planningDigits] <>
        " 位规划。已序列化节点不能补回精度；请按所需精度重新运行 FlintNDEPlanPath。",
      "Execution requests " <> ToString[workingDigits] <>
        " decimal digits, but this path was planned at " <>
        ToString[planningDigits] <>
        ". Serialized nodes cannot gain precision; rerun FlintNDEPlanPath at the requested precision."
    ];
    Print[notice];
    Return[Failure[
      "PlannedPathPrecisionInsufficient",
      <|
        "message" -> notice,
        "messageLanguage" -> language,
        "planningPrecisionDigits" -> planningDigits,
        "requestedPrecisionDigits" -> workingDigits
      |>
    ]]
  ];
  digits = OptionValue["OutputDigits"];
  request = <|
    "schema" -> $FlintNDERequestSchema,
    "action" -> "execute",
    "system" -> system,
    "initialVector" -> initialVector,
    "plannedResult" -> plan,
    "workingPrecisionDigits" -> workingDigits,
    "outputDigits" -> digits,
    "primaryOrder" -> OptionValue["PrimaryOrder"],
    "referenceOrder" -> OptionValue["ReferenceOrder"],
    "targetRelativeError" -> ToString[OptionValue["TargetRelativeError"]],
    "certificationMode" -> OptionValue["CertificationMode"],
    "radiusFraction" -> OptionValue["RadiusFraction"],
    "messageLanguage" -> language
  |>;
  result = flintNDEInvoke[
    request,
    OptionValue["Python"],
    OptionValue["WorkDirectory"]
  ];
  If[! AssociationQ[result], Return[result]];
  Print[If[
    language === "CN",
    "FlintNDEExecutePath：直接执行输入的已有计划；未再次规划路径。",
    "FlintNDEExecutePath: executed the supplied plan directly; no path replanning was performed."
  ]];
  flintNDEDecodeExecutionResult[result, digits]
];

FlintNDEExecutePath[___] := Failure[
  "InvalidExecuteArguments",
  <|"usage" -> "FlintNDEExecutePath[system, initialVector, plan]"|>
];


(* ::Chapter:: *)
(*固定 ep 任务的有界并行*)

Options[FlintNDEEvaluateEpBatch] = {
  ParallelTaskCount -> 12,
  "Python" -> Automatic,
  "WorkDirectory" -> Automatic,
  "WorkingPrecisionDigits" -> 200,
  "OutputDigits" -> 40,
  "PrimaryOrder" -> 40,
  "ReferenceOrder" -> 48,
  "TargetRelativeError" -> "1e-30",
  "CertificationMode" -> "embedded",
  MessageLanguage -> "EN",
  SingularityMode -> "Avoid",
  "RadiusFraction" -> 0.60,
  "MaxStepOverRadius" -> 0.45,
  "SingularityJumpThreshold" -> 0.5,
  "MatchFraction" -> 0.6,
  "MaxSingularityJumps" -> 16
};


flintNDEEpJobRequest[job_Association, options_Association] := Module[
  {required, missing, mode},
  required = {"ep", "system", "start", "points", "initialVector"};
  missing = Select[required, ! KeyExistsQ[job, #] &];
  If[missing =!= {},
    Return[Failure["InvalidEpJob", <|"missingKeys" -> missing, "job" -> job|>]]
  ];
  mode = flintNDENormalizeMode[options[SingularityMode]];
  If[Head[mode] === Failure, Return[mode]];
  <|
    "schema" -> $FlintNDERequestSchema,
    "action" -> "evaluate",
    "ep" -> job["ep"],
    "system" -> job["system"],
    "start" -> job["start"],
    "points" -> job["points"],
    "initialVector" -> job["initialVector"],
    "workingPrecisionDigits" -> options["WorkingPrecisionDigits"],
    "outputDigits" -> options["OutputDigits"],
    "primaryOrder" -> options["PrimaryOrder"],
    "referenceOrder" -> options["ReferenceOrder"],
    "targetRelativeError" -> ToString[options["TargetRelativeError"]],
    "certificationMode" -> options["CertificationMode"],
    "messageLanguage" -> options[MessageLanguage],
    "singularityMode" -> mode,
    "radiusFraction" -> options["RadiusFraction"],
    "maxStepOverRadius" -> options["MaxStepOverRadius"],
    "singularityJumpThreshold" -> options["SingularityJumpThreshold"],
    "matchFraction" -> options["MatchFraction"],
    "maxSingularityJumps" -> options["MaxSingularityJumps"]
  |>
];


FlintNDEEvaluateEpBatch[jobs_List, opts : OptionsPattern[]] := Module[
  {unknownOptions, parallelCount, language, optionValues, requests, failure,
   result, digits, decodedResults},
  unknownOptions = flintNDEUnknownOptionNames[{opts}, Options[FlintNDEEvaluateEpBatch]];
  If[unknownOptions =!= {},
    Return[Failure["UnknownOption", <|"function" -> "FlintNDEEvaluateEpBatch",
      "options" -> unknownOptions|>]]
  ];
  If[jobs === {}, Return[Failure["EpJobListEmpty", <||>]]];
  parallelCount = OptionValue[ParallelTaskCount];
  If[! IntegerQ[parallelCount] || parallelCount < 1,
    Return[Failure["ParallelTaskCountPositiveIntegerRequired",
      <|"value" -> parallelCount, "default" -> 12|>]]
  ];
  language = flintNDENormalizeLanguage[OptionValue[MessageLanguage]];
  If[Head[language] === Failure, Return[language]];
  optionValues = Association[Options[FlintNDEEvaluateEpBatch]];
  optionValues = Join[optionValues, Association[{opts}], <|MessageLanguage -> language|>];
  requests = flintNDEEpJobRequest[#, optionValues] & /@ jobs;
  failure = FirstCase[requests, _Failure, None];
  If[failure =!= None, Return[failure]];
  Print[If[language === "CN",
    "FlintNDE：不同 ep 任务的缺省并行数为 12；本次请求 " <>
      ToString[parallelCount] <> "，实际并行数为 " <>
      ToString[Min[parallelCount, Length[jobs]]] <> "。任务完成后自动续交队列。",
    "FlintNDE: default ep-task parallelism is 12; requested " <>
      ToString[parallelCount] <> ", effective " <>
      ToString[Min[parallelCount, Length[jobs]]] <>
      ". Queued jobs start automatically as workers finish."
  ]];
  result = flintNDEInvoke[
    <|"schema" -> $FlintNDERequestSchema, "action" -> "ep_batch",
      "requests" -> requests, "parallelTaskCount" -> parallelCount|>,
    OptionValue["Python"], OptionValue["WorkDirectory"]
  ];
  If[! AssociationQ[result], Return[result]];
  digits = OptionValue["OutputDigits"];
  decodedResults = MapThread[
    Function[{item, job},
      If[AssociationQ[Lookup[item, "execution", None]],
        Join[item, <|"ep" -> job["ep"],
          "execution" -> flintNDEDecodeExecutionResult[
            item["execution"], digits]|>],
        Join[item, <|"ep" -> job["ep"]|>]
      ]
    ],
    {Lookup[result, "results", {}], jobs}
  ];
  Join[result, <|"results" -> decodedResults|>]
];


FlintNDEEvaluateEpBatch[___] := Failure[
  "InvalidEpBatchArguments",
  <|"usage" -> "FlintNDEEvaluateEpBatch[{job1,...}, ParallelTaskCount->12]"|>
];


End[];


EndPackage[];
