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

FlintNDEBridgeError::error =
  "FlintNDE 的 Python 后端返回失败，后端信息为：`1`。请按该信息修正输入后重新调用同一函数。 \
The FlintNDE Python backend returned a failure with the message: `1`. Correct the input accordingly and call the same function again.";
FlintNDEBridgeError::pythonFlintMissing =
  "Python 已启动但无法导入 python-flint，后端日志为：`1`。请在该环境中安装 python-flint，或用 Python 选项指向已装好它的解释器后重试。 \
Python started but could not import python-flint; the backend log is: `1`. Install python-flint in that environment, or point the Python option at an interpreter that already has it, then try again.";
FlintNDEBridgeError::launchFailed =
  "无法启动 FlintNDE 的 Python 后端，日志为：`1`。请检查 Python 选项指向的解释器是否存在且可执行，并确认 WorkDirectory 可写后重试。 \
The FlintNDE Python backend could not be launched; the log is: `1`. Check that the interpreter given by the Python option exists and is executable, and that WorkDirectory is writable, then try again.";
FlintNDEBridgeError::outputMissing =
  "FlintNDE 的 Python 后端已退出但没有生成输出文件，日志为：`1`。请检查 WorkDirectory 是否可写、磁盘是否有余量，并按日志报错修正后重试。 \
The FlintNDE Python backend exited without creating its output file; the log is: `1`. Check that WorkDirectory is writable and has free disk space, fix what the log reports, then try again.";


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
  <|"value" -> other,
    "MessageTemplate" ->
      "Python 选项只接受 Automatic 或解释器命令字符串，当前收到的是其它类型的无效值；请把 Python 改为 Automatic 或有效的命令字符串后重试。 \
The Python option accepts only Automatic or an interpreter command string, and the given value is of another, invalid kind; set Python to Automatic or a valid command string and try again."|>
];


flintNDEResolveWorkDirectory[Automatic] := FileNameJoin[{Directory[], "results_temp"}];
flintNDEResolveWorkDirectory[path_String] := ExpandFileName[path];
flintNDEResolveWorkDirectory[other_] := Failure[
  "InvalidWorkDirectory",
  <|"value" -> other,
    "MessageTemplate" ->
      "WorkDirectory 选项只接受 Automatic 或目录路径字符串，当前收到的是其它类型的无效值；请把 WorkDirectory 改为 Automatic 或有效路径后重试。 \
The WorkDirectory option accepts only Automatic or a directory path string, and the given value is of another, invalid kind; set WorkDirectory to Automatic or a valid path and try again."|>
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
      "请通过 WorkDirectory 指定更短的临时运行目录；该错误发生在 Python 启动之前。",
    "MessageTemplate" ->
      "某个运行时文件路径超过 Windows 的 259 字符安全上限，该检查发生在启动 Python 之前；请通过 WorkDirectory 指定更短的运行目录后重试。 \
A runtime file path exceeds the Windows 259-character safe maximum, and this check runs before Python is launched; specify a shorter runtime directory via WorkDirectory and try again."
  |>]
];


flintNDENormalizeLanguage[value_] := If[
  MemberQ[{"EN", "CN"}, value],
  value,
  Failure[
    "InvalidMessageLanguage",
    <|"value" -> value, "acceptedValues" -> {"EN", "CN"},
      "MessageTemplate" ->
        "MessageLanguage 只接受 \"EN\" 或 \"CN\"，当前值不在其中；请把 MessageLanguage 改为 \"EN\" 或 \"CN\" 后重试。 \
MessageLanguage accepts only \"EN\" or \"CN\", and the current value is neither; set MessageLanguage to \"EN\" or \"CN\" and try again."|>
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
  <|"value" -> value, "acceptedValues" -> {"Avoid", "SingularityJump"},
    "MessageTemplate" ->
      "SingularityMode 只接受字符串 \"Avoid\" 或 \"SingularityJump\"，当前字符串不在其中；请改用这两个模式之一后重试。 \
SingularityMode accepts only the string \"Avoid\" or \"SingularityJump\", and the given string is neither; use one of these two modes and try again."|>
];
flintNDENormalizeMode[other_] := Failure[
  "InvalidSingularityMode",
  <|"value" -> other,
    "MessageTemplate" ->
      "SingularityMode 必须是字符串 \"Avoid\" 或 \"SingularityJump\"，当前收到的是非字符串值；请改用这两个模式之一后重试。 \
SingularityMode must be the string \"Avoid\" or \"SingularityJump\", and the given value is not a string; use one of these two modes and try again."|>
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
      <|"entry" -> expression, "variable" -> HoldForm[variable],
        "MessageTemplate" ->
          "某个矩阵元素化简后分子或分母不是变量的多项式，因此不是有理函数；请把该元素改为变量的有理函数后重试。 \
After simplification, one matrix entry has a numerator or denominator that is not a polynomial in the variable, so it is not a rational function; replace that entry with a rational function of the variable and try again."|>
    ]]
  ];
  coefficients = Join[
    CoefficientList[numerator, variable],
    CoefficientList[denominator, variable]
  ];
  If[! FreeQ[coefficients, _Real],
    Return[Failure[
      "InexactRationalMatrixEntry",
      <|"entry" -> expression,
        "MessageTemplate" ->
          "某个矩阵元素的系数中含有近似实数，而精确系统只接受精确有理数；请把该元素改写为精确有理系数后重试。 \
A matrix entry has approximate real numbers among its coefficients, while the exact system accepts only exact rationals; rewrite that entry with exact rational coefficients and try again."|>
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
    Return[Failure["SquareMatrixRequired", <|"dimensions" -> Dimensions[matrix],
      "MessageTemplate" ->
        "输入矩阵必须是方阵，当前矩阵为空或各行长度不一致；请传入 n 行 n 列的方阵后重试。 \
The input matrix must be square, but it is empty or its rows have different lengths; pass an n-by-n square matrix and try again."|>]]
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
  <|"usage" -> "FlintNDERationalSystem[matrix, variable]",
    "MessageTemplate" ->
      "FlintNDERationalSystem 收到的参数不符合调用形态，它需要一个由变量的有理函数组成的方阵和一个变量符号；请按 FlintNDERationalSystem[matrix, variable] 传入有理函数方阵与变量后重试。 \
FlintNDERationalSystem received arguments that do not match its calling form, which needs a square matrix of rational functions in the variable together with a variable symbol; call it as FlintNDERationalSystem[matrix, variable] with a rational-function square matrix and a variable, then try again."|>
];


flintNDEPolynomialMatricesQ[matrices_List] :=
  matrices =!= {} && AllTrue[matrices, MatrixQ];

flintNDEPartialFractionRecord[
  polynomialCoefficients_List,
  residues_List,
  poles_List
] := Module[{dimension, dimensions},
  If[! flintNDEPolynomialMatricesQ[polynomialCoefficients],
    Return[Failure["PolynomialMatrixListRequired", <|
      "MessageTemplate" ->
        "第一个参数必须是由至少一个矩阵组成的多项式系数列表 {P0,P1,...}；常数项请写成 {P0} 形式的矩阵列表后重试。 \
The first argument must be a nonempty list {P0,P1,...} of polynomial coefficient matrices; write a constant term as a matrix list such as {P0} and try again."|>]]
  ];
  dimension = Length[First[polynomialCoefficients]];
  dimensions = Dimensions /@ Join[polynomialCoefficients, residues];
  If[
    dimension == 0 ||
    Dimensions[First[polynomialCoefficients]] =!= {dimension, dimension} ||
    ! AllTrue[dimensions, # === {dimension, dimension} &],
    Return[Failure[
      "ConsistentSquareMatricesRequired",
      <|"dimensions" -> dimensions,
        "MessageTemplate" ->
          "多项式系数矩阵与留数矩阵必须全部是同一维数的方阵，当前维度不一致或不是方阵；请把所有矩阵统一为同一 n 行 n 列维数后重试。 \
The polynomial coefficient matrices and residue matrices must all be square with one common dimension, but the current dimensions differ or are not square; make every matrix the same n-by-n size and try again."|>
    ]]
  ];
  If[Length[residues] =!= Length[poles],
    Return[Failure[
      "PoleResidueLengthMismatch",
      <|"residueCount" -> Length[residues], "poleCount" -> Length[poles],
        "MessageTemplate" ->
          "留数矩阵个数与极点个数不相等，二者无法一一配对；请把 residues 与 poles 调整为相同长度后重试。 \
The number of residue matrices differs from the number of poles, so they cannot be paired one-to-one; adjust residues and poles to the same length and try again."|>
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
  <|"usage" -> "FlintNDEPartialFractionSystem[{P0,P1,...}, residues, poles]",
    "MessageTemplate" ->
      "FlintNDEPartialFractionSystem 收到的参数不符合调用形态，它需要多项式系数矩阵列表 {P0,P1,...}、留数列表和极点列表共三个参数；请按 FlintNDEPartialFractionSystem[{P0,P1,...}, residues, poles] 传入这三个列表后重试。 \
FlintNDEPartialFractionSystem received arguments that do not match its calling form, which needs three arguments: a list {P0,P1,...} of polynomial coefficient matrices, a list of residues and a list of poles; call it as FlintNDEPartialFractionSystem[{P0,P1,...}, residues, poles] with those three lists, then try again."|>
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
      Return[Failure["WorkDirectoryCreationFailed", <|"path" -> bridgeDirectory,
        "MessageTemplate" ->
          "无法创建 bridge 运行目录，可能是路径无效或没有写权限；请检查该路径的权限，或通过 WorkDirectory 指定一个可写目录后重试。 \
The bridge runtime directory could not be created, possibly because the path is invalid or not writable; check the permissions of that path or point WorkDirectory at a writable directory, then try again."|>]]
    ]
  ];
  requestWrite = Quiet@Check[Export[requestFile, flintNDEEncode[request], "JSON"], $Failed];
  If[requestWrite === $Failed || ! FileExistsQ[requestFile],
    Return[Failure["RuntimeInputWriteFailed", <|
      "path" -> requestFile, "pathLength" -> StringLength[requestFile],
      "MessageTemplate" ->
        "bridge 请求文件写入失败，可能是运行目录不可写或磁盘已满；请检查所列路径的写权限与磁盘剩余空间后重试。 \
Writing the bridge request file failed, possibly because the runtime directory is not writable or the disk is full; check write permission for the listed path and the free disk space, then try again."
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
    Return[Failure["BridgeLaunchFailed", <|"process" -> process,
      "MessageTemplate" ->
        "Python 后端进程根本未能启动，通常是 Python 选项指向的解释器不存在或不可执行；请核对解释器路径后重试。 \
The Python backend process could not be started at all, usually because the interpreter given by the Python option does not exist or is not executable; verify the interpreter path and try again."|>]]
  ];
  If[! FileExistsQ[outputFile],
    If[StringContainsQ[logText, Alternatives[
        "No module named 'flint'", "No module named \"flint\"",
        "ModuleNotFoundError: No module named 'flint'"
      ]],
      Message[FlintNDEBridgeError::pythonFlintMissing, logText];
      Return[Failure["PythonFlintUnavailable", <|
        "process" -> process, "requestFile" -> requestFile,
        "MessageTemplate" ->
          "Python 已启动但无法导入 python-flint；请在该解释器环境中安装 python-flint，或用 Python 选项指向已装好它的解释器后重试。 \
Python started but could not import python-flint; install python-flint in that interpreter environment, or point the Python option at an interpreter that already has it, then try again."
      |>]]
    ];
    If[exitCode =!= 0,
      Message[FlintNDEBridgeError::launchFailed, logText];
      Return[Failure["BridgeLaunchFailed", <|
        "process" -> process, "requestFile" -> requestFile,
        "MessageTemplate" ->
          "Python 后端以非零退出码结束且没有写出结果文件，具体原因见进程日志；请按日志修正 Python 环境或输入后重试。 \
The Python backend exited with a nonzero code and wrote no result file; see the process log for the cause, fix the Python environment or input accordingly, and try again."
      |>]]
    ];
    Message[FlintNDEBridgeError::outputMissing, logText];
    Return[Failure["BridgeOutputMissing", <|
      "process" -> process, "requestFile" -> requestFile,
      "MessageTemplate" ->
        "Python 后端已退出但没有生成输出文件；请检查 WorkDirectory 是否可写、磁盘是否有余量，并按日志报错修正后重试。 \
The Python backend exited without creating its output file; check that WorkDirectory is writable and has free disk space, fix what the log reports, then try again."
    |>]]
  ];
  Clear[Global`FlintNDEBridgeResult];
  outputLoad = Quiet@Check[Get[outputFile, CharacterEncoding -> "UTF-8"], $Failed];
  If[outputLoad === $Failed,
    Return[Failure["BridgeOutputInvalid", <|"path" -> outputFile, "process" -> process,
      "MessageTemplate" ->
        "后端输出文件存在但无法作为 Wolfram 表达式读取，内容可能已损坏；请检查该文件内容与进程日志后重试。 \
The backend output file exists but cannot be read as a Wolfram expression and may be corrupted; inspect the file contents and the process log, then try again."|>]]
  ];
  result = Global`FlintNDEBridgeResult;
  Quiet[DeleteFile /@ Select[{requestFile, outputFile, logFile}, FileExistsQ]];
  Clear[Global`FlintNDEBridgeResult];
  If[! AssociationQ[result],
    Return[Failure["InvalidBridgeResult", <|"result" -> result,
      "MessageTemplate" ->
        "后端输出可以读取，但没有定义符合 bridge schema 的 Association 结果；请核对后端模块与请求 schema 版本是否匹配后重试。 \
The backend output was readable but did not define a schema-conforming Association result; check that the backend module matches the request schema version, then try again."|>]]
  ];
  If[Lookup[result, "status", None] === "error",
    Message[FlintNDEBridgeError::error, Lookup[result, "message",
      "后端没有给出具体原因。 The backend did not report a specific cause."]];
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
      Lookup[result, "message",
        If[language === "CN",
          "FlintNDEPlanPath：规划出的线段会穿过奇点，已被拒绝。请显式设置 SingularityMode->\"SingularityJump\" 并自行确认分支，或调整输入点让路径避开奇点。",
          "FlintNDEPlanPath: a planned segment would cross a singularity, so the path was refused. Either set SingularityMode->\"SingularityJump\" and confirm the branch yourself, or move the input points so the path avoids the singularity."]],
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
      <|"function" -> "FlintNDEPlanPath", "options" -> unknownOptions,
        "MessageTemplate" ->
          "FlintNDEPlanPath 收到了其当前 Options 之外的未知选项；请删除多余规则或改用 FlintNDEPlanPath 的合法选项后重试。 \
FlintNDEPlanPath received option names that are not in its current Options; remove the extra rules or use only valid FlintNDEPlanPath options and try again."|>
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
  <|"usage" -> "FlintNDEPlanPath[system, start, {point1,...}]",
    "MessageTemplate" ->
      "FlintNDEPlanPath 收到的参数不符合调用形态，它需要一个系统、一个起点和一个目标点列表；请按 FlintNDEPlanPath[system, start, {point1,...}] 传入系统、起点与目标点列表后重试。 \
FlintNDEPlanPath received arguments that do not match its calling form, which needs a system, a start point and a list of target points; call it as FlintNDEPlanPath[system, start, {point1,...}] with a system, a start and a list of points, then try again."|>
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
      <|"function" -> "FlintNDEExecutePath", "options" -> unknownOptions,
        "MessageTemplate" ->
          "FlintNDEExecutePath 收到了其当前 Options 之外的未知选项；请删除多余规则或改用 FlintNDEExecutePath 的合法选项后重试。 \
FlintNDEExecutePath received option names that are not in its current Options; remove the extra rules or use only valid FlintNDEExecutePath options and try again."|>
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
      <|"reason" -> "expected the complete result returned by FlintNDEPlanPath",
        "MessageTemplate" ->
          "plan 参数不是 FlintNDEPlanPath 返回的完整规划结果，缺少要求的 schema、status、operation 或 plan 字段；请把 FlintNDEPlanPath 的返回值原样传入后重试。 \
The plan argument is not the complete result returned by FlintNDEPlanPath; the required schema, status, operation or plan field is missing. Pass the FlintNDEPlanPath result through unchanged and try again."|>
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
        "reason" -> "the plan does not record its planning precision; replan",
        "MessageTemplate" ->
          "该计划没有记录规划精度，无法与本次执行请求的精度比较；请重新运行 FlintNDEPlanPath 生成带精度记录的计划后重试。 \
The plan does not record its planning precision, so it cannot be compared with the precision requested for execution; rerun FlintNDEPlanPath to produce a plan that records its precision and try again."|>
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
  <|"usage" -> "FlintNDEExecutePath[system, initialVector, plan]",
    "MessageTemplate" ->
      "FlintNDEExecutePath 收到的参数不符合调用形态，它需要一个系统、一个初始向量和一份已有计划；请按 FlintNDEExecutePath[system, initialVector, plan] 传入系统、初始向量与计划后重试。 \
FlintNDEExecutePath received arguments that do not match its calling form, which needs a system, an initial vector and an existing plan; call it as FlintNDEExecutePath[system, initialVector, plan] with a system, an initial vector and a plan, then try again."|>
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
    Return[Failure["InvalidEpJob", <|"missingKeys" -> missing, "job" -> job,
      "MessageTemplate" ->
        "某个 ep 任务缺少必需字段（ep、system、start、points、initialVector 中的一项或多项）；请为每个任务补齐全部必需字段后重试。 \
An ep job is missing one or more required fields among ep, system, start, points and initialVector; complete every required field for each job and try again."|>]]
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
      "options" -> unknownOptions,
      "MessageTemplate" ->
        "FlintNDEEvaluateEpBatch 收到了其当前 Options 之外的未知选项；请删除多余规则或改用 FlintNDEEvaluateEpBatch 的合法选项后重试。 \
FlintNDEEvaluateEpBatch received option names that are not in its current Options; remove the extra rules or use only valid FlintNDEEvaluateEpBatch options and try again."|>]]
  ];
  If[jobs === {}, Return[Failure["EpJobListEmpty", <|
    "MessageTemplate" ->
      "FlintNDEEvaluateEpBatch 收到了空的任务列表；请至少提供一个包含 ep、system、start、points 和 initialVector 的任务后重试。 \
FlintNDEEvaluateEpBatch received an empty job list; provide at least one job containing ep, system, start, points and initialVector, then try again."|>]]];
  parallelCount = OptionValue[ParallelTaskCount];
  If[! IntegerQ[parallelCount] || parallelCount < 1,
    Return[Failure["ParallelTaskCountPositiveIntegerRequired",
      <|"value" -> parallelCount, "default" -> 12,
        "MessageTemplate" ->
          "ParallelTaskCount 必须是正整数（缺省为 12），当前值不是；请把 ParallelTaskCount 改为不小于 1 的整数后重试。 \
ParallelTaskCount must be a positive integer (default 12), and the current value is not; set ParallelTaskCount to an integer of at least 1 and try again."|>]]
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
  <|"usage" -> "FlintNDEEvaluateEpBatch[{job1,...}, ParallelTaskCount->12]",
    "MessageTemplate" ->
      "FlintNDEEvaluateEpBatch 收到的参数不符合调用形态，它需要一个 ep 任务列表，并可选地用 ParallelTaskCount 指定并行上限；请按 FlintNDEEvaluateEpBatch[{job1,...}, ParallelTaskCount->12] 传入任务列表后重试。 \
FlintNDEEvaluateEpBatch received arguments that do not match its calling form, which needs a list of ep jobs and optionally a ParallelTaskCount setting for the parallel limit; call it as FlintNDEEvaluateEpBatch[{job1,...}, ParallelTaskCount->12] with a job list, then try again."|>
];


End[];


EndPackage[];
