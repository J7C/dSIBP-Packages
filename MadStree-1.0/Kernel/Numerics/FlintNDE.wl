(* ::Package:: *)

(***
文件：FlintNDE.wl
用途：为单阶段 MSEvaluatePath 提供 exact 序列化、运行目录、缓存和后端进程调用。
范围：本模块不选择输运节点、不构造路径计划，也不保存或恢复计划对象。
***)

(* ::Chapter:: *)
(* Exact Q(i)(s) serialization with automatic decimal fallback *)

(* Exact algebraic constants like (-1)^(1/3) appear in dlog residue matrices;
   RootReduce turns them into Root objects that ComplexExpand resolves into
   explicit a+b I components, so the serialized strings stay parseable Q(i)
   literals on the Python side. *)
msGaussianRationalParts[value_] := Module[{real, imaginary},
  real = Together[ComplexExpand[Re[RootReduce[value]]]];
  imaginary = Together[ComplexExpand[Im[RootReduce[value]]]];
  (* RationalQ 不是内置符号且包内从未定义；此处谓词必须真正可求值，否则 Q(i) 门禁是空转。 *)
  If[! MatchQ[real, _Integer | _Rational] || ! MatchQ[imaginary, _Integer | _Rational], Return[$Failed]];
  {real, imaginary}
];

(* Q(i) 优先：整数/有理/有理复数与 RootReduce 可归一的精确值照旧输出精确 Q(i) 串；
   根式、超越常数等非 Q(i) 精确值按 digits 自动换算小数（后端按精确有理解析，扰动
   ~10^-digits 远低于目标误差）；仍非数值（未解析符号等）才 fail-closed。 *)
msGaussianRationalString[value_, digits_Integer] := Module[{parts, numeric},
  (* 分支条件必须恒为布尔值：Which 遇到非 True/False 的条件会整条不求值。 *)
  parts = Which[
    MatchQ[value, _Integer | _Rational], {value, 0},
    MatchQ[value, z_Complex /; MatchQ[Re[z], _Integer | _Rational] && MatchQ[Im[z], _Integer | _Rational]], {Re[value], Im[value]},
    True, msGaussianRationalParts[value]
  ];
  If[parts =!= $Failed,
    Return[<|
      "real" -> ToString[parts[[1]], InputForm],
      "imag" -> ToString[parts[[2]], InputForm]
    |>]
  ];
  numeric = Quiet[N[value, digits]];
  If[! NumericQ[numeric], Return[$Failed]];
  <|
    "real" -> msDecimalString[Re[numeric], digits],
    "imag" -> msDecimalString[Im[numeric], digits]
  |>
];

msPolynomialCoefficientStrings[polynomial_, parameter_Symbol, digits_Integer] := Module[{degree, coefficients, strings},
  If[! PolynomialQ[polynomial, parameter], Return[$Failed]];
  degree = Exponent[polynomial, parameter];
  coefficients = If[polynomial === 0, {0}, CoefficientList[polynomial, parameter]];
  strings = msGaussianRationalString[#, digits] & /@ coefficients;
  If[MemberQ[strings, $Failed], $Failed, strings]
];

msRationalFunctionRecord[expression_, parameter_Symbol, digits_Integer] := Module[
  {rational = Cancel[Together[expression]], numerator, denominator, numeratorStrings, denominatorStrings},
  numerator = Numerator[rational];
  denominator = Denominator[rational];
  numeratorStrings = msPolynomialCoefficientStrings[numerator, parameter, digits];
  denominatorStrings = msPolynomialCoefficientStrings[denominator, parameter, digits];
  If[MemberQ[{numeratorStrings, denominatorStrings}, $Failed],
    $Failed,
    <|"numerator" -> numeratorStrings, "denominator" -> denominatorStrings|>
  ]
];

msRationalMatrixRecords[matrix_?MatrixQ, parameter_Symbol, digits_Integer] := Module[{records},
  records = Map[msRationalFunctionRecord[#, parameter, digits] &, matrix, {2}];
  If[! FreeQ[records, $Failed], $Failed, records]
];

msDecimalString[value_, digits_Integer] := StringReplace[
  ToString[N[value, digits], InputForm],
  {RegularExpression["`[0-9.]*"] -> "", "*^" -> "e"}
];

msParseDecimalString[text_String] := ToExpression[StringReplace[text, "e" -> "*^"]];

msComplexDecimalRecord[value_?NumericQ, digits_Integer] := <|
  "real" -> msDecimalString[Re[N[value, digits]], digits],
  "imag" -> msDecimalString[Im[N[value, digits]], digits]
|>;

(* ::Chapter:: *)
(* Build the one-dimensional system along the affine kinematic path *)

msRuleValue[symbol_, rules_List] := symbol /. rules;


(* ::Chapter:: *)
(* Structured tensor letter serialization *)

(* 贡献记录只搬运张量因子：sector 对角向量、0/1 指示符、preVector 和稀疏中矩阵。
   2^n 维留数矩阵永远不进入请求；后端用局部基标签逐轴重建 U/U^{-1}。
   代入后的元绝大多数已是规范 Q(i) 字面量，msGaussianRationalString 对它们直接取分量；
   根式等未归一值自动按 digits 小数化，非数值才整体放弃。 *)
msTensorGaussianVectorRecord[vector_List, constantRules_List, digits_Integer] := Module[{records},
  records = msGaussianRationalString[#, digits] & /@ (vector /. constantRules);
  If[MemberQ[records, $Failed], $Failed, records]
];


msTensorSparseMatrixRecord[matrix_?MatrixQ, constantRules_List, digits_Integer] := Module[
  {substituted, positions, rows = {}, columns = {}, entries = {}, value, record},
  substituted = matrix /. constantRules;
  positions = SparseArray[substituted]["NonzeroPositions"];
  Do[
    value = substituted[[Sequence @@ position]];
    If[! TrueQ[PossibleZeroQ[value]],
      record = msGaussianRationalString[value, digits];
      If[record === $Failed, Return[$Failed]];
      AppendTo[rows, position[[1]]];
      AppendTo[columns, position[[2]]];
      AppendTo[entries, record]
    ],
    {position, positions}
  ];
  <|"rows" -> rows, "columns" -> columns, "entries" -> entries|>
];


(* JSON 没有 Missing；None 字段在线格式里写成 null。 *)
msTensorContributionRecord[contribution_Association, constantRules_List, digits_Integer] := Module[
  {vector, preVector, middle, postFactor},
  Switch[contribution["type"],
    "diagonal",
      vector = msTensorGaussianVectorRecord[contribution["vector"], constantRules, digits];
      If[vector === $Failed,
        $Failed,
        <|
          "kind" -> "diagonal",
          "scope" -> contribution["scope"],
          "sectorKey" -> contribution["sectorKey"],
          "vector" -> vector
        |>
      ],
    "conjugated",
      preVector = If[
        contribution["preVector"] === None,
        Null,
        msTensorGaussianVectorRecord[contribution["preVector"], constantRules, digits]
      ];
      If[preVector === $Failed, Return[$Failed]];
      middle = If[
        contribution["middleMatrix"] === None,
        Null,
        msTensorSparseMatrixRecord[contribution["middleMatrix"], constantRules, digits]
      ];
      If[middle === $Failed, Return[$Failed]];
      postFactor = msGaussianRationalString[contribution["postFactor"] /. constantRules, digits];
      If[postFactor === $Failed, Return[$Failed]];
      <|
        "kind" -> "conjugated",
        "scope" -> contribution["scope"],
        "sectorKey" -> contribution["sectorKey"],
        "preVector" -> preVector,
        "letterIndicator" -> contribution["letterIndicator"],
        "middle" -> middle,
        "postFactor" -> postFactor
      |>,
    _,
      $Failed
  ]
];


(* Along the affine path every certified dlog letter is alpha+beta s, so
   d Log[letter]/ds is a simple pole at -alpha/beta. The residue is the same
   structured tensor list for all segments; the backend merges coincident poles.
   Letters with beta==0 contribute nothing and are skipped. *)
msAffineTensorLetterRecord[
  letter_,
  contributions_List,
  pathRules_List,
  constantRules_List,
  parameter_Symbol,
  digits_Integer
] := Module[
  {along, alpha, beta, alphaRecord, betaRecord, contributionRecords},
  along = Simplify[letter /. pathRules /. constantRules];
  If[! PolynomialQ[along, parameter] || Exponent[along, parameter] > 1, Return[$Failed]];
  alpha = Together[along /. parameter -> 0];
  beta = Together[Coefficient[along, parameter, 1]];
  If[alpha === 0, Return[$Failed]];
  If[beta === 0, Return[<|"skipped" -> True|>]];
  alpha = RootReduce[alpha];
  beta = RootReduce[beta];
  alphaRecord = msGaussianRationalString[alpha, digits];
  betaRecord = msGaussianRationalString[beta, digits];
  contributionRecords = Table[
    msTensorContributionRecord[contribution, constantRules, digits],
    {contribution, contributions}
  ];
  If[
    alphaRecord === $Failed || betaRecord === $Failed ||
      MemberQ[contributionRecords, $Failed],
    $Failed,
    <|
      "alpha" -> alphaRecord,
      "beta" -> betaRecord,
      "contributions" -> contributionRecords
    |>
  ]
];


msAffineTensorLetterData[
  de_Association,
  anchorRules_List,
  targetRules_List,
  parameter_Symbol,
  digits_Integer
] := Module[{
  symbols, unresolved, pathRules, constantRules, records, record
},
  symbols = de["kinematicSymbols"];
  unresolved = Select[
    symbols,
    ! NumericQ[N[msRuleValue[#, anchorRules]]] || ! NumericQ[N[msRuleValue[#, targetRules]]] &
  ];
  If[unresolved =!= {}, Return[msFailure["IncompleteKinematicPath", <|"symbols" -> unresolved|>]]];
  pathRules = Map[
    Function[symbol,
      symbol -> (msRuleValue[symbol, anchorRules] +
        parameter (msRuleValue[symbol, targetRules] - msRuleValue[symbol, anchorRules]))
    ],
    symbols
  ];
  constantRules = Select[targetRules, FreeQ[symbols, First[#]] &];
  records = {};
  Do[
    record = msAffineTensorLetterRecord[
      letter, Lookup[de["tensorResidues"], letter, {}],
      pathRules, constantRules, parameter, digits
    ];
    If[record === $Failed,
      Return[msFailure[
        "FlintNDEExactPathRequired",
        <|"reason" -> "dlog letter is not affine with Q(i) tensor data along the path",
          "letter" -> letter|>
      ]]
    ];
    If[! TrueQ[Lookup[record, "skipped", False]], AppendTo[records, record]],
    {letter, de["letters"]}
  ];
  <|
    "parameter" -> parameter,
    "pathRules" -> pathRules,
    "constantRules" -> constantRules,
    "tensorSectors" -> Values[de["tensorSectors"]],
    "letterRecords" -> records
  |>
];


(* 沿 blow-up 曲线的 letter 拉回是 t 的有理函数 N(t)/D(t)；本记录只下发互素
   分子/分母系数与张量贡献，后端按 d Log[N/D] 把根装配成认证极点。次数不必
   限制在一次，t^2、t^3 因子由后端数值根隔离处理。与 t 无关的拉回 dlog 恒为
   零，直接跳过；其数值可含 Pi 等超越常数（如 pinch normalization），不进入
   Q(i) 序列化。 *)
msBoundaryTensorLetterRecord[
  along_,
  contributions_List,
  constantRules_List,
  parameter_Symbol,
  digits_Integer
] := Module[
  {record, numerator, denominator, contributionRecords},
  If[FreeQ[along, parameter], Return[<|"skipped" -> True|>]];
  record = msRationalFunctionRecord[along, parameter, digits];
  If[record === $Failed, Return[$Failed]];
  {numerator, denominator} = Lookup[record, {"numerator", "denominator"}];
  If[Length[numerator] === 1 && Length[denominator] === 1,
    Return[<|"skipped" -> True|>]
  ];
  contributionRecords = Table[
    msTensorContributionRecord[contribution, constantRules, digits],
    {contribution, contributions}
  ];
  If[MemberQ[contributionRecords, $Failed], Return[$Failed]];
  <|
    "numerator" -> numerator,
    "denominator" -> denominator,
    "contributions" -> contributionRecords
  |>
];


(* 逐 letter 序列化边界曲线拉回；系数与张量贡献按 Q(i) 优先、非 Q(i) 精确值自动
   小数化下发。仍出现非数值（未解析符号）时整体放弃并回退后端现行 Cauchy 路线，
   因此这里返回 Missing 而不是 Failure。 *)
msBoundaryTensorLetterData[boundary_Association, de_Association, digits_Integer] := Module[
  {parameter, constantRules, records},
  parameter = boundary["singularParameter"];
  constantRules = Lookup[boundary, "singularRetainedRules", {}];
  records = Select[
    Table[
      msBoundaryTensorLetterRecord[
        entry["along"],
        Lookup[de["tensorResidues"], entry["letter"], {}],
        constantRules, parameter, digits
      ],
      {entry, Lookup[boundary, "singularLetterPullbacks", {}]}
    ],
    ! TrueQ[Lookup[#, "skipped", False]] &
  ];
  If[MemberQ[records, $Failed], Missing["NotExact"], records]
];


(* ::Chapter:: *)
(* Python backend invocation *)

(* 缺省临时目录只追加一次 results_temp；显式 MSRuntimeDirectory 则直接表示临时目录本身。 *)
msDefaultRuntimeDirectory[] := FileNameJoin[{msRuntimeDirectory[], "results_temp"}];


msResolveRuntimeDirectory[Automatic] := msDefaultRuntimeDirectory[];


(* Absolute paths are expanded directly; relative paths are always resolved against the calling script directory to avoid depending on the current process directory. *)
msAbsoluteRuntimePathQ[path_String] := msAbsolutePathQ[path];


msResolveRuntimeDirectory[path_String] := Module[{base = msDefaultRuntimeDirectory[]},
  ExpandFileName[If[msAbsoluteRuntimePathQ[path], path, FileNameJoin[{base, path}]]]
];
msResolveRuntimeDirectory[other_] := msFailure[
  "RuntimeDirectoryRequired",
  <|"value" -> HoldForm[other]|>
];


(* 为参与缓存的源码生成固定宽度 SHA-256；文件缺失也进入身份记录，
   避免失败配置与某个既有成功缓存发生碰撞。 *)
msFlintNDESourceDigestRecord[label_String, file_String] := <|
  "file" -> label,
  "sha256" -> If[
    FileExistsQ[file],
    IntegerString[FileHash[file, "SHA256"], 16, 64],
    "missing"
  ]
|>;


(* 缓存必须同时绑定 MadStree adapter 与实际 FlintNDE 源码。这样即使版本号
   不变，原位修复也会自动生成新缓存目录，不会复用与当前源码身份不匹配的序列化计划。 *)
msFlintNDEBackendSourceIdentity[inputData_Association] := Module[{
  adapterFile, backendRoot, packageFile, moduleDirectory, moduleFiles
},
  adapterFile = FileNameJoin[{
    $MadStreePackageDirectory, "Backend", "flintnde_transport.py"
  }];
  backendRoot = Lookup[inputData, "backendPackagePath", Missing["Absent"]];
  If[! StringQ[backendRoot],
    Return[{
      msFlintNDESourceDigestRecord["Backend/flintnde_transport.py", adapterFile],
      <|"file" -> "FlintNDE", "sha256" -> "missing-backend-path"|>
    }]
  ];
  backendRoot = ExpandFileName[backendRoot];
  packageFile = FileNameJoin[{backendRoot, "pyproject.toml"}];
  moduleDirectory = FileNameJoin[{backendRoot, "flintnde"}];
  moduleFiles = If[
    DirectoryQ[moduleDirectory],
    Sort[FileNames["*.py", moduleDirectory]],
    {}
  ];
  Join[
    {
      msFlintNDESourceDigestRecord["Backend/flintnde_transport.py", adapterFile],
      msFlintNDESourceDigestRecord["FlintNDE/pyproject.toml", packageFile]
    },
    msFlintNDESourceDigestRecord[
      "FlintNDE/flintnde/" <> FileNameTake[#], #
    ] & /@ moduleFiles
  ]
];


(* 后端成功输出缓存在临时目录的 cache/<digest>.json。缓存键包含
   去除保存目录后的请求和后端源码身份；失败输出永不作为缓存命中。 *)
msFlintNDECacheFile[runtimeDirectory_String, inputData_Association] := Module[{
  stablePayload, sourceIdentity, digest
},
  sourceIdentity = msFlintNDEBackendSourceIdentity[inputData];
  stablePayload = Append[
    KeyDrop[inputData, {"backendPackagePath"}],
    "backendSourceIdentity" -> sourceIdentity
  ];
  digest = IntegerString[
    Hash[ExportString[stablePayload, "RawJSON"], "SHA256"], 16, 64
  ];
  FileNameJoin[{runtimeDirectory, "cache", digest <> ".json"}]
];


(* The backend interpreter is resolved without any filesystem probing: an
   explicit PythonExecutable string is used verbatim, Automatic honours the
   MADSTREE_PYTHON environment variable when it is set and non-empty, and
   otherwise the bare "python" command from PATH is used. Installing a
   suitable interpreter and keeping it reachable is the user's responsibility. *)
msResolvePythonExecutable[executable_String] := executable;


msResolvePythonExecutable[Automatic] := With[{override = Quiet[Environment["MADSTREE_PYTHON"]]},
  If[StringQ[override] && StringLength[StringTrim[override]] > 0,
    StringTrim[override],
    "python"
  ]
];


(* 短 token 只负责运行文件去重；独立 helper 允许回归测试固定写入目标。 *)
msFlintNDECreateToken[] := StringTake[StringDelete[CreateUUID[], "-"], 16];


msExecuteFlintNDEAdapter[inputData_Association, pythonExecutable_, runtimeDirectory_String] := Module[
  {cacheDirectory, cacheFile, cached, transportDirectory, identifier, inputFile,
   adapterFile, logFile, command, processResult, exitCode, standardOutput,
   standardError, logText, logTail, process, imported, launchFailure,
   pathFailure, inputWrite},
  cacheFile = msFlintNDECacheFile[runtimeDirectory, inputData];
  cacheDirectory = DirectoryName[cacheFile];
  transportDirectory = FileNameJoin[{runtimeDirectory, "nde"}];
  identifier = msFlintNDECreateToken[];
  inputFile = FileNameJoin[{transportDirectory, "i-" <> identifier <> ".json"}];
  (* Backend stdout/stderr is captured so a failed launch reports the actual
     Python error instead of leaving it on the console. *)
  logFile = FileNameJoin[{transportDirectory, "l-" <> identifier <> ".txt"}];
  adapterFile = FileNameJoin[{$MadStreePackageDirectory, "Backend", "flintnde_transport.py"}];
  pathFailure = msRuntimePathLengthFailure[{
    runtimeDirectory, cacheDirectory, cacheFile, transportDirectory, inputFile, logFile
  }];
  If[Head[pathFailure] === Failure, Return[pathFailure]];
  If[FileExistsQ[cacheFile],
    cached = Quiet@Check[Import[cacheFile, "RawJSON"], $Failed];
    If[AssociationQ[cached] && Lookup[cached, "status", None] === "success",
      Return[Join[cached, <|"cacheHit" -> True, "cacheDirectory" -> cacheDirectory|>]]
    ]
  ];
  If[msEnsureDirectory[cacheDirectory] === $Failed,
    Return[msFailure["RuntimeDirectoryCreationFailed", <|"path" -> cacheDirectory|>]]
  ];
  If[msEnsureDirectory[transportDirectory] === $Failed,
    Return[msFailure["RuntimeDirectoryCreationFailed", <|"path" -> transportDirectory|>]]
  ];
  (* 唯一 evaluate schema 的字段集合由 Python 端严格校验。 *)
  inputWrite = Quiet@Check[Export[inputFile, inputData, "RawJSON"], $Failed];
  If[inputWrite === $Failed || ! FileExistsQ[inputFile],
    Return[msFailure["RuntimeInputWriteFailed", <|
      "path" -> inputFile, "pathLength" -> StringLength[inputFile]
    |>]]
  ];
  command = {
    pythonExecutable,
    adapterFile,
    inputFile,
    cacheFile
  };
  (* 参数列表启动不经过 shell；stdout/stderr 由 Wolfram 捕获后写入 UTF-8 日志。 *)
  launchFailure = False;
  processResult = Quiet@Check[
    RunProcess[command, All, ProcessDirectory -> $MadStreePackageDirectory],
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
  logText = StringRiffle[
    Select[{standardOutput, standardError}, StringLength[#] > 0 &],
    "\n"
  ];
  Quiet@Check[Export[logFile, logText, "Text", CharacterEncoding -> "UTF-8"], Null];
  logTail = If[
    StringLength[logText] > 2000,
    "..." <> StringTake[logText, -2000],
    If[StringLength[logText] > 0, logText, "<no backend output captured>"]
  ];
  process = <|
    "ExitCode" -> exitCode,
    "Command" -> command,
    "WorkingDirectory" -> $MadStreePackageDirectory,
    "StandardOutput" -> standardOutput,
    "StandardError" -> standardError,
    "runtimeDirectory" -> runtimeDirectory,
    "inputFile" -> inputFile,
    "cacheFile" -> cacheFile,
    "logFile" -> logFile
  |>;
  If[launchFailure,
    Message[MSEvaluatePath::backendLaunchFailed, pythonExecutable, logTail];
    Return[msFailure["FlintNDELaunchFailed", <|"process" -> process, "stderr" -> logTail|>]]
  ];
  If[! FileExistsQ[cacheFile],
    If[StringContainsQ[logTail, Alternatives[
        "No module named 'flint'", "No module named \"flint\"", "ModuleNotFoundError: No module named 'flint'"
      ]],
      Message[MSEvaluatePath::pythonFlintMissing, pythonExecutable, logTail];
      Return[msFailure["PythonFlintUnavailable", <|"process" -> process, "stderr" -> logTail|>]]
    ];
    If[process["ExitCode"] =!= 0,
      Message[MSEvaluatePath::backendLaunchFailed, pythonExecutable, logTail];
      Return[msFailure["FlintNDELaunchFailed", <|"process" -> process, "stderr" -> logTail|>]]
    ];
    Message[MSEvaluatePath::backendOutputMissing, pythonExecutable, logTail];
    Return[msFailure["FlintNDEOutputMissing", <|"process" -> process, "stderr" -> logTail|>]]
  ];
  imported = Quiet@Check[Import[cacheFile, "RawJSON"], $Failed];
  If[! AssociationQ[imported],
    Return[msFailure["FlintNDEOutputInvalid", <|
      "process" -> process, "outputFile" -> cacheFile, "stderr" -> logTail
    |>]]
  ];
  (* A structured refusal (no-singularity mode found a pole on the polyline,
     or the leading-order local basis cannot be built, e.g. resonance) is an
     outcome, not a backend failure: it is returned verbatim with cacheHit
     False and never cached as a success. *)
  If[MemberQ[{"singularPathRefused", "leadingOrderRefused"}, Lookup[imported, "status", None]],
    Quiet[DeleteFile /@ Select[{inputFile, logFile}, FileExistsQ]];
    Return[Join[imported, <|"cacheHit" -> False, "cacheDirectory" -> cacheDirectory|>]]
  ];
  If[process["ExitCode"] =!= 0 || Lookup[imported, "status", None] =!= "success",
    (* 同时打印后端序列化的 traceback：只看外层 "error" 字段会掩蔽
       链式真因（如 exact_gaussian 包装的 int<->str 位数上限）。 *)
    Message[MSEvaluatePath::backendRunFailed, pythonExecutable,
      Lookup[imported, "error", logTail] <>
        If[StringQ[Lookup[imported, "traceback", ""]],
          "\n" <> Lookup[imported, "traceback", ""], ""]];
    Return[msFailure["FlintNDEProcessFailed", <|"process" -> process, "backend" -> imported, "stderr" -> logTail|>]]
  ];
  Quiet[DeleteFile /@ Select[{inputFile, logFile}, FileExistsQ]];
  Join[imported, <|"cacheHit" -> False, "cacheDirectory" -> cacheDirectory,
    "process" -> KeyDrop[process, {"inputFile", "cacheFile"}]|>]
];


(* Writes one exact-or-decimal normalized leading branch as the FlintNDE {a,b,C}
   contract: a 与 C 按 Q(i) 优先、非 Q(i) 精确值自动小数化；b 是整数 log 次数。 *)
msFrobeniusBranchRecord[branch_Association, digits_Integer] := Module[{exponent, vector},
  exponent = msGaussianRationalString[branch["frobeniusExponent"], digits];
  vector = msGaussianRationalString[#, digits] & /@ branch["normalizedLeadingVector"];
  If[exponent === $Failed || MemberQ[vector, $Failed], Return[$Failed]];
  <|
    "a" -> exponent,
    "b" -> branch["logPower"],
    "C" -> vector
  |>
];


(* A singular target may carry a deliberate text value such as Infinity. *)
msParseFlintVector[records_List] := Map[
  Function[record,
    If[KeyExistsQ[record, "text"], record["text"],
      msParseDecimalString[record["real"]] + I msParseDecimalString[record["imag"]]]],
  records
];
