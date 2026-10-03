(* ::Package:: *)

(* ::Chapter:: *)
(*tree vertex-family Private 公式内核*)

(* vertex-family 公式在冻结核心中；本模块只增加 DSInit context 适配，不复制递推公式。 *)

(* 从 root topology 统一建立 sector family，并给每个 family 保存同一个 root 引用。
   direct seed、tagged 迭代和 raw 迭代共用这一个构造，避免回退到旧 loop 投影。 *)
dsTreeFamilyContextFromRootTopology[rootTopology_Association] := Module[{familyContext, families},
   familyContext = makeTreeSectorFamilies[rootTopology];
   If[familyContext === $Failed, Return[$Failed]];
   families = Join[#, <|"rootTopology" -> rootTopology|>] & /@ familyContext["families"];
   Join[familyContext, <|"families" -> families, "topFamily" -> First[families]|>]
   ];


dsTreeFamilyContext[context_Association] /; dsContextQ[context] :=
   dsTreeFamilyContextFromRootTopology[context["topology"]];


(* raw J[vertexPacks] 的 source-aware 首步复用 sector-tagged direct seed，再仅在公开 raw
   返回边界去掉 sector token。这样传播子动量模长与共同-theta contact 都保持同源。 *)
treeSourceAwareStepFromTopology[
   int_J, vertexIndex_Integer, endpoint_Integer, data_Association
   ] /; AssociationQ[Lookup[data, "rootTopology", Missing["NoRootTopology"]]] := Module[
   {familyContext, token, taggedResult},
   familyContext = dsTreeFamilyContextFromRootTopology[data["rootTopology"]];
   If[familyContext === $Failed, Return[$Failed]];
   token = dsTreeToken[data["sector"], int];
   taggedResult = dsTreeTaggedSourceAwareStep[token, vertexIndex, endpoint, data, familyContext];
   If[taggedResult === $Failed, Return[$Failed]];
   taggedResult /. dsTreeToken[_, item_J] :> item
   ];


(* 公开 tree overload 统一在 PublicBoundary.wl 定义。这里仅保留 Private
   vertex-basis 内核，避免用户绕过三参数 J 的 sector 身份与表示审计。 *)


(* ::Chapter:: *)
(*多 sector dlog 数据汇总*)

(* 每个 sector 的对角块仍由 vertex-family 公式构造；非对角块从 loop dtau 的 tagged contact source 提取，
   因而共同 theta、多线 simultaneous contact 和 mixed-sign 禁用规则只在既有 loop 边界层出现一次。 *)

(* p=0 的 terminal contact vertex 没有 massive leg；空能量表必须是 {}，不能让 Lookup 产生 Missing。 *)
treeVertexDLogData[vertex_Association] := Module[
   {p, states, momentumMagnitudes, k0, omega0, omegaEx, tp, tpInv, m1, omega, letters, coeffs},
   p = vertex["p"];
   states = treeBinaryStates[p];
   momentumMagnitudes = If[vertex["massiveLegs"] === {}, {}, Lookup[vertex["massiveLegs"], "momentumMagnitude"]];
   k0 = vertex["signedExternalLegEnergy"];
   omega0 = -I DiagonalMatrix[Table[
      Log[k0 + Sum[(2 states[[row, i]] - 1) momentumMagnitudes[[i]], {i, p}]],
      {row, Length[states]}
      ]];
   omegaEx = DiagonalMatrix[Table[
     -Sum[states[[row, i]] (2 vertex["massiveLegs"][[i, "nu"]] + 1) Log[momentumMagnitudes[[i]]], {i, p}],
     {row, Length[states]}
     ]];
   tp = treeTp[vertex];
   tpInv = treeTpInverse[vertex];
   m1 = treeM1[vertex, vertex["nu0"] + 1];
   omega = Expand[omegaEx - I tpInv . omega0 . tp . m1];
   letters = DeleteDuplicates@Join[momentumMagnitudes, Cases[omega0, Log[arg_] :> arg, Infinity]];
   coeffs = Association@Table[letter -> Map[Coefficient[#, Log[letter]] &, omega, {2}], {letter, letters}];
   <|"vertex" -> vertex["id"], "states" -> states, "omega" -> omega, "letters" -> letters, "letterMatrices" -> coeffs|>
   ];


dsTreeBlockDiagonal[matrices_List] := Module[{dimensions = Length /@ matrices},
   If[matrices === {}, Return[{}]];
   ArrayFlatten@Table[
     If[i === j, matrices[[i]], ConstantArray[0, {dimensions[[i]], dimensions[[j]]}]],
     {i, Length[matrices]}, {j, Length[matrices]}
     ]
   ];


dsTreeMasterSectorOffsets[blocks_List] := Module[{counts, starts, ends},
   counts = Lookup[blocks, "masterCount"];
   starts = 1 + Most[Accumulate[Prepend[counts, 0]]];
   ends = Accumulate[counts];
   MapThread[
    <|"sectorKey" -> #1, "start" -> #2, "end" -> #3, "count" -> #4|> &,
    {Lookup[blocks, "sector"], starts, ends, counts}
    ]
   ];


(* ::Section::Closed:: *)
(*Contact selector 与 sector normalization*)

(* 对 parent sector 的 a=0 master 逐顶点调用 time-IBP；同 sector 项属于 M1/M0，只有 lower-sector 项进入 R。 *)
dsTreeContactRows[family_Association, masters_List, context_Association] := Module[
   {sourceSector = family["sector"], familyContext, rowsByVertex, vertexIndex, seedMaster, record, linearData, terms,
    sourceData, reducedSource},
   familyContext = dsTreeFamilyContext[context];
   If[familyContext === $Failed, Return[<|"status" -> "error", "reason" -> "treeFamilyInitializationFailed"|>]];
   rowsByVertex = Association@Table[
      vertexIndex = First@FirstPosition[family["vertexOrder"], vertexId];
      vertexId -> Table[
        (* dlog source 来自 f^(1) 的约化，即论文 Eq. (3.67) 中 R^(1)；binary state 不变。 *)
        seedMaster = J[ReplacePart[
           First[masters[[row]]],
           vertexIndex -> ReplacePart[First[masters[[row]]][[vertexIndex]], 1 -> 1]
           ]];
        record = dsDirectTreeSeedRecord[vertexId, seedMaster, family, familyContext];
        If[! AssociationQ[record] || Lookup[record, "status", "error"] =!= "generated",
         Return[<|"status" -> "error", "reason" -> "contactSeedGenerationFailed",
           "sectorKey" -> sourceSector, "vertex" -> vertexId, "row" -> row, "record" -> record|>]
         ];
        linearData = Lookup[record, "treeLinearData", <||>];
        If[! dsTreeLinearDataQ[linearData],
         Return[<|"status" -> "error", "reason" -> "contactSeedTaggingFailed",
           "sectorKey" -> sourceSector, "vertex" -> vertexId, "row" -> row|>]
         ];
        terms = Select[linearData["terms"], Lookup[#, "sectorKey", sourceSector] =!= sourceSector &];
        sourceData = <|"status" -> "generated", "terms" -> terms,
          "termCount" -> Length[terms], "sectorKeys" -> DeleteDuplicates[Lookup[terms, "sectorKey", {}]],
          "expression" -> Total[Lookup[#, "coefficient"] Lookup[#, "integral"] & /@ terms],
          "coefficientConvention" -> "complete physical powers: a+a0 and b+b0 or bS+bS0"|>;
        reducedSource = If[terms === {},
          Join[sourceData, <|"status" -> "reduced"|>],
          dsRepIterativeTreeLinearData[sourceData, Automatic, context]
          ];
        If[Lookup[reducedSource, "status", "error"] =!= "reduced",
         Return[<|"status" -> "error", "reason" -> "contactSourceReductionFailed",
           "sectorKey" -> sourceSector, "vertex" -> vertexId, "row" -> row,
           "sourceData" -> sourceData, "reducedSource" -> reducedSource|>]
         ];
        <|"row" -> row, "master" -> masters[[row]], "seedIntegral" -> seedMaster, "rawTerms" -> terms,
          "terms" -> Lookup[reducedSource, "terms", {}],
          "reductionSteps" -> Lookup[reducedSource, "steps", 0]|>,
        {row, Length[masters]}
        ],
      {vertexId, family["vertexOrder"]}
      ];
   <|"status" -> "generated", "sectorKey" -> sourceSector, "rowsByVertex" -> rowsByVertex|>
   ];


(* 每个 target sector 单独按其 master order 抽取矩阵，避免相同裸 J shape 跨 sector 混淆。 *)
dsTreeContactMatrices[
   rowData_Association,
   sectorOrder_List,
   masterLists_Association
   ] := Module[{sourceSector, targetSectors, unresolved = {}, matricesByVertex, positions, matrix},
   If[Lookup[rowData, "status", "error"] =!= "generated", Return[rowData]];
   sourceSector = rowData["sectorKey"];
   targetSectors = Drop[sectorOrder, First@FirstPosition[sectorOrder, sourceSector]];
   matricesByVertex = Association@KeyValueMap[
      Function[{vertexId, rows},
       vertexId -> Association@Table[
          positions = AssociationThread[masterLists[targetSector] -> Range[Length[masterLists[targetSector]]]];
          matrix = ConstantArray[0, {Length[rows], Length[masterLists[targetSector]]}];
          Do[
           Do[
            If[Lookup[term, "sectorKey", None] === targetSector,
             If[KeyExistsQ[positions, term["integral"]],
              matrix[[row, positions[term["integral"]]]] += term["coefficient"],
              AppendTo[unresolved, <|"sourceSector" -> sourceSector, "targetSector" -> targetSector,
                "vertex" -> vertexId, "row" -> row, "term" -> term|>]
              ]
             ],
            {term, rows[[row, "terms"]]}
            ],
           {row, Length[rows]}
           ];
          targetSector -> Expand[matrix],
          {targetSector, targetSectors}
          ]
       ],
      rowData["rowsByVertex"]
      ];
   If[unresolved =!= {},
    <|"status" -> "error", "reason" -> "contactSourceNotInTargetMasterBasis", "unresolved" -> unresolved|>,
    <|"status" -> "generated", "sectorKey" -> sourceSector,
      "rowsByVertex" -> rowData["rowsByVertex"], "matricesByVertex" -> matricesByVertex|>
    ]
   ];


dsTreeContactTransitions[contactData_List] := Flatten[Map[
    Function[sourceData,
     KeyValueMap[
      Function[{vertexId, targetMatrices},
       Flatten@KeyValueMap[
         Function[{targetSector, matrix},
          Cases[
           Flatten[MapIndexed[
             Function[{value, position},
              <|"sourceSector" -> sourceData["sectorKey"], "targetSector" -> targetSector,
                "vertex" -> vertexId, "row" -> position[[1]], "column" -> position[[2]],
                "coefficient" -> value|>
              ],
             matrix,
             {2}
             ]],
           item_Association /; ! TrueQ[Together[item["coefficient"]] === 0]
           ]
          ],
         targetMatrices
         ]
       ],
      sourceData["matricesByVertex"]
      ]
     ],
    contactData
    ], 2];


(* contact selector 把各 lower sector 的入边化为常矩阵；物理 sector prefactor 另由
   sector metadata 保存。两层不能混名，否则容易在 DE 中漏掉或重复连续零点幂。 *)
dsTreeDLogSelectorCoefficients[sectorOrder_List, contactData_List, families_List] := Module[
   {transitions, selectorCoefficients = <|First[sectorOrder] -> 1|>, audits = {}, energyExpressions, energyVariables,
    child, candidates, chosen, selectorCoefficient, ratios, nonconstant},
   transitions = dsTreeContactTransitions[contactData];
   energyExpressions = Flatten[Cases[
      families,
      vertex_Association /; KeyExistsQ[vertex, "signedExternalLegEnergy"] :>
       Join[{vertex["signedExternalLegEnergy"]}, Lookup[vertex["massiveLegs"], "momentumMagnitude", {}]],
      Infinity
      ]];
   energyVariables = DeleteDuplicates[Variables[energyExpressions]];
   Do[
    child = sectorOrder[[index]];
    candidates = Select[transitions,
      Lookup[#, "targetSector", None] === child &&
        KeyExistsQ[selectorCoefficients, Lookup[#, "sourceSector", None]] &];
    If[candidates === {},
     Return[<|"status" -> "error", "reason" -> "selectorSourceMissing", "sectorKey" -> child|>]
     ];
    chosen = First[candidates];
    selectorCoefficient = Expand[
      selectorCoefficients[chosen["sourceSector"]] chosen["coefficient"]
      ];
    AssociateTo[selectorCoefficients, child -> selectorCoefficient];
    ratios = Map[
      Function[item,
       Join[item, <|"normalizedCoefficient" -> Together[
           selectorCoefficients[item["sourceSector"]] item["coefficient"]/selectorCoefficient
           ]|>]
       ],
      candidates
      ];
    nonconstant = If[energyVariables === {}, {},
      Select[ratios, ! FreeQ[Lookup[#, "normalizedCoefficient"], Alternatives @@ energyVariables] &]
      ];
    AppendTo[audits, <|"sectorKey" -> child, "selectorCoefficient" -> selectorCoefficient,
      "chosenTransition" -> chosen, "incomingTransitions" -> ratios,
      "energyIndependentRatiosQ" -> (nonconstant === {}), "nonconstantRatios" -> nonconstant|>];
    If[nonconstant =!= {},
     Return[<|"status" -> "error", "reason" -> "nonconstantContactSelector",
       "sectorKey" -> child, "selectorCoefficients" -> selectorCoefficients, "audits" -> audits|>]
     ],
    {index, 2, Length[sectorOrder]}
    ];
   <|"status" -> "generated", "selectorCoefficients" -> selectorCoefficients, "audits" -> audits,
     "energyVariables" -> energyVariables, "transitions" -> transitions|>
   ];


(* 每个 tree master 的完整 normalization 定义为
   selectorCoefficient * physicalSectorPrefactor。前者化简跨 sector contact，后者保存
   Wronskian 连续幂；所有导数路线都从本记录读取二者，禁止再次从 contact 反推。 *)
dsTreeMasterNormalizationRecords[
   sectorOrder_List,
   selectorData_Association,
   context_Association
   ] := Module[{selectorCoefficients, records, physicalRecord, selectorCoefficient},
   selectorCoefficients = Lookup[selectorData, "selectorCoefficients", <||>];
   records = Association@Table[
      physicalRecord = sectorPrefactorRecordForKey018[context, sectorKey];
      If[! AssociationQ[physicalRecord],
       Return[<|"status" -> "error", "reason" -> "physicalSectorPrefactorMissing",
         "sectorKey" -> sectorKey, "record" -> physicalRecord|>]
       ];
      selectorCoefficient = Lookup[
        selectorCoefficients,
        sectorKey,
        Missing["NoSelectorCoefficient", sectorKey]
        ];
      If[Head[selectorCoefficient] === Missing,
       Return[<|"status" -> "error", "reason" -> "selectorCoefficientMissing",
         "sectorKey" -> sectorKey|>]
       ];
      sectorKey -> Join[
        <|
         "sectorKey" -> sectorKey,
         "selectorCoefficient" -> selectorCoefficient
         |>,
        physicalRecord,
        <|
         "completeNormalization" -> Expand[
           selectorCoefficient physicalRecord["physicalSectorPrefactor"]
           ],
         "normalizationConvention" ->
          "master = selectorCoefficient * normalized J; normalized J = physicalSectorPrefactor * bare integral"
         |>
        ],
      {sectorKey, sectorOrder}
      ];
   <|"status" -> "generated", "records" -> records|>
   ];


(* ::Section::Closed:: *)
(*Block-triangular primitive connection*)

dsTreeVertexSourcePrimitive[vertex_Association] := Module[{states, momentumMagnitudes, cuts},
   states = treeBinaryStates[vertex["p"]];
   momentumMagnitudes = If[vertex["massiveLegs"] === {}, {}, Lookup[vertex["massiveLegs"], "momentumMagnitude"]];
   cuts = Table[
     vertex["signedExternalLegEnergy"] + Sum[(2 states[[row, i]] - 1) momentumMagnitudes[[i]], {i, vertex["p"]}],
     {row, Length[states]}
     ];
   treeTpInverse[vertex] . (-I DiagonalMatrix[Log /@ cuts]) . treeTp[vertex]
   ];


dsTreeAssembleConnection[
   blocks_List,
   families_List,
   sectorOrder_List,
   contactData_List,
   normalizationRecords_Association
   ] := Module[{dimensions, sectorPositions, selectorCoefficients, completeNormalizations,
    omegaBlocks, sourceIndex, targetIndex,
    sourceFamily, vertexIndex, sourcePrimitive, rawMatrix, normalizedMatrix, normalizedContactData},
   dimensions = Lookup[blocks, "masterCount"];
   sectorPositions = AssociationThread[sectorOrder -> Range[Length[sectorOrder]]];
   selectorCoefficients = AssociationMap[
     Lookup[normalizationRecords[#], "selectorCoefficient"] &,
     sectorOrder
     ];
   completeNormalizations = AssociationMap[
     Lookup[normalizationRecords[#], "completeNormalization"] &,
     sectorOrder
     ];
   omegaBlocks = Table[
     If[i === j,
      Expand[
       blocks[[i, "omega"]] +
        Log[completeNormalizations[sectorOrder[[i]]]] IdentityMatrix[dimensions[[i]]]
       ],
      ConstantArray[0, {dimensions[[i]], dimensions[[j]]}]
      ],
     {i, Length[blocks]}, {j, Length[blocks]}
     ];
   normalizedContactData = Map[
      Function[sourceData,
      sourceIndex = sectorPositions[sourceData["sectorKey"]];
      sourceFamily = families[[sourceIndex]];
      Join[sourceData, <|"normalizedMatricesByVertex" -> Association@KeyValueMap[
          Function[{vertexId, targetMatrices},
           vertexIndex = First@FirstPosition[sourceFamily["vertexOrder"], vertexId];
           sourcePrimitive = treeEmbedVertexMatrix[
             dsTreeVertexSourcePrimitive[sourceFamily["vertices"][[vertexIndex]]],
             vertexIndex,
             2^Lookup[sourceFamily["vertices"], "p"]
             ];
           vertexId -> Association@KeyValueMap[
              Function[{targetSector, matrix},
                targetIndex = sectorPositions[targetSector];
                rawMatrix = matrix;
                normalizedMatrix = Expand[
                  selectorCoefficients[sourceData["sectorKey"]]/
                   selectorCoefficients[targetSector] rawMatrix
                  ];
               omegaBlocks[[sourceIndex, targetIndex]] = Expand[
                 omegaBlocks[[sourceIndex, targetIndex]] - I sourcePrimitive . normalizedMatrix
                 ];
               targetSector -> normalizedMatrix
               ],
              targetMatrices
              ]
           ],
          sourceData["matricesByVertex"]
          ]|>]
      ],
     contactData
     ];
   <|"omega" -> ArrayFlatten[omegaBlocks], "omegaBlocks" -> omegaBlocks,
     "normalizedContactData" -> normalizedContactData|>
   ];


dsTreeMultiSectorDLog[context_Association, seedData_: Automatic, auditLevel_: "standard"] := Module[
   {familyContext, families, sectorOrder, blocks, dimensions, masterLists, contactRows, contactData,
    selectorData, normalizationData, normalizationRecords, assembled, completeNormalizations,
    letters, omega, letterMatrices, offsets, taggedMasters,
    bareMasters, dlogResidual, sourceEquations},
   familyContext = dsTreeFamilyContext[context];
   If[familyContext === $Failed, Return[<|"status" -> "error", "reason" -> "treeFamilyInitializationFailed"|>]];
   families = familyContext["families"];
   sectorOrder = familyContext["sectorOrder"];
   blocks = dsTreeDLogBlock018 /@ families;
   If[AnyTrue[blocks, Lookup[#, "status", "error"] =!= "generated" &],
    Return[<|"status" -> "error", "reason" -> "sectorDLogFailed", "sectorBlocks" -> blocks|>]
    ];
   dimensions = Lookup[blocks, "masterCount"];
   masterLists = AssociationThread[sectorOrder -> Lookup[blocks, "masters"]];
   contactRows = MapThread[dsTreeContactRows[#1, #2, context] &, {families, Lookup[blocks, "masters"]}];
   If[AnyTrue[contactRows, Lookup[#, "status", "error"] =!= "generated" &],
    Return[<|"status" -> "error", "reason" -> "contactRowGenerationFailed", "contactRows" -> contactRows|>]
    ];
   contactData = dsTreeContactMatrices[#, sectorOrder, masterLists] & /@ contactRows;
   If[AnyTrue[contactData, Lookup[#, "status", "error"] =!= "generated" &],
    Return[<|"status" -> "error", "reason" -> "contactMatrixGenerationFailed", "contactData" -> contactData|>]
    ];
   selectorData = dsTreeDLogSelectorCoefficients[sectorOrder, contactData, families];
   If[Lookup[selectorData, "status", "error"] =!= "generated",
    Return[Join[<|"status" -> "error", "reason" -> "dlogSelectorFailed"|>, selectorData]]
     ];
   normalizationData = dsTreeMasterNormalizationRecords[sectorOrder, selectorData, context];
   If[Lookup[normalizationData, "status", "error"] =!= "generated",
    Return[Join[<|"status" -> "error", "reason" -> "masterNormalizationFailed"|>, normalizationData]]
    ];
   normalizationRecords = normalizationData["records"];
   assembled = dsTreeAssembleConnection[
     blocks, families, sectorOrder, contactData, normalizationRecords
     ];
   completeNormalizations = AssociationMap[
     Lookup[normalizationRecords[#], "completeNormalization"] &,
     sectorOrder
     ];
   omega = assembled["omega"];
   letters = DeleteDuplicates@Join[
      Flatten[Lookup[blocks, "letters"]],
      DeleteCases[Lookup[completeNormalizations, sectorOrder], 1]
      ];
   letterMatrices = Association@Table[
      letter -> Map[Coefficient[#, Log[letter]] &, omega, {2}],
      {letter, letters}
      ];
   dlogResidual = If[
     auditLevel === "full",
     Expand[omega - Total[MapThread[Log[#1] #2 &, {letters, Values[letterMatrices]}]]],
     Missing["NotRunAtStandardAuditLevel"]
     ];
   offsets = dsTreeMasterSectorOffsets[blocks];
    taggedMasters = Flatten[Map[
      Function[block,
       treeTaggedIntegral[
          block["sector"],
          #,
          normalizationRecords[block["sector"]]["selectorCoefficient"]
          ] & /@ block["masters"]
       ],
      blocks
      ]];
   bareMasters = Lookup[taggedMasters, "integral"];
   sourceEquations = Lookup[contactData, "rowsByVertex"];
   <|
    "status" -> "generated",
    "connectionStructure" -> "sectorDAGBlockTriangular",
    "offDiagonalSourceStatus" -> "assembledFromDirectCompiledTheta",
    "sectorOrder" -> sectorOrder,
    "sectorBlocks" -> blocks,
    "masterSectorOffsets" -> offsets,
    "masters" -> taggedMasters,
    "bareMasters" -> bareMasters,
    "masterCount" -> Length[taggedMasters],
    "omega" -> omega,
    "letters" -> letters,
    "letterMatrices" -> letterMatrices,
    "matrixDimension" -> Dimensions[omega],
    "auditLevel" -> auditLevel,
    "dlogResidual" -> dlogResidual,
    "dlogQ" -> If[Head[dlogResidual] === Missing, dlogResidual,
      TrueQ[dlogResidual === ConstantArray[0, Dimensions[omega]]]],
    "masterNormalizationRecords" -> normalizationRecords,
    "dlogSelectorCoefficients" -> selectorData["selectorCoefficients"],
    "selectorAudits" -> selectorData["audits"],
    "contactMaps" -> assembled["normalizedContactData"],
    "omegaBlocks" -> assembled["omegaBlocks"],
    "sourceEquations" -> sourceEquations,
    "inputSourceData" -> seedData,
    "sourceConvention" -> "Eq. (3.66)-(3.68): -I T^-1 Omega0 T times direct compiled-WT contact selector"
    |>
   ];


(* DSTreeDLogDE 的公开定义由 PublicBoundary.wl 负责把本模块结果映射回统一 J。 *)
