(* ::Package:: *)

(***
File: Numerics.wl
Purpose: Turns a certified connection along a user-given path into a one-dimensional numerical initial-value problem.
Scope: This module does not guess physical boundary conditions; it returns only a Missing status when ordered boundaryValues are absent.
***)

(* ::Chapter:: *)
(* Numerical path system *)

MSNumericalSystem[de_Association, spec_Association : <||>] := Module[
  {parameter, pathRules, substitutions, interval, boundaryPoint, boundaryValues,
   omegaAlongPath, connection, dimension, functions, equations, initialConditions},
  If[Lookup[de, "status", None] =!= "generated" || ! MatrixQ[Lookup[de, "omegaPotential", None]],
    Return[msFailure["InvalidDLogData", <||>]]
  ];
  parameter = Lookup[spec, "pathParameter", t];
  pathRules = Lookup[spec, "pathRules", {}];
  substitutions = Lookup[spec, "substitutions", {}];
  interval = Lookup[spec, "interval", {0, 1}];
  boundaryPoint = Lookup[spec, "boundaryPoint", First[interval]];
  boundaryValues = Lookup[spec, "boundaryValues", Missing["BoundaryData"]];
  dimension = de["masterCount"];
  omegaAlongPath = Simplify[de["omegaPotential"] /. pathRules /. substitutions];
  connection = Simplify[D[omegaAlongPath, parameter]];
  functions = Array[msNumericalMaster, dimension];
  equations = Thread[D[Through[functions[parameter]], parameter] == connection.Through[functions[parameter]]];
  initialConditions = If[
    Head[boundaryValues] === Missing,
    Missing["BoundaryData"],
    If[! ListQ[boundaryValues] || Length[boundaryValues] =!= dimension,
      Return[msFailure["BoundaryVectorDimension", <|"expected" -> dimension, "actual" -> Length[boundaryValues]|>]],
      Thread[Through[functions[boundaryPoint]] == boundaryValues]
    ]
  ];
  <|
    "status" -> If[Head[initialConditions] === Missing, "boundaryDataRequired", "ready"],
    "masters" -> de["masters"],
    "masterDigest" -> de["masterDigest"],
    "pathParameter" -> parameter,
    "interval" -> interval,
    "omegaAlongPath" -> omegaAlongPath,
    "connectionAlongPath" -> connection,
    "equations" -> equations,
    "initialConditions" -> initialConditions,
    "unknownFunctions" -> functions
  |>
];
