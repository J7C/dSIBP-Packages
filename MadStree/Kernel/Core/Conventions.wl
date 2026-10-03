(* ::Package:: *)

(***
File: Conventions.wl
Purpose: Centralizes the Pauli atoms, the h-prefactor convention and the local H/h transformations.
Conventions: nu=|nu| is always enforced. The default is h=z^nu H_nu; the negative prefactor h=z^(-nu) H_nu corresponds to 2401.
***)

(* ::Chapter:: *)
(* Local two-state atoms *)

msIdentity2 = IdentityMatrix[2];
msSigma1 = {{0, 1}, {1, 0}};
msSigma2 = {{0, -I}, {I, 0}};
msSigma3 = {{1, 0}, {0, -1}};
msProjector1 = (msIdentity2 - msSigma3)/2;

msPaperT = 1/Sqrt[2] {{1, -I}, {-I, 1}};
msPaperTInverse = 1/Sqrt[2] {{1, I}, {I, 1}};
msHadamard = 1/Sqrt[2] {{1, 1}, {1, -1}};


(* ::Section:: *)
(* H/h forward and inverse transformations *)

(* Private constructors consume only already-resolved convention strings. *)
msHTohMatrix[nu_, z_, nuConvention_String] := Module[{power},
  power = If[nuConvention === "Positive", nu, -nu];
  {{z^power, 0}, {power z^(power - 1), z^power}}
];

msHToHMatrix[nu_, z_, nuConvention_String] := Module[{power},
  power = If[nuConvention === "Positive", nu, -nu];
  {{z^(-power), 0}, {-power z^(-power - 1), z^(-power)}}
];

(* Public local matrices must read the initialized context; no per-call convention option is provided. *)
MSHTohMatrix[nu_, z_, context_?MSContextQ] := msHTohMatrix[
  nu, z, context["convention"]["nuConvention"]
];

MShToHMatrix[nu_, z_, context_?MSContextQ] := msHToHMatrix[
  nu, z, context["convention"]["nuConvention"]
];

MSHTohMatrix[___] := msFailure["InitializedContextRequired", <|"function" -> "MSHTohMatrix"|>];
MShToHMatrix[___] := msFailure["InitializedContextRequired", <|"function" -> "MShToHMatrix"|>];

MSConvertBasis[
  vector_List,
  direction_Rule,
  nu_,
  z_,
  context_?MSContextQ
] /; Length[vector] === 2 := Module[
  {nuConvention = context["convention"]["nuConvention"]},
  Switch[
    direction,
    "H" -> "h", Simplify[msHTohMatrix[nu, z, nuConvention].vector],
    "h" -> "H", Simplify[msHToHMatrix[nu, z, nuConvention].vector],
    _, Message[MSConvertBasis::unsupported, direction]; msFailure["UnsupportedBasisConversion", <|"direction" -> direction|>]
  ]
];

MSConvertBasis[object_, direction_, ___] := (
  Message[MSConvertBasis::unsupported, <|"object" -> HoldForm[object], "direction" -> direction|>];
  msFailure["UnsupportedBasisConversion", <|"direction" -> direction|>]
);


(* ::Section:: *)
(* Function-system presets *)

msFunctionSystemPreset["h", nu_, variable_: z, nuConvention_: "Positive"] := Module[{formulaNu},
 formulaNu = If[nuConvention === "Positive", -nu, nu];
 <|
  "basis" -> "h",
  "nu" -> nu,
  "formulaNu" -> formulaNu,
  "nuConvention" -> nuConvention,
  "variable" -> variable,
  "P" -> (2 formulaNu + 1)/variable,
  "Q" -> 1,
  "hToH" -> msHToHMatrix[nu, variable, nuConvention],
  "HToh" -> msHTohMatrix[nu, variable, nuConvention]
 |>
];

msFunctionSystemPreset["H", nu_, variable_: z] := <|
  "basis" -> "H",
  "nu" -> nu,
  "variable" -> variable,
  "P" -> 1/variable,
  "Q" -> 1 - nu^2/variable^2,
  "hToH" -> msHToHMatrix[nu, variable, "Positive"],
  "HToh" -> msHTohMatrix[nu, variable, "Positive"]
|>;
