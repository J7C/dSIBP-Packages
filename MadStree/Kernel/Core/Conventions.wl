(* ::Package:: *)

(***
File: Conventions.wl
Purpose: Centralizes the Pauli atoms, the h-prefactor convention and the local H/h transformations.
Conventions: The package's baseline formulas are written in the +nu (s_nu=+1) convention: h=z^(s_nu nu) H_nu, where the input nu is used verbatim as the Hankel order (real, pure imaginary i mu or complex; no |nu| and no Re is ever taken). The Hankel order never changes with the prefactor convention; only the sign s_nu carried by formulaNu:=s_nu nu does. "Positive" (default) gives formulaNu=+nu and "Negative" gives formulaNu=-nu; the printed arXiv:2401.00129 / arXiv:2411.03088 shape is recovered from the baseline by the replacement nu_ref=-nu.
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
  (* ToString 是 HoldAll，这里不再额外包 HoldForm，避免打印出 HoldForm[...] 包装。 *)
  Message[MSConvertBasis::unsupported,
    ToString[direction, InputForm] <> "（对象 object: " <>
      ToString[object, InputForm] <> "）"];
  msFailure["UnsupportedBasisConversion", <|"direction" -> direction|>]
);


(* ::Section:: *)
(* Function-system presets *)

msFunctionSystemPreset["h", nu_, variable_: z, nuConvention_: "Positive"] := Module[{formulaNu},
 formulaNu = If[nuConvention === "Positive", nu, -nu];
 <|
  "basis" -> "h",
  "nu" -> nu,
  "formulaNu" -> formulaNu,
  "nuConvention" -> nuConvention,
  "variable" -> variable,
  (* +nu baseline EOM h''+(1-2 nu) h'/z+h=0: the P coefficient is (1-2 formulaNu)/variable. *)
  "P" -> (1 - 2 formulaNu)/variable,
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
