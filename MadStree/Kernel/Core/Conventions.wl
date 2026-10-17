(* ::Package:: *)

(***
File: Conventions.wl
Purpose: Centralizes the Pauli atoms and the h-prefactor convention.
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
