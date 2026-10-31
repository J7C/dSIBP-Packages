# MadStree changelog

## 1.01

Fixed the handoff between the automatic infinity boundary and fixed external energies supplied via `ParameterRules`: every external energy is now transported from its boundary matching value to the first user point before being held at the requested value on the scan; affected calls must be recomputed.

Updated the MadStree package-paper citation to arXiv:2610.28526 in the load reminder and bilingual README files.
