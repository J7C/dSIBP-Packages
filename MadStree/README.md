# MadStree

MadStree assembles the de Sitter tree and time-only formula in the Wolfram Language.

Given a graph — tree topologies, time-only looped incidence graphs, or a standalone
single-vertex function family — MadStree derives the contact-reachable sectors, the master
integrals in exact matrix order, the recurrence metadata that reduces any shifted integral, and
the analytic `dlog`-form differential equation. It then builds the automatic boundary from the
`k0 -> Infinity` Frobenius formula and evaluates the resulting initial-value problem to
arbitrary precision by handing the equation and the boundary vector to FlintNDE.

MadStree assembles the formula; it does not generate general IBP systems or reduction input.
That is dSIBP's job, and the two packages speak through `MSToDSIBPJ` and `MSFromDSIBPJ`.

## Loading

```wl
AppendTo[$Path, "/absolute/path/to/dSIBP-Packages/MadStree"];
Needs["MadStree`"];
```

```
MadStree/
  MadStree.m               top-level loader
  Kernel/                  module sources
  Backend/                 the Python transport entry point MadStree launches
  Vendor/FlintNDE/         bundled snapshot of the FlintNDE Python package
  Examples/                seven runnable cases
  Documentation/           handbook (.tex and .pdf, English and Chinese)
```

MadStree runs its numerical back end as a Python child process against the FlintNDE snapshot in
`Vendor/FlintNDE`, so no separate installation is required. `MSFlintNDEConfiguration[]` reports
the resolved backend directory and whether the bundled `flintnde` Python sources are present; it
does not start an interpreter, so a missing or unusable Python environment — including
`python-flint` — surfaces as an error from the numeric entry points at run time.
`MSSetFlintNDERelativePath[path]` points the package at a different copy. The result cache is
keyed by the content identity of the backend sources, so a changed vendor copy never silently
reuses cached numbers.

## From a graph to a number

```wl
context = MSInitTree[spec];                  (* or MSInitTimeGraph / MSInitVertexFamily *)

MSDLogDE[context]                            (* masters + block-triangular dlog connection *)
MSFormulaData[context]                       (* masters, recurrence metadata, full analytic DE *)
MSReduce[expr, context, MasterBasis -> basis]

boundary = MSBoundaryData[context, targetRules];

result = MSEvaluatePath[
  context,
  {{k0, k1}, {1, 1}, {11/10, 1}, {{6/5, 1}, "tmp"}},
  ParameterRules -> {q -> 1, nu -> 1/2, a -> 2},
  FlintNDEPathPlanning -> True
];
```

`MSEvaluatePath` is the only numerical entry point, for one point or many. The first row of
`pointSequence` is the ordered coordinate-symbol header, every later row is an equal-width point
of the same schema; a single point is just a header plus one value row. Parameters that are not
used as running coordinates are given once in `ParameterRules`, and the header plus
`ParameterRules` must together cover every symbol the analytic DE and the boundary require.
A row written as `{{values...}, "tmp"}` is a transient waypoint: it takes part in segment
recognition and transport but is not returned. If anything is missing or stays non-numeric after
substitution, the call returns a `Failure` naming the offending symbols and states that the
numerical transport never started.

`MSBoundaryData` and `MSEvaluatePath` both default to 200 digits of `WorkingPrecision`.
Before Python starts, every numerical entry writes or reuses the complete analytic formula
artifacts of the same context — masters, recurrence metadata, the analytic `dlog` DE and a
manifest — and returns their actual paths under `"formulaArtifacts"`; if that write fails, the
transport does not start. `MSWriteFormulaArtifacts[context]` does the same step on its own for
purely symbolic work, and `MSExportEvaluationData[evaluation]` writes saved ordinary-point
records to CSV and JSON.

## Path planning and singularities

MadStree identifies the longest consecutive complex-affine segments `x(s) = x0 + s v` in the
input order, pulls the `dlog` DE back once per segment, and passes all segments and boundary
data to a single Python process. It does not choose transport nodes itself.

- `FlintNDEPathPlanning -> True` (default): FlintNDE plans the nodes. User points covered by the
  same expansion node are evaluated as a bucket through fast multipoint evaluation.
- `FlintNDEPathPlanning -> False`: every supplied point is used as a node, in order, with no
  insertion, removal or silent replanning. The point list must stay inside the successive
  convergence discs and must not cross a singularity.

`SingularityMode -> "Automatic"` (default) honours singularities the user has stated
explicitly: one local basis covers both sides of the singularity inside its convergence domain,
a genuinely divergent value comes back as the text `Infinity`, and transport continues from the
outgoing ordinary point. A singular endpoint gets a hidden matching point when it is needed.
`"Avoid"` refuses any crossing, and `"SingularityJump"` chooses the crossing branch explicitly
and requires the multivalued ambiguity of the equivalent detour class to be accepted.
Non-collinear turns at a singularity, and interior singularities with planning switched off, fail
closed rather than being guessed.

## Regulator series

Limits carrying the common normalization parameter `ep` go through one entry point:

```wl
series = MSReconstructEpSeries[
  context, ep, pointSequence,
  ParameterRules -> {q -> 1, a1 -> 1 + ep, a2 -> 1 + ep},
  MaximumEpPower -> 0,
  EpGoalDigits -> 20,
  ParallelTaskCount -> 12
];
```

`MaximumEpPower -> 0` asks for coefficients through `ep^0`. Before any numerical transport starts,
the package certifies the lowest integer power from the actual symbolic boundary conditions and
the `dlog` DE, and hands that certificate to the sampling planner — no leading power is ever
guessed from a terminal numerical pilot. It then chooses production points, independent validation
points, working precision and transport orders automatically. You can instead supply
`EpSamplePoints` as an ordered redundant candidate pool and `EpValidationPoints` as fixed
independent checks, or give only `EpSampleAngleRange -> {thetaMin, thetaMax}` and let the package
pick at most three interior complex rays at magnitudes set by the target accuracy. Candidates are
consumed incrementally and old values are reused; the pool is never extended past the range you
gave, and if it runs out the current best coefficients are returned with an uncertified flag. The
working precision is never below 200 digits, and `ParallelTaskCount` (default 12) bounds the
number of concurrent `ep` processes, not threads inside one Python process.

## Where output goes

`MSRuntimeDirectory -> Automatic` makes the calling script's directory the run root and writes
scratch under its `results_temp/`; an explicit path is used as the run root itself. Transient
inputs and logs live in `nde/`, successful cache entries in `cache/`. Durable products are
written by the caller under `results/` — `results/madstree_formula/run-<UUID>/` for the analytic
artifacts, `results/madstree_evaluation/run-<UUID>/` for exported evaluation data. Nothing is
written into the package source tree.

## Examples

| Example | Covers |
| --- | --- |
| `01_massless_full_edge.wl` | A single theta-carrying massless full edge: initialization to masters, recurrence and the `dlog` DE, then single-point and multipoint numerical evaluation. |
| `02_single_vertex_family.wl` | A single-vertex function family in both compact and explicit form, analytic saving, masters, reduction, and one table-driven numerical run. |
| `03_time_only_cycle_chart.wl` | Time-only cycle initialization, the common-theta contact sector, the `dlog` DE and all strict time-rank blow-up chart certificates. |
| `04_three_vertex_tree.wl` | The minimal massless three-vertex, two-propagator tree (`+++` contour) through to a shared-anchor batch multipoint evaluation with CSV/JSON export. |
| `05_massive_three_vertex_tree.wl` | Three-vertex tree with two massive propagators: masters, recurrence, `dlog` DE, and two separately planned point lists sharing one context. |
| `06_massless_three_vertex_ep_regularization.wl` | Common time-power analytic regularization `a1 = a2 = a3 = 1 + ep`: certification of the leading power, exact-`ep` sampling and Laurent reconstruction. |
| `07_zero_external_leg_energy.wl` | Vertices without a physical external leg: omitted or explicit `0` energy, and transport to the physical zero point. |

```powershell
wolframscript -file Examples/01_massless_full_edge.wl
```

## Documentation

`Documentation/tree_formula_en_full.pdf` is the handbook in English; `tree_formula.pdf` is the
same handbook in Chinese. Both cover the conventions, the integral normalization, the formula
construction, the boundary certificates and the numerical interface, and both are built from the
`.tex` sources in the same directory with `references.bib`.
