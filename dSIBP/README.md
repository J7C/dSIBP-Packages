# dSIBP

dSIBP is a topology-driven Wolfram Language generator of integration-by-parts (IBP) relations
for correlators in de Sitter spacetime.

You describe a graph — any number of loops, any topology, with massive and massless internal
lines mixed — and dSIBP produces the complete relation system for it: time and loop IBP
equations, parameter-derivative seeds, sector metadata, and a backend-neutral `linearData`
representation that can be serialized for Kira or any other linear reducer. dSIBP generates and
exports those relations; it never runs a reduction itself.

## Loading

```wl
AppendTo[$Path, "/absolute/path/to/dSIBP-Packages/dSIBP"];
Needs["dSIBP`"];
```

The layout is the standard Wolfram module tree:

```
dSIBP/
  dSIBP.m          top-level loader
  Kernel/          module sources
  Examples/        runnable cases
  Documentation/   user manual (.tex and .pdf)
```

The first successful load prints a citation reminder. Source files are UTF-8; the loader passes
the encoding explicitly, so no option is needed.

## Public interface

`DSPublicAPI[]` is the authoritative listing: it returns every public symbol grouped by
workflow, together with the default options of each high-level entry point.

| Group | Symbols |
| --- | --- |
| Initialization | `DSInit`, `DSInfo`, `DSKinematics`, `DSParameterNotation`, `DSRedefineParameters` |
| Loop workflow | `DSSeeds`, `DSAllSeeds`, `DSSeedGroups`, `DSSeedGroupMetadata`, `DSMetaSeedRange`, `DSGenerateIBP`, `DSLinear`, `DSReorderIntegrals`, `DSUserMI`, `DSKiraPlan`, `DSKiraExport`, `DSKiraImport`, `DSDE`, `DSScaleCheck` |
| Pure-time workflow | `DSTreeSeeds`, `repIterative`, `DSTreeNaiveIBP`, `DSTreeNaiveDE`, `DSTreeDLogDE` |
| Atomic operations | `dtau`, `dqq`, `dqk`, `ds`, `rep2innerform`, `rep2outform`, `rep2Integrand`, `symmetry`, `repSymmetry0` |
| Messages | `DSMessagesOn`, `DSMessagesOff`, `DSMessagesQ` |
| Introspection | `DSPublicAPI` |

All sectors share one integral head, `J`. In full loop mode an integral is
`J[aList, linePacks, ispList]`; the time-only mode exposes the compact
`J[sectorKey, timeShifts, stateBits]`, where `sectorKey` is a fixed-length `0`/`1` string in
root-propagator order and the remaining slots carry the time powers and the discrete
building-block state. Loop momenta, ISPs and kinematics are declared explicitly — you supply
`loopExternalMomenta` and `independentExternalMomenta` separately, and dSIBP does not infer
either from symbol names.

## Reduction is external

`DSKiraPlan` and `DSKiraExport` write a Kira-ready input tree; `DSKiraImport` reads the
reduction table back and `DSDE` turns it into the differential equation. The reducer itself runs
outside this repository, and the examples that use it require the environment variable
`DSIBP_KIRA_WORKSPACE` to point at a writable workspace **outside** the package tree — they
refuse to run otherwise. Only the generating scripts and lightweight input/result summaries are
shipped here; no `init/`, `kira/`, database, save, log, reduction table, cache or DE run
directory is part of the package.

## Examples

`Examples/` contains six complementary cases, all loaded through
`Examples/load_current_package.wl`, which resolves the source tree next to it:

| Example | Covers |
| --- | --- |
| `01_mixed_bubble_workflow.wl` | Mixed massive/massless bubble through the public entry points: discrete templates, continuous sample points, `linearData`. Builds the reduction plan in memory only. |
| `02_function_system_hankel.wl` | Hankel / function-system input and the associated state transforms. |
| `03_single_massive_sunrise/` | Two-loop single-massive sunrise: general seeds and the `{ss11, kE}` parameter-derivative operators. |
| `04_pure_massive_bubble_closed_loop/` | Pure massive bubble end to end: formal Kira input, external reduction, readback, the 19-dimensional DE and the scaling check. |
| `05_tree_two_vertex_time_ibp/` | Two-vertex tree time IBP, the naive DE and the formula route. |
| `06_mix_bubble_tree/` | Mixed cycle/bridge topology: momentum roles, massless convention, contact contraction, and the 81-master Kira/DE record. |

```powershell
wolframscript -file Examples/01_mixed_bubble_workflow.wl
```

`03` deliberately stops at seeds and operators; it never enters sampled relations, Kira, DE or
scaling. `04` and `06` are the two cases that drive an external reduction.

## Documentation

`Documentation/dSIBP_user_manual.pdf` is the user manual: conventions, the whole public
surface, and a walk-through of each example. The editable `dSIBP_user_manual.tex` sits next to
it along with the `jheppub` class files it needs; compile it from inside `Documentation/`.

dSIBP is independent of MadStree and FlintNDE. The time-only `J` representation it publishes is
the same one MadStree reads and writes through `MSToDSIBPJ` and `MSFromDSIBPJ`.
