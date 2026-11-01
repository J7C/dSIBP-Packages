# dSIBP-Packages

Packages for integration-by-parts (IBP) reduction and differential equations in de Sitter
spacetime: **MadStree**, **dSIBP** and **FlintNDE**.

| Package | Language | Role |
| --- | --- | --- |
| [dSIBP](dSIBP-1.0/) | Wolfram Language | Topology-driven generation of dS IBP relations, parameter-derivative seeds, sector metadata and backend-neutral reduction input for graphs with any number of loops, any topology and mixed massive/massless propagators. |
| [MadStree](MadStree-1.01/) | Wolfram Language | Assembly of the dS tree and time-only formula: master integrals, recurrence metadata, the analytic `dlog`-form differential equation, and high-precision evaluation of the resulting boundary-value problem. |
| [FlintNDE](FlintNDE-1.01/) | Python, with a Wolfram Language front end and an optional C++ engine | Arbitrary-precision numerical transport of first-order matrix differential equations `Y'(x) = A(x) Y(x)`, plus reconstruction of half-analytic Laurent/power series in a regulator. An optional C++/FLINT engine runs the fixed-regulator batch contract with bit-identical outputs. |

## How the packages fit together

The dependency direction is fixed: **MadStree → FlintNDE**.

MadStree turns a graph into an analytic `dlog` differential equation and an exact starting
vector, then hands both to FlintNDE, which does the numerics. MadStree ships its own bundled copy
of FlintNDE under `MadStree-1.01/Vendor/FlintNDE`, so it works without a separate installation of
the Python package.

dSIBP calls neither of the other two. It generates and serializes relations; it never runs a
reduction and never reads the other packages. The three packages do talk to each other where
the representations overlap: `MSToDSIBPJ` and `MSFromDSIBPJ` convert the time-only integral
object losslessly between MadStree and dSIBP.

FlintNDE is a general single-variable matrix ODE solver. It knows nothing about de Sitter
graphs, master-integral ordering or normalization; turning a physics problem into a
`RationalMatrixSystem` is the caller's job, and MadStree is that job's reference caller.

## Repository layout

```
dSIBP-Packages/
  MadStree-1.01/   Kernel/, Backend/, Vendor/FlintNDE/, Examples/, Documentation/
  dSIBP-1.0/     Kernel/, Examples/, Documentation/
  FlintNDE-1.01/   flintnde/, Mathematica/, Kernel/, engine/, examples/, config/, Documentation/
```

## Current versions

| Package | Development version | Release version | Release directory |
| --- | --- | --- | --- |
| dSIBP | `030.0` | `1.0` | `dSIBP-1.0/` |
| MadStree | `v1.01` | `1.01` | `MadStree-1.01/` |
| FlintNDE | `v1.01` | `1.01` | `FlintNDE-1.01/` |

Each package is self-contained and is loaded from its own directory; nothing is read through a
shared parent path.

## Requirements

- **Wolfram Language** for dSIBP and MadStree, and for the Wolfram front end of FlintNDE.
  `wolframscript` is enough to run every example in this repository.
- **Python 3.10 or newer** with [`python-flint`](https://github.com/flintlib/python-flint)
  `>= 0.6` and `sympy >= 1.12` for the numerical transport used by FlintNDE and MadStree.
  Installing FlintNDE once pulls both in:

  ```powershell
  python -m pip install path/to/FlintNDE-1.01
  ```

  MadStree launches the Python back end as a child process from `MadStree-1.01/Backend/` against
  its bundled `Vendor/FlintNDE`, so `python -m pip install` is optional for MadStree users, but
  a working `python-flint` in the interpreter is not.
- The FlintNDE C++ engine is optional and not needed for any Python or Wolfram entry point.
  Building it requires an MSYS2 UCRT64 toolchain (`g++`, `make`) and the FLINT 3.x development
  library; on Windows the `libflint` runtime must be reachable through PATH or the
  `"EngineLibraryPath"` option. See `FlintNDE-1.01/README.md` for the build and usage.
- External IBP reduction (Kira) is not part of any package. dSIBP emits reduction input and
  reads results back; the solver itself runs outside this repository.

## Quick start

Load one package per code block; each snippet is independent.

```wl
(* dSIBP *)
AppendTo[$Path, "/absolute/path/to/dSIBP-Packages/dSIBP-1.0"];
Needs["dSIBP`"];
DSPublicAPI[]                      (* the full public surface, grouped by workflow *)

(* MadStree *)
AppendTo[$Path, "/absolute/path/to/dSIBP-Packages/MadStree-1.01"];
Needs["MadStree`"];
```

```python
# FlintNDE
import flintnde
```

MadStree reaches FlintNDE through its own bundled copy, so a MadStree session never needs
`Needs["FlintNDE`"]`. If you do load the standalone FlintNDE Wolfram front end in the same
kernel as MadStree, three option names — `ParallelTaskCount`, `MessageLanguage` and
`SingularityMode` — exist in both contexts and the kernel reports shadowing; the unqualified
name then resolves to whichever package was loaded last. Write `MadStree`ParallelTaskCount` (or
load MadStree last) when you pass an option to a MadStree entry point.

Each package prints a citation reminder on its first successful load. The package release
identifiers are `1.0` for dSIBP and `1.01` for MadStree and FlintNDE. At run time,
`dSIBP`$dSIBPVersion` reports `1.0`, while `MadStree`$MadStreeVersion` and
`FlintNDE`$FlintNDEVersion` report `1.01`.
Source and documentation are UTF-8, and every internal `Get` passes the encoding explicitly, so
no extra option is needed when loading a package.

## Documentation

User manuals and formula references are compiled PDFs next to their editable `.tex` sources:

- `dSIBP-1.0/Documentation/dSIBP_manual-cn.pdf` — the Chinese dSIBP manual: conventions, the public
  interface and a worked walk-through of every example.
- `dSIBP-1.0/Documentation/dSIBP_manual-en.pdf` — the English “dSIBP Handbook” covering the same
  material.
- `MadStree-1.01/Documentation/handbook.pdf` — the MadStree handbook in English;
  `用户手册.pdf` is the same handbook in Chinese.
- `FlintNDE-1.01/Documentation/FlintNDE_manual-en.pdf` — the FlintNDE paper: the numerical method,
  the public interface and its documented limits; `FlintNDE_manual-cn.pdf` is the Chinese edition.

## Examples

Every example resolves its own package directory from `$InputFileName` or from a single path
variable at the top of the file, so all of them run directly from a clone of this repository:

```powershell
wolframscript -file MadStree-1.01/Examples/01_massless_time_ordered_edge.wl
wolframscript -file dSIBP-1.0/Examples/01_mixed_bubble_workflow.wl
python FlintNDE-1.01/examples/qnm_2x2.py
```

See `MadStree-1.01/Examples/README.md`, `dSIBP-1.0/Examples/README.md` and `FlintNDE-1.01/examples/README.md`
for what each case covers. Two of the dSIBP cases (`04`, `06`) drive an external Kira reduction
and expect the environment variable `DSIBP_KIRA_WORKSPACE` to point at a writable workspace
outside this repository; they refuse to run without it, which is the intended guard, not a
failure.

Run output is written next to the calling script under `results/` (durable) and `results_temp/`
(scratch); nothing is written into a package source directory.

## Citing

Please cite the underlying papers when you use these packages:

- Jiaqi Chen and Bo Feng, *Towards Systematic Evaluation of de Sitter Correlators via
  Generalized Integration-by-Parts Relations*, [arXiv:2401.00129](https://arxiv.org/abs/2401.00129).
- Jiaqi Chen, Bo Feng and Yi-Xiao Tao, *Multivariate hypergeometric solutions of cosmological
  (dS) correlators by d log-form differential equations*,
  [arXiv:2411.03088](https://arxiv.org/abs/2411.03088).
- Jiaqi Chen, Bo Feng, Zhehan Qin and Yi-Xiao Tao, *Loop integrals in de Sitter spacetime: The
  parity-split IBP system and d log-form differential equations*,
  [arXiv:2604.14549](https://arxiv.org/abs/2604.14549).
- Jiaqi Chen, Yi-Xiao Tao, Xinyi Wang and Dong Zhang, *MadStree: automated evaluation of
  tree-level de Sitter correlators*, [arXiv:2610.28526](https://arxiv.org/abs/2610.28526).

The dSIBP and FlintNDE package papers are in preparation and their arXiv identifiers are not yet
assigned.
