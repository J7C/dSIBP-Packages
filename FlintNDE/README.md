# FlintNDE

FlintNDE is a Python package for arbitrary-precision numerical work with first-order matrix
differential equations

```
Y'(x) = A(x) Y(x)
```

built on [`python-flint`](https://github.com/flintlib/python-flint) ball arithmetic. Everything
is carried as `arb`/`acb` balls, so a computed value comes with a radius that is valid at the
precision you asked for rather than at machine precision.

Two capabilities sit at the top:

- **Numerical transport.** Given a system, a starting vector and a path, FlintNDE moves the
  solution between points to a requested relative error, choosing series orders and step sizes
  automatically.
- **Power-series reconstruction in a regulator.** When the system and its boundary depend on a
  regulator `ep`, FlintNDE solves at a set of exact `ep` values and rebuilds the half-analytic
  answer — analytic in the power of `ep`, numerical in its coefficients — as a Laurent/power
  series, including the adaptive sampling needed to certify the leading power.

Both rest on reusable machinery: exact `Q(i)(x)` matrix entry reduction and singularity
discovery, Cauchy–DFT Taylor reconstruction at ordinary points, matrix indicial equations with
Jordan and resonance gates, power-log local bases at regular singular points, exact Lee–Moser
projector balancing for higher-order poles, certified decoupled exponential sectors, and bounded
multiprocess execution over independent `ep` values.

FlintNDE is a general single-variable ODE solver. It has no notion of graphs, master integrals or
physical normalization; building `A(x)` and the starting vector is the caller's job. MadStree is
that job's reference caller and ships its own snapshot of this package.

## Install and requirements

Python 3.10 or newer, `python-flint >= 0.6`, `sympy >= 1.12`:

```powershell
python -m pip install path/to/FlintNDE
```

`pip` pulls the two dependencies in. After one install, any script anywhere can simply
`import flintnde`; the examples in this directory instead set a `PACKAGE_ROOT` variable at the
top so they run straight from a clone with no installation.

```
FlintNDE/
  pyproject.toml
  flintnde/            the Python package (import name: flintnde)
  FlintNDE.m           Wolfram loader
  Kernel/init.m
  Mathematica/         Wolfram front end
  examples/            seven runnable cases
  config/              configuration file used by the QNM example
  Documentation/       the paper (.tex and .pdf)
```

## Use from Python

```python
from flintnde import (
    RationalMatrixSystem, rational_function, column_vector,
    configure_working_precision, build_adaptive_path, transport_path_refined,
)

configure_working_precision(200)

system = RationalMatrixSystem(
    ((rational_function(("1/5",), (1, 1)),),),   # A(x) = 1/5/(1+x)
    variable_name="x",
)
path = build_adaptive_path(system, 0, 1)
result = transport_path_refined(system, column_vector([1]), path,
                                primary_order=48, reference_order=64,
                                target_relative_error="1e-18")
```

Default working precision is 200 decimal digits, and the internal bit count is
`ceil(WorkingPrecisionDigits * log2(10)) + 32`.

Exact input goes through `RationalMatrixSystem`: each matrix entry is a `rational_function` whose
numerator and denominator coefficients are ordered by increasing power of the variable, and each
exact complex number may be given as an integer, an `fmpq`, a fraction or decimal string, an
`"a+b*I"` string, an `(real, imag)` pair or a `{"real": ..., "imag": ...}` dict. The package
first reduces every entry exactly in `Q(i)[x]`, so poles that cancel inside one entry never
reach the singularity list, then square-free-factors the unique denominators and isolates all
complex roots with `acb` balls. `A(x) = P(x) + Σ R_j/(x - p_j)` with any finite-degree `P` is
recognized and takes the fast residue-recursion route automatically; you never pre-extract
singularities or supply `dlog` letters. Python `float`, `complex`, `arb` and `acb` values are not
accepted on the exact route — use the numerical system and state the reliable input precision.

## Paths, save points and singularities

`build_adaptive_path` plans, `transport_path_refined` executes. Points where you want the answer
are marked directly in the path, with no separate naming step:

```python
path = [(z0, "save"), z1, (z2, "save")]
result = transport_path_refined(system, boundary, path)
```

Each save point is written as soon as it is reached (`flintnde_save_001.json`, …), and
`flintnde_save_points.json` summarizes the run only after the whole path succeeds. Only the
primary chain is saved; the reference chain exists to verify the error estimate. Ordinary points
record the coordinate and the `acb` column vector, regular singular points record a reusable
`{a, b, C}` boundary, and certified strictly decoupled exponential singularities record
`{phi, a, b, C}` for `exp(phi(z)) z^a log(z)^b C`.

The convergence radius at an ordinary point is the distance to the nearest singularity; at a
singularity it is the distance to the nearest *other* singularity, and `max_step_over_radius`
tightens the step relative to that radius. By default `singularity_mode="avoid"` refuses any
segment crossing a singularity and returns a structured rejection naming the singularity and the
segment; supply `detour_points` to route around it yourself. Only an explicit
`singularity_mode="singularity_jump"` builds the two-sided matching points and the local-basis
bridge, and it requires you to accept the multivalued branch implied by the equivalent detour
class. A non-collinear turn at a singularity does not uniquely specify a branch and fails closed.

Call `build_adaptive_path_plan` before committing to a run: it executes the same local scheduler
on the start, the target and every interior singularity and returns `continuation_ready` with the
classification, the certified method and the rejection reason for each singularity. When it says
`False` — and when `build_adaptive_path` raises `LocalReductionError` — the current package
cannot complete this path; change the detour points or supply external local/Stokes data. The
package never quietly degrades to a plain Taylor route after an exact gate, a Lee–Moser step or a
formal gate fails.

## Boundaries at singular points

A path that starts at a regular singular point takes its boundary as
`frobenius_boundary([{"a": a, "b": b, "C": [...]}], ...)` describing the leading
`z**a log(z)**b C` term, with `z = s - s0` at a finite point and `z = sinv = 1/s` at infinity.
The package verifies `a`, `b` and the full coefficient vector `C` exactly against the indicial
and generalized-eigenspace structure, generates the lower log completions a Jordan branch needs,
and initializes transport at the first ordinary matching point. Malformed input, a `C` outside
the required eigenspace, or a finite ordinary value passed off as a singular boundary is refused
before any computation. Exponential singularities accept only the explicit
`exponential_boundary([{"phi": [{"power": -1, "coefficient": -k}], "a": a, "b": b, "C": [...]}])`
format, where `phi` contains negative integer powers only and `phi=[]` means `exp(phi) = 1` for
that sector.

`transport_frobenius_boundaries_refined` moves several Frobenius initial vectors of the same
exact system at once: each chain builds its local basis once and each ordinary segment its Taylor
matrix once, with the initial values carried as columns. Per-column
`relative_differences_inf` and `target_relative_error_met` are still reported. The batch entry
requires every transition after the singular start to be ordinary Taylor and does not handle
`save` labels.

## Regulator series

`reconstruct_series_solution` is the high-level entry. You hand it the DE matrix, the boundary,
the path, the target `maximum_power` and the symbolically certified integer `leading_power` with
its certificate. There is no numerical pilot for the leading power, because a bare callable
carries no symbolic information about its Laurent support — the upstream symbolic stage is
expected to certify it. Production and independent validation points both reuse the same
transport routine, and production samples are fitted through FLINT `acb` square-matrix
interpolation, never least squares or a pseudo-inverse. Sample counts, parameter magnitudes and
working precision follow published empirical formulas from AMFlow 2.0; inside the package the
feature is called power-series reconstruction and no AMFlow-named interface is exposed.

## Output layout and parallel work

`initialize_output_layout(__file__, run_name=...)` puts everything it manages under
`results/<run_name>/` next to the calling script, split into `configuration/`, `singularities/`,
`frobenius/`, `transport/`, `regularization/` and `summary/`. The package never writes into its
own installation directory or into an unpredictable current working directory.

`run_ep_tasks(ep_values, solve_ep, parallel_task_count=12)` runs independent fixed-`ep` jobs in a
bounded process pool; the effective worker count is `min(parallel_task_count, number of tasks)`
and queued tasks start automatically as workers finish. This is task-level parallelism, distinct
from the `ctx.threads` setting of a single python-flint process. Task functions must live at
module top level so Windows spawn can pickle them.

## Use from Wolfram Language

```wl
AppendTo[$Path, "/absolute/path/to/dSIBP-Packages/FlintNDE"];
Needs["FlintNDE`"];
```

`FlintNDERationalSystem`, `FlintNDEPlanPath` and `FlintNDEExecutePath` mirror the Python plan-then-
execute split, and `FlintNDEEvaluateEpBatch` covers the fixed-`ep` pool. `"WorkDirectory"` is the
temporary run root itself; `Automatic` uses `results_temp/` in the calling directory and the
interface appends only a `bridge/` segment there. The bridge launches Python through an argument
list `RunProcess` — no shell string, no quoting helper, no redirection, no retry fallback — and
keeps the real failure boundaries apart: `RuntimePathTooLong`, `RuntimeInputWriteFailed`,
`PythonFlintUnavailable`, `BridgeLaunchFailed`, `BridgeOutputMissing`. `MessageLanguage` selects
`"EN"` (default) or `"CN"` runtime messages.

## Current limits

Power-log transport is the foundation and is complete for what it claims: exact
`Q(i)(x)` input that Lee–Moser projector balances reduce to simple poles, plus systems whose
higher-order Laurent coefficients are diagonal in the same exact constant basis and fully
decoupled across distinct exponential signatures, reducing to at most a simple pole once the
exponential factor is pulled out. For `A_{-2} = k` the factor is `exp(-k/z)`; `k = 0` means only
that a sector has no such scalar exponent, it does not relax the nilpotent or off-diagonal
higher-order constraints. Pulling out a certified scalar higher part and re-checking the residual
system means an `exp(-k/z)` result is not automatically asymptotic: when the residual is at most
a simple pole, the convergent Frobenius power-log series is used, with its radius still set by
the nearest other singularity. A residual with higher non-diagonal terms, or one that needs an
infinite formal gauge to eliminate order by order, is Gevrey/asymptotic and Stokes-sector
dependent; only the starting-point recursion for the simple double-pole case is implemented, and
`step/R` limits the step only when it is smaller than the generated scale. General Katz /
Levelt–Turrittin formal gauges, ramification and Stokes matching are not implemented and fail
closed. Duplicate or defective formal blocks, and exact Frobenius over general algebraic number
fields, likewise remain unimplemented.

`five_term_tail_diagnostic(terms, N, threshold=...)` is the public scalar check on a caller's own
term sequence: it stores numerator, denominator, ratio, threshold and the strict-less-than gate
for `abs(sum(T[N+1:N+6])) / sum(abs(T[N-4:N+1]))`. It does not replace the vector-block ratio used
inside the matrix formal branches, and it fails closed when the denominator is zero or the
interval cannot prove the bound.

## Examples

`examples/` keeps seven runnable cases; each one resolves the package through a path variable at
the top of the file.

| Example | Shows |
| --- | --- |
| `qnm_2x2.py` | An exact 2x2 first-order system in a unified `u`: two runs start from `{a, b, C}` boundaries at the horizon and at infinity, transport to a common matching point, and check the components that must vanish on the other side. |
| `regular_singular_save.py` | Transport from an `{a, b, C}` start at a scalar regular singular point to an ordinary point, saving both boundaries. |
| `exponential_boundary_save.py` | A strictly decoupled exponential singularity: transport from a `{phi, a, b, C}` start, save that boundary, read it back in. |
| `ep_series_reconstruction.py` | Adaptive reconstruction of the `1/ep` and finite parts from a `leading_power = -1` certificate, with a redundant candidate pool, independent validation points, order growth on failure, and an open complex-angle range variant. |
| `ep_parallel.py` | Bounded multiprocess transport over different fixed `ep` values. |
| `mathematica_interface_example.wl` | The Wolfram route: exact rational system, explicit path plan, execution, comparison against the closed form. |
| `ep_parallel_mathematica.wl` | The same fixed-`ep` parallel contract through `FlintNDEEvaluateEpBatch`. |

```powershell
python qnm_2x2.py
wolframscript -file examples/mathematica_interface_example.wl
```

`examples/README.md` records the details of each run.

## Documentation

`Documentation/FlintNDE.pdf` is the package paper: the numerical method, the exact input
protocol, the local-basis scheduling, the public interface and the limits above, with the
worked numerical checks. `FlintNDE.tex` and `references.bib` next to it are the sources.
