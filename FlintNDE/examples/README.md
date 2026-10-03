# Examples

Every example here loads the package through the `PACKAGE_ROOT` (Python) or
`versionDirectory` (Wolfram) variable at the top of the script, which points at the
`FlintNDE` directory next to `examples/`. Nothing is copied out of the package and no other
project code is read.

For normal use, install the package once (`python -m pip install path/to/FlintNDE`) and then
`import flintnde` from any working directory. The path variable only serves examples that run
directly from this source tree without an installation.

## Typical cases

- `qnm_2x2.py` — an exact 2x2 first-order system in a unified `u`. Two runs start from the
  `{a, b, C}` boundaries at the literal horizon `1` and at infinity `"inf"`, transport to a
  common matching point, and check the components that must vanish on the other side. A
  slightly shifted frequency runs as a counterfactual.
- `regular_singular_save.py` — a scalar regular singular point. Transports from an
  `{a, b, C}` start to an ordinary point and saves both the reusable singular boundary and
  the ordinary endpoint.
- `exponential_boundary_save.py` — a strictly decoupled exponential singularity. Transports
  from a `{phi, a, b, C}` start, saves that boundary, and reads the same boundary back in.
- `ep_series_reconstruction.py` — reconstructs `1/ep` and the finite part adaptively from the
  `leading_power = -1` certificate supplied by the upstream symbolic stage. It consumes only 4
  of 7 explicit production candidates in the first round, adds 2 more after a failed
  validation while reusing every cached value, and leaves 1 candidate unused; outer
  parallelism defaults to 12 processes. Explicit validation points never enter the fit, and no
  leading power is ever guessed from numerical samples. The same script shows an open angular
  range `(-pi/3, pi/3)`, where at most three interior complex rays are chosen and the
  magnitudes remain determined by the target accuracy.
- `ep_parallel.py` — bounded multiprocess transport over different fixed `ep` values.
  `parallel_task_count` defaults to 12, and the actual concurrency is
  `min(parallel_task_count, number of tasks)`.
- `mathematica_interface_example.wl` — the Wolfram route: builds a system from an exact
  rational matrix, plans an ordinary-point path explicitly, executes it, compares against the
  closed form, and writes its results under `results/` in this directory.
- `ep_parallel_mathematica.wl` — the same fixed-`ep` parallel contract through
  `FlintNDEEvaluateEpBatch[..., ParallelTaskCount -> 12]`.

## Running

```powershell
python qnm_2x2.py
# or select a convention explicitly
python qnm_2x2.py --config ../config/qnm_u_unified_it0_3_it1_minus1.json
python regular_singular_save.py
python exponential_boundary_save.py
python ep_parallel.py
python ep_series_reconstruction.py
wolframscript -file ep_parallel_mathematica.wl
```

Both Wolfram examples set `"WorkDirectory"` to the single `results_temp/` inside `examples/`;
the package only appends a short `bridge/` segment and a 16-character token file name. If you
pass a directory yourself, that path is likewise taken as the runtime root itself.

Durable output goes through the public layout API, so it follows the calling script rather than
the current working directory: the `qnm_2x2.py` summary lands in
`results/qnm_2x2/summary/qnm_2x2_summary.json` and its layout description in
`results/qnm_2x2/configuration/output_layout.json`. `it0/it1` are read from the actual config
file and written back into the summary. Infinity automatically uses the start-only
`formal_exponential_asymptotic` route; the nearest exponential-root separation and the default
minimum-term rule generate the first matching point, and the result always keeps the requested
`N` orders together with the five-term block ratio, the five-order relative refinement and the
auxiliary minimum-term diagnostics. A block ratio of at least 1 raises a warning but still
saves the result. Automatic points and nearby points are cross-checked against an independent
second-order scalar recurrence.
