# Examples

Every example here resolves the package from `$InputFileName`, so the whole file runs directly
from this tree:

```powershell
wolframscript -file 01_massless_full_edge.wl
```

Each one is also written to be executed section by section in the Mathematica front end; the
chapter cells are ordered so that a section only reads what the sections above it produced.

| Example | Covers |
| --- | --- |
| `01_massless_full_edge.wl` | A single theta-carrying massless full propagator: topology initialization, masters, recurrence and the `dlog` DE, then the single-point and multipoint numerical workflow on one table. |
| `02_single_vertex_family.wl` | A single-vertex function family defined both compactly and explicitly. The two initialization schemas build the same context; `ParameterRules` is given once, and the single-point and multipoint runs differ only in the number of value rows. |
| `03_time_only_cycle_chart.wl` | Time-only cycle initialization, the common-theta contact sector, the `dlog` DE, and the strict time-rank blow-up chart certificates. |
| `04_three_vertex_tree.wl` | The minimal massless three-vertex, two-propagator tree (`+++` contour vertices) through to a batch multipoint evaluation sharing one finite anchor, exported to CSV and JSON. |
| `05_massive_three_vertex_tree.wl` | A three-vertex tree with two massive propagators (`nu12 = 3/4`, `nu23 = 1/3`, so no half-integer representation degeneracy). The original point list and the list with middle external-leg energy `k2 = 0` share one context but are planned, evaluated and exported separately. |
| `06_massless_three_vertex_ep_regularization.wl` | Common time-power analytic regularization `a1 = a2 = a3 = 1 + ep`. The leading integer power is certified from the symbolic boundary and the `dlog` DE before any numerical transport, exact-`ep` points are chosen automatically, and the Laurent coefficients are reconstructed through `ep^0`. The default run uses the open angle range `{-Pi/3, Pi/3}`; pass a positive integer as the last command-line argument to override the requested `ParallelTaskCount` for a batch benchmark. |
| `07_zero_external_leg_energy.wl` | Vertices with no physical external leg: omitting `externalLegEnergy` and passing `0` build the same analytic `dlog`, the numerical interface never exposes the private auxiliary energy, and transport continues to the physical zero point. |

## Output

Run products are written next to this directory, not into the package sources:

- `results_temp/` — transient bridge input, logs and the successful result cache.
- `results/` — durable formula artifacts and exported evaluation data.

Both follow the directory contract described in `../README.md`, and both are regenerated on
demand; nothing in them is required as input by any example.
