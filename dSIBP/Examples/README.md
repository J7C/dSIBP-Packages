# Examples

This directory contains runnable Wolfram Language sources for the public dSIBP interface.
Only source files and the lightweight input/result records of the two Kira cases live here.
External reduction is never executed inside these directories, so no `init/`, `kira/`,
database, save, log, reduction table, cache or DE run directory is part of the package.

| Example | Coverage |
| --- | --- |
| [`01_mixed_bubble_workflow.wl`](01_mixed_bubble_workflow.wl) | Mixed massive/massless bubble through the unified public entry points: discrete templates, continuous sample points and `linearData`. Builds the reduction plan in memory only; no Kira export or run. |
| [`02_function_system_hankel.wl`](02_function_system_hankel.wl) | Hankel/function-system input and the associated state transforms. |
| [`03_single_massive_sunrise/`](03_single_massive_sunrise/README.md) | Two-loop single-massive sunrise: general seeds and the `{ss11, kE}` parameter-derivative operators. |
| [`04_pure_massive_bubble_closed_loop/`](04_pure_massive_bubble_closed_loop/README.md) | Pure massive bubble end to end: formal Kira input, external reduction, readback, the 19-dimensional DE and the scaling check. |
| [`05_tree_two_vertex_time_ibp/`](05_tree_two_vertex_time_ibp/main.wl) | Two-vertex tree time IBP, the naive DE and the formula route. |
| [`06_mix_bubble_tree/`](06_mix_bubble_tree/README.md) | Mixed cycle/bridge topology: momentum roles, massless convention, contact contraction, and the 81-master Kira/DE record. |

Every example loads the package through
[`load_current_package.wl`](load_current_package.wl), which resolves the source tree that
sits next to this directory. Run a single-file example with:

```powershell
wolframscript -file 01_mixed_bubble_workflow.wl
```

The two Kira cases generate their input in a workspace outside this repository
(`DSIBP_KIRA_WORKSPACE`) and read the reduction results back from the same workspace; their
input-scale and final-result records are the `kira_input_summary.wl` and
`kira_result_summary.wl` files inside each case directory. Those records describe the case;
they are not importer input.
