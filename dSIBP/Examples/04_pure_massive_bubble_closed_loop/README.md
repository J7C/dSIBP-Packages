# Pure massive bubble closed loop

Run with `wolframscript -file main.wl`.

The case fixes the `--` contour branch, even parity, and equal external energies. Package and
reference energies satisfy `P_pkg = P0 = -P_ref`. A generic symbolic `analyticRegulator` is
attached to the zero point of both vertices and only set to `analyticRegulator -> 0` after the
symbolic seeds are generated, so the time-endpoint gate can distinguish a regulated family
from an undefined divergent IBP. The remaining flow runs from initialization to formal Kira
input, external reduction, result readback, the 19-dimensional differential equation, and the
Eq. (51)/(64) scaling checks. Vertex-exchange symmetry is enabled only under the equal-energy
constraint of this case; a family with independent `P1/P2` must not reuse it.

`reference_user_mi_basis.wl` stores only the 19 reference candidate linear combinations (all
active, with no auxiliary relations), `activeIndices = Range[19]`, the physical `ks = ss11`,
the derivative variables, and the 19 uniform scaling degrees. `main.wl` hands this data
straight to the package entry `DSUserMI`; the
linear rank, the invertible `J/userMI` map, backend ids, derivative closure and the manifest
are all produced by the package, so the example implements no basis adapter of its own.
`dlog_basis.wl` keeps only the reference-readable notation and is not a formal basis on its
own.

Kira is not run inside this directory. Input generation happens in a Wolfram session outside
the repository and the reduction itself runs in an external Kira installation:

1. Set `DSIBP_KIRA_WORKSPACE` to a directory outside the repository, then run
   `wolframscript -file main.wl`. The script writes `init/` and `kira/` under its
   `pure_massive_bubble/` subdirectory and does not start Kira; the first clean state is
   `awaitingExternalKira`.
2. From a shell that can reach the same external workspace, enter
   `pure_massive_bubble/kira/` and run `kira --parallel=10 jobs.yaml`. Full results stay in
   that workspace and are not copied back into the example directory.
3. Run `main.wl` again, or run only `post_kira_check.wl`. The latter regenerates nothing and
   executes `DSKiraImport -> DSDE[{ss11, P0}] -> DSScaleCheck[{1, 1}]`, writing into the same
   external workspace.

The final DE and its ordered master list are written to `results/dlogDE/` of the external
workspace, and the closed-loop summary to `results/post_kira_summary.wl`. If the family,
active basis, package or Kira input changes, the reduction must be regenerated and rerun; the
scripts refuse to read results older than the current export manifest. The two summaries kept
here record the input scale and the final outcome of this case, and are not importer input.
