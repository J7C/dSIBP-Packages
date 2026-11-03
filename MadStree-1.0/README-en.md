# MadStree

Chinese version: [README-cn.md](README-cn.md)。

Changelog: per-version incremental updates live in `更新日志.md` at the package root
(development repository); the release carries no changelog.

`MadStree` is a Wolfram Language package that works directly with the time-integral tensor formula. The input can be an ordered dS tree topology or a pure time-only incidence graph that does not integrate loop momenta; the output includes contact-reachable sectors, ordered master integrals, stepwise or full iterative reduction, a block-triangular dlog connection and automatic boundary certificates. The main algorithm does not generate general IBP systems and does not call Kira.

The default conventions are

```text
time power: (-tau)^A
nu: the Hankel order, used verbatim, never an absolute value
h(nu,0;z): z^(sNu nu) H_nu(z)
formulaNu: sNu nu, with sNu=+1 by default
NuConvention: "Positive"
```

Only the function name `h` is used, with `h(nu,0;z)=z^(sNu nu) H_nu(z)`: the input `nu` is the Hankel order verbatim, never an absolute value. The baseline formulas are written in the `+nu` (sNu=+1) shape, and every position where the order appears reads `formulaNu := sNu nu`. `NuConvention -> "Negative"` selects the negative-prefactor convention of 2401.00129, where `formulaNu=-nu`; the two formula sets differ only by the sign of `formulaNu`, while the Hankel order itself never changes.

`NuConvention` is chosen only at `MSInitTree` initialization. Once a context is built, the formula matrices, recurrence and dlog read this value fixed; there is no per-call convention override.

## C++ engine lane and precision-adaptive serialization

This version (release re-numbered `1.0`) differs from `v0.20` in two areas: connecting the C++ certified dual-chain engine inside the bundled backend, and deriving all vendor serialization widths and tolerances from the requested working digits. Public entry points, master definitions, normalizations and numerical-result definitions are unchanged from `v0.20`.

- `Vendor/FlintNDE-dS/engine/` is a one-time frozen copy of the C++ engine from package-FlintNDE `v1.02` (certified dual-chain transport, `auto`/`sparse`/`dense` residue kernels, streaming batch parsing), extended to this system's request contract: exact C++ tensor reconstruction from structured letters (`tensorSectors` plus per-letter `{alpha,beta,contributions}`, no dense residue matrices shipped), `certificationMode "embedded"` (a single reference chain with prefix snapshots when a segment has exactly one point, otherwise forced back to certified), and per-segment reporting of the effective certification mode.
- Adapter engine lane: a request runs entirely in the engine when it satisfies `pathPlanning=False`, `SingularityMode "Avoid"`, a finite boundary, and no segment anchor or user point on a letter pole; a missing engine or the environment variable `MADSTREE_FLINTNDE_ENGINE=python` falls back to the Python pool, runtime errors are surfaced fail-closed; `FLINTNDE_ENGINE_LIBPATH` locates the libflint runtime directory; each segment reports the lane actually taken in `planReport.engine`. Requests with a regular-singular boundary or automatic planning enabled always use the Python lane by contract.
- Cross-validation: on real tensor requests the engine's structured letters agree digit-for-digit with the dense form and the engine lane agrees with the Python lane at tolerance level; engine-lane suite 11/11, adapter 23/23, bundled backend 195/195, check_core 54/54, version assets 23/23.
- Three-vertex zd n100 single-round retest: production 128.71 s (0.62x versus the `v0.20` pre-C7 207.28 s, the gain coming from the C7 boundary letter fast path taking effect at n100 for the first time); this production geometry does not use the engine lane by contract; registered as Step 7 in `bench/eff_bench.md`.

## Coarse-block transport and known-roots boundary certificate

The coarse-block DE tensor transport has been measured on the two- and five-vertex cases. Identity factors are extracted exactly from sector/source rectangles in the actual master basis, and the complete solution is packed as active rows by (spectator indices times solution columns) for FLINT matrix multiplication. No solution factorization or removal of subsector couplings is required. Contributions with the same pole and output layout are merged; each recurrence order shares input packing, single-contribution rows are referenced directly and multi-contribution rows use one small C-level matrix sum, without per-2x2-slot scheduling.

Ordinary affine segments use compiled blocks directly. Ordinary boundary segments use the same recurrence when the exact rational system is certified as a polynomial plus simple Q(i) poles; other systems retain the general analytic route. Frobenius/known-root initialization keeps its existing algorithm, and explicit-matrix interfaces and singular-point operations materialize residues on demand.

Measured on the ordinary affine segment with a finite initial vector (same request, median of 3 rounds, against the v0.19 dense backend): two vertices 0.02326 s vs 0.03275 s (0.710x); five vertices 16.154 s vs 64.773 s (0.2494x, 4.01x faster), primary chain 6.797 s vs 26.691 s and reference chain 9.243 s vs 37.770 s. All 6/162 output component midpoints agree, with identical node counts and refinement differences. For the three-vertex case at the zd target points (N=1/10), isolating the target segment gives primary 2.315/2.628 s vs 2.371/2.726 s and reference 2.754/3.149 s vs 3.092/3.534 s; per-branch primary new/old is 0.90--1.03 and all components agree (maximum relative difference 0). At dimensions 16--25 with a non-tiled layout, block packing and assembly offset the savings. Full chain including the exact boundary: two vertices 0.6408 s vs 0.6441 s; five vertices 931.9 s vs 968.7 s (0.962x) with all 162 components agreeing. About 84% of the five-vertex wall time is the shared boundary: all 149 poles are simple, but a non-Q(i) exact root (finite_001) makes simple-pole certification fail, so the boundary keeps the general Cauchy route; the block gain therefore appears only in the ordinary segment. The adapter suite passes 17/17 and the bundled backend passes 191 Python tests. Benchmark scripts are `test/bench_03_block_transport.wls/.py` and `test/bench_04_three_vertex_zd.*` in this version.

Relative to `v0.19`, `v0.20` changed only this directory and its bundled `Vendor/FlintNDE-dS`, carrying two MadStree optimizations with measured gains: per-letter substitution, differentiation and residue extraction in the boundary curve and low-level user path, and direct back substitution of nonresonant leading vectors in reverse sector-DAG order, with the full constrained solve retained for zero denominators or nonzero residuals.

The ordinary affine segment and the ordinary boundary segment were also prototyped with Python block tensor-residue operators. Same-request measurement showed a slowdown from 10.857 s to 16.521 s (52.2%) in total wall time and from 4.304 s to 6.799 s (58.0%) in the primary chain, so they have been reverted to the dense FLINT route of `v0.19`. The old implementation multiplies dense matrices in C, while per-contribution Python scheduling cost exceeds the saved scalar arithmetic. The associated `residue_actions`, `TensorResidueOperator`, `TensorBoundarySystem` and dedicated tests have been removed from this version.

The retained boundary optimization reduces the three-run median from 6.405 s to 4.511 s, a 29.6% speedup. Of the 81 five-vertex boundary branches, 17 use direct substitution and 64 use the constrained fallback, all agreeing with the original result; all 81 final components have identical serialized midpoints at 60 and 200 digits. The handbook keeps the per-letter boundary-residue and sector back-substitution formulas, but no longer presents the reverted ordinary tensor recurrence as the current implementation.

The singular boundary of `MSBoundaryData` now carries a power-log certificate: the exact value, algebraic and geometric multiplicity, Jordan defect and positive-integer differences of every indicial root, all derived from the existing leading-vector data plus exact nullspace/rank analysis where needed, without solving a characteristic polynomial. The certificate travels with the evaluation request, and the bundled backend builds a verification-only manifest from it: ascending kernel chains check that each claimed root has generalized-eigenspace dimension equal to its claimed multiplicity and that the multiplicities sum to the system dimension (together these pin down the characteristic polynomial uniquely), so `charpoly()` and root finding are no longer called. When logs or mixed-defective Jordan structure appear, the power-log recurrence continues inside the known spectrum (spectral projectors plus resonant block systems; the first log order emerges naturally from incompatibility-driven budget raises) instead of falling back to a spectrum-recomputing entry point. Measured on the same five-vertex request: backend manifest construction drops from 3.344 s to 1.464 s (−56%), power-log basis construction is unchanged between routes, and Wolfram-side boundary construction pays +0.103 s (+2.8%) for the certificate.

## Tensor-product engine

This version keeps reduction, dlog assembly and the NDE input in tensor-product form over local 2x2 matrices, and never multiplies out 2^n matrices. Formulas, public interface, boundary algorithms and numerical results are the same as `v0.18`; only the representation used for computation and transmission changes.

### What is stored

Each sector stores three kinds of data and no 2^n matrix:

- the local base tag of every slot: `massiveT` (sigma2 diagonalization), `masslessH` (sigma1 diagonalization), `identity` (contracted massless quotient slot);
- the quotient maps: the 0/1 table `canonicalOfRaw` from raw states to master columns, and the representative raw row `representativeRaw` of each master;
- the M1 diagonal vector of each component and the per-slot local matrices of M0.

Hence `U = P (⊗ 2x2) S` and `U^{-1} = P (⊗ inverse 2x2) S`, where P and S are sparse projections and embeddings that only select rows and sum over preimages. All three local bases are integer Gaussian matrices times one common factor 1/sqrt(2), which cancels in `U^{-1} D U`; residues stay exact in Q(i) with no algebraic-number approximation.

### The three operations

- Reduction: one recurrence step is pure tensor action. Raising toward zero is `M1^{-1} M0`; lowering toward zero is `U^{-1} diag(1/kappa) U M1`. M1 is exactly diagonal on masters and is read from the representative raw row; M0 acts slot by slot, costing O(n 2^n) and never a 2^n matrix-vector product.
- dlog: each letter carries contribution records instead of a residue matrix. A record is either `diagonal` (a state-diagonal vector, used by Omega_ex and the normalization gauge) or `conjugated` (a `preVector`, a 0/1 `letterIndicator`, a sparse middle matrix and a `postFactor`, used by `U^{-1} diag(Log lambda) U` and by contact blocks).
- Public queries: `MSDLogDE` and `MSFormulaMatrices` still return dense matrices of the same shape as `v0.18`, but materialize them from the structured data on demand; the numerical transport path never touches those dense matrices.
- Construction certificates: checking `U M0_c U^{-1} = diag(I kappa_c)` per component costs one `U^{-1} e_c` and one M0 action per column, after which the U rows come from the transpose accessor and the whole matrix is tested once; sending `U^{-1} e_c` back through all U columns would grow the leaf count like 4^n.
- Letter-list simplify-then-substitute: substituted contribution entries are already canonical Q(i) literals (integers, rationals, rational complex numbers), so real and imaginary parts are stringified directly instead of running `Simplify[Together[...]]` per entry; only other values go through RootReduce normalization. The Q(i) gate predicate is also replaced by an evaluable `MatchQ[_, _Integer | _Rational]` (the previous `RationalQ` spelling is undefined and the gate never fired), so non-Q(i) values now genuinely fail closed with `$Failed`.

### Wire format v2

Requests to FlintNDE use the schema `madstree_flintnde_evaluate_v2` and carry a new `tensorSectors` field; each letter is written as `{alpha, beta, contributions}` with the structured records above. The backend `Vendor/FlintNDE-dS/flintnde/tensor_dlog.py` rebuilds U and U^{-1} by applying local 2x2 matrices along each slot axis and then forms `U^{-1} diag(indicator) U` at cost O(n 2^n) per sector, so no 2^n residue matrix ever enters the request.

### Where the optimization lives

| File | Change |
| --- | --- |
| `Kernel/Formula/TensorProduct.wl` (new) | local bases and inverses, slot metadata, tensor actions of U/U^{-1}/U^T/M0/M1, sector tensor records and certificates |
| `Kernel/Formula/Recurrence.wl` | one recurrence step uses the tensor actions above |
| `Kernel/DE/DLog.wl` | Omega and dlog blocks become structured contribution lists; emits `tensorSectors`; dense matrices only at public queries |
| `Kernel/Formula/TensorAtoms.wl` | `MSFormulaMatrices` materializes from the certified local diagonalizers |
| `Kernel/Numerics/FlintNDE.wl` | letters are serialized as `{alpha, beta, contributions}`; letter-list simplify-then-substitute and a Q(i) gate that actually evaluates |
| `Kernel/Numerics/PathEvaluation.wl` | request schema becomes v2 and carries `tensorSectors` |
| `Backend/flintnde_transport.py` | v2 request validation and tensor residue entry points |
| `Vendor/FlintNDE-dS/flintnde/tensor_dlog.py` (new) | slot-wise U/U^{-1} reconstruction and exact Q(i) residue accumulation |

### Agreement and efficiency

On the paper's 5-vertex massless chain (81 masters, 41 letters) the two builds have the same 5313 nonzero residues with maximum exact difference 0, and final values agree to a worst relative difference of 1.242e-210. Per-letter serialization drops from 50.40 s to 0.05 s (the median over 100 `ep` values, range 0.047--0.080 s); the letter payload drops from 12.5 MB to 221--227 KB; the complete v2 request is 298.8 KB against 16.84 MB for `v0.18`. Backend transport wall time goes from 11.23 s to 10.76 s, and the primary chain from 4.377 s to 4.285 s. Initialization and DE assembly stay in the same range (1.09 s / 1.86 s versus 1.18 s / 1.99 s); both stages spend their time in local bases and certificates, so the difference on this graph is measurement noise and the gain is concentrated in letter assembly and transmission.

## User points at finite singularities

`MSEvaluatePath` sends each maximal continuous complex-affine segment to FlintNDE. Candidate
singularities are determined only from the exact rational matrix or letters obtained by pulling back
the complete analytic dlog DE on that segment. Variable names, a zero coordinate value, graph names,
and vertex counts do not enter the test. When a user point really is singular, the local power--log
solution classifies each component as removable, a true pole, or logarithmically divergent. Removable
components return Acb values; divergent components return the text value `Infinity`. The top-level
`singularityClassifications` table records the evidence. No production branch is specialized to
`k0`, `k2=0`, or a particular topology.

## Loading

```wl
packageRoot = "/absolute/path/to/MadStree";
AppendTo[$Path, packageRoot];
Needs["MadStree`"];
```

The first successful load prints one citation reminder. arXiv identifiers are clickable
in notebooks, while a headless kernel prints full URLs. MadStree lists only
[arXiv:2401.00129](https://arxiv.org/abs/2401.00129),
[arXiv:2411.03088](https://arxiv.org/abs/2411.03088), and Jiaqi Chen, Yi-Xiao Tao, Xinyi Wang
and Dong Zhang, [*MadStree: automated evaluation of tree-level de Sitter correlators*](https://arxiv.org/abs/2610.28526),
arXiv:2610.28526.

## Minimal workflow

```wl
spec = <|
  "vertices" -> {
    <|"id" -> 1, "externalLegEnergy" -> k1, "timePower" -> a1, "vertexType" -> "+"|>,
    <|"id" -> 2, "externalLegEnergy" -> k2, "timePower" -> a2, "vertexType" -> "+"|>
  },
  "lines" -> {
    <|"type" -> "massless", "endpoints" -> {1, 2},
      "momentum" -> q, "nu" -> 1/2|>
  }
|>;

context = MSInitTree[spec];
masters = MSMasterIntegrals[context];
masterDefinitions = KeyTake[#, {"integral", "normalization", "bareIntegral", "definition"}] & /@ masters;
topKey = First[context["sectorOrder"]];
matrices = MSFormulaMatrices[context, topKey];
de = MSDLogDE[context];

formulaData = MSFormulaData[context];
written = MSWriteFormulaArtifacts[
  context,
  TimePowerRules -> Automatic
];
written["outputDirectory"]
written["files", "dlogDE"]

targetRules = {k1 -> 9 I, k2 -> 3 I, q -> 1, a1 -> 1, a2 -> 1};
chart = MSBoundaryChartCertificate[context, targetRules];
```

`MSIntegral[sectorKey,timeShifts,stateBits]` is the normalized master
$J_s(\mathbf n;\mathbf a)$; the same three entries uniquely identify the bare integral
$I_s(\mathbf n;\mathbf a)$. `MSMasterIntegrals` directly reports the exact `normalization`
and $J_s=\text{normalization}\,I_s$ for every master. Use
`MSIntegralDefinition[MSIntegral[sectorKey,timeShifts,stateBits],context]` for one shifted integral.
For a pinched massive (anti-)time-ordered propagator, this normalization contains its one contour factor: after
suppressing the common momentum power, a `++` subsector is `-4 I/Pi` and a `--` subsector is `+4 I/Pi`.
This sign is not multiplied again into the dlog contact block.
A two-state factor is either the two-dimensional state of one vertex of a massive propagator, or the quotient state shared by a whole massless (anti-)time-ordered propagator; `MSSlotRegistry` gives their bit order. The $++/--$ h combinations are defined once in the manual and are not expanded in every record. `MSSectors[context]` returns the deterministic order of contact-reachable sectors with their keys, pinched propagators and master counts, consistent with `context["sectorOrder"]`; `MSContextQ[obj]` tests whether an object is a valid MadStree context and is the shared input gate of every public entry.

`"vertexType" -> "+"|"-"` is the vertex's Schwinger--Keldysh contour branch and is
required explicitly on every vertex. A propagator no longer accepts `skType`, `sigma`, or
vertex-sign input. Its T/W/External class, contour label, overall sign, the signs of
its two vertices, and default Hankel branches are derived from the vertices listed in `"endpoints"`.
The public `"type"` of a propagator accepts only `"massive"` or `"massless"`. Propagators have no
public ID field; MadStree always assigns `1,2,...` internally in the exact `"lines"` input order.
Only an explicit `thetaBundles` declaration refers to propagators through these position numbers.

Here T denotes the same-contour components: the time-ordered (Feynman) propagator
`G++` and the anti-time-ordered (anti-Feynman) propagator `G--`, collectively called
(anti-)time-ordered propagators. W denotes the Wightman propagators `G+-` and `G-+`.
The derived internal types are `massiveT`, `masslessT`, `massiveW` and `masslessW`;
they are not public input types.

`"externalLegEnergy"` does not mean that a vertex itself carries an energy. It is the
parameter of the theta-free external-leg exponential attached to the vertex.
`vertexType="+"` gives `Exp[-I externalLegEnergy tau]`, while `"-"` gives
`Exp[+I externalLegEnergy tau]`. The SK contour sign is also derived from `vertexType`,
but is stored separately from the exponent sign internally. For a positive damping
coordinate `K`, use `externalLegEnergy=+I K` on a `+` vertex and `-I K` on a `-` vertex.

The massless graph above can generate a complete multi-sector Frobenius boundary and transport to one or more user points through one entry point:

```wl
result = MSEvaluatePath[
  context,
  {{k1, k2}, {9 I, 3 I}},
  ParameterRules -> {q -> 1, a1 -> 1, a2 -> 1},
  FlintNDEPathPlanning -> True,
  BoundaryScale -> 4,
  RankOrder -> {v1, v2},
  PythonExecutable -> "...",
  WorkingPrecision -> 200,
  TransportOrder -> 80,
  ReferenceTransportOrder -> 104,
  TargetRelativeError -> "1e-20",
  MessageLanguage -> "EN"
];
```

The second argument is always named `pointSequence`. Its first row is the ordered header of
running coordinates and every later row is one equal-width coordinate value row. A single point
has one value row; multiple points only append rows. Fixed, non-running parameters are supplied
once through `ParameterRules`. The header and `ParameterRules` must be disjoint and together cover
every symbol required by the analytic DE and boundary. Missing or still nonnumeric values return a
structured `Failure` listing the relevant symbols before a runtime directory is created or
Python/FlintNDE is started. A transient row is `{{values...},"tmp"}`.

External energies that remain constant on the user scan may be supplied through
`ParameterRules`. Since the automatic infinity boundary moves all external energies,
MadStree adds those energies to the internal path coordinates, transports from the actual
boundary matching point to the first user point, and keeps their supplied values at all user
points. The returned `coordinateSymbols` and point coordinates can therefore include these
additional energies; the returned `parameterRules` no longer contains them. Listing every
external energy explicitly in the header, with repeated constants in each row, also remains
valid. Internal momenta, masses, and time powers are unchanged.

`TransportOrder -> 80` is the production order: returned nodes, dense user points, segment
endpoints, and regulator-fit samples come from the lower-order primary chain.
`ReferenceTransportOrder -> 104` runs only the higher-order reference chain used by
`TargetRelativeError`. Across complex-affine segments, both chains propagate their own endpoint;
the reference value neither replaces user output nor becomes the next primary initial state.

For a Laurent limit in a shared regulator, use

```wl
epSeries = MSReconstructEpSeries[
  context, ep,
  {{k1, k2, k3}, {9 I, 3 I, 5 I}},
  ParameterRules -> {q12 -> 1, q23 -> 2,
    a1 -> 1 + ep, a2 -> 1 + ep, a3 -> 1 + ep},
  MaximumEpPower -> 0,
  EpGoalDigits -> 20,
  ParallelTaskCount -> 12
];
```

`MaximumEpPower` is the highest requested power; its default `0` asks for poles, if present,
and the finite part. By default users do not supply `ep` samples. Before any numerical NDE solve, the actual
symbolic boundary condition and dlog DE certify the lowest integer power. A negative DE Laurent
power, a non-Laurent boundary, or an unproved leading coefficient fails closed. The fit initially
includes two powers above the user maximum; failed independent validation adds two powers per round,
solves only new production points, and reuses separate production/validation caches. Three failed
rounds return the current best coefficients without relaxing the tolerance. Such a result has
`status -> "computed_with_warning"` and `precisionTargetMet -> False`; it is not precision-certified.
`ParallelTaskCount` defaults to 12 and controls outer
process parallelism, not python-flint's in-process `ctx.threads` setting.

To restrict every regulator value to a user-approved range, set `EpSamplePoints -> {...}` and
`EpValidationPoints -> {...}` explicitly. The former is an ordered production candidate pool and
may contain surplus points. With `EpInitialInternalMaximumPower -> q`, the first round consumes only
the prefix required to fit through `ep^q`; failed validation incrementally consumes later candidates
while reusing cached values. Validation points never enter the fit and must be disjoint from the
whole production pool. Exhaustion never generates an out-of-range point and is reported as
`candidate_pool_exhausted` while retaining the current coefficients.

Alternatively, `EpSampleAngleRange -> {thetaMin,thetaMax}` selects an open angular interval in
radians for automatic complex samples. Up to three uniformly spaced interior rays are used; the
magnitudes remain controlled by the precision, pole depth, and fit order, with no user radius cap.
Validation samples keep the corresponding angles and reduce only the magnitudes. MadStree converts
the high-precision generated points to exact Gaussian rationals before the exact dlog pullback.
These options default to `Automatic`, so the default call and backend schema are unchanged.

`WorkingPrecision` defaults to 200 decimal digits for both `MSBoundaryData` and `MSEvaluatePath`.
An explicit positive integer overrides the default; the backend uses
`ceil(WorkingPrecision log2(10))+32` bits.

For multiple points, replace the second argument with `{pointP1, pointP2, pointP3}`. MadStree partitions the input order into maximal consecutive complex-affine one-variable segments `x(s)=x0+s v`, pulls each segment back once, and sends all exact complex parameters to FlintNDE. MadStree does not plan nodes or detours. Bare coordinates are returned; `{coord,"tmp"}` is transient, and every other string tag is rejected.

Any multipoint transport through intermediate nodes is a jump; only a crossing of a singularity through a local-basis connection is a singularity jump. `SingularityMode -> "Automatic"` is the default and follows explicitly supplied singular user points instead of silently detouring around them. FlintNDE removes the exact singularity from the ordinary-node chain, uses one local basis for user points on both sides inside its convergence disk, returns textual `Infinity` for divergent components, and resumes from an outgoing ordinary point. A terminal singular target receives a hidden in-disk match point when required. `"Avoid"` explicitly refuses a crossing; `"SingularityJump"` explicitly selects the local continuation branch, whose multivalued choice is equivalent to a detour class and must be confirmed by the user. A non-collinear turn at a singularity and an intermediate singularity with planning disabled remain fail closed. `MessageLanguage -> "EN"|"CN"` is case-sensitive and defaults to English.

`FlintNDEPathPlanning -> True` lets FlintNDE plan nodes inside each segment. User points covered by one node form an evaluation bucket and use that node's stored vector series: buckets of at least eight points use subproduct/remainder-tree fast multipoint evaluation, while small buckets use the iterative algorithm. `False` preserves every user point as a node and never calls the planner. Different affine segments never share local coefficients.

`de["masters"]` is the unique authority for DE rows and boundary-vector order. `MSEvaluatePath` sends exact `{a,b,C}` branches, all affine segments, and the selected options in one UTF-8 request; FlintNDE initializes the singular boundary and transports every finite segment in one process. Production code never evaluates defining integrals or calls `NIntegrate`; these are independent-validation oracles only.

Generic entries never dispatch on graph ids, vertex counts, master counts or paper numerical points. They fail closed in a structured way only when the dlog is not closed, the rank/chart does not pass the normal-crossing certificate, the late-time exponents do not decay, the pullback system is not an exact regular singularity, or the target/anchor lies on a DE letter.

When both vertices of a propagator lie on the same contour branch, a public `"massless"`
propagator is internally derived as `masslessT` and uniquely contributes one quotient
two-state slot shared by the whole propagator. The product of the four vertex states is used
only to derive the quotient in the manual; it is neither a program input nor a second
production representation. The current interface has no representation selection option.

The `sectorKey` of every sector is a fixed-length string in the root `lines` order: `0` means the propagator is pinched and `1` means it is not; propagators that cannot be pinched keep bit `1`. The top is the all-`1` string. Leading zeros are part of the identity and the key must never be treated as an integer. `context["sectorKeySchema"]` gives `rootLineOrder`, the width and the bit semantics; for example, pinching `{e1,e2,e4}` in `{e1,e2,e3,e4}` gives `"0010"`. The full identity of a master is `MSIntegral[sectorKey,timeShifts,stateBits]`, so different subsectors are never the same integral even if their local shifts/bits coincide. The `sectorIdentityCertificate` of `MSInitTree` checks at once the pinch sets, sector keys, full masters, global indices and an SHA-256 digest over the full master order; collisions reject initialization. The massless top-to-sub dlog contact atoms multiply by `(-1)^N` for the actual number of massless propagators pinched in the event; simultaneous contacts use no graph-specific symbol table.

`MSBoundaryChartCertificate[context,targetRules]` builds the nested blow-up `1/K[sigma[j]]=Product[x[r],{r,j,V}]` before the boundary computation, fixes all thetas sector by sector, and checks that the complete dlog letters, normalization, coordinate Jacobian and the new denominators from shifted contacts are coordinate monomials times a nonzero boundary unit. `RankOrder -> All` checks all root-time strict charts; `MSBoundaryData` returns `BoundaryChartNotCertified` when the certificate fails. Local quantities in certificates and boundary records use two public symbols: `MSBlowupCoordinate[i]` is the `i`-th local coordinate of the nested blow-up, and `MSDampingEnergy[v]` is the positive damping energy `K_v` of vertex `v`.

## Common theta and time-only cycles

(anti-)time-ordered propagators on the same current component pair automatically form a common-theta bundle. The program generates only nonempty odd-subset events with coefficient `2^(1-Length[selected])`; one event deletes the selected propagators and merges the vertices exactly once. Sectors come from an event BFS, not the power set of (anti-)time-ordered propagators. After merging, an unselected massless (anti-)time-ordered propagator keeps only the even state, a massive (anti-)time-ordered propagator uses the `10 -> 01` equal-time canonical form, and the raw tensor matrix is projected onto the true master space with the stored projection/embedding.

Pure time-only cycles use the dedicated entry:

```wl
cycleContext = MSInitTimeGraph[<|
  "vertices" -> {...},
  "lines" -> {...}
|>];
```

It allows only time integrals over fixed propagator momenta and neither reads nor generates loop momenta, ISPs, Landau or threshold data. Active self-loop propagators no longer trigger contacts; multiple paths to the same vertex partition/pinched set are canonicalized to one sector.

## Single-vertex function families

A single-vertex family needs no fake topology. The compact input follows the `ki/nui` convention of the reference code:

```wl
vertexContext = MSInitVertexFamily[<|
  "ki" -> {k0, k1, k2},
  "nui" -> {a0, nu1, nu2},
  "hankelBranches" -> {1, 2}
|>];

vertexDE = MSDLogDE[vertexContext];
singlePointSequence = {{k0}, {3 I}};
vertexValue = MSEvaluatePath[
  vertexContext,
  singlePointSequence,
  ParameterRules -> {k1 -> 1, k2 -> 2, a0 -> 2, nu1 -> 1/2, nu2 -> 3/2}
];
```

Here `First[ki]` is the theta-free external-leg exponent parameter and `First[nui]` is the base power of `(-tau)`; the remaining entries give one by one the momentum and the Hankel order `nu`, used verbatim, of each h block. The explicit model uses
`externalLegEnergy/timePower/hBlocks/exponentialBlocks/vertexType/normalization`; each h block has
`id/momentum/nu/hankelBranch`, and each extra exponential block has
`id/momentum/exponentType`. This initialization schema differs from tree `vertices/lines`
because it defines a function product at one vertex rather than an incidence graph. Once initialized,
both context kinds use the same masters, recurrence, DE, boundary, and `MSEvaluatePath` interface.
A pure exponential block is `1x1` and adds no state bits; every h block is `2x2`. The `NuConvention` chosen at initialization stays fixed and cannot be overridden in recurrence, DE or numerics.

## Direct reduction

```wl
reduction = MSReduce[
  2 MSIntegral[First[vertexContext["sectorOrder"]], {1}, {0, 0}] -
  3 MSIntegral[First[vertexContext["sectorOrder"]], {-1}, {1, 0}],
  vertexContext,
  MasterBasis -> Automatic
];
```

`MSReduce` handles only valid `MSIntegral[sectorKey,timeShifts,stateBits]` objects in a fixed initialized context: `sectorKey` must belong to the contact-reachable sector DAG, `timeShifts` are integer shifts per component of that sector, and `stateBits` traverse all two-dimensional slots of the sector in the binary order of 2401.00129. It also accepts `List`, regular `Table`, empty lists and ragged nested lists; `result`, per-element `coefficientVector/masterRules`, residuals and remaining shifts preserve the input shape and share one recursion memo. A failed element reports its exact position. It does not support new propagator powers, new discrete indices, sectors outside the context, or another function family. Reduction recurses along the shift/contact DAG to the zero-shift masters and caches repeated subproblems. For scalar input, the returned `masterBasis`, `coefficientVector` and `masterRules` are strictly ordered; `result` is an explicit linear combination of master integrals; `nonMasterResidual` and `remainingShiftedIntegrals` must both be empty for a complete reduction. `MasterBasis` may be any permutation of all context masters, but missing, duplicated or extra objects are rejected. `singularLayers` records stepwise the direction, component and denominators: raising toward zero reports only the `M0` energy letters, and lowering toward zero reports only the `M1` eigenvalues. The reduction runs on the vertex time powers fixed at initialization: to substitute numerical or regulated powers (e.g. `a1=-1+ep`), re-initialize the context with these values as `timePower` and reduce there (Example 06 demonstrates this), or substitute into the symbolic result afterwards; `MSReduce` itself carries no substitution option.

The single-step entry is `MSRecurrenceStep[integral, component, context]`: it moves one nonzero shift of the chosen component one step toward zero, takes the first nonzero shift when the component is omitted, and returns `reducedOneStep`, `nearerShift` and the denominators actually inverted in that step; `MSReduce` is its memoized complete entry. The contact injection structure is available separately through `MSContactMaps[context, sectorKey]`: one record per contact event, with the per-component injection matrices (`matricesByComponent`), the target sector and the zero/raised target time shifts (`R0TargetShift`, `R1TargetShiftsByComponent`), the normalization ratio and the list of pinched propagators; omitting the sector key, or passing `All`, returns every sector keyed by `context["sectorOrder"]`.

Initialization uses the known two leading small-time powers of h to store a safe time-IBP shift
lower bound for every reachable component. A zero-shift master whose Gamma argument is already a
non-positive integer is rejected immediately. A positive-integer starting argument is accepted,
but a reduction is rejected only when it crosses the safe lower bound; the failure lists the
original vertices, sector, branch, Gamma argument and requested shift. Adding a symbolic analytic
regulator leaves that parameter symbolic and permits the family to proceed.

If the top-to-sub `R^(1)` lands on a subsector with a nonzero shift, `MSDLogDE` calls the same formula recurrence column by column, builds the reduction matrix from the shifted subsector to the global ordered masters, and right-composes it with the contact block. Each contact block keeps `shiftReductionRecords`, the residual, remaining shifts and singular layers; if any column does not close, the dlog status is `contactShiftReductionFailed` and no falsely closed connection is emitted.

The formula layer never constructs a generic large matrix inverse. The fixed `sigma2` of the two-dimensional state at one vertex of a massive propagator uses the paper's `2x2` transform and its explicit inverse, while the fixed `sigma1` of the quotient slot shared by a whole massless (anti-)time-ordered propagator uses the self-inverse Hadamard; an all-sector matrix is only the Kronecker product of these local matrices. `M1` takes reciprocals of diagonal entries in the state-bit basis and `M0` takes reciprocals letter by letter in the corresponding local tensor basis. If a future slot contains multiple components that cannot be diagonalized by the same local transform, this fast path must fail closed.

The FlintNDE numerical backend ships as a bundled copy inside this version directory, independent
of later standalone releases; its default location is defined in a single relative-path variable:

```wl
MSFlintNDEConfiguration[]
(* relativePath -> Vendor/FlintNDE-dS *)

MSSetFlintNDERelativePath[FileNameJoin[{"new-name", "code", "package"}]]
```

Paths are always relative to the current MadStree version directory; if Vendor is renamed or moved later, only this item changes. `MSNumericalSystem` remains the low-level entry for users who supply their own boundary vectors.

The default runtime root for automatic MMA calls is the single `results_temp/` next to the calling script:

```text
results_temp/
  nde/    # current request input and log
  cache/  # successful JSON cache keyed by request and backend source identity
```

An explicit `MSRuntimeDirectory -> path` denotes the runtime root itself. Absolute paths are used directly and relative paths are resolved against the calling script directory; MadStree does not append another `results_temp`. On Windows, complete cache/input/log paths are checked before directory creation or Python launch. An overlong path returns `RuntimePathTooLong` with its actual length. `RuntimeInputWriteFailed`, `PythonFlintUnavailable`, `FlintNDELaunchFailed`, `FlintNDEOutputMissing`, and `FlintNDEOutputInvalid` retain distinct failure boundaries. A successful-result cache key includes both the request and SHA-256 identities of `Backend/flintnde_transport.py`, the Vendor `pyproject.toml`, and sorted `flintnde/*.py`; an in-place source fix therefore cannot reuse a plan with a different source identity. Apart from caller-side caches, the package source directory never receives run artifacts.
Both the MadStree adapter and the Vendor Wolfram bridge launch Python with argument-list `RunProcess`; shell `Run`, quoting helpers, redirection, and retry fallbacks do not exist. Current numerical gates and the Vendor `26/26` interface regression cover repeated Windows DLL initialization.

The analytic `de = MSDLogDE[context]` is itself a durable user result. A symbolic-only workflow may
call `MSWriteFormulaArtifacts[context]` explicitly; its two dedicated options are `TimePowerRules`,
which specializes only the written recurrence metadata and never the analytic DE, and
`MSOutputDirectory`, which redirects the artifact directory. `MSEvaluatePath` and `MSReconstructEpSeries`
automatically write or reuse the complete analytic artifacts before any Python/FlintNDE task starts;
a write failure prevents numerical NDE, and the successful result reports all paths under
`"formulaArtifacts"`. The default directory is `results/madstree_formula/run-<UUID>/`, containing
`masters.wl`, `recurrence_metadata.wl`, `dlog_de.wl`, and `manifest.wl`. The dlog file always contains
the complete analytic DE; `ParameterRules` only creates a function-local evaluation copy and never
overwrites it. Repeated evaluations of the same context in one Wolfram session reuse the files while
they remain complete; deleted files are regenerated. The manifest records version, case, convention,
sector/master digests, and paths.

Ordinary user points are returned under the `"saved"` records of `MSEvaluatePath`; `"tmp"` points participate only in transport. Export saved points with

```wl
MSExportEvaluationData[result, MSOutputDirectory -> outputDirectory]
```

to produce CSV/JSON. The default directory is `results/madstree_evaluation/run-<UUID>/`, with
fixed names `evaluation_data.csv` and `evaluation_data.json`. If any requested format fails, the
function returns `EvaluationExportFailed` and never reports the whole export as `"written"`; an
overlong target remains a separate `RuntimePathTooLong`. `results/` contains durable user outputs,
`results_temp/` contains requests, logs, and cache, and `test/results_test/` contains development
test artifacts. Package source directories receive no runtime products.
The full formulas, the massless `4 -> 2` quotient, contact shifts and the top-to-sub dlog derivation are in [Documentation/handbook.pdf](Documentation/handbook.pdf); the one-line table of every public symbol is in its "Public interface reference" section.

## Examples

- [01_massless_time_ordered_edge.wl](Examples/01_massless_time_ordered_edge.wl): massless quotient, master integrals, recurrence, dlog and the automatic boundary/numerical entry (including a user custom finite boundary and batch multi-point evaluation with CSV/JSON export).
- [02_single_vertex_family.wl](Examples/02_single_vertex_family.wl): compact and explicit single-vertex-family definitions, masters, reduction, and a single/multipoint cross-check through the same `pointSequence` interface.
- [03_time_only_cycle_chart.wl](Examples/03_time_only_cycle_chart.wl): time-only cycle, common theta, contact sector and all strict-rank chart certificates.
- [04_three_vertex_tree.wl](Examples/04_three_vertex_tree.wl): three-vertex massless tree with the `+++` vertex structure, followed by a batch multi-point evaluation with CSV/JSON export.
- [05_massive_three_vertex_tree.wl](Examples/05_massive_three_vertex_tree.wl): three-vertex massive
  tree with the `+++` vertex structure and non-half-integer nu, with a normalized/bare definition
  self-check: per master, the integral/bareIntegral/normalization fields of `MSIntegralDefinition`
  must match the `MSMasterIntegrals` records, and the shifted integral must read the normalization
  of its sector. The same context first evaluates and
  exports the original nonzero point sequence, then uses a second independent `MSEvaluatePath` call
  for four points with middle external-leg energy `k2=0`; planning and results are not mixed.
- [06_massless_three_vertex_ep_regularization.wl](Examples/06_massless_three_vertex_ep_regularization.wl): a massless three-vertex chain with `a1=-1+ep, a2=a3=1+ep`, whose `ep=0` time integrals genuinely diverge. Before the numerical stage it runs two symbolic self-checks: `MSReduce` acts once on the sum of two integrals whose powers lie one unit below the master powers (with status/residual/remaining/singularLayers verified); the same reduction is repeated on a context re-initialized with the regulated powers substituted into `timePower`, so the coefficients carry explicit `ep`. The symbolic boundary and DE then certify the leading power `-1` before numerical NDE work, after which the program selects production/validation points and extracts the `ep^-1` pole coefficient together with the finite part.
- [07_zero_external_leg_energy.wl](Examples/07_zero_external_leg_energy.wl): a vertex without a physical external leg may omit `externalLegEnergy` or set it to zero; a private auxiliary energy remains in the analytic dlog and is transported automatically to the physical zero.

All seven examples pass fresh runs from an empty output directory, printing `Example PASSED` and exiting with code `0`. Example 06 passes `17/17` checks with four production points and two independent validation points; the default request of 12 workers is automatically capped at the six actual ep tasks.

## Statement on AI use

The notation, the functionality and the public interfaces of this package were designed by the
authors, who also specified the derivation recorded in
[Documentation/handbook.pdf](Documentation/handbook.pdf) and audited every formula against
it. The implementation was produced under their direction by a number of different agent programs
driven by large language models, and the authors designed the verification tasks that establish the
reliability of what those programs produced. Only the tasks that cross-check computed results are
listed here, together with the source of the reference data in each case.

- **Single-vertex families and the massive two-vertex correlator**: the exact differential equations,
  the boundary coefficients and the numerically transported results are compared with the known
  results of [arXiv:2411.03088](https://arxiv.org/abs/2411.03088), which are computed by an
  independent multivariate hypergeometric construction.
- **The massive three-vertex tree graph**: the reference data were produced by an agent program that
  generalized the previously hand-written two-vertex code to three vertices along an independent
  route. Its output was first verified against the published eight-branch result of
  [arXiv:2309.10849](https://arxiv.org/abs/2309.10849) and only then used to check the package, so
  neither the reference nor the package shares a derivation.
- **IBP relations and differential equations**: these are compared with the `dSIBP` package (in
  preparation, released in the same repository as this one,
  [github.com/J7C/dSIBP-Packages](https://github.com/J7C/dSIBP-Packages)). `dSIBP` samples naive IBP
  seeds and solves the resulting linear reduction system, whereas the present package applies the
  general formulas derived in the handbook directly and never re-solves an IBP system, so the
  agreement tests the closed formulas themselves.

## Current boundaries

- Common-theta odd-subset simultaneous contacts, coincident massless/massive quotients, event-reachable sector BFS and pure time-only cycles are supported; this does not mean loop-momentum integration is supported.
- Production boundaries never dispatch on graph names or master counts. The single-vertex massiveExternal case follows the same general route (the explicit 2411.03088 Sec.3.3 series is kept only as a test reference baseline); closed tree/time-only contexts build the nested curve, complete dlog pullback residue and ancestor-sector leading system from the sector DAG, component/slot metadata, normalization and strict time rank. When the formula/dlog/chart is not closed, the late-time exponents do not decay, or the pullback system is not an exact regular singularity, the package fails closed and never falls back to finite-point defining integrals.
- FlintNDE transport requires the pulled-back connection to be in exact `Q(i)(s)` or `Q(i)(t)`.
  This fresh development run covers single-vertex, pure-massless, massive (anti-)time-ordered, and mixed
  three-vertex systems through the current generic boundary producer.
- FlintNDE retains exact Lee--Moser, high-pole, and singularity-jump capabilities, and adds node-bucket fast multipoint evaluation plus the public strict user-node route. Uncertified high poles and internal/final points requiring ramification, algebraic extensions, or general Stokes connections continue to fail closed.

The current fresh regression ran 25 Wolfram development scripts, all exiting with code `0`:
core `54/54`, formula artifacts `23/23`, point sequences/planning `29/29`, integral
definitions `20/20`, propagator classes `52/52`, the five boundary suites `41/41`, `14/14`,
`10/10`, `6/6` and `6/6`, 900-point complex-plane grouping `8/8`, single-vertex NDE `13/13`,
mixed three-vertex `8/8`, runtime/export `10/10`, auxiliary external-leg energy `6/6`,
time-IBP domain `15/15`, simultaneous chart `22/22`, vertex-family reduction `19/19`, and
Examples 01--07 `7/7`. The Python adapter `17/17`, the bundled FlintNDE-dS's 191 Python tests
and its `26/26` Wolfram interface checks also pass.
