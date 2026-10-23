# FlintNDE Changelog

## 1.01

Added a standalone C++/FLINT engine under `engine/` (built as `flintnde_engine` with `make` in an MSYS2 UCRT64 shell) that performs the prefix certified dual-chain transport under the same file-protocol contract: both chains share one precision and differ only in order, so the first primaryOrder+1 coefficients of the reference series are exactly the primary series. Engine outputs are bit-identical to the Python side; the public interface of the Python package `flintnde` is unchanged.

Added the residue kernel option `--residues auto|sparse|dense`, with the default `auto` selecting the sparse skip-multiply kernel whenever the nonzero density of the residue matrices is at most 0.25; all three kernels share one precision and arithmetic sequence and produce bit-identical outputs.

Batch jobs now run on a single-process `std::thread` pool with per-task streaming parsing (`--workers` falls back to the request's `parallelTaskCount` by default); the parsing wall time of large request JSON files drops sharply (median 56.8 s down to 11.8 s on the 100-task benchmark).

Added the Wolfram entry `FlintNDEEngineEvaluateEpBatch`, which hands fixed-regulator jobs directly to the engine without the Python bridge; on success `finalValues` and the per-point values are decoded into Complex.

Dependency changes: the Python dependencies are unchanged, still `python-flint>=0.6` and `sympy>=1.12`; the engine is an optional component whose build adds a dependency on the MSYS2 UCRT64 toolchain and the FLINT 3.x development library, and whose Windows runtime needs the `libflint` runtime reachable through PATH or explicitly given by `"EngineLibraryPath"`.

From this version on FlintNDE evolves independently and is no longer synced file by file with the MadStree Vendor; MadStree keeps using the FlintNDE v1.0 frozen inside its Vendor.
