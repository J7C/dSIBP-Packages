# Vendored FlintNDE

This directory is a synchronized copy of the FlintNDE Python package, shipped so that MadStree
can transport dlog systems without asking the user to install anything separately. MadStree
resolves it through a single relative path:

```wl
MSFlintNDERelativePath  (* defaults to "Vendor/FlintNDE", relative to the MadStree root *)
```

Python import name is `flintnde`; the Wolfram wrapper in `Mathematica/FlintNDE.wl` is loaded
with `Needs["FlintNDE`"]` after adding this directory to `$Path`.

Edit the standalone FlintNDE package and resynchronize this copy. Do not develop features here:
the copy exists only so that a released MadStree keeps a fixed numerical backend, and MadStree's
result cache is keyed by the content identity of `flintnde/*.py` and `pyproject.toml`.

See the top-level `FlintNDE/` directory of this repository for the full interface
documentation and examples.
