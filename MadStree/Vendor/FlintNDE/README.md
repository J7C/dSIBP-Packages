# Vendored FlintNDE

This directory is a bundled copy of the FlintNDE Python package, fixed at the backend version
MadStree was released with. It does not follow later standalone FlintNDE releases. MadStree
resolves it through a single relative path:

```wl
MSFlintNDERelativePath  (* defaults to "Vendor/FlintNDE", relative to the MadStree root *)
```

Python import name is `flintnde`; the Wolfram wrapper in `Mathematica/FlintNDE.wl` is loaded
with `Needs["FlintNDE`"]` after adding this directory to `$Path`.

Do not develop features here: the copy exists only so that a released MadStree keeps a fixed
numerical backend, and MadStree's result cache is keyed by the content identity of
`flintnde/*.py` and `pyproject.toml`.

See the top-level `FlintNDE/` directory of this repository for the documentation of the
standalone package, which evolves independently of this copy.
