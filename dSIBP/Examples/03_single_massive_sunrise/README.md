# Single-massive sunrise

Run with `wolframscript -file main.wl`.

The case is fixed to two vertices, three parallel propagators and two loop momenta: the
first propagator is a massive `h` block, the other two are massless exponentials, and both
independent scalar products are given explicitly.

Both vertices share one external-leg energy `kE`, while the independent loop external
momentum `kL` is kept separate, so the input carries two scales. The two ISPs are chosen as
a pair that is exchanged when the two massless parallel propagators are swapped, and
`symmetryRules` implements the vertex exchange together with that massless propagator/ISP exchange.

The script runs `DSKinematics`, `DSInit` and `DSSeeds`/`DSAllSeeds`, building the general IBP
templates for every contact-reachable sector while all continuous indices stay symbolic. It
also reads the `{ss11, kE}` general parameter-derivative operators produced by initialization
and shows the two public `ds` results on one general top integral. `ss11` is realized through
the chain rule `2 ss11 d/d sp[kL,kL]`, and `kE` acts on the shared energy of both vertices.

This example does not call `DSMetaSeedRange`, `DSGenerateIBP`, `DSLinear`, Kira, `DSDE` or
`DSScaleCheck`. It selects no numerical point, target or master, and writes no run artifacts.

Scope is deliberately limited to seeds and operators: no external reduction is registered for
the sunrise, and the two closed-loop cases in this directory are the only examples that go
through Kira. This is the only sunrise example in the package; other examples must not
duplicate it under another name.
