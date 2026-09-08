# Notice

`LICENSE` is the MIT license and nothing else, so that GitHub and other tooling
identify it as MIT rather than "Other". This file carries the scope notes that
used to sit above the grant.

## Scope

The license covers the code in this repository written by the copyright holder.
That includes the archived packages under `replications/*/legacy/`, whose
original file headers carry earlier and more restrictive wording ("This code is
free to use for academic purposes only, provided that the paper is cited as:
..."). Those headers are preserved unaltered as part of the verbatim archive and
are superseded by the repository license; citing the corresponding paper remains
expected scholarly practice rather than a licensing condition.

## Not covered

Third-party code retains its original license and authorship: everything under
`third_party/`, and the third-party files bundled inside some legacy packages.
That includes the whole of
`replications/chan_matthes_yu2026_qe_svarsign/legacy/auxFunctions/` — 25 files
implementing the algorithm of Read (2022), copied from the code that paper
provides — as well as:

- `gigrnd.m` (Makalic and Schmidt, 2015)
- `EvalFore.m` (Roque Montero, 2016)
- `tmult.m`, from the AD_OptHyper package's `MatCode/`
- `heatmap.m` (MathWorks File Exchange), bundled inside the OISV package and
  vendored here as `third_party/heatmap_fx.m` so that it does not shadow
  MATLAB's built-in `heatmap`

See those files' own headers.
