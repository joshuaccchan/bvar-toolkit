# Golden runs (step 2 of the consolidation plan)

Before any refactoring, each legacy package is run once on current MATLAB and its printed
output and saved results are captured as *golden files*. Later, the functionized replication
drivers must reproduce these numbers within Monte Carlo tolerance.

## Protocol

1. `run_golden.ps1 -Slug <slug> -Entry <script.m>` does, in order:
   - copies `replications/<slug>/legacy/` to `%LOCALAPPDATA%\bvar-toolkit\golden_runs\<slug>\` (fresh copy;
     deliberately outside the repo - Dropbox sync locks freshly written files and should not sync MCMC scratch),
   - overlays `tests/golden_runs/patches/<slug>/` if that folder exists (minimal run-enablement
     patches only - e.g. a deprecated `xlsread` call or a legacy `rand('seed',...)` line; every
     patch file must carry a header comment stating exactly what was changed and why),
   - locates the entry script inside the copy, `cd`s to its folder, and runs it via
     `matlab -batch` with a diary,
   - copies the diary and every file created or modified by the run into
     `tests/golden/<slug>/<entry>_<date>/`.
2. **Legacy folders are never edited.** All patching happens on the `build/` copy.
3. Runs are stochastic and the original results were produced under clock-seeded legacy RNG,
   so goldens match published tables in distribution, not bit-for-bit. The diary records the
   MATLAB version; comparisons use Monte Carlo tolerances.
4. `manifest.md` lists, per package: entry scripts, toolbox requirements, expected runtime,
   and known blockers. Capture the quick wins first; schedule the multi-day jobs separately
   (and for OISV forecasting, use the shipped `results_mat/*.mat` - the drivers are cluster
   fragments and are not runnable as shipped).
