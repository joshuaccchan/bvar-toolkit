# bvar-toolkit

MATLAB code for large Bayesian VARs by [Joshua Chan](https://joshuachan.org) — the packages
distributed at [joshuachan.org/code.html](https://joshuachan.org/code.html), consolidated into a
documented library with the original code preserved verbatim alongside it.

```matlab
run setup.m                 % adds core/ and third_party/ to the path
cd examples
ex03_minnesota_bvar         % a small BVAR, start to finish
```

Requirements: MATLAB with the Statistics and Machine Learning Toolbox. A few replication
drivers also want Optimization (`fminunc`) or System Identification; `setup.m` warns about
what is missing.

For the methods behind the code, see the book *Bayesian Macroeconometrics: Methods and
Applications* (Chapman & Hall/CRC, forthcoming) —
[sample chapters](https://joshuachan.org/papers/BayesMacroBook_sample.pdf) and
[its own code repository](https://github.com/joshuaccchan/bayesian-macroeconometrics), with
MATLAB, R and Python for all fourteen chapters. Each example in `examples/` names the chapter
that develops it.

## Two ways to use this repo

**Reproduce a paper.** Every package is here exactly as published, never edited, under
`replications/<paper>/legacy/`, with a permanent `as-published/<paper>` git tag and the source
zip's md5 recorded in `provenance.md`. Run those files as you would the original download.

**Build on the code.** The samplers, priors, and forecasting machinery are factored into the
`bvar` package under `core/`, each function tested to reproduce its legacy counterpart
draw-for-draw under a fixed seed. Call the blocks directly, or copy the nearest `run_all.m` as
a template.

## Which model do I want?

| If you want | Paper | Folder | Driver |
|---|---|---|---|
| Shrinkage priors for a large BVAR (the default choice) | Chan (2021, IJF) | `chan2021_ijf_mahp` | `run_all('MNG',…)` |
| A VAR-SV that does not depend on variable ordering | Chan, Koop & Yu (2024, JBES) | `chan_koop_yu2024_jbes_oisv` | `run_all('OI',…)` |
| Non-Gaussian / serially dependent errors, and marginal likelihoods | Chan (2020, JBES) | `chan2020_jbes_kronecker` | `run_all`, `run_ml` |
| Asymmetric conjugate prior, closed-form ML, sign restrictions | Chan (2022, QE) | `chan2022_qe_acp` | legacy only |
| Which SV specification for a large VAR? | Chan (2023, JoE) | `chan2023_joe_mlvarsv` | legacy only |
| Time-varying parameters, decided per equation | Chan (2023, JBES) | `chan2023_jbes_hybtvp` | legacy only |
| Forecast comparison across priors and volatility models | Chan (2020, Springer) | `chan2020_springer_largebvar` | legacy only |
| The precision sampler for state space models | Chan & Jeliazkov (2009) | `chan_jeliazkov2009_statespace` | legacy only |
| Prior sensitivity by automatic differentiation | Chan, Jacobi & Zhu (2019/2020/2022) | `cjz2018_ad_var`, `cjz2019_ad_opthyper`, `cjz2021_jae_ad_ml` | legacy only |

"Legacy only" means the package has not been functionized yet — the original code is there and
runs; a `run_all.m` will follow. Full citations are in `provenance.md`.

## The `bvar` library

A Bayesian VAR is estimated by Markov chain Monte Carlo. Once the prior has been
constructed, each sweep draws the VAR coefficients, the log-volatility path and the
shrinkage hyperparameters in turn, each conditional on the rest; forecasts and marginal
likelihoods are computed afterwards from the stored draws. Across the twelve packages
those steps were written out again and again — the auxiliary mixture sampler that draws
the log-volatility path appears in eight of them, under four names. `bvar` is those steps
factored into one function each.

They are not rewrites. Each function's body is taken from a specific published package,
and a unit test runs the original code alongside it and requires identical output — draw
for draw, bitwise, under a fixed seed. Calling `bvar.sv.ksc_rw_h0` runs the computation
the paper ran.

| Namespace | What it is for |
|---|---|
| `bvar.priors` | Building priors. `resid_var_ar4`, `minnesota_C` and `vtheta` compute the Minnesota scaling every prior here rests on; `minn`, `niw` and `acp_stru`/`acp_redu` are the prior constructors themselves — Minnesota, natural conjugate, and the asymmetric conjugate prior of Chan (2022) whose marginal likelihood is available in closed form. |
| `bvar.sv` | Drawing stochastic volatility. The `ksc_*` functions are the Kim–Shephard–Chib auxiliary-mixture sampler, one per state equation (random walk with a known initial value, random walk with a diffuse one, stationary AR(1)); `csv_armh` draws a single common volatility factor; `sv_params` and `nu_studentt` draw the parameters governing them. |
| `bvar.samplers` | Drawing everything else in the Gibbs loop: VAR coefficients equation by equation (`eq_gauss` for the structural form, `eq_svar_oi` for the order-invariant one, `eq_tri_cs` and `alp_tri_cs` for the Cholesky benchmark) and the hierarchical shrinkage blocks (`gig_shrinkage`, `horseshoe_kappa_psi`, `nu_psi_ng`). |
| `bvar.forecast` | Producing forecasts from a chain. `iterate` runs one draw forward and scores it, `tables` accumulates RMSFEs and log predictive likelihoods, `realtime_loaddata` assembles a real-time data vintage. |
| `bvar.structural` | Contemporaneous structure: `construct_Sigt` builds the time-varying covariance from the impact matrix, `b0_row_sampler` draws that matrix row by row for the order-invariant model. |
| `bvar.ml` | Marginal likelihoods, for model comparison. Chib's method for the VARs with non-Gaussian, heteroscedastic and serially dependent innovations of Chan (2020), plus the integrated-likelihood evaluators and log densities it needs. |
| `bvar.util` | The small shared pieces: `build_lags` (the lag matrix, intercept first), `diffmat` (the state-equation difference matrix that makes the precision samplers banded), `surform`/`surform2` (two different sparse expansions — see their headers), `logsumexp`, `igrnd`, and a few one-liners. |

Where two legacy versions of a step turned out to differ numerically, both survive under
separate names rather than being merged: `ksc_rw_h0` and `ksc_rw_diffuse` are the same
sampler under different initial conditions, `resid_var_ar4` and `resid_var_allvars_ridge`
compute the same scaling from different regressions. `tests/variant_map.md` records for
every function which legacy copies it stands in for, how that was checked, and a
never-merge list of the pairs that must stay apart.

## Examples

`examples/` holds short scripts, each runnable in seconds, from the precision sampler up to a
BVAR with stochastic volatility assembled from core blocks. See `examples/README.md`.

## Verification

Every core function is covered by a unit test that runs the corresponding legacy code and
requires exact agreement — bitwise, draw-for-draw under a fixed seed for the stochastic ones.
The functionized drivers are tested the same way against the original scripts in full.

```matlab
run tests/unit/run_unit_tests.m
```

`tests/golden/` holds captured output from the original packages (the log marginal likelihoods,
forecast metric tables and figures they print), with `tests/golden_runs/manifest.md` recording
what was run, how long it took, and which scripts do not run as shipped.

Three marginal-likelihood scripts in the Chan (2020, JBES) package evaluate an ordinate at a leftover
chain value where the posterior mean is intended. The core functions fix this by default and
reproduce the published computation under `'bugcompat', true`; the corrections change the
reported values by at most 2.45 log points and do not affect the paper's model ranking. The
audit and the full comparison are in `tests/variant_map.md`.

## Citation

Cite the paper whose code you use — full references in `provenance.md`. For the toolkit itself:

> Chan, J. C. C. *bvar-toolkit: MATLAB code for large Bayesian VARs*. https://github.com/joshuaccchan/bvar-toolkit

## License

MIT — see `LICENSE`. This relicenses the archived packages too: their original headers say
"free to use for academic purposes only", wording preserved unaltered as part of the verbatim
archive and superseded by the repository license. Citing the paper you use is expected
scholarly practice, not a licensing condition. Third-party code keeps its own license — the
files under `third_party/`, and the third-party files bundled inside some legacy packages
(`gigrnd.m`, `EvalFore.m`, `heatmap.m`).
