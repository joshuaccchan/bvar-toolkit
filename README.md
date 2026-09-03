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
`bvt` package under `core/`, each function tested to reproduce its legacy counterpart
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

## What is in `core/+bvt`

| Namespace | Contents |
|---|---|
| `bvt.util` | `build_lags`, `diffmat` (state-equation difference matrix), `surform`/`surform2` (different operators — see their headers), `igrnd`, `logsumexp`, `vec`, `vech`, `ldet`, `mgammaln`, `tnormrnd`, `anormrnd` |
| `bvt.priors` | Minnesota scaling (`resid_var_ar4`, `minnesota_C`, `vtheta`), constructors `minn`, `niw`, `acp_stru`/`acp_redu`, `impact_B0` |
| `bvt.sv` | KSC auxiliary-mixture samplers (`ksc_rw_h0`, `ksc_rw_diffuse`, `ksc_ar1_mean`), common SV (`csv_armh`), SV parameters (`sv_params`, `sv0_params`), t degrees of freedom (`nu_studentt`) |
| `bvt.samplers` | `eq_gauss`, `gig_shrinkage`, `nu_psi_ng`, `eq_svar_oi`, `eq_tri_cs`, `alp_tri_cs`, `horseshoe_kappa_psi` |
| `bvt.forecast` | `iterate` (per-draw forecasts and predictive likelihoods), `tables` (RMSFE / ALPL), `realtime_loaddata` |
| `bvt.structural` | `construct_Sigt`, `b0_row_sampler` |
| `bvt.ml` | Marginal likelihoods for the Kronecker model family, integrated-likelihood evaluators, log-density utilities |

Variants that look interchangeable but are not — two `SVRW` samplers with different initial
conditions, two `prior_NCP` signatures, `get_resid_var` versus its all-variables variant — are
kept as separately named functions. `tests/variant_map.md` records, for every core function,
which legacy copies it canonicalizes and how that was verified, plus a never-merge list.

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

Three marginal-likelihood scripts in the Kronecker package evaluate an ordinate at a leftover
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
