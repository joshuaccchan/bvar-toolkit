# bvar-toolkit

MATLAB code for large Bayesian VARs by [Joshua Chan](https://joshuachan.org) — the packages
distributed at [joshuachan.org/code.html](https://joshuachan.org/code.html), consolidated into a
documented library with the original code preserved verbatim alongside it.

```matlab
run setup.m                 % adds core/ and third_party/ to the path
cd examples
ex03_minnesota_bvar         % a small BVAR, start to finish
```

Requirements: MATLAB with the Statistics and Machine Learning Toolbox, which nearly every
sampler needs. The Optimization Toolbox is needed by the marginal-likelihood functions for
the MA-error models in `bvar.ml`, and by two replication drivers; the System Identification
Toolbox by one replication driver only. `setup.m` checks all three and names what each
missing one will break.

For the methods behind the code, see the book *Bayesian Macroeconometrics: Methods and
Applications* (Chapman & Hall/CRC, forthcoming) —
[sample chapters](https://joshuachan.org/papers/BayesMacroBook_sample.pdf) and
[its own code repository](https://github.com/joshuaccchan/bayesian-macroeconometrics), with
MATLAB, R and Python for all fourteen chapters. `examples/README.md` maps each example to the
chapter that develops it.

## Three ways to use this repo

**Learn the methods.** `examples/` holds seven short scripts, each runnable in seconds and
each printing its reasoning as it goes: the precision sampler, stochastic volatility by
auxiliary mixture, a small BVAR end to end, a BVAR with stochastic volatility assembled from
core blocks, marginal likelihoods and model comparison, and identification by sign
restrictions, both on impact and over a stretch of the impulse response. Each works on data
small enough to check the numbers against the truth, and they are meant to be read in
order. `examples/README.md` lists what each one teaches and which core functions it calls.

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
| Asymmetric conjugate prior, closed-form ML, sign restrictions | Chan (2022, QE) | `chan2022_qe_acp` | `run_all`, `run_jointden` |
| Which SV specification for a large VAR? | Chan (2023, JoE) | `chan2023_joe_mlvarsv` | `run_all('VAR-SV',…)`, `run_ml` |
| Time-varying parameters, equation by equation | Chan (2023, JBES) | `chan2023_jbes_hybtvp` | `run_all` |
| Forecast comparison across priors and volatility models | Chan (2020, Springer) | `chan2020_springer_largebvar` | legacy only |
| The precision sampler for state space models | Chan & Jeliazkov (2009) | `chan_jeliazkov2009_statespace` | legacy only |
| Sign and ranking restrictions in a large structural VAR | Chan, Matthes & Yu (2026, QE) | `chan_matthes_yu2026_qe_svarsign` | legacy; algorithm in `bvar.structural.sign_assign` |
| Prior sensitivity by automatic differentiation | Chan, Jacobi & Zhu (2019/2020/2022) | `cjz2018_ad_var`, `cjz2019_ad_opthyper`, `cjz2021_jae_ad_ml` | legacy only |

"Legacy only" means the package has not been functionized: the original code is there, and a
`run_all.m` may follow. Read it with two qualifications. The SVAR-sign package is not slated
for one, because its main programs depend on third-party code and produce figures rather than
reusable computation, so only its algorithm was extracted. And a few legacy scripts do not run
as shipped; `tests/golden_runs/manifest.md` names them and says why. Full citations are in
`provenance.md`.

## The `bvar` library

A Bayesian VAR is estimated by Markov chain Monte Carlo. Once the prior has been
constructed, each sweep draws the VAR coefficients, the log-volatility path and the
shrinkage hyperparameters in turn, each conditional on the rest; forecasts and marginal
likelihoods are computed afterwards from the stored draws. Across the thirteen packages
those steps were written out again and again — the auxiliary mixture sampler that draws
the log-volatility path appears in seven of them, under three names. `bvar` is those steps
factored into one function each.

They are not rewrites. Each function's body is taken from a specific published package,
and a unit test runs the original code alongside it and requires identical output — draw
for draw, bitwise, under a fixed seed. Calling `bvar.sv.ksc_rw_h0` runs the computation
the paper ran. A small number of functions are new code rather than extractions — they are
marked as such in their headers and listed in `tests/variant_map.md`, and are pinned
instead to the inline spelling they generalize.

| Namespace | What it is for |
|---|---|
| `bvar.priors` | Building priors. `resid_var_ar4`, `minnesota_C` and `vtheta` compute the Minnesota scaling every prior here rests on; `minn`, `niw` and `acp_stru`/`acp_redu` are the prior constructors themselves — Minnesota, natural conjugate, and the asymmetric conjugate prior of Chan (2022), whose marginal likelihood is available in closed form; `acp_opt_kappa` uses that to choose the shrinkage hyperparameters by maximizing it. |
| `bvar.sv` | Drawing stochastic volatility. The `ksc_*` functions are the Kim–Shephard–Chib auxiliary-mixture sampler, one per state equation (random walk with a known initial value, random walk with a diffuse one, stationary AR(1)); `csv_armh` draws a single common volatility factor; `sv_params` and `nu_studentt` draw the parameters governing them. |
| `bvar.samplers` | Drawing everything else in the Gibbs loop: VAR coefficients equation by equation (`eq_gauss` for the structural form, `eq_var_redu_tri` and `eq_svar_oi` for the reduced form, `eq_var_oi` for the same order-invariant conditional as `eq_svar_oi` at `O(T k^2 + k^3)` per equation instead of `O(T n k^2 + k^3)`, `eq_tri_cs` for the Cholesky benchmark), the factor blocks (`factor_fsv`, `eq_fsv_load`), `eq_hyb_tvp` for the hybrid TVP-VAR, where each equation's coefficients are drawn jointly with the indicators for whether they vary at all, `acp_theta_sig` for the asymmetric conjugate prior, whose conjugacy means it returns every draw in one call rather than a chain, and the hierarchical shrinkage blocks (`gig_shrinkage`, `horseshoe_kappa_psi`, `nu_psi_ng`). |
| `bvar.forecast` | Producing forecasts from a chain. `iterate` runs one draw forward and scores it, `tables` accumulates RMSFEs and log predictive likelihoods, `realtime_loaddata` assembles a real-time data vintage. |
| `bvar.structural` | Contemporaneous structure and identification. `construct_Sigt` builds the time-varying covariance from the impact matrix and `b0_row_sampler` draws that matrix row by row for the order-invariant model; `reduced_form` maps structural draws to their reduced form, and `qr_sign`, `sign_restrict` and `irf_redu` are the three steps of a sign-restricted SVAR - draw a rotation, test it against the sign and inequality restrictions, and compute the impulse responses of the draws that survive. `sign_assign` replaces the middle step with the search of Chan, Matthes and Yu (2026), which accepts a rotation whenever every shock has some admissible column rather than requiring the columns to arrive in order. |
| `bvar.ml` | Marginal likelihoods, for model comparison. Chib's method for the VARs with non-Gaussian, heteroscedastic and serially dependent innovations of Chan (2020), adaptive importance sampling for the stochastic volatility specifications of Chan (2023), and `acp`, which is closed form - the property that motivates the asymmetric conjugate prior of Chan (2022), and the reason selecting its hyperparameters is an optimization rather than a second round of estimation. Plus the integrated-likelihood evaluators and log densities the simulation-based ones share. |
| `bvar.util` | The small shared pieces: `build_lags` (the lag matrix, intercept first), `diffmat` (the state-equation difference matrix that makes the precision samplers banded), `surform`/`surform2` (two different sparse expansions — see their headers), `logsumexp`, `igrnd`, and a few one-liners. |

Where two legacy versions of a step turned out to differ numerically, both survive under
separate names rather than being merged: `ksc_rw_h0` and `ksc_rw_diffuse` are the same
sampler under different initial conditions, `resid_var_ar4` and `resid_var_allvars_ridge`
compute the same scaling from different regressions. Where the two differ by a single
setting rather than by the computation, one function takes an option instead, with the
legacy behaviour as the default: `bvar.ml.acp` has a `ridge` argument because one package
adds a jitter to the posterior precision that another does not. `tests/variant_map.md`
records for every function which legacy copies it stands in for, how that was checked, and
a never-merge list of the pairs that must stay apart.

## Verification

Every core function extracted from a published package is covered by a unit test that runs
the corresponding legacy code and requires exact agreement — bitwise, draw-for-draw under a
fixed seed for the stochastic ones. The functionized drivers are tested the same way against
the original scripts in full. The few functions that are new code, listed under new
functions in `tests/variant_map.md`, have no legacy counterpart to run: they are tested
against the inline expression they replace, to floating-point rather than bitwise agreement,
and `eq_var_oi` additionally against the published block it computes more cheaply.

```matlab
run tests/unit/run_unit_tests.m
```

`tests/golden/` holds captured output from the original packages (the log marginal likelihoods,
forecast metric tables and figures they print), with `tests/golden_runs/manifest.md` recording
what was run, how long it took, and which scripts do not run as shipped.

Two of the marginal-likelihood scripts in the Chan (2020, JBES) package — for models 4 and 8 —
contain three places between them where a value left over from the estimation loop is used
where the computation calls for a freshly evaluated one. The core functions compute the
intended quantity by default and reproduce the published computation bitwise under
`'bugcompat', true`.

The VAR-SVO marginal likelihood in the Chan (2023, JoE) package has three defects of its own,
all in its outlier block. `bvar.ml.mlvarsv_arsvo_redu` takes the same approach — corrected by
default, `'bugcompat', true` to reproduce the original — and `run_ml` prints which computation
produced the number.

Both audits are in `tests/variant_map.md`, with the full comparison for the 2020 package.

## Citation

Cite the paper whose code you use — full references in `provenance.md`. For the toolkit
itself, `CITATION.cff` holds the machine-readable record, which GitHub's "Cite this
repository" button reads:

> Chan, J. C. C. *bvar-toolkit: MATLAB code for large Bayesian VARs*. https://github.com/joshuaccchan/bvar-toolkit

## License

MIT — see `LICENSE`, with the scope recorded separately in `NOTICE.md` so that `LICENSE`
stays the unmodified MIT text and tooling can identify it. This relicenses the archived
packages too: their original headers say "free to use for academic purposes only", wording
preserved unaltered as part of the verbatim archive and superseded by the repository
license. Citing the paper you use is expected scholarly practice, not a licensing condition.
Third-party code keeps its own license — the files under `third_party/`, the third-party
files bundled inside some legacy packages, and the whole of the SVAR-sign package's
`auxFunctions/`. `NOTICE.md` lists them.
