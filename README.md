# bvar-toolkit <img src=".github/logo.png" align="right" height="139" alt="bvar-toolkit logo" />

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22804032.svg)](https://doi.org/10.5281/zenodo.22804032) [![View on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/184885-bvar-toolkit-large-bayesian-vars)

A MATLAB library for large Bayesian VARs: samplers, shrinkage priors, stochastic volatility,
marginal likelihoods and forecasting. The repository archives fourteen replication packages from
[joshuachan.org](https://joshuachan.org/code.html) exactly as published and factors their shared
machinery into the library, where unit tests check it against the published code. The
[how-to guides](guides/), [examples](examples/) and [tutorials](tutorials/) show how to use it.

A Python version, [bvar-toolkit-python](https://github.com/joshuaccchan/bvar-toolkit-python), has
the same functions under the same names, each tested against its MATLAB twin, and the tutorials
and how-to guides with Python code.

## Installation

Download the repository as a ZIP file (the Code button, then Download ZIP) and unzip it, or clone
it with git. The latest release is also on
[File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/184885-bvar-toolkit-large-bayesian-vars).
In MATLAB, go to the folder that contains `setup.m` and run it:

```matlab
run setup.m                 % adds core/ and third_party/ to the path
```

The code needs MATLAB with the Statistics and Machine Learning Toolbox. `setup.m` also checks
for the Optimization and System Identification Toolboxes, which a few functions and drivers
use, and names what each missing one would break. The library and the examples are tested on
MATLAB R2025b, on Windows and, in the repository's continuous integration, on Linux.

## A First Example

The asymmetric conjugate prior of Chan (2022) shrinks the coefficients on a variable's own lags
and on other variables' lags by different amounts, and its marginal likelihood is available in
closed form. For a VAR of GDP growth, PCE inflation and the 3-month Treasury bill rate with
four lags, the two hyperparameters that maximize it take one call:

```matlab
data = readmatrix('replications/chan2019wp_acp/legacy/macrodata_Q_2018Q4.csv');
data = data(1:238, [1 95 145]);
p = 4;
Y0 = data(1:8, :);                          % initial conditions
Y = data(9:end, :);                         % the estimation sample
[~, Z] = bvar.util.build_lags([Y0(end-p+1:end, :); Y], p);
[lml, kappa] = bvar.priors.acp_opt_kappa(Y0, Y, Z, p, [.04 .04], 'stru');
```

Here `kappa(1)` and `kappa(2)` are the own-lag and other-lag hyperparameters, 0.542 and
0.0128, and `lml` is the maximized log marginal likelihood, -1237.5. The
[shrinkage tutorial](tutorials/shrinkage/) takes the same steps for 21 variables. For a small
BVAR from start to finish, run [`ex03_minnesota_bvar`](examples/ex03_minnesota_bvar.m) in
`examples/`.

## Learn the Methods

The [how-to guides](guides/) show how to use the core functions, one step at a time, on a small
example with generated data. The first series, Getting Started with Bayesian VARs, stacks the
equations of a VAR, draws its coefficients under the Minnesota prior, estimates it by Gibbs
sampling, and chooses its shrinkage by the marginal likelihood:

- [How to Stack the Equations of a VAR](guides/var_stacking/)
- [How to Draw the VAR Coefficients under the Minnesota Prior](guides/coefficient_draw/)
- [How to Estimate a VAR by Gibbs Sampling](guides/gibbs_sampler/)
- [How to Choose the Shrinkage of a VAR by the Marginal Likelihood](guides/natural_conjugate_ml/)

The second series, Forecasting with a VAR, works with data simulated from a VAR fitted to US GDP
growth, inflation and unemployment, so that every forecast can be compared with the forecast made
with the true parameters. It simulates the data, forecasts with a BVAR, evaluates the forecasts
out of sample, and chooses the lag length and the shrinkage:

- [How to Simulate Data from an Estimated VAR](guides/var_simulation/)
- [How to Forecast with a Bayesian VAR](guides/var_forecast/)
- [How to Evaluate Forecasts Out of Sample](guides/var_forecast_evaluation/)
- [How to Choose the Lag Length and the Shrinkage for Forecasting](guides/var_lag_shrinkage/)

The third series, VARs for a Small Open Economy, works with data simulated from a VAR fitted to
US and Australian data, in which the Australian variables do not enter the US equations. It
estimates the VAR block by block, tests this restriction, traces a US shock through the Australian
economy, and forecasts Australia under a US scenario:

- [How to Estimate a VAR for a Small Open Economy](guides/soe_var/)
- [How to Test Block Exogeneity by the Marginal Likelihood](guides/soe_exogeneity_test/)
- [How to Trace a Foreign Shock through a Small Open Economy](guides/soe_foreign_shock/)
- [How to Forecast a Small Open Economy under a Foreign Scenario](guides/soe_conditional_forecast/)

The [tutorials](tutorials/) each answer one research question with the settings of a paper, on
the data of its replication package, of the book or of FRED, and can be read without MATLAB.
Each comes with a script that runs the same analysis on your own data.

- [Does the Order of the Variables Change My VAR Results?](tutorials/variable_ordering/) Under
  the Cholesky model of stochastic volatility it does: reversing the order of a 20-variable VAR
  moves one estimated correlation from about 0.82 to about 0.19, while the order-invariant model
  gives the same estimates in both orders.
- [How Should I Choose the Shrinkage Hyperparameters of My BVAR?](tutorials/shrinkage/) By
  maximizing the marginal likelihood under the asymmetric conjugate prior of Chan (2022), which
  in a 21-variable VAR improves the one-quarter-ahead point forecasts of 20 of the 21 variables
  relative to the best symmetric prior.
- [Which Stochastic Volatility Specification Should My VAR Use?](tutorials/sv_specification/)
  The marginal likelihoods of Chan (2023) rank five of them; in a five-variable quarterly VAR,
  the Cholesky specification with an outlier component fits best.
- [Does Modeling the Volatility Improve My Forecasts?](tutorials/forecasting/) Over 140
  recursive forecasts from 1990, a common volatility factor and a volatility process per
  equation both improve on a homoskedastic VAR in point and density forecasts, with gains in
  calm quarters as well as in the pandemic.
- [Can I Use Many Sign Restrictions in My Large VAR?](tutorials/sign_restrictions/) Yes, with
  the algorithm of Chan, Matthes and Yu (2026), which imposes 98 sign and 7 ranking restrictions
  on eight shocks in a 35-variable VAR, where the conventional accept-reject algorithm is not
  computationally feasible.
- [How Do I Handle Missing and Mixed-Frequency Data in My VAR?](tutorials/mixed_frequency/) By
  using data augmentation and drawing all the missing values in one block, as Chan, Poon and Zhu
  (2023) do. In a monthly VAR with quarterly GDP, this puts the fall in monthly GDP in April 2020
  at 7.9 percent.

The fourteen scripts in [`examples/`](examples/) each run in under a minute and are best read in
order: the building blocks (ex01–ex02, the precision sampler and stochastic volatility), VAR
specifications (ex03–ex08), model comparison and forecasting with an estimated VAR
(ex09–ex10), structural identification by sign restrictions (ex11–ex12), diagnostics for
the output of a sampler (ex13), and missing data and mixed frequencies (ex14).
[`examples/README.md`](examples/README.md) lists what each one teaches and which library
functions it calls.

The methods are developed in the book *Bayesian Macroeconometrics: Methods and Applications*
(Chapman & Hall/CRC, 2027): see the
[sample chapters](https://joshuachan.org/papers/BayesMacroBook_sample.pdf) and
[its code repository](https://github.com/joshuaccchan/bayesian-macroeconometrics), with
MATLAB, R and Python for all fourteen chapters. [`RELATED.md`](RELATED.md) maps the models here
to the book's R and Python scripts and to the R package `bvars`.

## Reproduce a Paper

Every package is here exactly as published, never edited, under
`replications/<paper>/legacy/`, with a permanent `as-published/<paper>` git tag and the source
zip's md5 recorded in [`provenance.md`](provenance.md). Run those files as you would the
original download. Seven packages also have a driver, `run_all.m`, that runs the same computation
with the library functions.

### Which Model Do I Want?

| If you want | Paper | Folder | Driver |
|---|---|---|---|
| Shrinkage priors for a large BVAR (the default choice) | Chan (2021, IJF) | `chan2021_ijf_mahp` | `run_all('MNG',…)` |
| A VAR-SV that does not depend on variable ordering | Chan, Koop & Yu (2024, JBES) | `chan_koop_yu2024_jbes_oisv` | `run_all('OI',…)` |
| Non-Gaussian / serially dependent errors, and marginal likelihoods | Chan (2020, JBES) | `chan2020_jbes_kronecker` | `run_all`, `run_ml` |
| Asymmetric conjugate prior, closed-form ML, sign restrictions | Chan (2022, QE) | `chan2022_qe_acp` | `run_all`, `run_jointden` |
| Which SV specification for a large VAR? | Chan (2023, JoE) | `chan2023_joe_mlvarsv` | `run_all('VAR-SV',…)`, `run_ml` |
| Time-varying parameters, equation by equation | Chan (2023, JBES) | `chan2023_jbes_hybtvp` | `run_all` |
| Which time-varying parameter VAR, by marginal likelihood or DIC | Chan & Eisenstat (2018, JAE) | `chan_eisenstat2018_jae_mltvpsv` | `run_all('TVP-SV',…)`, `run_ml` |
| Forecast comparison across priors and volatility models | Chan (2020, Springer) | `chan2020_springer_largebvar` | legacy only |
| The precision sampler for state space models | Chan & Jeliazkov (2009) | `chan_jeliazkov2009_statespace` | legacy only |
| Sign and ranking restrictions in a large structural VAR | Chan, Matthes & Yu (2026, QE) | `chan_matthes_yu2026_qe_svarsign` | legacy; algorithm in `bvar.structural.sign_assign` |
| Prior sensitivity by automatic differentiation | Chan, Jacobi & Zhu (2019/2020/2022) | `cjz2018_ad_var`, `cjz2019_ad_opthyper`, `cjz2021_jae_ad_ml` | legacy only |

"Legacy only" means the package has no driver yet. The SVAR-sign package will not get one,
because its main programs depend on third-party code and produce figures, so only its algorithm
was extracted. The precision sampler of `chan_jeliazkov2009_statespace` is the function
`ssm.simulate_states` in [statespace-toolkit](https://github.com/joshuaccchan/statespace-toolkit),
the companion library for Bayesian state space models, whose unit tests check it against that
package's `UC.m`. [`tests/golden_runs/manifest.md`](tests/golden_runs/manifest.md) names the
few legacy scripts that do not run as shipped. What to cite for each method, with BibTeX, is in
[`CITING.md`](CITING.md).

## Build on the Code

Across the fourteen packages the same steps of the Gibbs sampler were written out again and
again: the auxiliary mixture sampler for the log-volatility path appears in eight of them, under
three names. The `bvar` package under `core/` has each step once. Call the blocks directly,
start from `bvar.models.var_sv`, the complete sampler for a VAR with stochastic volatility, or
copy the nearest `run_all.m` as a template.

| Namespace | Contents |
|---|---|
| `bvar.priors` | Minnesota, natural conjugate and asymmetric conjugate priors (`minn`, `niw`, `acp_stru`, `acp_redu`), the Minnesota scaling they share (`resid_var_ar4`, `minnesota_C`, `vtheta`), and `acp_opt_kappa`, which chooses the shrinkage hyperparameters by maximizing the closed-form marginal likelihood. |
| `bvar.sv` | Stochastic volatility: the auxiliary mixture sampler for three state equations (`ksc_rw_h0`, `ksc_rw_diffuse`, `ksc_ar1_mean`), a common volatility factor (`csv_armh`, and `csv_armh_block`, which updates it one block at a time), the state-equation parameters (`sv_params`, `sv0_params`) and the Student-t degrees of freedom (`nu_studentt`). |
| `bvar.samplers` | The other Gibbs blocks: VAR coefficients equation by equation (`eq_gauss`, `eq_var_redu_tri`, `eq_svar_oi`, `eq_var_oi`, `eq_tri_cs`, `alp_tri_cs`), factor stochastic volatility (`factor_fsv`, `eq_fsv_load`), the hybrid TVP-VAR (`eq_hyb_tvp`), direct draws under the asymmetric conjugate prior (`acp_theta_sig`), hierarchical shrinkage (`gig_shrinkage`, `horseshoe_kappa_psi`, `nu_psi_ng`), a joint draw of all VAR coefficients under an independent normal prior (`var_coef_joint`) and the missing values of a VAR in one block, under hard or soft linear restrictions such as a quarterly aggregation (`missing_var`). |
| `bvar.structural` | Impact matrices and identification: the order-invariant impact matrix (`b0_row_sampler`, `construct_Sigt`), the map to the reduced form (`reduced_form`) and sign restrictions (`qr_sign`, `sign_restrict`, `sign_assign`, `irf_redu`), and two summaries of the identified shocks: the forecast error variance decomposition (`fevd`) and the historical decomposition (`hist_decomp`). |
| `bvar.ml` | Marginal likelihoods: Chib's method for the models of Chan (2020, JBES) (`kron_bvar*`), adaptive importance sampling for those of Chan (2023, JoE) (`mlvarsv_*`), importance sampling for the time-varying, constant-coefficient and regime-switching VARs of Chan and Eisenstat (2018, JAE) (`mltvpsv_*`), the closed forms under the natural conjugate prior (`ncp`) and the asymmetric conjugate prior (`acp`), and the integrated likelihoods and log densities they share. |
| `bvar.forecast` | Forecasts from a chain: `iterate` runs one draw forward and scores it, `predictive` returns the mean and standard deviation of the h-step predictive distribution of every draw and `simulate` does the same by simulation when the error covariance varies over time, `tables` accumulates RMSFEs and log predictive likelihoods, and `realtime_loaddata` assembles a real-time data vintage. |
| `bvar.diag` | Diagnostics for MCMC output: inefficiency factors (`inefficiency_factor`), Monte Carlo standard errors (`mcse`) and Geweke's convergence diagnostic (`geweke`), each from the long-run variance that `specvar0` estimates. |
| `bvar.models` | Complete samplers: `var_sv`, a VAR with Cholesky or order-invariant stochastic volatility under the prior of Chan, Koop and Yu (2024), `var_csv`, a VAR with one common volatility factor, and `mfvar_csv`, the same model with missing data and mixed frequencies under the priors of Chan, Poon and Zhu (2023). Each can also return its draws. `var_ncp` takes independent draws of a VAR under the natural conjugate prior, with the shrinkage chosen by the marginal likelihood, and `var_soe` those of a VAR for a small open economy with block exogeneity, one block at a time. |
| `bvar.util` | Shared pieces: the lag matrix (`build_lags`), data simulated from a VAR (`simulate_var`), the reduced form of a VAR with block exogeneity (`soe_reduced_form`), the difference matrix of the precision samplers (`diffmat`), sparse expansions (`surform`, `surform2`), credible bands for figures (`shaded_band`), the selection matrices of a missing data pattern (`select_obs`), the restrictions that tie quarterly growth rates to monthly ones (`mm_constraint`) and growth rates across a gap in a series to the change over it (`dlog_gaps`), the csv and mat files a script writes at the end (`report`) and small numerical helpers. |

Where two legacy versions of a step differ numerically, both survive under separate names.
[`tests/variant_map.md`](tests/variant_map.md) records which legacy copies each function
stands in for and which pairs must never be merged.

## Verification

Most library functions are extracted from a published package, and a unit test runs the
original code alongside the function and requires identical output, draw for draw under a fixed
seed. Where the library factors a matrix once and the original code factors it twice, the test
first applies the substitutions declared in `tests/unit/private/one_factor_patch.m` to the
original code. Where a function has since been improved, for example with an option in place
of a hard-coded constant, its default still reproduces the published computation and the
improvement has a test of its own. New functions, such as `bvar.models.var_sv`, are tested
against the code they replace. The marginal likelihood code of two packages, Chan (2020, JBES)
and Chan (2023, JoE), contains defects; the library corrects them by default and reproduces the
published computation with `'bugcompat', true`. On every push, CI checks the archived packages
against their `as-published` tags and runs the unit suite, `setup.m`, the examples, each
tutorial's `your_data.m` and each guide's `guide.m`. To run the tests locally:

```matlab
run tests/unit/run_unit_tests.m
run tests/run_examples.m
```

[`tests/variant_map.md`](tests/variant_map.md) has the full record, including the audits of the
two packages and every deviation from the legacy code, and
[`tests/golden_runs/manifest.md`](tests/golden_runs/manifest.md) lists the captured output of
the original packages in `tests/golden/`.

## Citation

Cite the paper behind each method you use: [`CITING.md`](CITING.md) maps every replication
package, library function and example to its paper, with BibTeX. To cite the toolkit itself:

> Chan, J. C. C. (2026). *bvar-toolkit: MATLAB code for large Bayesian VARs*. Zenodo. https://doi.org/10.5281/zenodo.22804032

The DOI resolves to the latest release. `CITATION.cff` holds the machine-readable record that
GitHub's "Cite this repository" button reads.

## License

MIT; see [`LICENSE`](LICENSE). The license covers the archived packages too: their original
headers, which say "free to use for academic purposes only", are kept verbatim and superseded
by it. Third-party files keep their own licenses, and [`NOTICE.md`](NOTICE.md) lists them.
