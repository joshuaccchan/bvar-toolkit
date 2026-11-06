# Examples

Fourteen short scripts in six groups: the two computational building blocks the toolkit rests
on, VAR specifications of increasing flexibility, model comparison and forecasting with
an estimated VAR, structural identification by sign restrictions, diagnostics for the
output of a sampler, and missing data and mixed frequencies. Each prints its
reasoning as it goes. Seven of them (ex01, ex02, ex04, ex05, ex07, ex08 and ex14) run on
simulated data, so their estimates can be checked against the truth. The other seven use data
from the replication packages. For reproducing a published table, use `replications/<paper>/`
instead.

Every script puts the toolkit on the path itself, so any of them runs from a clean session:

```matlab
cd examples
ex01_precision_sampler
```

Read them in order; each group builds on the ones before it.

### Building blocks

| | Script | What it teaches | Book | Runs in |
|---|---|---|---|---|
| 1 | `ex01_precision_sampler.m` | Drawing an entire state path in one block, with no filtering recursion — the Chan–Jeliazkov (2009) precision sampler. Derives it from the banded precision matrix, checks the estimated path against the simulated truth, and shows the sparse structure that makes it linear in *T*. | Ch. 9 | 11 s |
| 2 | `ex02_sv_ksc.m` | Stochastic volatility by the Kim–Shephard–Chib auxiliary mixture: how squaring and logging the data turns a nonlinear model into the linear Gaussian one ex01 already solves, and what the seven-component mixture is for. | Ch. 10 | 5 s |

### VAR specifications

| | Script | What it teaches | Book | Runs in |
|---|---|---|---|---|
| 3 | `ex03_minnesota_bvar.m` | A small BVAR end to end with a Minnesota / natural-conjugate prior — how the prior is built, what the shrinkage hyperparameter does, and why the natural-conjugate restriction yields an analytic posterior. No MCMC: samples are directly drawn from the posterior. | Ch. 12 | 3 s |
| 4 | `ex04_bvar_sv_blocks.m` | Assembling a reduced-form BVAR with stochastic volatility from core blocks, drawn equation by equation — the sampler of `VAR_ARSV_redu.m` from Chan (2023, JoE), on simulated data with the truth known. | Ch. 14 | 4 s |
| 5 | `ex05_var_csv.m` | A BVAR whose whole covariance matrix is scaled by one common volatility path: the common stochastic volatility of Carriero, Clark and Marcellino (2016, JBES), implemented in Chan (2023, JoE), with the accept-reject Metropolis step of Chan (2020, JBES). Given that path the model is a weighted VAR under a natural-conjugate prior, so every coefficient and the covariance come from a single draw, against the equation-by-equation loop of ex04. On simulated data the path comes back with correlation 0.93. The script also runs the accept-reject Metropolis step at three values of its envelope constant, and estimates the model on data whose four equations have volatility paths of their own, where one path can only follow their average. | Ch. 14 | 3 s |
| 6 | `ex06_variable_ordering_sv.m` | Whether the estimates depend on the order of the variables: the Cholesky SV model of ex04 against the order-invariant model of Chan, Koop and Yu (2024), each run in the published order and reversed, with a second seed as a Monte Carlo yardstick. On four FRED-MD series the ordering moves the Cholesky correlation paths by up to 0.11 and the order-invariant ones by about as much as a change of seed. | Ch. 13, 14 | 34 s |
| 7 | `ex07_factor_sv.m` | A VAR whose errors load on two latent factors, with stochastic volatility on each factor and each idiosyncratic error, so the covariance of eight variables moves with 10 volatility paths and 13 free loadings instead of 36 free elements at every date. On simulated data the factors come back with correlations near 0.97, the loadings track their true values, and the 68% bands cover the truth in about two thirds of the element-dates. | Ch. 14 | 8 s |
| 8 | `ex08_hybrid_tvp_var.m` | A hybrid TVP-VAR in the structural form of Chan (2023, JBES): each equation has two switches, for whether its VAR coefficients and its impact-matrix elements vary over time, drawn with the state paths integrated out. On data simulated from a different configuration in each equation, the posterior puts 1.00, 0.85 and 1.00 on the true configurations. The prior on the size of the time variation is wider than the paper's preset, as the header states. | Ch. 13 | 14 s |

### Using an estimated VAR

| | Script | What it teaches | Book | Runs in |
|---|---|---|---|---|
| 9 | `ex09_marginal_likelihood.m` | Marginal likelihoods and model comparison: the three pieces of Chib's identity, why the posterior ordinate is *subtracted* (the Ockham factor), and a three-model comparison from Chan (2020, JBES). | Ch. 5 | 1 s |
| 10 | `ex10_forecast_evaluation.m` | A recursive forecasting exercise end to end: the structural VAR with stochastic volatility and the Minnesota-type adaptive hierarchical prior of Chan (2021), estimated at each vintage, forecast one and four quarters ahead, and scored by RMSFE and log predictive likelihood. The only example that uses `bvar.forecast`. | Ch. 14 | 15 s |

### Structural identification

| | Script | What it teaches | Book | Runs in |
|---|---|---|---|---|
| 11 | `ex11_sign_restrictions.m` | Identifying a structural VAR by sign restrictions, and the cost of the search: two acceptance rules run over the same posterior draws and the same rotations, one requiring each shock in its own column and one searching over assignments. | &mdash; | 5 s |
| 12 | `ex12_dynamic_sign_restrictions.m` | Sign restrictions imposed at several horizons of the impulse response — Uhlig's (2005) agnostic identification of a monetary shock. Collects 1000 accepted draws under each acceptance rule and compares the resulting credible bands: the two rules have different acceptance rates and produce the same bands. | &mdash; | 13 s |

### Checking a sampler

| | Script | What it teaches | Book | Runs in |
|---|---|---|---|---|
| 13 | `ex13_mcmc_diagnostics.m` | Inefficiency factors, Monte Carlo standard errors and Geweke's convergence diagnostic, from `bvar.diag`, for the parameter draws of the two samplers of ex06. The VAR coefficients mix well under both models, while the shrinkage parameters, the volatility parameters and the diagonal of the order-invariant impact matrix need tens to over a hundred draws for each independent one. The inefficiency factor rises with the truncation lag until the lag passes the autocorrelations, and Geweke's statistic, computed with two lags, separates slow mixing from a chain that has not settled. | Ch. 6 | 35 s |

### Missing data and mixed frequencies

| | Script | What it teaches | Book | Runs in |
|---|---|---|---|---|
| 14 | `ex14_mixed_frequency.m` | A monthly VAR with common stochastic volatility in which one series is observed only as a quarterly growth rate: its monthly values are missing data, drawn in one block from their banded conditional distribution under the aggregation of Mariano and Murasawa (2003), the approach of Chan, Poon and Zhu (2023). On simulated data the monthly values come back with correlation 0.94, against 0.88 for a path that sets each month to a third of its quarter's growth, and 91 percent of them lie inside the 90 percent bands. The tutorial [How Do I Handle Missing and Mixed-Frequency Data in My VAR?](../tutorials/mixed_frequency/) applies it to US data. | &mdash; | 11 s |

Timings are from one warm R2025b session on a desktop machine; treat them as orders of
magnitude. All but ex09 and ex11 draw figures as well as printing. The chapter column
refers to *Bayesian Macroeconometrics: Methods and Applications* (Chan, Chapman &
Hall/CRC, 2027), whose [sample chapters](https://joshuachan.org/papers/BayesMacroBook_sample.pdf)
and [code repository](https://github.com/joshuaccchan/bayesian-macroeconometrics) are
online: 9 Linear Gaussian State Space Models, 10 Stochastic Volatility Models, 12
Vector Autoregressions, 13 Time-Varying Vector Autoregressions, 14 Large VARs with
Stochastic Volatility, 5 Bayesian Model Comparison, 6 Foundations of Bayesian Computation. The auxiliary mixture ex02 uses is developed in Chapter 4, Mixture Models.
Each script repeats its chapter in the header. Five scripts draw on Chapter 14: ex04 for
the sampler, ex05 for the common-volatility model, ex06 and ex07 for the order-invariant
and factor models, and ex10 for the forecasting exercise. ex11 and ex12 have no entry because the book identifies structural
VARs recursively, in Chapter 12, and cites sign restrictions only as further reading. ex14 has none
because the book cites Chan, Poon and Zhu (2023) only as further reading, in Chapters 11 and 14.

## What each one exercises

Useful if you are looking for a worked call of a particular core function.

| Script | Core functions called | Data |
|---|---|---|
| ex01 | `bvar.util.surform`, `bvar.util.shaded_band` — the precision sampler itself is written out line by line, since deriving it is the point | simulated |
| ex02 | `bvar.sv.ksc_rw_h0`, `bvar.util.shaded_band` | simulated |
| ex03 | `bvar.priors.minn`, `bvar.priors.niw`, `bvar.priors.resid_var_ar4`, `bvar.util.build_lags` | `replications/chan2020_jbes_kronecker/legacy/data_Q.csv`, read-only |
| ex04 | `bvar.priors.minn`, `bvar.priors.impact_B0`, `bvar.samplers.alp_tri_cs`, `bvar.sv.ksc_ar1_mean`, `bvar.sv.sv_params`, `bvar.sv.init_approx1N`, `bvar.util.build_lags`, `bvar.util.vec`, `bvar.util.shaded_band` | simulated |
| ex05 | `bvar.sv.csv_armh`, `bvar.sv.sv0_params`, `bvar.priors.niw`, `bvar.util.build_lags`, `bvar.util.shaded_band` — the only example that calls the common-volatility sampler | simulated |
| ex06 | `bvar.models.var_sv`, whose body calls `bvar.structural.b0_row_sampler`, `bvar.structural.construct_Sigt`, `bvar.samplers.eq_var_oi`, `bvar.samplers.eq_tri_cs`, `bvar.samplers.alp_tri_cs`, `bvar.samplers.horseshoe_kappa_psi`, `bvar.sv.ksc_ar1_mean`, `bvar.sv.sv0_params`, `bvar.sv.sv_params`, `bvar.priors.resid_var_ar4`, `bvar.priors.minnesota_C`, `bvar.priors.vtheta` and `bvar.util.build_lags` | `replications/chan_koop_yu2024_jbes_oisv/legacy/FRED_MD_20vars.csv`, read-only |
| ex07 | `bvar.samplers.factor_fsv`, `bvar.samplers.eq_fsv_load`, `bvar.sv.ksc_ar1_mean`, `bvar.sv.sv_params`, `bvar.sv.init_approx1N`, `bvar.priors.minn`, `bvar.util.build_lags`, `bvar.util.shaded_band` | simulated |
| ex08 | `bvar.samplers.eq_hyb_tvp`, `bvar.sv.ksc_rw_h0`, `bvar.priors.resid_var_allvars_ridge`, `bvar.priors.minnesota_C`, `bvar.priors.vtheta`, `bvar.util.shaded_band` | simulated |
| ex09 | `replications/chan2020_jbes_kronecker/run_ml.m`, which calls the `bvar.ml.*` evaluators | that package's `data_Q.csv` |
| ex10 | `bvar.forecast.iterate`, `bvar.forecast.tables`, `bvar.samplers.eq_gauss`, `bvar.samplers.gig_shrinkage`, `bvar.samplers.nu_psi_ng`, `bvar.priors.minnesota_C`, `bvar.priors.vtheta`, `bvar.priors.resid_var_ar4`, `bvar.sv.ksc_rw_h0`, `bvar.util.build_lags` | `replications/chan2021_ijf_mahp/legacy/macrodata_Q_2018Q4.csv`, read-only |
| ex11 | `bvar.priors.resid_var_ar4`, `bvar.priors.acp_redu`, `bvar.samplers.acp_theta_sig`, `bvar.structural.reduced_form`, `bvar.structural.qr_sign`, `bvar.structural.sign_restrict`, `bvar.structural.sign_assign` | `replications/chan_matthes_yu2026_qe_svarsign/legacy/data/database_2019Q4.csv`, read-only |
| ex12 | the same seven, plus `bvar.structural.irf_redu` — the only example that computes an impulse response — and `bvar.ml.acp`, `bvar.priors.acp_opt_kappa` and `bvar.util.build_lags` for the marginal likelihood, and `bvar.util.shaded_band` for the bands | that package's `data/Uhlig_monthly.csv`, read-only |
| ex13 | `bvar.models.var_sv` with `'draws'`, `bvar.diag.inefficiency_factor`, `bvar.diag.mcse`, `bvar.diag.geweke` | `replications/chan_koop_yu2024_jbes_oisv/legacy/FRED_MD_20vars.csv`, read-only |
| ex14 | `bvar.models.mfvar_csv`, whose body calls `bvar.samplers.missing_var` (which calls `bvar.util.select_obs`), `bvar.samplers.var_coef_joint`, `bvar.sv.csv_armh_block`, `bvar.sv.sv0_params` and `bvar.util.build_lags`; `bvar.util.mm_constraint` for the aggregation and `bvar.util.shaded_band` for the band | simulated |

We note two points about reading these scripts. First, ex01 and ex04 spell out inline what a
core function would otherwise do in one call: the precision-sampler draw in ex01, and the
equation-by-equation coefficient block in ex04. The construction is what these two scripts
teach, and their headers name the packaged version to use in practice. Second, the settings in
ex09 are far smaller than those of the published run, a few hundred draws against 30,000. The
ranking of the three models is still informative, but the values are not those reported in
the paper. The closing lines of that script give the settings for a full-length run.
