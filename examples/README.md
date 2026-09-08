# Examples

Six short scripts that build up from the one computational idea the toolkit rests on to
model comparison and structural identification. They are teaching material, not
replications: each works on data small enough that you can check the numbers against the
truth, and each prints its reasoning as it goes. For reproducing a published table, use
`replications/<paper>/` instead.

Every script puts the toolkit on the path itself, so any of them runs from a clean session:

```matlab
cd examples
ex01_precision_sampler
```

Read them in order. Each one uses what the last one built.

| | Script | What it teaches | Book | Runs in |
|---|---|---|---|---|
| 1 | `ex01_precision_sampler.m` | Drawing an entire state path in one block, with no filtering recursion — the Chan–Jeliazkov (2009) precision sampler. Derives it from the banded precision matrix, checks the draws against the Kalman smoother, and shows the sparse structure that makes it linear in *T*. | Ch. 9 | 11 s |
| 2 | `ex02_sv_ksc.m` | Stochastic volatility by the Kim–Shephard–Chib auxiliary mixture: how squaring and logging the data turns a nonlinear model into the linear Gaussian one ex01 already solves, and what the seven-component mixture is for. | Ch. 10 | 5 s |
| 3 | `ex03_minnesota_bvar.m` | A small BVAR end to end with a Minnesota / natural-conjugate prior — how the prior is built, what the shrinkage hyperparameter does, and why the natural-conjugate restriction yields an analytic posterior. No MCMC: samples are directly drawn from the posterior. | Ch. 12 | 3 s |
| 4 | `ex04_bvar_sv_blocks.m` | Assembling a reduced-form BVAR with stochastic volatility from core blocks, drawn equation by equation — the sampler of `VAR_ARSV_redu.m` from Chan (2023, JoE), on simulated data with the truth known. | Ch. 14 | 4 s |
| 5 | `ex05_marginal_likelihood.m` | Marginal likelihoods and model comparison: the three pieces of Chib's identity, why the posterior ordinate is *subtracted* (the Ockham factor), and a three-model comparison from Chan (2020, JBES). | Ch. 5 | 1 s |
| 6 | `ex06_sign_restrictions.m` | Identifying a structural VAR by sign restrictions, and the cost of the search: two acceptance rules run over the same posterior draws and the same rotations, one requiring each shock in its own column and one searching over assignments. | &mdash; | 5 s |

Timings are from one warm R2025b session on a desktop machine; treat them as orders of
magnitude. The first four each draw figures as well as printing. The chapter column
refers to *Bayesian Macroeconometrics: Methods and Applications* (Chan, Chapman &
Hall/CRC, forthcoming), whose [sample chapters](https://joshuachan.org/papers/BayesMacroBook_sample.pdf)
and [code repository](https://github.com/joshuaccchan/bayesian-macroeconometrics) are
online: 9 Linear Gaussian State Space Models, 10 Stochastic Volatility Models, 12
Vector Autoregressions, 14 Large VARs with Stochastic Volatility, 5 Bayesian Model
Comparison. The auxiliary mixture ex02 uses is developed in Chapter 4, Mixture Models.
Each script repeats its chapter in the header. ex06 has no entry because the book has
no chapter on structural VARs.

## What each one exercises

Useful if you are looking for a worked call of a particular core function.

| Script | Core functions called | Data |
|---|---|---|
| ex01 | `bvar.util.surform` — the precision sampler itself is written out line by line, since deriving it is the point | simulated |
| ex02 | `bvar.sv.ksc_rw_h0` | simulated |
| ex03 | `bvar.priors.minn`, `bvar.priors.niw`, `bvar.priors.resid_var_ar4`, `bvar.util.build_lags` | `replications/chan2020_jbes_kronecker/legacy/data_Q.csv`, read-only |
| ex04 | `bvar.priors.minn`, `bvar.priors.impact_B0`, `bvar.samplers.alp_tri_cs`, `bvar.sv.ksc_ar1_mean`, `bvar.sv.sv_params`, `bvar.sv.init_approx1N`, `bvar.util.build_lags`, `bvar.util.vec` | simulated |
| ex05 | `replications/chan2020_jbes_kronecker/run_ml.m`, which calls the `bvar.ml.*` evaluators | that package's `data_Q.csv` |
| ex06 | `bvar.priors.resid_var_ar4`, `bvar.priors.acp_redu`, `bvar.samplers.acp_theta_sig`, `bvar.structural.reduced_form`, `bvar.structural.qr_sign`, `bvar.structural.sign_restrict`, `bvar.structural.sign_assign` | `replications/chan_matthes_yu2026_qe_svarsign/legacy/data/database_2019Q4.csv`, read-only |

We note two points about reading these scripts. First, ex01 and ex04 write out by hand what a
core function would otherwise do in one call: the precision-sampler draw in ex01, and the
equation-by-equation coefficient block in ex04. The construction is what these two scripts
teach, and their headers name the packaged version to use in practice. Second, the settings in
ex05 are far smaller than those of the published run, a few hundred draws against 30,000. The
ranking of the three models is still informative, but the values are not those reported in
the paper. The closing lines of that script give the settings for a full-length run.
