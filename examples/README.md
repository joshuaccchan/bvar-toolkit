# Examples

Five self-contained scripts, in reading order. Each one adds `core/` and
`third_party/` to the path by running `setup.m` at the repo root, so you can run
any of them from anywhere:

```matlab
run('examples/ex01_precision_sampler.m')      % or: cd examples; ex01_precision_sampler
```

They are written to be read, not to produce results: every run finishes in
seconds, uses a fixed seed, prints numbers you can check against a simulated
truth, and draws a non-blocking figure. Data are either simulated or read
**read-only** from a `replications/*/legacy/` folder.

| # | Script | What it teaches | Book chapter |
|---|---|---|---|
| 1 | `ex01_precision_sampler.m` | The Chan-Jeliazkov precision sampler on a local-level model: build a banded precision matrix, draw the whole state path in one block, and see the O(T) cost against a dense O(T^3) factorization. Ends with `bvt.util.surform` for a multi-dimensional state. | 9. Linear Gaussian State Space Models |
| 2 | `ex02_sv_ksc.m` | Univariate stochastic volatility via the Kim-Shephard-Chib 7-component auxiliary mixture (`bvt.sv.ksc_rw_h0`): recover a simulated log-volatility path, and see why the `log(y^2 + c)` offset makes the sampler sensitive to the units your data are in. | 10. Stochastic Volatility Models (the auxiliary mixture itself: 4. Mixture Models) |
| 3 | `ex03_minnesota_bvar.m` | A small BVAR on real quarterly US data: `bvt.util.build_lags`, `bvt.priors.resid_var_ar4`, and the two prior constructors `bvt.priors.minn` and `bvt.priors.niw` side by side; analytic posterior draws and a scored one-step-ahead forecast. | 12. Vector Autoregressions |
| 4 | `ex04_bvar_sv_blocks.m` | Equation-by-equation estimation of a **reduced-form** BVAR with stochastic volatility, following Chan (2023, JoE): the coefficient matrix `A` drawn one equation at a time (equation `ii` is *column* `ii` of `A`, since `A` is `k x n` with `beta = vec(A)`), then `bvt.samplers.alp_tri_cs` for the triangular impact matrix, `bvt.sv.ksc_ar1_mean` for the volatility paths and `bvt.sv.sv_params` for their AR(1) parameters, with `bvt.priors.minn`/`bvt.priors.impact_B0` supplying the priors. The point of the one-equation-at-a-time draw is that a lower-triangular `B0` truncates equation `ii`'s stacked system to `(n-ii+1)*T` rows - the exact saving that makes reduced-form estimation feasible at scale, and the one index range that separates this block from `bvt.samplers.eq_svar_oi`. Illustrates (does not bitwise reproduce) `replications/chan2023_joe_mlvarsv/legacy/VAR_ARSV_redu.m`; checked against a simulated truth. | 14. Large VARs with Stochastic Volatility |
| 5 | `ex05_marginal_likelihood.m` | Marginal likelihoods by Chib's method through `replications/chan2020_jbes_kronecker/run_ml.m` and `bvt.ml.*`: the three pieces of the identity printed separately, then Gaussian vs Student-t vs common-SV errors ranked by log Bayes factor. | 5. Bayesian Model Comparison |

Examples 1 and 2 are the computational foundation; 3 and 4 are the modelling
layer built on top of it; 5 is model comparison. After 4, read
`replications/chan2023_joe_mlvarsv/legacy/VAR_ARSV_redu.m` for the published
version of that sampler (that package is not functionized yet, so it has no
`run_all.m`) and `replications/chan2021_ijf_mahp/run_all.m` for the
structural-form flagship with the hierarchical shrinkage block switched on; then
`tests/variant_map.md` for what each core function canonicalizes and how that
was verified.

The chapter column refers to *Bayesian Macroeconometrics: Methods and
Applications* (Chapman & Hall/CRC, forthcoming), which develops the theory
these scripts implement. The
[sample](https://joshuachan.org/papers/BayesMacroBook_sample.pdf) contains the
preface and Chapters 1-4; the book's own code, in MATLAB, R and Python, is at
[joshuaccchan/bayesian-macroeconometrics](https://github.com/joshuaccchan/bayesian-macroeconometrics).

Requirements: MATLAB with the Statistics and Machine Learning Toolbox
(`gamrnd`, `iwishrnd`, `normpdf`, `quantile`). The models example 5 runs need
nothing further; the MA models in that package additionally use `fminunc` from
the Optimization Toolbox.
