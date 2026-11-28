# Citing bvar-toolkit

The code here implements methods from the papers below. If you use it in your work, please
cite the paper behind each method you use. The table maps each replication package, library
function and example to that paper, and the BibTeX for every entry follows. A study usually
needs more than one row: one for the model it estimates and one for each sampler inside it.
To credit the software as well, cite the toolkit record at the end. [`RELATED.md`](RELATED.md)
lists implementations of the same models in R and Python.

## What to cite

| Method | In this repository | Cite |
|---|---|---|
| The precision sampler for state space models | `replications/chan_jeliazkov2009_statespace`, `bvar.util.diffmat`; ex01 | [Chan and Jeliazkov (2009)](#chan-and-jeliazkov-2009) |
| Stochastic volatility by the auxiliary mixture sampler | `bvar.sv.ksc_rw_h0`, `ksc_rw_diffuse`, `ksc_ar1_mean`; ex02 | [Kim, Shephard and Chib (1998)](#kim-shephard-and-chib-1998) for the mixture approximation; [Chan and Jeliazkov (2009)](#chan-and-jeliazkov-2009) for drawing the log-volatility path in one block from its banded precision |
| The natural conjugate prior and its closed-form marginal likelihood | `bvar.ml.ncp`, `bvar.models.var_ncp`; guides/natural_conjugate_ml | [Chan (2027)](#chan-2027) |
| The Minnesota prior and its scaling | `bvar.priors.minn`, `minnesota_C`, `resid_var_ar4`, `vtheta`; ex03 | [Doan, Litterman and Sims (1984)](#doan-litterman-and-sims-1984); [Litterman (1986)](#litterman-1986) |
| Minnesota-type adaptive hierarchical priors | `replications/chan2021_ijf_mahp`, `bvar.samplers.gig_shrinkage`, `nu_psi_ng`, `eq_gauss`; ex10 | [Chan (2021)](#chan-2021) |
| The asymmetric conjugate prior, its closed-form marginal likelihood and optimal hyperparameters | `replications/chan2022_qe_acp`, `replications/chan2019wp_acp`, `bvar.priors.acp_stru`, `acp_redu`, `acp_opt_kappa`, `bvar.samplers.acp_theta_sig`, `bvar.ml.acp`; ex11, ex12 | [Chan (2022)](#chan-2022) |
| Non-Gaussian, heteroscedastic or serially dependent errors, and their marginal likelihoods | `replications/chan2020_jbes_kronecker`, `bvar.ml.kron_*`, `intlike_*`, `llike_*`; ex09 | [Chan (2020a)](#chan-2020a) |
| Common stochastic volatility | `bvar.models.var_csv`, `bvar.sv.csv_armh`, `csv_armh_block`; ex05, ex14, tutorials/forecasting, tutorials/mixed_frequency | [Carriero, Clark and Marcellino (2016)](#carriero-clark-and-marcellino-2016) for the model; [Chan (2020a)](#chan-2020a) for the estimation algorithm, whose volatility step is the accept-reject Metropolis-Hastings of [Chan (2017)](#chan-2017) |
| Cholesky stochastic volatility | `bvar.samplers.eq_tri_cs`, `alp_tri_cs`, `bvar.models.var_sv` with `'model'` set to `'CS'`; ex04, ex06 | [Cogley and Sargent (2005)](#cogley-and-sargent-2005) for the model; [Carriero, Chan, Clark and Marcellino (2022)](#carriero-chan-clark-and-marcellino-2022) for the corrected equation-by-equation algorithm in `eq_tri_cs`, and [Carriero, Clark and Marcellino (2019)](#carriero-clark-and-marcellino-2019) for the original |
| Choosing a stochastic volatility specification by marginal likelihood; factor stochastic volatility | `replications/chan2023_joe_mlvarsv`, `bvar.ml.mlvarsv_*`, `bvar.samplers.factor_fsv`, `eq_fsv_load`, `eq_var_redu_tri`, `bvar.priors.impact_B0`, `bvar.sv.init_approx1N`, `svo_outlier`; ex04, ex07 | [Chan (2023a)](#chan-2023a) |
| Order-invariant stochastic volatility | `replications/chan_koop_yu2024_jbes_oisv`, `bvar.models.var_sv`, `bvar.samplers.eq_svar_oi`, `eq_var_oi`, `horseshoe_kappa_psi`, `bvar.structural.b0_row_sampler`, `construct_Sigt`, `bvar.sv.sv_params`, `sv0_params`; ex06 | [Chan, Koop and Yu (2024)](#chan-koop-and-yu-2024) for the model; [Waggoner and Zha (2003)](#waggoner-and-zha-2003) and [Villani (2009)](#villani-2009) for the row-by-row draw of the impact matrix in `b0_row_sampler` |
| Hybrid time-varying parameter VARs | `replications/chan2023_jbes_hybtvp`, `bvar.samplers.eq_hyb_tvp`, `bvar.priors.resid_var_allvars_ridge`; ex08 | [Chan (2023b)](#chan-2023b) |
| Comparing time-varying parameter VARs by marginal likelihood and DIC | `replications/chan_eisenstat2018_jae_mltvpsv`, `bvar.ml.mltvpsv_*`, `intlike_tvp`, `intlike_tvpsv`, `intlike_cvarsv`, `intlike_rsvar`, `dirifit`, `ldiripdf`, `bvar.util.dirirnd` | [Chan and Eisenstat (2018)](#chan-and-eisenstat-2018) |
| Forecasting with large BVARs across priors and volatility models | `replications/chan2020_springer_largebvar`, `bvar.forecast.iterate`, `tables`, `realtime_loaddata`, `bvar.priors.niw`, `bvar.sv.nu_studentt` | [Chan (2020b)](#chan-2020b) |
| Sign restrictions | `bvar.structural.qr_sign`, `sign_restrict`, `irf_redu`, `reduced_form`; ex11, ex12 | [Rubio-Ramírez, Waggoner and Zha (2010)](#rubio-ramírez-waggoner-and-zha-2010) |
| Sign and ranking restrictions in large structural VARs | `replications/chan_matthes_yu2026_qe_svarsign`, `bvar.structural.sign_assign`; ex11, ex12, tutorials/sign_restrictions | [Chan, Matthes and Yu (2026)](#chan-matthes-and-yu-2026) |
| Sign restrictions over several horizons | ex12 | [Uhlig (2005)](#uhlig-2005) |
| Block exogeneity for a small open economy, estimated block by block | `bvar.models.var_soe`, `bvar.util.soe_reduced_form`; guides/soe_var, guides/soe_exogeneity_test, guides/soe_foreign_shock, guides/soe_conditional_forecast | [Cushman and Zha (1997)](#cushman-and-zha-1997) for the model; [Zha (1999)](#zha-1999) for the estimation block by block |
| A shock to the US excess bond premium, identified recursively | guides/soe_foreign_shock | [Gilchrist and Zakrajšek (2012)](#gilchrist-and-zakrajšek-2012) |
| Forecasts conditional on a path for some of the variables | guides/soe_conditional_forecast | [Waggoner and Zha (1999)](#waggoner-and-zha-1999) |
| Missing data and mixed frequencies: the missing values of a VAR drawn in one block, under restrictions that tie them to observed quarterly values | `bvar.models.mfvar_csv`, `bvar.samplers.missing_var`, `var_coef_joint`, `bvar.util.select_obs`, `mm_constraint`, `dlog_gaps`; ex14, tutorials/mixed_frequency | [Chan, Poon and Zhu (2023)](#chan-poon-and-zhu-2023); [Schorfheide and Song (2015)](#schorfheide-and-song-2015) for the mixed-frequency VAR; [Mariano and Murasawa (2003)](#mariano-and-murasawa-2003) for the aggregation of monthly growth rates into quarterly ones; [Chan (2020a)](#chan-2020a) for the draws of the common stochastic volatility and its parameters |
| Sensitivity of forecasts to prior hyperparameters | `replications/cjz2018_ad_var` | [Chan, Jacobi and Zhu (2019)](#chan-jacobi-and-zhu-2019) |
| Hyperparameter selection by automatic differentiation | `replications/cjz2019_ad_opthyper`, `bvar.priors.niw(..., 'opthyper_ncp')` | [Chan, Jacobi and Zhu (2020)](#chan-jacobi-and-zhu-2020) |
| Prior robustness of marginal likelihoods | `replications/cjz2021_jae_ad_ml` | [Chan, Jacobi and Zhu (2022)](#chan-jacobi-and-zhu-2022) |
| Diagnostics for MCMC output: inefficiency factors, Monte Carlo standard errors and the convergence diagnostic | `bvar.diag.inefficiency_factor`, `mcse`, `geweke`, `specvar0`; ex13 | [Geweke (1992)](#geweke-1992); [Newey and West (1987)](#newey-and-west-1987) for the long-run variance and [Newey and West (1994)](#newey-and-west-1994) for the default lag of `geweke` |
| FRED-MD and FRED-QD data | ex06, `replications/chan_koop_yu2024_jbes_oisv`, tutorials/mixed_frequency | [McCracken and Ng (2016)](#mccracken-and-ng-2016); [McCracken and Ng (2021)](#mccracken-and-ng-2021) |

## References

### Brave, Butters and Kelley (2019)

Brave, S. A., Butters, R. A. and Kelley, D. (2019). A New 'Big Data' Index of U.S. Economic
Activity. *Economic Perspectives, Federal Reserve Bank of Chicago* 43(1).
[Journal version](https://doi.org/10.21033/ep-2019-1)

The collapsed dynamic factor model behind the monthly GDP series that the mixed-frequency
tutorial compares its estimates with. The model combines a large panel of monthly indicators with
quarterly real GDP growth.

```bibtex
@article{BBK19,
  author  = {Brave, S. A. and Butters, R. A. and Kelley, D.},
  title   = {A New `Big Data' Index of {U.S.} Economic Activity},
  journal = {Economic Perspectives, Federal Reserve Bank of Chicago},
  year    = {2019},
  volume  = {43},
  number  = {1},
  doi     = {10.21033/ep-2019-1}
}
```

### Carriero, Chan, Clark and Marcellino (2022)

Carriero, A., Chan, J. C. C., Clark, T. E. and Marcellino, M. (2022). Corrigendum to: Large
Bayesian Vector Autoregressions with Stochastic Volatility and Non-Conjugate Priors. *Journal
of Econometrics* 227(2): 506-512.
[Journal version](https://doi.org/10.1016/j.jeconom.2021.11.010) ·
[Working paper](https://joshuachan.org/papers/CCCM.pdf)

Corrects the equation-by-equation algorithm of Carriero, Clark and Marcellino (2019), whose
conditional distributions for the VAR coefficients omitted part of the information, at the
same computational cost.

```bibtex
@article{CCCM22,
  author  = {Carriero, A. and Chan, J. C. C. and Clark, T. E. and Marcellino, M.},
  title   = {Corrigendum to: Large {B}ayesian Vector Autoregressions with Stochastic Volatility and Non-Conjugate Priors},
  journal = {Journal of Econometrics},
  year    = {2022},
  volume  = {227},
  number  = {2},
  pages   = {506--512},
  doi     = {10.1016/j.jeconom.2021.11.010}
}
```

### Carriero, Clark and Marcellino (2016)

Carriero, A., Clark, T. E. and Marcellino, M. (2016). Common Drifting Volatility in Large
Bayesian VARs. *Journal of Business and Economic Statistics* 34(3): 375-390.
[Journal version](https://doi.org/10.1080/07350015.2015.1040116)

The common stochastic volatility model that `bvar.sv.csv_armh` and ex05 estimate.

```bibtex
@article{CCM16,
  author  = {Carriero, A. and Clark, T. E. and Marcellino, M.},
  title   = {Common Drifting Volatility in Large {B}ayesian {VAR}s},
  journal = {Journal of Business and Economic Statistics},
  year    = {2016},
  volume  = {34},
  number  = {3},
  pages   = {375--390},
  doi     = {10.1080/07350015.2015.1040116}
}
```

### Carriero, Clark and Marcellino (2019)

Carriero, A., Clark, T. E. and Marcellino, M. (2019). Large Bayesian Vector Autoregressions
with Stochastic Volatility and Non-Conjugate Priors. *Journal of Econometrics* 212(1): 137-154.
[Journal version](https://doi.org/10.1016/j.jeconom.2019.04.024)

The original equation-by-equation algorithm for large VARs with Cholesky stochastic
volatility, corrected in Carriero, Chan, Clark and Marcellino (2022).

```bibtex
@article{CCM19,
  author  = {Carriero, A. and Clark, T. E. and Marcellino, M.},
  title   = {Large {B}ayesian Vector Autoregressions with Stochastic Volatility and Non-Conjugate Priors},
  journal = {Journal of Econometrics},
  year    = {2019},
  volume  = {212},
  number  = {1},
  pages   = {137--154},
  doi     = {10.1016/j.jeconom.2019.04.024}
}
```

### Chan (2017)

Chan, J. C. C. (2017). The Stochastic Volatility in Mean Model with Time-Varying Parameters: An
Application to Inflation Modeling. *Journal of Business and Economic Statistics* 35(1): 17-28.
[Journal version](https://doi.org/10.1080/07350015.2015.1052459) ·
[Working paper](https://joshuachan.org/papers/SVM.pdf)

A stochastic volatility in mean model whose coefficients evolve over time, estimated by an
accept-reject Metropolis-Hastings step that draws the whole log-volatility path in one block from
a Gaussian proposal at the mode of its conditional density. That step is what `bvar.sv.csv_armh`
adapts to the common volatility factor.

```bibtex
@article{chan17jbes,
  author  = {Chan, J. C. C.},
  title   = {The Stochastic Volatility in Mean Model with Time-Varying Parameters: An Application to Inflation Modeling},
  journal = {Journal of Business and Economic Statistics},
  year    = {2017},
  volume  = {35},
  number  = {1},
  pages   = {17--28},
  doi     = {10.1080/07350015.2015.1052459}
}
```

### Chan (2020a)

Chan, J. C. C. (2020). Large Bayesian VARs: A Flexible Kronecker Error Covariance Structure.
*Journal of Business and Economic Statistics* 38(1): 68-79.
[Journal version](https://doi.org/10.1080/07350015.2018.1451336) ·
[Working paper](https://joshuachan.org/papers/BVAR.pdf) ·
Code: [`replications/chan2020_jbes_kronecker`](replications/chan2020_jbes_kronecker)

A class of large BVARs with non-Gaussian, heteroscedastic and serially dependent innovations,
made tractable by a Kronecker structure of the likelihood and estimated by MCMC. In an
application with 20 macroeconomic variables these outperform the standard BVAR in both
in-sample fit and forecasts.

```bibtex
@article{chan20jbes,
  author  = {Chan, J. C. C.},
  title   = {Large {B}ayesian {VAR}s: A Flexible {K}ronecker Error Covariance Structure},
  journal = {Journal of Business and Economic Statistics},
  year    = {2020},
  volume  = {38},
  number  = {1},
  pages   = {68--79},
  doi     = {10.1080/07350015.2018.1451336}
}
```

### Chan (2020b)

Chan, J. C. C. (2020). Large Bayesian Vector Autoregressions. In: P. Fuleky (Ed.),
*Macroeconomic Forecasting in the Era of Big Data*, 95-125. Springer, Cham.
[Chapter](https://doi.org/10.1007/978-3-030-31150-6_4) ·
[Working paper](https://joshuachan.org/papers/large_BVAR.pdf) ·
Code: [`replications/chan2020_springer_largebvar`](replications/chan2020_springer_largebvar)

A review of shrinkage priors for large BVARs, and of extensions with stochastic volatility and
non-Gaussian and serially correlated errors, illustrated by a real-time forecasting exercise.

```bibtex
@incollection{chan20springer,
  author    = {Chan, J. C. C.},
  title     = {Large {B}ayesian Vector Autoregressions},
  booktitle = {Macroeconomic Forecasting in the Era of Big Data},
  editor    = {Fuleky, P.},
  series    = {Advanced Studies in Theoretical and Applied Econometrics},
  volume    = {52},
  publisher = {Springer},
  address   = {Cham},
  year      = {2020},
  pages     = {95--125},
  doi       = {10.1007/978-3-030-31150-6_4}
}
```

### Chan (2021)

Chan, J. C. C. (2021). Minnesota-Type Adaptive Hierarchical Priors for Large Bayesian VARs.
*International Journal of Forecasting* 37(3): 1212-1226.
[Journal version](https://doi.org/10.1016/j.ijforecast.2021.01.002) ·
[Working paper](https://joshuachan.org/papers/BVAR-MAHP.pdf) ·
Code: [`replications/chan2021_ijf_mahp`](replications/chan2021_ijf_mahp)

A family of shrinkage priors that combines adaptive hierarchical priors, which shrink only
small coefficients strongly toward zero, with the cross-variable shrinkage and lag decay of the
Minnesota prior, together with a fast posterior sampler. In a forecasting exercise the new
priors outperform both parent classes.

```bibtex
@article{chan21ijf,
  author  = {Chan, J. C. C.},
  title   = {{M}innesota-Type Adaptive Hierarchical Priors for Large {B}ayesian {VAR}s},
  journal = {International Journal of Forecasting},
  year    = {2021},
  volume  = {37},
  number  = {3},
  pages   = {1212--1226},
  doi     = {10.1016/j.ijforecast.2021.01.002}
}
```

### Chan (2022)

Chan, J. C. C. (2022). Asymmetric Conjugate Priors for Large Bayesian VARs. *Quantitative
Economics* 13(3): 1145-1169.
[Journal version](https://doi.org/10.3982/QE1381) ·
[Working paper](https://joshuachan.org/papers/BVAR-ACP-R1.pdf) ·
Code: [`replications/chan2022_qe_acp`](replications/chan2022_qe_acp)

A prior that allows cross-variable shrinkage while keeping the analytical results of the
natural conjugate prior, including a closed-form marginal likelihood, with fast posterior
simulation. Illustrated by a 15-variable VAR with sign restrictions identifying five shocks.

```bibtex
@article{chan22qe,
  author  = {Chan, J. C. C.},
  title   = {Asymmetric Conjugate Priors for Large {B}ayesian {VAR}s},
  journal = {Quantitative Economics},
  year    = {2022},
  volume  = {13},
  number  = {3},
  pages   = {1145--1169},
  doi     = {10.3982/QE1381}
}
```

### Chan (2023a)

Chan, J. C. C. (2023). Comparing Stochastic Volatility Specifications for Large Bayesian VARs.
*Journal of Econometrics* 235(2): 1419-1446.
[Journal version](https://doi.org/10.1016/j.jeconom.2022.11.003) ·
[Working paper](https://joshuachan.org/papers/ML-LargeVARSV.pdf) ·
Code: [`replications/chan2023_joe_mlvarsv`](replications/chan2023_joe_mlvarsv)

Marginal likelihood estimators, combining conditional Monte Carlo with adaptive importance
sampling, for choosing among stochastic volatility specifications and shrinkage priors in large
BVARs.

```bibtex
@article{chan23joe,
  author  = {Chan, J. C. C.},
  title   = {Comparing Stochastic Volatility Specifications for Large {B}ayesian {VAR}s},
  journal = {Journal of Econometrics},
  year    = {2023},
  volume  = {235},
  number  = {2},
  pages   = {1419--1446},
  doi     = {10.1016/j.jeconom.2022.11.003}
}
```

### Chan (2023b)

Chan, J. C. C. (2023). Large Hybrid Time-Varying Parameter VARs. *Journal of Business and
Economic Statistics* 41(3): 890-905.
[Journal version](https://doi.org/10.1080/07350015.2022.2080683) ·
[Working paper](https://joshuachan.org/papers/HYB-TVPVAR.pdf) ·
Code: [`replications/chan2023_jbes_hybtvp`](replications/chan2023_jbes_hybtvp)

Hybrid TVP-VARs, in which each equation's coefficients and contemporaneous relations are
either constant or time-varying, with an efficient sparsification method that decides which,
equation by equation.

```bibtex
@article{chan23jbes,
  author  = {Chan, J. C. C.},
  title   = {Large Hybrid Time-Varying Parameter {VAR}s},
  journal = {Journal of Business and Economic Statistics},
  year    = {2023},
  volume  = {41},
  number  = {3},
  pages   = {890--905},
  doi     = {10.1080/07350015.2022.2080683}
}
```

### Chan (2027)

Chan, J. C. C. (2027). *Bayesian Macroeconometrics: Methods and Applications*. Chapman &
Hall/CRC.
[Sample chapters](https://joshuachan.org/papers/BayesMacroBook_sample.pdf) ·
Code: [bayesian-macroeconometrics](https://github.com/joshuaccchan/bayesian-macroeconometrics)

The textbook treatment of the models and samplers in this repository, and the reference for the
how-to guides.

```bibtex
@book{chan-bmar,
  author    = {Chan, Joshua C. C.},
  title     = {Bayesian Macroeconometrics: Methods and Applications},
  publisher = {Chapman \& Hall/CRC},
  year      = {2027}
}
```

### Chan and Eisenstat (2018)

Chan, J. C. C. and Eisenstat, E. (2018). Bayesian Model Comparison for Time-Varying Parameter
VARs with Stochastic Volatility. *Journal of Applied Econometrics* 33(4): 509-532.
[Journal version](https://doi.org/10.1002/jae.2617) ·
[Working paper](https://joshuachan.org/papers/ml-tvpsv.pdf) ·
Code: [`replications/chan_eisenstat2018_jae_mltvpsv`](replications/chan_eisenstat2018_jae_mltvpsv)

Importance sampling estimators of the marginal likelihood and the DIC for TVP-VARs with
stochastic volatility, based on the integrated likelihood.

```bibtex
@article{CE18,
  author  = {Chan, J. C. C. and Eisenstat, E.},
  title   = {{B}ayesian Model Comparison for Time-Varying Parameter {VAR}s with Stochastic Volatility},
  journal = {Journal of Applied Econometrics},
  year    = {2018},
  volume  = {33},
  number  = {4},
  pages   = {509--532},
  doi     = {10.1002/jae.2617}
}
```

### Chan, Jacobi and Zhu (2019)

Chan, J. C. C., Jacobi, L. and Zhu, D. (2019). How Sensitive Are VAR Forecasts to Prior
Hyperparameters? An Automated Sensitivity Analysis. *Advances in Econometrics* 40A: 229-248.
[Chapter](https://doi.org/10.1108/S0731-90532019000040A010) ·
[Working paper](https://joshuachan.org/papers/AD_VAR.pdf) ·
Code: [`replications/cjz2018_ad_var`](replications/cjz2018_ad_var)

A method based on automatic differentiation for computing the sensitivities of point and
interval forecasts from a VAR to any prior hyperparameter.

```bibtex
@incollection{CJZ19,
  author    = {Chan, J. C. C. and Jacobi, L. and Zhu, D.},
  title     = {How Sensitive Are {VAR} Forecasts to Prior Hyperparameters? {A}n Automated Sensitivity Analysis},
  booktitle = {Topics in Identification, Limited Dependent Variables, Partial Observability, Experimentation, and Flexible Modeling: Part A},
  editor    = {Jeliazkov, I. and Tobias, J. L.},
  series    = {Advances in Econometrics},
  volume    = {40A},
  publisher = {Emerald Publishing Limited},
  year      = {2019},
  pages     = {229--248},
  doi       = {10.1108/S0731-90532019000040A010}
}
```

### Chan, Jacobi and Zhu (2020)

Chan, J. C. C., Jacobi, L. and Zhu, D. (2020). Efficient Selection of Hyperparameters in Large
Bayesian VARs Using Automatic Differentiation. *Journal of Forecasting* 39(6): 934-943.
[Journal version](https://doi.org/10.1002/for.2660) ·
[Working paper](https://joshuachan.org/papers/AD_OptHyper.pdf) ·
Code: [`replications/cjz2019_ad_opthyper`](replications/cjz2019_ad_opthyper)

Optimal prior hyperparameters for large BVARs with the natural conjugate prior, computed by
automatic differentiation, much faster than a grid search.

```bibtex
@article{CJZ20,
  author  = {Chan, J. C. C. and Jacobi, L. and Zhu, D.},
  title   = {Efficient Selection of Hyperparameters in Large {B}ayesian {VAR}s Using Automatic Differentiation},
  journal = {Journal of Forecasting},
  year    = {2020},
  volume  = {39},
  number  = {6},
  pages   = {934--943},
  doi     = {10.1002/for.2660}
}
```

### Chan, Jacobi and Zhu (2022)

Chan, J. C. C., Jacobi, L. and Zhu, D. (2022). An Automated Prior Robustness Analysis in
Bayesian Model Comparison. *Journal of Applied Econometrics* 37(3): 583-602.
[Journal version](https://doi.org/10.1002/jae.2889) ·
[Working paper](https://joshuachan.org/papers/AD_ML.pdf) ·
Code: [`replications/cjz2021_jae_ad_ml`](replications/cjz2021_jae_ad_ml)

The sensitivities of simulation-based marginal likelihood estimates to any prior
hyperparameter, computed by automatic differentiation alongside the MCMC algorithm.

```bibtex
@article{CJZ22,
  author  = {Chan, J. C. C. and Jacobi, L. and Zhu, D.},
  title   = {An Automated Prior Robustness Analysis in {B}ayesian Model Comparison},
  journal = {Journal of Applied Econometrics},
  year    = {2022},
  volume  = {37},
  number  = {3},
  pages   = {583--602},
  doi     = {10.1002/jae.2889}
}
```

### Chan and Jeliazkov (2009)

Chan, J. C. C. and Jeliazkov, I. (2009). Efficient Simulation and Integrated Likelihood
Estimation in State Space Models. *International Journal of Mathematical Modelling and
Numerical Optimisation* 1(1/2): 101-120.
[Journal version](https://doi.org/10.1504/IJMMNO.2009.030090) ·
[Working paper](https://joshuachan.org/papers/statespace1.pdf) ·
Code: [`replications/chan_jeliazkov2009_statespace`](replications/chan_jeliazkov2009_statespace)

A derivation of the posterior distribution of the states that leads to a modular, scalable
precision-based simulation algorithm for state space models, and a simple way to evaluate the
integrated likelihood.

```bibtex
@article{CJ09,
  author  = {Chan, J. C. C. and Jeliazkov, I.},
  title   = {Efficient Simulation and Integrated Likelihood Estimation in State Space Models},
  journal = {International Journal of Mathematical Modelling and Numerical Optimisation},
  year    = {2009},
  volume  = {1},
  number  = {1/2},
  pages   = {101--120},
  doi     = {10.1504/IJMMNO.2009.030090}
}
```

### Chan, Koop and Yu (2024)

Chan, J. C. C., Koop, G. and Yu, X. (2024). Large Order-Invariant Bayesian VARs with
Stochastic Volatility. *Journal of Business and Economic Statistics* 42(2): 825-837.
[Journal version](https://doi.org/10.1080/07350015.2023.2252039) ·
[Working paper](https://joshuachan.org/papers/OISV.pdf) ·
Code: [`replications/chan_koop_yu2024_jbes_oisv`](replications/chan_koop_yu2024_jbes_oisv)

A multivariate stochastic volatility specification that avoids the lower triangular
parameterization of the error covariance matrix, is identified through the stochastic
volatility and is invariant to the ordering of the variables, with an MCMC algorithm for
estimation and prediction.

```bibtex
@article{CKY24,
  author  = {Chan, J. C. C. and Koop, G. and Yu, X.},
  title   = {Large Order-Invariant {B}ayesian {VAR}s with Stochastic Volatility},
  journal = {Journal of Business and Economic Statistics},
  year    = {2024},
  volume  = {42},
  number  = {2},
  pages   = {825--837},
  doi     = {10.1080/07350015.2023.2252039}
}
```

### Chan, Matthes and Yu (2026)

Chan, J. C. C., Matthes, C. and Yu, X. (2026). Large Structural VARs with Multiple Sign and
Ranking Restrictions. *Quantitative Economics* 17(3): 709-740.
[Journal version](https://doi.org/10.3982/QE2529) ·
[Working paper](https://joshuachan.org/papers/SVAR-sign.pdf) ·
Code: [`replications/chan_matthes_yu2026_qe_svarsign`](replications/chan_matthes_yu2026_qe_svarsign)

A numerically efficient algorithm for large structural VARs identified by many sign and
ranking restrictions on impulse responses, illustrated with a 35-variable VAR in which over 100
restrictions identify 8 shocks.

```bibtex
@article{CMY26,
  author  = {Chan, J. C. C. and Matthes, C. and Yu, X.},
  title   = {Large Structural {VAR}s with Multiple Sign and Ranking Restrictions},
  journal = {Quantitative Economics},
  year    = {2026},
  volume  = {17},
  number  = {3},
  pages   = {709--740},
  doi     = {10.3982/QE2529}
}
```

### Chan, Poon and Zhu (2023)

Chan, J. C. C., Poon, A. and Zhu, D. (2023). High-Dimensional Conditionally Gaussian State
Space Models with Missing Data. *Journal of Econometrics* 236(1): 105468.
[Journal version](https://doi.org/10.1016/j.jeconom.2023.05.005) ·
[Working paper](https://joshuachan.org/papers/BVAR-MF-R1.pdf)

An efficient approach to sampling the missing values of a conditionally Gaussian state space
model in one block, from a conditional distribution whose precision matrix is banded for common
missing data patterns, with restrictions that tie missing high-frequency values to observed
low-frequency ones. Its applications are a weekly mixed-frequency VAR with common stochastic
volatility, which produces weekly GDP estimates, and a dynamic factor model with stochastic
volatility on unbalanced FRED-MD vintages.

```bibtex
@article{CPZ23,
  author  = {Chan, J. C. C. and Poon, A. and Zhu, D.},
  title   = {High-Dimensional Conditionally {G}aussian State Space Models with Missing Data},
  journal = {Journal of Econometrics},
  year    = {2023},
  volume  = {236},
  number  = {1},
  pages   = {105468},
  doi     = {10.1016/j.jeconom.2023.05.005}
}
```

### Cogley and Sargent (2005)

Cogley, T. and Sargent, T. J. (2005). Drifts and Volatilities: Monetary Policies and Outcomes
in the Post WWII US. *Review of Economic Dynamics* 8(2): 262-302.
[Journal version](https://doi.org/10.1016/j.red.2004.10.009)

The VAR with stochastic volatility in the Cholesky parameterization, the benchmark that ex04
and ex06 estimate.

```bibtex
@article{CS05,
  author  = {Cogley, T. and Sargent, T. J.},
  title   = {Drifts and Volatilities: {M}onetary Policies and Outcomes in the Post {WWII} {US}},
  journal = {Review of Economic Dynamics},
  year    = {2005},
  volume  = {8},
  number  = {2},
  pages   = {262--302},
  doi     = {10.1016/j.red.2004.10.009}
}
```

### Cushman and Zha (1997)

Cushman, D. O. and Zha, T. (1997). Identifying Monetary Policy in a Small Open Economy under
Flexible Exchange Rates. *Journal of Monetary Economics* 39(3): 433-448.
[Journal version](https://doi.org/10.1016/S0304-3932(97)00029-9)

A structural VAR for Canada in which the Canadian variables do not enter the equations of the US
and world variables, within the period or with a lag: block exogeneity.

```bibtex
@article{CZ97,
  author  = {Cushman, D. O. and Zha, T.},
  title   = {Identifying Monetary Policy in a Small Open Economy under Flexible Exchange Rates},
  journal = {Journal of Monetary Economics},
  year    = {1997},
  volume  = {39},
  number  = {3},
  pages   = {433--448},
  doi     = {10.1016/S0304-3932(97)00029-9}
}
```

### Doan, Litterman and Sims (1984)

Doan, T., Litterman, R. B. and Sims, C. A. (1984). Forecasting and Conditional Projection Using
Realistic Prior Distributions. *Econometric Reviews* 3(1): 1-100.
[Journal version](https://doi.org/10.1080/07474938408800053)

With Litterman (1986), the origin of the Minnesota prior.

```bibtex
@article{DLS84,
  author  = {Doan, T. and Litterman, R. B. and Sims, C. A.},
  title   = {Forecasting and Conditional Projection Using Realistic Prior Distributions},
  journal = {Econometric Reviews},
  year    = {1984},
  volume  = {3},
  number  = {1},
  pages   = {1--100},
  doi     = {10.1080/07474938408800053}
}
```

### Geweke (1992)

Geweke, J. (1992). Evaluating the Accuracy of Sampling-Based Approaches to the Calculation of
Posterior Moments. In J. M. Bernardo, J. O. Berger, A. P. Dawid and A. F. M. Smith (Eds),
*Bayesian Statistics 4*, 169-193. Oxford University Press.
[Publisher version](https://doi.org/10.1093/oso/9780198522669.003.0010) ·
[Working paper](https://doi.org/10.21034/sr.148)

Numerical standard errors and relative numerical efficiency from the spectral density at
frequency zero, and the convergence diagnostic that compares the means of an early and a late
segment of a chain.

```bibtex
@incollection{Geweke92,
  author    = {Geweke, J.},
  title     = {Evaluating the Accuracy of Sampling-Based Approaches to the Calculation of Posterior Moments},
  booktitle = {{B}ayesian Statistics 4},
  editor    = {Bernardo, J. M. and Berger, J. O. and Dawid, A. P. and Smith, A. F. M.},
  publisher = {Oxford University Press},
  year      = {1992},
  pages     = {169--193},
  doi       = {10.1093/oso/9780198522669.003.0010}
}
```

### Gilchrist and Zakrajšek (2012)

Gilchrist, S. and Zakrajšek, E. (2012). Credit Spreads and Business Cycle Fluctuations. *American
Economic Review* 102(4): 1692-1720.
[Journal version](https://doi.org/10.1257/aer.102.4.1692)

The excess bond premium: the part of a corporate bond spread that expected defaults do not
explain. A shock to it is identified recursively, with the premium after activity and inflation
and before the interest rates.

```bibtex
@article{GZ12,
  author  = {Gilchrist, S. and Zakraj{\v{s}}ek, E.},
  title   = {Credit Spreads and Business Cycle Fluctuations},
  journal = {American Economic Review},
  year    = {2012},
  volume  = {102},
  number  = {4},
  pages   = {1692--1720},
  doi     = {10.1257/aer.102.4.1692}
}
```

### Kim, Shephard and Chib (1998)

Kim, S., Shephard, N. and Chib, S. (1998). Stochastic Volatility: Likelihood Inference and
Comparison with ARCH Models. *Review of Economic Studies* 65(3): 361-393.
[Journal version](https://doi.org/10.1111/1467-937X.00050)

The seven-component normal mixture approximation that the auxiliary mixture samplers in
`bvar.sv` use.

```bibtex
@article{KSC98,
  author  = {Kim, S. and Shephard, N. and Chib, S.},
  title   = {Stochastic Volatility: Likelihood Inference and Comparison with {ARCH} Models},
  journal = {Review of Economic Studies},
  year    = {1998},
  volume  = {65},
  number  = {3},
  pages   = {361--393},
  doi     = {10.1111/1467-937X.00050}
}
```

### Litterman (1986)

Litterman, R. B. (1986). Forecasting with Bayesian Vector Autoregressions: Five Years of
Experience. *Journal of Business and Economic Statistics* 4(1): 25-38.
[Journal version](https://doi.org/10.1080/07350015.1986.10509491)

With Doan, Litterman and Sims (1984), the origin of the Minnesota prior.

```bibtex
@article{litterman86,
  author  = {Litterman, R. B.},
  title   = {Forecasting with {B}ayesian Vector Autoregressions: Five Years of Experience},
  journal = {Journal of Business and Economic Statistics},
  year    = {1986},
  volume  = {4},
  number  = {1},
  pages   = {25--38},
  doi     = {10.1080/07350015.1986.10509491}
}
```

### Mariano and Murasawa (2003)

Mariano, R. S. and Murasawa, Y. (2003). A New Coincident Index of Business Cycles Based on
Monthly and Quarterly Series. *Journal of Applied Econometrics* 18(4): 427-443.
[Journal version](https://doi.org/10.1002/jae.695)

A coincident index of business cycles from monthly and quarterly series, with the log-linear
approximation that ties a quarterly growth rate to the monthly growth rates of the same series.

```bibtex
@article{MM03,
  author  = {Mariano, R. S. and Murasawa, Y.},
  title   = {A New Coincident Index of Business Cycles Based on Monthly and Quarterly Series},
  journal = {Journal of Applied Econometrics},
  year    = {2003},
  volume  = {18},
  number  = {4},
  pages   = {427--443},
  doi     = {10.1002/jae.695}
}
```

### McCracken and Ng (2016)

McCracken, M. W. and Ng, S. (2016). FRED-MD: A Monthly Database for Macroeconomic Research.
*Journal of Business and Economic Statistics* 34(4): 574-589.
[Journal version](https://doi.org/10.1080/07350015.2015.1086655)

The monthly database that ex06, the order-invariant package and the mixed-frequency tutorial
draw their data from.

```bibtex
@article{MN16,
  author  = {McCracken, M. W. and Ng, S.},
  title   = {{FRED-MD}: A Monthly Database for Macroeconomic Research},
  journal = {Journal of Business and Economic Statistics},
  year    = {2016},
  volume  = {34},
  number  = {4},
  pages   = {574--589},
  doi     = {10.1080/07350015.2015.1086655}
}
```

### McCracken and Ng (2021)

McCracken, M. W. and Ng, S. (2021). FRED-QD: A Quarterly Database for Macroeconomic Research.
*Federal Reserve Bank of St. Louis Review* 103(1): 1-44.
[Journal version](https://doi.org/10.20955/r.103.1-44)

The quarterly database that the mixed-frequency tutorial takes real GDP and real private fixed
investment from.

```bibtex
@article{MN21,
  author  = {McCracken, M. W. and Ng, S.},
  title   = {{FRED-QD}: A Quarterly Database for Macroeconomic Research},
  journal = {Federal Reserve Bank of St. Louis Review},
  year    = {2021},
  volume  = {103},
  number  = {1},
  pages   = {1--44},
  doi     = {10.20955/r.103.1-44}
}
```

### Newey and West (1987)

Newey, W. K. and West, K. D. (1987). A Simple, Positive Semi-Definite, Heteroskedasticity and
Autocorrelation Consistent Covariance Matrix. *Econometrica* 55(3): 703-708.
[Journal version](https://doi.org/10.2307/1913610)

The Bartlett-window estimator of the long-run variance that `bvar.diag.specvar0` computes.

```bibtex
@article{NW87,
  author  = {Newey, W. K. and West, K. D.},
  title   = {A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix},
  journal = {Econometrica},
  year    = {1987},
  volume  = {55},
  number  = {3},
  pages   = {703--708},
  doi     = {10.2307/1913610}
}
```

### Newey and West (1994)

Newey, W. K. and West, K. D. (1994). Automatic Lag Selection in Covariance Matrix Estimation.
*Review of Economic Studies* 61(4): 631-653.
[Journal version](https://doi.org/10.2307/2297912)

The rule of thumb for the truncation lag that `bvar.diag.geweke` uses by default.

```bibtex
@article{NW94,
  author  = {Newey, W. K. and West, K. D.},
  title   = {Automatic Lag Selection in Covariance Matrix Estimation},
  journal = {Review of Economic Studies},
  year    = {1994},
  volume  = {61},
  number  = {4},
  pages   = {631--653},
  doi     = {10.2307/2297912}
}
```

### Rubio-Ramírez, Waggoner and Zha (2010)

Rubio-Ramírez, J. F., Waggoner, D. F. and Zha, T. (2010). Structural Vector Autoregressions:
Theory of Identification and Algorithms for Inference. *Review of Economic Studies* 77(2):
665-696.
[Journal version](https://doi.org/10.1111/j.1467-937X.2009.00578.x)

The accept-reject algorithm for sign restrictions: draw rotations uniformly and keep those
whose impact responses have the required signs.

```bibtex
@article{RWZ10,
  author  = {Rubio-Ram{\'\i}rez, J. F. and Waggoner, D. F. and Zha, T.},
  title   = {Structural Vector Autoregressions: Theory of Identification and Algorithms for Inference},
  journal = {Review of Economic Studies},
  year    = {2010},
  volume  = {77},
  number  = {2},
  pages   = {665--696},
  doi     = {10.1111/j.1467-937X.2009.00578.x}
}
```

### Schorfheide and Song (2015)

Schorfheide, F. and Song, D. (2015). Real-Time Forecasting with a Mixed-Frequency VAR. *Journal of
Business and Economic Statistics* 33(3): 366-380.
[Journal version](https://doi.org/10.1080/07350015.2014.954707)

A VAR for monthly and quarterly series, specified at the monthly frequency and cast in state space
form, in which the monthly values of the quarterly series are unobserved, estimated by Bayesian
methods under a Minnesota-style prior, with real-time forecasts compared with those of a quarterly
VAR and of MIDAS regressions.

```bibtex
@article{SS15,
  author  = {Schorfheide, F. and Song, D.},
  title   = {Real-Time Forecasting with a Mixed-Frequency {VAR}},
  journal = {Journal of Business and Economic Statistics},
  year    = {2015},
  volume  = {33},
  number  = {3},
  pages   = {366--380},
  doi     = {10.1080/07350015.2014.954707}
}
```

### Uhlig (2005)

Uhlig, H. (2005). What Are the Effects of Monetary Policy on Output? Results from an Agnostic
Identification Procedure. *Journal of Monetary Economics* 52(2): 381-419.
[Journal version](https://doi.org/10.1016/j.jmoneco.2004.05.007)

Identification of a monetary policy shock by sign restrictions over several horizons, the
application of ex12.

```bibtex
@article{uhlig05,
  author  = {Uhlig, H.},
  title   = {What Are the Effects of Monetary Policy on Output? {R}esults from an Agnostic Identification Procedure},
  journal = {Journal of Monetary Economics},
  year    = {2005},
  volume  = {52},
  number  = {2},
  pages   = {381--419},
  doi     = {10.1016/j.jmoneco.2004.05.007}
}
```

### Villani (2009)

Villani, M. (2009). Steady-State Priors for Vector Autoregressions. *Journal of Applied
Econometrics* 24(4): 630-650.
[Journal version](https://doi.org/10.1002/jae.1065)

The extension of the Waggoner and Zha (2003) algorithm to priors with nonzero means, and the
two-component normal approximation to the absolute-normal density, both used by
`bvar.structural.b0_row_sampler`.

```bibtex
@article{villani09,
  author  = {Villani, M.},
  title   = {Steady-State Priors for Vector Autoregressions},
  journal = {Journal of Applied Econometrics},
  year    = {2009},
  volume  = {24},
  number  = {4},
  pages   = {630--650},
  doi     = {10.1002/jae.1065}
}
```

### Waggoner and Zha (1999)

Waggoner, D. F. and Zha, T. (1999). Conditional Forecasts in Dynamic Multivariate Models.
*Review of Economics and Statistics* 81(4): 639-651.
[Journal version](https://doi.org/10.1162/003465399558508)

The distribution of forecasts conditional on a path for some of the variables of a VAR. The
parameters are drawn given both the data and the path. guides/soe_conditional_forecast computes
this distribution under block exogeneity, given a path for all the foreign variables, where draws
given the data alone suffice.

```bibtex
@article{WZ99,
  author  = {Waggoner, D. F. and Zha, T.},
  title   = {Conditional Forecasts in Dynamic Multivariate Models},
  journal = {Review of Economics and Statistics},
  year    = {1999},
  volume  = {81},
  number  = {4},
  pages   = {639--651},
  doi     = {10.1162/003465399558508}
}
```

### Waggoner and Zha (2003)

Waggoner, D. F. and Zha, T. (2003). A Gibbs Sampler for Structural Vector Autoregressions.
*Journal of Economic Dynamics and Control* 28(2): 349-366.
[Journal version](https://doi.org/10.1016/S0165-1889(02)00168-9)

The algorithm that draws the rows of an unrestricted impact matrix one at a time, which
`bvar.structural.b0_row_sampler` uses for the order-invariant model.

```bibtex
@article{WZ03,
  author  = {Waggoner, D. F. and Zha, T.},
  title   = {A {G}ibbs Sampler for Structural Vector Autoregressions},
  journal = {Journal of Economic Dynamics and Control},
  year    = {2003},
  volume  = {28},
  number  = {2},
  pages   = {349--366},
  doi     = {10.1016/S0165-1889(02)00168-9}
}
```

### Zha (1999)

Zha, T. (1999). Block Recursion and Structural Vector Autoregressions. *Journal of Econometrics*
90(2): 291-316.
[Journal version](https://doi.org/10.1016/S0304-4076(98)00045-1)

Inference for structural VARs whose contemporaneous coefficient matrix is block recursive, with
restrictions on the lags such as block exogeneity: the likelihood splits into blocks, and the
posterior is drawn block by block.

```bibtex
@article{Zha99,
  author  = {Zha, T.},
  title   = {Block Recursion and Structural Vector Autoregressions},
  journal = {Journal of Econometrics},
  year    = {1999},
  volume  = {90},
  number  = {2},
  pages   = {291--316},
  doi     = {10.1016/S0304-4076(98)00045-1}
}
```

## The toolkit

Chan, J. C. C. (2026). *bvar-toolkit: MATLAB code for large Bayesian VARs*. Zenodo.
https://doi.org/10.5281/zenodo.22804032

The DOI resolves to the latest release. `CITATION.cff` holds the same record in machine-readable
form, which GitHub's "Cite this repository" button reads.

```bibtex
@misc{chan26bvartoolkit,
  author       = {Chan, J. C. C.},
  title        = {bvar-toolkit: {MATLAB} Code for Large {B}ayesian {VAR}s},
  year         = {2026},
  howpublished = {Zenodo},
  doi          = {10.5281/zenodo.22804032},
  url          = {https://github.com/joshuaccchan/bvar-toolkit}
}
```
