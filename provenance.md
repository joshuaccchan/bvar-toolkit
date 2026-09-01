# Provenance

Verbatim imports from joshuachan.org, one commit and one `as-published/<slug>` tag per package.
Import date: 2026-09-01. `.gitattributes` sets `* -text` so git stores every byte unmodified.

| Package (slug) | Original zip | MD5 of zip | Zip date | Files | Paper |
|---|---|---|---|---|---|
| `chan_jeliazkov2009_statespace` | [sp_code.zip](https://joshuachan.org/code/sp_code.zip) | `be21a9ec76126bb5af36c6df0a5938ee` | 2015-01-27 | 9 | Chan, J.C.C. and Jeliazkov, I. (2009). "Efficient Simulation and Integrated Likelihood Estimation in State Space Models," International Journal of Mathematical Modelling and Numerical Optimisation, 1, 101-120. |
| `cjz2018_ad_var` | [AD_VAR_code.zip](https://joshuachan.org/code/AD_VAR_code.zip) | `bca355e95f95ec3fdf235efae797aebd` | 2018-05-30 | 18 | Chan, J.C.C., Jacobi, L. and Zhu, D. (2019). "How Sensitive Are VAR Forecasts to Prior Hyperparameters? An Automated Sensitivity Analysis," Advances in Econometrics, 40A: 229-248. (Code headers cite CAMA Working Paper 25/2018.) |
| `cjz2019_ad_opthyper` | [AD_OptHyper_code.zip](https://joshuachan.org/code/AD_OptHyper_code.zip) | `de1b9d3451c42f664ed1b6483df0d964` | 2019-06-28 | 47 | Chan, J.C.C., Jacobi, L. and Zhu, D. (2020). "Efficient Selection of Hyperparameters in Large Bayesian VARs Using Automatic Differentiation," Journal of Forecasting, 39(6): 934-943. (Code headers cite CAMA Working Paper 46/2019.) |
| `chan2019wp_acp` | [BVAR_ACP_code.zip](https://joshuachan.org/code/BVAR_ACP_code.zip) | `8a35858d5b34843f82c53147feb97398` | 2019-10-12 | 10 | Chan, J.C.C. (2019). "Asymmetric Conjugate Priors for Large Bayesian VARs," working-paper version with forecasting application; published as Chan (2022), Quantitative Economics, 13(3): 1145-1169. |
| `chan2020_springer_largebvar` | [large_BVAR_code.zip](https://joshuachan.org/code/large_BVAR_code.zip) | `c3146389152b3fd9206312631a108ab6` | 2019-12-05 | 39 | Chan, J.C.C. (2020). "Large Bayesian Vector Autoregressions." In: P. Fuleky (Ed.), Macroeconomic Forecasting in the Era of Big Data, 95-125, Springer, Cham. |
| `chan2020_jbes_kronecker` | [BVAR_code.zip](https://joshuachan.org/code/BVAR_code.zip) | `ca0cf486afb9f6bac9de64c7d7c4f4e8` | 2022-05-13 | 68 | Chan, J.C.C. (2020). "Large Bayesian VARs: A Flexible Kronecker Error Covariance Structure," Journal of Business and Economic Statistics, 38(1): 68-79. |
| `chan2021_ijf_mahp` | [BVAR_MAHP_code.zip](https://joshuachan.org/code/BVAR_MAHP_code.zip) | `105e7542e7208fe82267301db6406831` | 2021-01-06 | 16 | Chan, J.C.C. (2021). "Minnesota-Type Adaptive Hierarchical Priors for Large Bayesian VARs," International Journal of Forecasting, 37(3): 1212-1226. |
| `cjz2021_jae_ad_ml` | [AD_ML_code.zip](https://joshuachan.org/code/AD_ML_code.zip) | `87b55bd3f4747b29232597419579a7fc` | 2021-10-08 | 91 | Chan, J.C.C., Jacobi, L. and Zhu, D. (2022). "An Automated Prior Robustness Analysis in Bayesian Model Comparison," Journal of Applied Econometrics, 37(3): 583-602. (Code README cites the then-forthcoming version.) |
| `chan2022_qe_acp` | [BVAR_ACP_R1_code.zip](https://joshuachan.org/code/BVAR_ACP_R1_code.zip) | `ff689656c4952c734c311650d6f8f4ec` | 2026-08-27 | 15 | Chan, J.C.C. (2022). "Asymmetric Conjugate Priors for Large Bayesian VARs," Quantitative Economics, 13(3): 1145-1169. |
| `chan2023_jbes_hybtvp` | [HYB_TVPVAR_code.zip](https://joshuachan.org/code/HYB_TVPVAR_code.zip) | `84415eb7fbbe6ad207b49ecc1f6a6f3b` | 2022-05-12 | 12 | Chan, J.C.C. (2023). "Large Hybrid Time-Varying Parameter VARs," Journal of Business and Economic Statistics, 41(3): 890-905. |
| `chan2023_joe_mlvarsv` | [ml_varsv_code.zip](https://joshuachan.org/code/ml_varsv_code.zip) | `53a77cc17cd42ec153b74f463c7b4a11` | 2023-01-13 | 33 | Chan, J.C.C. (2023). "Comparing Stochastic Volatility Specifications for Large Bayesian VARs," Journal of Econometrics, 235(2): 1419-1446. |
| `chan_koop_yu2024_jbes_oisv` | [OISV_code.zip](https://joshuachan.org/code/OISV_code.zip) | `4072a2a1deec84cfe4c3b83343ac2de0` | 2026-09-01 | 41 | Chan, J.C.C., Koop, G. and Yu, X. (2024). "Large Order-Invariant Bayesian VARs with Stochastic Volatility," Journal of Business and Economic Statistics, 42(2): 825-837. |

## Notes

- **chan2022_qe_acp**: Zip updated 2026-08-27; includes the IRredu h=1 fix applied on that date.
- **Git LFS decision**: not used. Largest file is 6.2 MB (`replications/chan_koop_yu2024_jbes_oisv/legacy/results_mat/forecastingOI2-cluster.mat`), total working tree 55.3 MB - both well under GitHub limits (100 MB/file). Revisit only if a future artifact exceeds 50 MB.
