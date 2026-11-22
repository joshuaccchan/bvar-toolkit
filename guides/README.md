# How-To Guides

Each guide shows how to use the toolkit's core functions on one small example with generated data.
A guide is a page, `README.md`, and a script, `guide.m`, that prints the numbers and draws the
figures on the page. The guides come in series, and within a series each guide builds on the one
before.

## Getting Started with Bayesian VARs

The guides use the same generated VAR(2) of three variables.

| Guide | What it shows | Core functions |
|---|---|---|
| [How to Stack the Equations of a VAR](var_stacking/) | Stack the equations in matrix form and as one long regression | `bvar.util.build_lags`, `bvar.util.surform2` |
| [How to Draw the VAR Coefficients under the Minnesota Prior](coefficient_draw/) | Get posterior draws of the VAR coefficients | `bvar.priors.minn`, `bvar.samplers.var_coef_joint` |
| [How to Estimate a VAR by Gibbs Sampling](gibbs_sampler/) | Draw the coefficients and the error covariance matrix in turn, and check the chain | `bvar.samplers.var_coef_joint`, `bvar.diag.inefficiency_factor`, `bvar.diag.mcse`, `bvar.diag.geweke` |
| [How to Choose the Shrinkage of a VAR by the Marginal Likelihood](natural_conjugate_ml/) | Choose the shrinkage and the lag length by the closed-form marginal likelihood | `bvar.ml.ncp`, `bvar.ml.kron_bvar` |

## Forecasting with a VAR

The guides use data simulated from a VAR(4) fitted to US real GDP growth, PCE inflation and the
unemployment rate, 1960Q1 to 2019Q4.

| Guide | What it shows | Core functions |
|---|---|---|
| [How to Simulate Data from an Estimated VAR](var_simulation/) | Simulate data from a VAR fitted to US GDP growth, inflation and unemployment | `bvar.util.build_lags`, `bvar.util.simulate_var` |
| [How to Forecast with a Bayesian VAR](var_forecast/) | Forecast from posterior draws, and compare with the true model's forecasts | `bvar.ml.ncp`, `bvar.util.simulate_var`, `bvar.forecast.predictive`, `bvar.forecast.mixquantile` |
| [How to Evaluate Forecasts Out of Sample](var_forecast_evaluation/) | Score forecasts made at each quarter with the data of the time | `bvar.ml.ncp`, `bvar.forecast.predictive`, `bvar.util.logsumexp` |
| [How to Choose the Lag Length and the Shrinkage for Forecasting](var_lag_shrinkage/) | Choose the lag length and the shrinkage by the marginal likelihood at each forecast origin | `bvar.ml.ncp`, `bvar.util.build_lags` |

## VARs for a Small Open Economy

The guides use data simulated from a VAR(3) fitted to US and Australian data, 1984Q1 to 2019Q4,
in which the Australian variables do not enter the US equations.

| Guide | What it shows | Core functions |
|---|---|---|
| [How to Estimate a VAR for a Small Open Economy](soe_var/) | Estimate a VAR with block exogeneity one block at a time | `bvar.models.var_soe`, `bvar.ml.ncp`, `bvar.util.soe_reduced_form` |
| [How to Test Block Exogeneity by the Marginal Likelihood](soe_exogeneity_test/) | Compare the VAR with and without block exogeneity by their marginal likelihoods | `bvar.models.var_soe`, `bvar.ml.ncp`, `bvar.ml.lniwpdf` |
| [How to Trace a Foreign Shock through a Small Open Economy](soe_foreign_shock/) | Compute the responses to a foreign shock and their credible bands | `bvar.models.var_soe`, `bvar.structural.irf_redu`, `bvar.structural.fevd`, `bvar.util.shaded_band` |
| [How to Forecast a Small Open Economy under a Foreign Scenario](soe_conditional_forecast/) | Forecast the Australian variables given a path for the US variables | `bvar.models.var_soe`, `bvar.forecast.predictive`, `bvar.forecast.mixquantile` |

The [examples](../examples/) run the core functions on additional models, with real and generated
data, and the [tutorials](../tutorials/) answer empirical questions with the settings of published
papers.
