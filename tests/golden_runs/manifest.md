# Golden-run manifest

One row per entry script worth capturing. Runtime classes: **min** (< 30 min), **hours**, **days**.
Toolboxes beyond base MATLAB are listed per package. "Blockers" are known issues found in the
2026-09-01 code audit; patch on the `build/` copy only.

## Quick wins (capture first)

| Package | Entry | Runtime | Notes |
|---|---|---|---|
| `chan_jeliazkov2009_statespace` | `TVPVAR.m` | min | 4-variable TVP-VAR via precision sampler. Stats Tbx (gamrnd, wishrnd). Scripts sit in nested `sp_code/` subfolder; must run from that folder (bare `load USdata.csv`). CAPTURED 2026-09-01: runs clean (10k sweeps, ~37 s) but prints progress only - results go to figures, so this golden certifies run-health, not numbers; add a save-posterior-means patch if numeric goldens are needed. |
| `chan_jeliazkov2009_statespace` | `UC.m` | min | UC-SV on US CPI inflation. CAPTURED 2026-09-01: runs clean (10k sweeps, ~2 s); same figures-only caveat as TVPVAR.m. |
| `chan2019wp_acp` | `main_BVAR_ACP.m` | min | Exact sampling + closed-form ML; no MCMC chains to wait on. CAPTURED 2026-09-01: optimal-kappa table (symmetric vs asymmetric, log-ML -9436 vs -9201), ~2 s. |
| `chan2022_qe_acp` | `main_ACP_jointden.m` | min-hours | Analytic log-ML on a (kappa1,kappa2) meshgrid. Uses deprecated `xlsread` (still works in R2025b, emits warning). CAPTURED 2026-09-01: runs clean, ~41 s, but output is the contour figure only (lost in -batch) - patch a savefig/save if a numeric golden is needed. |
| `chan_koop_yu2024_jbes_oisv` | `Table3_forecasting.m` | min | Builds published Table 3 (RMSFE/ALPL) from the SHIPPED `results_mat/*.mat` cluster outputs - ideal golden, fully deterministic. CAPTURED 2026-09-01: full table with DM stars + `Table3_LaTex_cluster.xls`, ~2 s. |
| `chan_koop_yu2024_jbes_oisv` | `Fig89_plot_varcov_first2last2.m`, `Fig10_plot_B0_fullsample.m` | min | Figure inputs from shipped fullsample .mat files. Bundled `heatmap.m` shadows the R2017a+ builtin - fine inside the package folder. CAPTURED 2026-09-01: Fig 8/9 .eps pair, Fig 10 heatmap .eps, and `Table2_kappa.csv` (deterministic kappa estimates). |

## Medium jobs

| Package | Entry | Runtime | Notes |
|---|---|---|---|
| `cjz2018_ad_var` | `main_AD_VAR.m` | min-hours | Needs System Identification Tbx (`ar()` in Min_Prior) + Stats Tbx + deprecated `xlsread`. 15 globals; run standalone only. |
| `cjz2019_ad_opthyper` | `main_AD_OptHpyer.m` (sic) | hours | fmincon ML-maximization over kappa. Optimization + Stats Tbx. `main_AD_forecasting.m` calls `addpath(genpath('..'))` - climbs out of package root; run from an isolated build copy only (the harness does this). |
| `cjz2021_jae_ad_ml` | `main_ADML_VAR.m` | hours | Stats Tbx. Legacy seed syntax verified working on R2025b - no patch expected. |
| `cjz2021_jae_ad_ml` | `main_ADML_FactorModel.m` | hours | Uses `readmatrix` (R2019a+), `rng('default')` - reproducible as shipped. |
| `chan2021_ijf_mahp` | `main_BVAR.m` | min | 23-variable MCMC, three priors (MNG/NG/Minn selected in-script; as-shipped default model=1 MNG, nsim=10000+1000). Legacy `randn('seed',clock)` syntax VERIFIED WORKING on R2025b - no patch; runs are clock-seeded as shipped. CAPTURED 2026-09-01 (~2.5 min): as-shipped run-health, plus a savegolden variant (patch appends print/save of the computed posterior means): kappa1 0.0343, kappa2 0.000611, nu_psi 0.215, and golden_posterior_means.mat (beta_hat, alp_hat, h_hat, kappa_hat) - the regression anchor for the step-5 functionization. |
| `chan2023_joe_mlvarsv` | `main_varsv.m` | min | As shipped: n=15, model=2. CAPTURED 2026-09-01, all 5 models via model-selector patches (log-ML, n=15): VAR-NCP -6918.8, VAR-CSV -6618.1, VAR-SV -6443.7, VAR-FSV -6475.5, VAR-SVO -6386.0. Runtimes 0.2 s - 30 min, far below the audit's "hours" guess. KNOWN SUSPECT: `ml_var_arsvo_redu` (ngrid 31/32 mismatch, o_hat linear-indexing) ran clean; -6386.0 is the documented as-shipped value pending adjudication against the published JoE table. |
| `chan2023_jbes_hybtvp` | `main_HYB_TVPSV.m` | hours | The audit flagged the `xlsread`-with-cell-range-on-csv load as broken, but it was VERIFIED WORKING on R2025b (2026-09-01): returns exactly 238x248, matching the csv - no patch needed. |
| `chan2020_jbes_kronecker` | `main_BVAR.m` | hours | Full-sample, 8 models, 30000 draws. Legacy `randn('seed',...)` syntax verified working on R2025b - no patch expected. `addpath('./realtime_forecasts')` creates path-order shadowing between root and realtime copies of sample_h/llike_CSV_MA - run models one at a time and record which copy resolved. |
| `chan2020_jbes_kronecker` | `ml_*` scripts | hours | Marginal likelihoods. KNOWN BUGS to capture as-is (they are behind the published numbers): `ml_BVAR_MA.m` stale `psi`; `ml_BVAR_CSV_t_MA.m` stale `Sig` + leftover-workspace reliance. Do NOT fix in the golden run. |

## Multi-day jobs (schedule separately)

| Package | Entry | Runtime | Notes |
|---|---|---|---|
| `chan2020_springer_largebvar` | `main_forecasting.m` | days | Recursive 1975Q1-2015Q4 x 8 models x 5000 draws. Run per-model via the in-script `model` selector (1-8); each model saves `<model_name>.mat` and prints RMSFE/ALPL tables. Stats + Optimization Tbx; deprecated `xlsread` on 20 vintage files. |
| `chan2020_jbes_kronecker` | `main_forecasting.m` | days | Recursive real-time exercise, 9 forecast scripts. |
| `chan2022_qe_acp` | `main_ACP_apps.m` (15-var) | days | Sign-restriction accept-reject documented as taking days. The 6-variable variant is hours - capture that first. |
| `cjz2019_ad_opthyper` | `main_AD_forecasting.m` | days | Recursive 1984Q1-2018Q4; kappa warm-starts persist across vintages (workspace state - a functionization hazard the golden must pin down). |
| `chan_koop_yu2024_jbes_oisv` | `submain_forecasting_*.m` | n/a | NOT runnable as shipped: cluster-array fragments (loop and save commented out, scalar `t` supplied externally). Golden = the shipped `results_mat/*.mat` + `Table3_forecasting.m` output. |
