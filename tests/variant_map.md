# Variant map

For every core function: which legacy copies it canonicalizes, how identity was verified,
and — in the **never-merge** section — same-name files that must NOT be unified because they
are numerically different. Update this file with every extraction.

## Canonicalized in step 3 (zero-risk extractions, 2026-09-01)

Verification: "diff" = byte/comment-stripped diff (differences only in comments, whitespace,
variable spelling, or a closing `end`); "unit" = exact-equivalence unit test against the
legacy copy in `tests/unit/` (stochastic functions compared draw-for-draw under a fixed seed).

| Core function | Canonical source (legacy) | Also canonicalizes | Verified |
|---|---|---|---|
| `bvt.util.surform` | chan2023_jbes_hybtvp `utility/SURform.m` | chan_jeliazkov2009_statespace `sp_code/SURform.m` | diff + unit |
| `bvt.util.surform2` | chan2023_joe_mlvarsv `utility/SURform2.m` | chan2020_springer_largebvar, chan2020_jbes_kronecker `realtime_forecasts/`, chan_koop_yu2024_jbes_oisv (dead there) | diff + unit |
| `bvt.util.vec` | chan2023_joe_mlvarsv `utility/vec.m` | chan_koop_yu2024_jbes_oisv (byte-identical) | md5 + unit |
| `bvt.util.vech` | chan2023_joe_mlvarsv `utility/vech.m` | (single copy) | unit |
| `bvt.util.ldet` | chan2023_joe_mlvarsv `utility/ldet.m` | cjz2019_ad_opthyper | diff + unit |
| `bvt.util.mgammaln` | chan2023_joe_mlvarsv `utility/mgammaln.m` | cjz2019_ad_opthyper | diff + unit |
| `bvt.util.tnormrnd` | chan2023_joe_mlvarsv `utility/tnormrnd.m` | (single copy) | unit (seeded draws) |
| `bvt.sv.init_approx1N` | chan2023_joe_mlvarsv `utility/getARh_approx1N.m` (live) | chan_koop_yu2024_jbes_oisv (dead there) | diff + unit |
| `bvt.forecast.realtime_loaddata` | chan2020_springer_largebvar `loaddata.m` | chan2020_jbes_kronecker `realtime_forecasts/loaddata.m` | md5 (byte-identical); no unit test yet - needs vintage structs, covered by future forecasting regressions |
| `third_party/gigrnd.m` | chan2021_ijf_mahp `gigrnd.m` (Makalic-Schmidt 2015 / Devroye 2014) | chan2023_jbes_hybtvp, chan2023_joe_mlvarsv (all md5-identical) | md5 + unit (seeded draws) |
| `third_party/EvalFore.m` | chan_koop_yu2024_jbes_oisv (Roque Montero 2016) | (single copy) | copy |
| `third_party/tmult.m` | cjz2019_ad_opthyper `MatCode/tmult.m` | (single copy; needed by AD_ML's GetA) | copy |
| `third_party/heatmap_fx.m` | chan_koop_yu2024_jbes_oisv `heatmap.m` (MathWorks FX) | RENAMED: the original shadows MATLAB's built-in `heatmap` (R2017a+) | copy + rename |

New functions with no legacy counterpart (behavior fixed by unit tests only):
`bvt.util.build_lags` (codifies the inline `Z=[1, lags]` construction repeated in every
package - test reproduces the inline pattern exactly), `bvt.util.logsumexp`.

Edits made during extraction, in full: provenance header prepended; function renamed where
the table says so (surform, surform2, init_approx1N, realtime_loaddata, heatmap_fx). Bodies
are otherwise verbatim from the canonical source.

## NEVER MERGE - same name, numerically different

A future deduplication must not unify any of these; doing so silently changes published results.

- **`SVRW.m`**: sp_code's variant uses a DIFFUSE initial condition h_1 ~ N(0,Vh), lower-Cholesky,
  returns `[h S]`; the large_BVAR/BVAR_code/MAHP variant takes a KNOWN h0, upper-Cholesky.
  Different model, same name.
- **`sample_SV` / `sample_SVRW`**: OISV/ml_varsv stationary AR(1)-with-mean vs HYB random-walk -
  different state equations.
- **`llike_CSV_MA.m`**: the BVAR_code ROOT copy includes the `-n/2*sum(h)` term; the large_BVAR
  and realtime_forecasts copies omit it. Interchangeable inside the psi-MH at fixed h, WRONG to
  swap for marginal-likelihood ordinates (ml_BVAR_CSV_MA depends on the root version).
- **`prior_NCP.m`**: two incompatible signatures - AD_OptHyper's 5-kappa version with lag-decay
  exponent `l^kappa2` vs ml_varsv's `(p,c1,c2,Y0,Y)` with fixed `l^2` decay.
- **`prior_Minn.m`**: large_BVAR uses `Y0(end-p+1:end,:)`, ml_varsv uses `Y0(end-4+1:end,:)` -
  different data enter the AR(4) fits when p differs from 4; ml_varsv adds a U_hat output.
- **`get_resid_var.m` vs `get_resid_var_v2.m`**: the `_v2` (HYB) regresses on 4 lags of ALL
  variables plus a 1e-4 ridge - numerically different sig2, hence different Minnesota scalings.
- **`getVtheta.m`**: HYB copy hard-codes kappa3=.2, kappa4=1 inside the body; MAHP copy takes
  them from the kappa vector.
- **`sample_SVpara.m`**: ml_varsv copy has the `if mu ~= 0` vectorized gate (skips the whole mu
  block if ANY element is exactly zero); OISV splits into sample_SVpara/sample_SV0para with
  different phi truncation bounds (.998/.999/.99).
- **`macrodata_Q_2018Q4.csv`**: byte-identical between MAHP and HYB but a DIFFERENT file in
  BVAR_ACP (md5-verified). Never key a shared data folder by this filename.
- **`Min_Prior.m`** (AD packages): dual-number `{.v,.d}` builders; AD_VAR fits AR via the System
  Identification Toolbox `ar()` (different residual variances than OLS AR(4) in AD_ML). Not
  interchangeable with each other or with the plain Minnesota builders.
