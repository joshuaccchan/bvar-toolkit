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
| `bvt.forecast.realtime_loaddata` | chan2020_springer_largebvar `loaddata.m` | chan2020_jbes_kronecker `realtime_forecasts/loaddata.m` | md5 (byte-identical); real-vintage isequaln verification added in step 6 (see below) |
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

## Canonicalized in step 4 (SV/prior core, 2026-09-01)

| Core function | Canonical source (legacy) | Also canonicalizes | Verified |
|---|---|---|---|
| `bvt.sv.sv_params` | chan_koop_yu2024_jbes_oisv `utility/sample_SVpara.m` | (nothing else; the chan2023_joe_mlvarsv `sample_SVpara.m` is NOT canonicalized - never-merge, see below) | unit (seeded draws, 3 cases: r=0, r>0, mu gated at zero) |
| `bvt.sv.sv0_params` | chan_koop_yu2024_jbes_oisv `utility/sample_SV0para.m` | (single copy) | unit (seeded draws) |

Edits made during extraction, in full: provenance header prepended; functions renamed
(`sample_SVpara` -> `sv_params`, `sample_SV0para` -> `sv0_params`); the hard-coded phi
MH truncation bound promoted to an optional trailing argument `phi_bnd` DEFAULTING to
the OISV value (`.999` in `sv_params`, `.99` in `sv0_params`), so default calls
reproduce the OISV copies exactly (same draws, same rand/randn/gamrnd sequence under
the same seed). Bodies otherwise verbatim, including the `if mu~=0` gate in `sv_params`.
Passing `phi_bnd = .998` to `sv_params` does NOT reproduce the ml_varsv copy - that
copy also differs structurally (no n+r column split; mu demeans ALL columns of h).

## Canonicalized in step 5 (MAHP flagship functionization, 2026-09-01)

| Core function | Canonical source (legacy) | Also canonicalizes | Verified |
|---|---|---|---|
| `bvt.samplers.eq_gauss` | chan2021_ijf_mahp `BVAR_MNG.m` lines 40-59 (inline "sample alp and beta" block) | the verbatim inline copies in `BVAR_NG.m` 38-57, `BVAR_Minn.m` 31-50, and `forecast_BVAR_MNG/_NG/_Minn.m` (Yt/Zt/Tt renaming only). Valp/Vbeta pre-scaling (the MNG/forecast `*2`) stays with the CALLER. | unit (`test_mahp_equivalence`: draw-for-draw isequal on all stores + terminal rng state, all three estimation models) |
| `bvt.samplers.gig_shrinkage` | variant `'mng'` = chan2021_ijf_mahp `BVAR_MNG.m` 68-81; `'ng'` = `BVAR_NG.m` 66-78; `'minn'` = `BVAR_Minn.m` 59-62 | `'mng'` with `psi_floor=1e-16` also reproduces `forecast_BVAR_MNG.m` 70-83; `'minn'` also reproduces `forecast_BVAR_Minn.m` 65-68. `forecast_BVAR_NG.m` is NOT canonicalized (never-merge, below). The three variants are numerically different - never unify. | unit (same test) |
| `bvt.samplers.nu_psi_ng` | chan2021_ijf_mahp `sample_nu_psi.m` (only copy; renamed) | all four call sites (BVAR_MNG/NG and forecast_BVAR_NG two-output - the forecast flag is captured but never accumulated - and forecast_BVAR_MNG one-output) | unit (same test) |

New replication drivers (not core): `replications/chan2021_ijf_mahp/run_all.m`
(functionized estimation pipeline, models MNG/NG/Minn) and `preset.m` (every
hard-coded legacy constant with per-line citations, plus the documented
estimation-vs-forecast divergences). The legacy clock-seed line
(`randn('seed',sum(clock*100)); rand('seed',sum(clock*1000))`) is deliberately
NOT reproduced in `run_all` - it is irreproducible by construction and switches
MATLAB to the legacy v4/v5 generators; the equivalence test removes exactly that
line (its sole patch) from tempdir copies of the legacy scripts and runs both
pipelines from `rng(seed,'twister')`.

Edits made during extraction, in full: provenance headers prepended; `eq_gauss`
wrapped as a function with `np = size(Z,2)-1` replacing the literal `n*p`
(identical integers) and fresh-zeros `beta`/`alp` (the legacy scripts fully
overwrite both every sweep); `gig_shrinkage` wrapped with kappa/psi state passed
in/out, the hard-coded `1e-10` psi floor promoted to the `psi_floor` argument
(estimation `1e-10`, forecast-MNG `1e-16`), and the `Psi(idx)` reassembly
(BVAR_MNG 82-83) left with the caller; `sample_nu_psi` renamed `nu_psi_ng`, body
verbatim including its dead `count` guard. The h0 and Sigh draws (identical
3-and-2-line inline blocks in all three scripts) stay inline in `run_all`.

## Canonicalized in step 6 (forecast engine, 2026-09-01)

One entry point `bvt.forecast.iterate(branch, draw, cfg)` - called once per posterior
draw - with VERBATIM named branches, plus `bvt.forecast.tables` for the accumulation /
RMSFE / ALPL table tails. Equivalence tests run the legacy forecast scripts wholesale
from tempdir copies (byte-verbatim - unlike the MAHP estimation scripts, NO forecast
script in either package carries a clock-seed line, so the step-5 sole-patch is not
needed here; the tests assert that premise) and compare draw-for-draw against the
functionized pipeline (core priors/samplers/sv blocks + iterate), including the
terminal rng state.

| Core function (branch) | Canonical source (legacy) | Also canonicalizes | Verified |
|---|---|---|---|
| `bvt.forecast.iterate('mahp_sv')` | chan2021_ijf_mahp `forecast_BVAR_MNG.m` forecast-loop body (lines 112-147) | the textually identical tails of `forecast_BVAR_NG.m` (113-148) and `forecast_BVAR_Minn.m` (92-128) | unit (`test_forecast_iterate_mahp`: full MNG pipeline at vintages t=91 and t=T-2, isequal on tmpyhat1/tmpyhat4/all stores/kappa_hat/kappaCI + terminal rng state; t=T-2 exercises the tt==4 guard-off zeros path) |
| `bvt.forecast.iterate('springer_gauss')` | chan2020_springer_largebvar `forecast_BVAR_Minn.m` lines 36-57 | `forecast_BVAR_small.m` 41-62 (caller passes `data_tpk(:,var_small)`), `forecast_BVAR_NCP.m` 40-61, `forecast_BVAR_IP.m` 45-66, `forecast_BVAR_SSVS.m` 52-73 - given caller-supplied `A`, `CSig` IN THE LEGACY STORAGE CLASS (sparse diag Minn/small, dense chol NCP/IP/SSVS) and `dSig` (`Sig_hat'` vs `diag(Sig)'`) | unit (`test_forecast_iterate_springer`, model 2 at vintages t=41 complete and t=129 missing-latest: isequal tmpyhat0/tmpyhat1 + rng state). small/NCP/IP/SSVS callers verified textually identical, not yet run end-to-end (springer family pass) |
| `bvt.forecast.iterate('springer_csv')` | `forecast_BVAR_CSV.m` lines 67-94 | (single copy) | unit (same test, model 6 at t=129) |
| `bvt.forecast.iterate('springer_csv_t')` | `forecast_BVAR_CSV_t.m` lines 77-106 | (single copy) | unit (same test, model 7 at t=129) |
| `bvt.forecast.iterate('springer_csv_t_ma')` | `forecast_BVAR_CSV_t_MA.m` lines 109-142 | (single copy) | unit (same test, model 8 at t=129, incl. the fminunc/fminbnd psi-MH estimation stage) |
| `bvt.forecast.tables('accum_row')` | chan2020_springer_largebvar `main_forecasting.m` lines 147-154 | chan2021_ijf_mahp `main_forecasting.m` lines 98-105 (same formula); storage GUARDS (t<=T-1 / t<=T-4) stay with the caller | unit (`test_forecast_tables`: legacy lines sliced from the frozen files at test time and dispatched on synthetic arrays, incl. a complex-typed zero-imag vintage locking the magnitude-max behavior) |
| `bvt.forecast.tables('springer')` | `main_forecasting.m` lines 160-172 (both the model==1 all-variable form and the model~=1 var_core form) | (single copy) | unit (same test, both forms) |
| `bvt.forecast.tables('mahp')` | chan2021_ijf_mahp `main_forecasting.m` lines 111-116 | (single copy) | unit (same test) |

Also verified in step 6: `bvt.forecast.realtime_loaddata` now has its first REAL-vintage
verification - `test_forecast_iterate_springer` loads the 20 legacy vintage files
(cached to a tempdir .mat) and asserts isequaln against the legacy `loaddata.m` copy for
every vintage t = 41..129 (the step-3 row's "no unit test yet" caveat is retired).

Deferred to family passes (NOT canonicalized by iterate): the 9
chan2020_jbes_kronecker `realtime_forecasts/forecast_*.m` blocks (tt=1:5 skeleton
evaluating tt=1,2,3,5 into tmpyhat0/1/2/4; homoskedastic-t / MA / t-CSV / t-MA /
CSV-MA / t-CSV-MA / n-variate-SV-small error families; subsetted
`[n+var_small, 2*n+1]` accumulation in its main_forecasting tail),
cjz2019_ad_opthyper `forecast_BVAR_NCP.m` (Gaussian tt=1:4 evaluating tt=1,4 on
final-vintage data, no is_last_miss step, fmincon kappa warm-starts across
vintages), and the two chan_koop_yu2024_jbes_oisv cluster fragments
`forecast_CS_MH.m` / `forecast_SVARSV_MH.m` (12-horizon monthly design,
2n+3-column rows, `get_frcst_lkhd`/EvalFore evaluation).

Edits made during extraction, in full: provenance headers; the per-draw blocks
wrapped as branch subfunctions with `tmpyhat*(isave,:)` renamed `fc(1,:)`/`fc(2,:)`
(rows returned to the caller, which stores them - unevaluated rows return
zeros(1,2n+1), matching the legacy untouched-prealloc convention); `A_id`/`A = eye(n)`
(computed once outside the legacy MAHP loop, fully overwritten every draw) moved
inside the branch; the springer-gauss `dSig = Sig_hat'` / `diag(Sig)'` assignment
replaced by `dSig = draw.dSig;` at the same loop position (the caller computes the
legacy expression once - it is constant across tt, bit-identical); everything else
byte-verbatim, including the springer `t<=T-tt` guard that skips the h=1 evaluation
at t = T-1 and the complex-typed `sum(diag(log(CS)))` joint-density path. The
`tables` actions are the legacy lines verbatim; the legacy fprintf display blocks
are not reproduced.

Step-6 verification notes:
- Teeth check passed: a 1e-7 relative perturbation of one constant in the
  `mahp_sv` branch makes `test_forecast_iterate_mahp` FAIL on tmpyhat1.
- Path nuance (same as step 5): inside both forecast equivalence tests, the
  unqualified `gigrnd` (via `bvt.samplers.gig_shrinkage`) and `llike_CSV_MA`
  (via the model-8 psi-MH kept inline in the test pipeline) resolve to the
  tempdir LEGACY copies; the springer/realtime `llike_CSV_MA` omits the
  `-n/2*sum(h)` term of the BVAR_code root copy (see never-merge) - the
  springer family pass must keep them separate.
- The springer estimation-stage conditional draws with no core counterpart yet
  (iwishrnd Sig/A, lam/sigh2 gamrnd, rho MH with bound .999 (CSV) vs .99
  (CSV-t/CSV-t-MA), the psi MH) live VERBATIM in the test's run_core_*
  pipelines; functionizing them belongs to the springer family pass.

## Canonicalized in step 7 (OISV family pass, 2026-09-02)

Full-sample estimation pipeline of chan_koop_yu2024_jbes_oisv (Chan, Koop and Yu 2024, JBES
42(2): 825-837): main_SVAR_fullsample.m -> func_main_SVAR_v2.m -> SVARSV_MH.m ('OI') /
CS_MH.m ('CS'). Equivalence test: `tests/unit/test_oisv_equivalence.m` runs the legacy
scripts from tempdir copies at small nsim - SVARSV_MH with its ACTIVE clock-seed line
(line 24) removed as the sole patch (asserted exactly-one-occurrence, not commented);
CS_MH BYTE-VERBATIM (its clock-seed line 49 ships commented out; asserted) - OI at the
default ordering, CS at the reversed ordering, and asserts isequal on all stores, the
script-tail summaries, the six func_main outputs, and the terminal rng state.

| Core function | Canonical source (legacy) | Also canonicalizes | Verified |
|---|---|---|---|
| `bvt.util.anormrnd` | chan_koop_yu2024_jbes_oisv `utility/anormrnd.m` (single copy) | both call sites (SVARSV_MH.m line 63, forecast_SVARSV_MH.m line 57) via `b0_row_sampler`. NEVER fold into `bvt.util.tnormrnd` - different density, different rng sequence. | unit (`test_anormrnd` seeded draws + `test_oisv_equivalence`) |
| `bvt.structural.construct_Sigt` | chan_koop_yu2024_jbes_oisv `utility/construct_Sigt.m` | the private subfunction copy inside `func_main_SVAR_v2.m` lines 67-73 (comment-stripped identical, diff 2026-09-02; that copy is what the legacy func resolves at runtime) | diff + unit (`test_construct_sigt`, and end-to-end through the Sig_mean assertion in `test_oisv_equivalence`) |
| `bvt.structural.b0_row_sampler` | SVARSV_MH.m lines 49-72 (inline row-wise "sammple B0" rotation loop; caller supplies U = Y-X*A) | forecast_SVARSV_MH.m lines 43-66 (textually identical modulo Y/X/T -> Yt/Xt/Tt) | unit (`test_oisv_equivalence`) |
| `bvt.samplers.eq_svar_oi` | SVARSV_MH.m lines 76-87 (inline "sample alpha" block; caller computes tmpdV and keeps alpha = A(:)) | NOTHING else - the forecast fragment REWRITES this step (never-merge, below) | unit (same test) |
| `bvt.samplers.eq_tri_cs` | CS_MH.m lines 54-72 (inline "sample B" block; caller computes tmpdV and keeps beta = reshape(B',k_beta,1); the dead `zi` assignment kept verbatim) | forecast_CS_MH.m lines 45-63 (identical modulo Yt/Xt/Tt) | unit (same test) |
| `bvt.samplers.alp_tri_cs` | CS_MH.m lines 77-87 (inline "sample alp" count_alp loop; caller computes E = Y-XB and keeps A(A_id) = alp) | forecast_CS_MH.m lines 68-78 (identical modulo Tt) | unit (same test) |
| `bvt.samplers.horseshoe_kappa_psi` | SVARSV_MH.m lines 102-120 (psi -> z_psi -> kappa(1:2) -> z_kappa ladder, theta = alpha) | CS_MH.m lines 102-120 (theta = beta), forecast_SVARSV_MH.m lines 96-114, forecast_CS_MH.m lines 93-111 - all four textually identical modulo the coefficient vector's name. NEVER merge with `bvt.samplers.gig_shrinkage` (MAHP normal-gamma GIG ladder - different prior family). | unit (same test) |
| `bvt.priors.vtheta` (Vbeta output; row added, no new function) | - | chan_koop_yu2024_jbes_oisv `utility/getVbeta.m`: exactly vtheta's three Vbeta assignment lines on the same inputs; OISV callers use `[~,Vbeta] = bvt.priors.vtheta(...)` and discard Valp (NaN under the OI kappa(3) = NaN, never read) | diff + unit (same test) |

Reused as-is (extracted in steps 3-4, headers already list the OISV copies): `bvt.priors.resid_var_ar4`
(legacy get_resid_var), `bvt.priors.minnesota_C` (get_C), `bvt.sv.ksc_ar1_mean` (sample_SV),
`bvt.sv.sv0_params` (sample_SV0para, phi bound default .99), `bvt.sv.sv_params` (sample_SVpara,
phi bound default .999), `bvt.util.vec`, `bvt.util.build_lags` (the func_main inline lag loop),
`third_party/EvalFore.m`, `third_party/heatmap_fx.m`. `bvt.priors.impact_B0` is NOT used by OISV
(its inline Hyper.B0 = eye(n)/VB0 = ones(n) prior is a different object - see impact_B0's header).
The dead OISV utilities SURform2.m / getARh_approx1N.m stay covered by their step-3 rows
(`bvt.util.surform2`, `bvt.sv.init_approx1N`); genSV.m / gendata_SVARSV.m (simulation-study
helpers, no callers in the shipped pipeline) are left un-canonicalized in legacy/.

New replication drivers (not core): `replications/chan_koop_yu2024_jbes_oisv/run_all.m`
(functionized full-sample pipeline, `run_all(model, flip, nsim, burnin, seed)` covering the
four main_SVAR_fullsample configurations OI/CS x default/flipped) and `preset.m` (every
hard-coded legacy constant with per-line citations, plus the documented forecast-fragment
divergences under `pr.forecast`). The OI clock-seed line is deliberately NOT reproduced in
run_all (same rationale and mechanics as the step-5 MAHP note above); CS needs no patch.

Edits made during extraction, in full: provenance headers prepended; blocks wrapped as
functions with sizes recomputed from arguments ([T,n] = size(U)/size(Y)/size(E),
k = size(X,2), np = numel(idx_kappa1), nnp = numel(idx_kappa2) - identical integers);
Hyper.B0/VB0/beta0/Valp passed as explicit arguments; alpha/beta unified as `theta` in
`horseshoe_kappa_psi` (the four legacy copies differ only in that name) with the
rng-neutral Psi reassembly left to the caller (legacy position: between the psi and z_psi
draws; Psi is not read inside the ladder); `alp_tri_cs` returns alp preallocated
zeros(1,k_alp) instead of the legacy dynamic growth into the same 1 x k_alp row (fully
overwritten every sweep before any read); unqualified anormrnd/vec calls now
bvt.util.anormrnd/bvt.util.vec (code-identical). Everything else byte-verbatim, including
the OI sign fix B0(ii,:) = phii*sign(phii(ii)) and CS's dead `zi` line. The h loops
(3 lines per model), the SV-parameter calls, and the chain-init draws stay inline in
run_all. The legacy wall-clock timing displays are not reproduced.

Deferred (NOT functionized in step 7): the two forecast cluster fragments
`forecast_SVARSV_MH.m` / `forecast_CS_MH.m` and their four submain_forecasting_* drivers -
the vintage loop `for t = T0:T-1` ships commented out (submain line 52) and `t` is a
per-cluster-job input (README.txt), so the fragments are not runnable as shipped; their
goldens are the shipped `legacy/results_mat/forecasting{OI1,OI2,CS1,CS2}-cluster.mat`.
Also deferred: the evaluation/plot tails Table3_forecasting.m + get_frcst_lkhd.m +
EvalFore.m and Fig89/Fig10 plot scripts (read results_mat only). getISden_ARSS: absent
from this package (glob-verified 2026-09-02).

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
- **`sample_SVpara.m`**: the ml_varsv and OISV copies are numerically different and must not be
  unified (record corrected 2026-09-01 after full reads: BOTH copies carry the `if mu ~= 0`
  vectorized gate that skips the whole mu block if ANY element of mu is exactly zero - the gate
  is not what separates them). Real differences: (i) phi truncation bound .998 (ml_varsv) vs
  .999 (OISV); (ii) OISV takes h with n+r columns and demeans only the first n (mu applies to
  n series; r extra zero-mean columns share the phi/sig2 draws), ml_varsv demeans ALL columns;
  (iii) OISV's mu block indexes phi(1:n)/sig2(1:n), ml_varsv uses full vectors. OISV
  additionally splits the zero-mean case into `sample_SV0para.m` with bound .99. Step 4
  canonicalized the OISV pair as `bvt.sv.sv_params` / `bvt.sv.sv0_params` (phi bound exposed
  as an argument, defaults reproduce OISV exactly); the ml_varsv copy stays un-canonicalized
  on this list.
- **`macrodata_Q_2018Q4.csv`**: byte-identical between MAHP and HYB but a DIFFERENT file in
  BVAR_ACP (md5-verified). Never key a shared data folder by this filename.
- **MAHP `forecast_BVAR_NG.m` kappa/psi block**: NOT reproduced by
  `bvt.samplers.gig_shrinkage('ng',...)` at any psi_floor - its conditionals carry an extra
  factor 2 (`tmpc_j = sum(beta_j.^2./(2*psi_j))`, `tmpv_j = beta_j.^2/(2*kappa)`, lines
  72-78) pairing with its doubled Valp/Vbeta (line 43), which estimation `BVAR_NG.m` does
  NOT double. A different parameterization of the NG prior - functionize separately if the
  MAHP forecast pipeline is ever consolidated. Further estimation-vs-forecast divergences
  (psi floors 1e-10 vs 1e-16, halved psi_kappa1 init scale, Minn kappa init [.04,.04,1,100]
  vs [.4,.001,1,100]) are recorded field-by-field in
  `replications/chan2021_ijf_mahp/preset.m` under `pr.forecast`.
- **`Min_Prior.m`** (AD packages): dual-number `{.v,.d}` builders; AD_VAR fits AR via the System
  Identification Toolbox `ar()` (different residual variances than OLS AR(4) in AD_ML). Not
  interchangeable with each other or with the plain Minnesota builders.
- **OISV OI alpha step, estimation vs forecast**: `bvt.samplers.eq_svar_oi` reproduces ONLY
  SVARSV_MH.m lines 76-87 (yi = vec((Y-X*A)*B0')./Lambda, Wi = kron(B0(:,ii),X)./Lambda -
  equation-stacked, row-scaled). forecast_SVARSV_MH.m lines 69-81 REWRITE the same
  conditional: zi = reshape(B0*(Yt-[XA cols, zeroed ii])',Tt*n,1), Wi = kron(Xt,B0(:,ii))
  (time-interleaved stacking) with an explicit exp(-reshape(h',Tt*n,1)) weighting matrix.
  Same posterior, different floating-point path and stacking order - functionize separately
  if the OISV forecast fragments are ever consolidated (see
  replications/chan_koop_yu2024_jbes_oisv/preset.m `pr.forecast.oi_alpha_rewritten`).
- **OISV CS h step, estimation vs forecast**: CS_MH.m line 94 passes `mu(ii)` to sample_SV;
  forecast_CS_MH.m line 85 passes `0` - the forecast h path is drawn zero-mean while
  sample_SVpara keeps updating mu and the predictive recursion uses muh (preset
  `pr.forecast.cs_h_mu_zero`). A forecast functionization must NOT reuse the estimation
  call as-is.
- **`anormrnd.m` vs `tnormrnd.m`**: anormrnd is the OISV bimodal two-component draw for the
  first B0 rotation coordinate (one rand + one randn); tnormrnd is an inverse-cdf truncated
  normal. Same "restricted normal draw" vibe, entirely different densities and rng
  sequences - never unify.
- **`horseshoe_kappa_psi` vs `gig_shrinkage`**: the OISV horseshoe ladder (inverse-gamma /
  auxiliary z draws, kappa(1:2) global scales) and the MAHP normal-gamma (GIG) ladder are
  different prior families with different draw sequences - never unify.

## Verification notes (step 5 adversarial review, 2026-09-01)

- The perturbation "teeth" check passed: altering a single preset constant (sv_offset) in a
  scratch mirror makes test_mahp_equivalence FAIL on the stored draws - the equivalence test
  detects one-constant deviations.
- Path nuance in test_mahp_equivalence: inside the test, the unqualified gigrnd call in
  bvt.samplers.gig_shrinkage resolves to the tempdir copy of the LEGACY gigrnd.m, not
  third_party/gigrnd.m. The two differ only by the provenance header (test_gigrnd covers the
  third_party copy separately), but the guarantee silently depends on those files staying
  code-identical - never edit one without the other.
- The step-2 golden .mat for MAHP stores alp_hat/beta_hat/h_hat/kappa_hat only; the golden
  nu_psi mean (0.214989) lives in the golden run log, not the .mat.

## Verification notes (step 6 adversarial review, 2026-09-01)

- Verdict EQUIVALENT. Teeth checks: a horizon-index perturbation and a 1e-7 constant
  perturbation in scratch mirrors both fail the forecast tests; suite green on revert.
- Of the 11 canonicalized blocks, 6 rest on independent byte-diffs rather than end-to-end
  tests: the MAHP NG/Minn forecast tails (byte-identical to the tested MNG canonical, so
  effectively covered) and the springer small/NCP/IP/SSVS bodies (verbatim modulo the
  declared caller-supplied pieces). The CALLER contract for those four springer models
  (A reshape, CSig storage class - sparse diag for Minn/small vs dense chol for NCP/IP/SSVS -
  dSig expression, model-1 var_small subsetting) is untested until the springer family pass.
- Byte-level audit expectation: forecast_BVAR_IP.m lines 56/61 and forecast_BVAR_SSVS.m
  lines 63/68 write "-.5*(...)" where the canonical Minn body has "- .5*(...)" - same
  parse, numerically identical.
- The springer guard-off vintage path (t = T-1/T skipping the h=1 evaluation) is not
  exercised end-to-end; the guard expression is byte-verbatim, the analogous MAHP
  guard-off path IS tested end-to-end (t = T-2), and the accumulation-side guard is
  covered on synthetic arrays including t = T.
- Complex-typing nuance (headers corrected): on R2025b diag() demotes the zero-imag
  complex diagonal to real, so real-data rows stay real on both sides; older MATLABs may
  retain the complex attribute - equivalence holds either way (byte-identical expression).

## Verification notes (step 7 adversarial review, 2026-09-02)

- Teeth check passed: a 1e-7 relative perturbation of one ladder constant in
  `bvt.samplers.horseshoe_kappa_psi` makes `test_oisv_equivalence` FAIL on store_kappa;
  suite green on revert.
- The six "textually identical modulo renaming" claims of the step-7 table (B0 block,
  CS B block, CS alp block, and the ladder across all four scripts) were verified
  mechanically: comment-stripped, whitespace-normalized diffs with Yt/Xt/Tt -> Y/X/T and
  alpha/beta -> theta come back empty; the OI alpha estimation-vs-forecast diff is
  NON-empty, confirming the never-merge entry's direction.
- Record corrected in `bvt.sv.sv_params`'s header: the step-4 note "in the OISV CS_MH run
  mu is initialized at zero, so mu is never updated there" was WRONG - CS_MH.m lines 34-37
  initialize mu(ii) = mean(log(s2i)) from the data (almost surely all-nonzero), so the
  `if mu~=0` gate passes and mu IS updated every sweep. The equivalence test's store_hpara
  comparison (columns 1:n are the mu draws) covers it draw-for-draw.
- Path nuance in test_oisv_equivalence (same shape as steps 5-6): the legacy side resolves
  the tempdir copies of anormrnd/vec/sample_SV*/get_*/getVbeta/construct_Sigt, the core
  side resolves only qualified bvt.* names - no unqualified name in any step-7 core
  function can be shadowed by the tempdir. Both packages define a `run_all`; the test
  pins resolution to the OISV package with a `which` assertion.
- CS shape quirks reproduced by construction: legacy `alp` is a 1 x k_alp ROW (dynamic
  growth) and `alp_tri_cs` returns the same row shape; `z_kappa` enters sweep 1 as a 2x1
  column and leaves every ladder call as a 1x2 row - kept verbatim. The equivalence
  test's isequal locks the STORES and terminal rng state (values); these internal
  orientations are draw-irrelevant (orientation-agnostic indexing) and are not
  independently asserted. Also: the test pairs OI-with-default and CS-with-flipped
  ordering; the cross pairings are unexercised end-to-end (flip only reverses var_id
  upstream of the model switch - negligible risk, swap the pairings once if paranoid).
