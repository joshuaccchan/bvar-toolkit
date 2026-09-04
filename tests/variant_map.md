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
| `bvar.util.surform` | chan2023_jbes_hybtvp `utility/SURform.m` | chan_jeliazkov2009_statespace `sp_code/SURform.m` | diff + unit |
| `bvar.util.surform2` | chan2023_joe_mlvarsv `utility/SURform2.m` | chan2020_springer_largebvar, chan2020_jbes_kronecker `realtime_forecasts/`, chan_koop_yu2024_jbes_oisv (dead there) | diff + unit |
| `bvar.util.vec` | chan2023_joe_mlvarsv `utility/vec.m` | chan_koop_yu2024_jbes_oisv (byte-identical) | md5 + unit |
| `bvar.util.vech` | chan2023_joe_mlvarsv `utility/vech.m` | (single copy) | unit |
| `bvar.util.ldet` | chan2023_joe_mlvarsv `utility/ldet.m` | cjz2019_ad_opthyper | diff + unit |
| `bvar.util.mgammaln` | chan2023_joe_mlvarsv `utility/mgammaln.m` | cjz2019_ad_opthyper | diff + unit |
| `bvar.util.tnormrnd` | chan2023_joe_mlvarsv `utility/tnormrnd.m` | (single copy) | unit (seeded draws) |
| `bvar.sv.init_approx1N` | chan2023_joe_mlvarsv `utility/getARh_approx1N.m` (live) | chan_koop_yu2024_jbes_oisv (dead there) | diff + unit |
| `bvar.forecast.realtime_loaddata` | chan2020_springer_largebvar `loaddata.m` | chan2020_jbes_kronecker `realtime_forecasts/loaddata.m` | md5 (byte-identical); real-vintage isequaln verification added in step 6 (see below) |
| `third_party/gigrnd.m` | chan2021_ijf_mahp `gigrnd.m` (Makalic-Schmidt 2015 / Devroye 2014) | chan2023_jbes_hybtvp, chan2023_joe_mlvarsv (all md5-identical) | md5 + unit (seeded draws) |
| `third_party/EvalFore.m` | chan_koop_yu2024_jbes_oisv (Roque Montero 2016) | (single copy) | copy |
| `third_party/tmult.m` | cjz2019_ad_opthyper `MatCode/tmult.m` | (single copy; needed by AD_ML's GetA) | copy |
| `third_party/heatmap_fx.m` | chan_koop_yu2024_jbes_oisv `heatmap.m` (MathWorks FX) | RENAMED: the original shadows MATLAB's built-in `heatmap` (R2017a+) | copy + rename |

New functions with no legacy counterpart (behavior fixed by unit tests only):
`bvar.util.build_lags` (codifies the inline `Z=[1, lags]` construction repeated in every
package - test reproduces the inline pattern exactly), `bvar.util.logsumexp`.

Edits made during extraction, in full: provenance header prepended; function renamed where
the table says so (surform, surform2, init_approx1N, realtime_loaddata, heatmap_fx). Bodies
are otherwise verbatim from the canonical source.

## Canonicalized in step 4 (SV/prior core, 2026-09-01)

| Core function | Canonical source (legacy) | Also canonicalizes | Verified |
|---|---|---|---|
| `bvar.sv.sv_params` | chan_koop_yu2024_jbes_oisv `utility/sample_SVpara.m` | chan2023_joe_mlvarsv `utility/sample_SVpara.m` at `phi_bnd = .998` (established in step 9 - the original "never-merge" verdict here was WRONG, see the step-9 correction below) | unit (seeded draws, 3 cases: r=0, r>0, mu gated at zero; ml_varsv in `test_sv_params_mlvarsv`) |
| `bvar.sv.sv0_params` | chan_koop_yu2024_jbes_oisv `utility/sample_SV0para.m` | (single copy) | unit (seeded draws) |

Edits made during extraction, in full: provenance header prepended; functions renamed
(`sample_SVpara` -> `sv_params`, `sample_SV0para` -> `sv0_params`); the hard-coded phi
MH truncation bound promoted to an optional trailing argument `phi_bnd` DEFAULTING to
the OISV value (`.999` in `sv_params`, `.99` in `sv0_params`), so default calls
reproduce the OISV copies exactly (same draws, same rand/randn/gamrnd sequence under
the same seed). Bodies otherwise verbatim, including the `if mu~=0` gate in `sv_params`.

**Superseded (step 9, 2026-09-03).** The step-4 claim "passing `phi_bnd = .998` does NOT
reproduce the ml_varsv copy - it also differs structurally (no n+r column split; mu
demeans ALL columns of h)" is refuted. Those structural differences are no-ops whenever
`numel(mu) == size(h,2)`, i.e. r = 0, and EVERY ml_varsv call site satisfies that
(`VAR_CSV.m` 61 h is T x 1 with mu = 0; `VAR_ARSV_redu.m` 84 and `VAR_ARSVO_redu.m` 91
h is T x n with mu n x 1; `VAR_FSV.m` 82 h is T x (n+r) with mu (n+r) x 1 - so r = 0
inside `sv_params` there too, the factor columns arriving as extra "series" with their
own mu). The phi bound is the only difference that can bite, and it does bite:
`test_sv_params_mlvarsv` forces a candidate into [.998,.999) and shows .999 accepting
where .998 and the legacy reject, and a scratch-mirror teeth check that hard-codes .999
inside `sv_params` makes `test_mlvarsv_equivalence` fail on model 3. r > 0 remains
uncovered by the claim: there the two bodies genuinely differ.

## Canonicalized in step 5 (MAHP flagship functionization, 2026-09-01)

| Core function | Canonical source (legacy) | Also canonicalizes | Verified |
|---|---|---|---|
| `bvar.samplers.eq_gauss` | chan2021_ijf_mahp `BVAR_MNG.m` lines 40-59 (inline "sample alp and beta" block) | the verbatim inline copies in `BVAR_NG.m` 38-57, `BVAR_Minn.m` 31-50, and `forecast_BVAR_MNG/_NG/_Minn.m` (Yt/Zt/Tt renaming only). Valp/Vbeta pre-scaling (the MNG/forecast `*2`) stays with the CALLER. | unit (`test_mahp_equivalence`: draw-for-draw isequal on all stores + terminal rng state, all three estimation models) |
| `bvar.samplers.gig_shrinkage` | variant `'mng'` = chan2021_ijf_mahp `BVAR_MNG.m` 68-81; `'ng'` = `BVAR_NG.m` 66-78; `'minn'` = `BVAR_Minn.m` 59-62 | `'mng'` with `psi_floor=1e-16` also reproduces `forecast_BVAR_MNG.m` 70-83; `'minn'` also reproduces `forecast_BVAR_Minn.m` 65-68. `forecast_BVAR_NG.m` is NOT canonicalized (never-merge, below). The three variants are numerically different - never unify. | unit (same test) |
| `bvar.samplers.nu_psi_ng` | chan2021_ijf_mahp `sample_nu_psi.m` (only copy; renamed) | all four call sites (BVAR_MNG/NG and forecast_BVAR_NG two-output - the forecast flag is captured but never accumulated - and forecast_BVAR_MNG one-output) | unit (same test) |

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

One entry point `bvar.forecast.iterate(branch, draw, cfg)` - called once per posterior
draw - with VERBATIM named branches, plus `bvar.forecast.tables` for the accumulation /
RMSFE / ALPL table tails. Equivalence tests run the legacy forecast scripts wholesale
from tempdir copies (byte-verbatim - unlike the MAHP estimation scripts, NO forecast
script in either package carries a clock-seed line, so the step-5 sole-patch is not
needed here; the tests assert that premise) and compare draw-for-draw against the
functionized pipeline (core priors/samplers/sv blocks + iterate), including the
terminal rng state.

| Core function (branch) | Canonical source (legacy) | Also canonicalizes | Verified |
|---|---|---|---|
| `bvar.forecast.iterate('mahp_sv')` | chan2021_ijf_mahp `forecast_BVAR_MNG.m` forecast-loop body (lines 112-147) | the textually identical tails of `forecast_BVAR_NG.m` (113-148) and `forecast_BVAR_Minn.m` (92-128) | unit (`test_forecast_iterate_mahp`: full MNG pipeline at vintages t=91 and t=T-2, isequal on tmpyhat1/tmpyhat4/all stores/kappa_hat/kappaCI + terminal rng state; t=T-2 exercises the tt==4 guard-off zeros path) |
| `bvar.forecast.iterate('springer_gauss')` | chan2020_springer_largebvar `forecast_BVAR_Minn.m` lines 36-57 | `forecast_BVAR_small.m` 41-62 (caller passes `data_tpk(:,var_small)`), `forecast_BVAR_NCP.m` 40-61, `forecast_BVAR_IP.m` 45-66, `forecast_BVAR_SSVS.m` 52-73 - given caller-supplied `A`, `CSig` IN THE LEGACY STORAGE CLASS (sparse diag Minn/small, dense chol NCP/IP/SSVS) and `dSig` (`Sig_hat'` vs `diag(Sig)'`) | unit (`test_forecast_iterate_springer`, model 2 at vintages t=41 complete and t=129 missing-latest: isequal tmpyhat0/tmpyhat1 + rng state). small/NCP/IP/SSVS callers verified textually identical, not yet run end-to-end (springer family pass) |
| `bvar.forecast.iterate('springer_csv')` | `forecast_BVAR_CSV.m` lines 67-94 | (single copy) | unit (same test, model 6 at t=129) |
| `bvar.forecast.iterate('springer_csv_t')` | `forecast_BVAR_CSV_t.m` lines 77-106 | (single copy) | unit (same test, model 7 at t=129) |
| `bvar.forecast.iterate('springer_csv_t_ma')` | `forecast_BVAR_CSV_t_MA.m` lines 109-142 | (single copy) | unit (same test, model 8 at t=129, incl. the fminunc/fminbnd psi-MH estimation stage) |
| `bvar.forecast.tables('accum_row')` | chan2020_springer_largebvar `main_forecasting.m` lines 147-154 | chan2021_ijf_mahp `main_forecasting.m` lines 98-105 (same formula); storage GUARDS (t<=T-1 / t<=T-4) stay with the caller | unit (`test_forecast_tables`: legacy lines sliced from the frozen files at test time and dispatched on synthetic arrays, incl. a complex-typed zero-imag vintage locking the magnitude-max behavior) |
| `bvar.forecast.tables('springer')` | `main_forecasting.m` lines 160-172 (both the model==1 all-variable form and the model~=1 var_core form) | (single copy) | unit (same test, both forms) |
| `bvar.forecast.tables('mahp')` | chan2021_ijf_mahp `main_forecasting.m` lines 111-116 | (single copy) | unit (same test) |

Also verified in step 6: `bvar.forecast.realtime_loaddata` now has its first REAL-vintage
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
  unqualified `gigrnd` (via `bvar.samplers.gig_shrinkage`) and `llike_CSV_MA`
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
| `bvar.util.anormrnd` | chan_koop_yu2024_jbes_oisv `utility/anormrnd.m` (single copy) | both call sites (SVARSV_MH.m line 63, forecast_SVARSV_MH.m line 57) via `b0_row_sampler`. NEVER fold into `bvar.util.tnormrnd` - different density, different rng sequence. | unit (`test_anormrnd` seeded draws + `test_oisv_equivalence`) |
| `bvar.structural.construct_Sigt` | chan_koop_yu2024_jbes_oisv `utility/construct_Sigt.m` | the private subfunction copy inside `func_main_SVAR_v2.m` lines 67-73 (comment-stripped identical, diff 2026-09-02; that copy is what the legacy func resolves at runtime) | diff + unit (`test_construct_sigt`, and end-to-end through the Sig_mean assertion in `test_oisv_equivalence`) |
| `bvar.structural.b0_row_sampler` | SVARSV_MH.m lines 49-72 (inline row-wise "sammple B0" rotation loop; caller supplies U = Y-X*A) | forecast_SVARSV_MH.m lines 43-66 (textually identical modulo Y/X/T -> Yt/Xt/Tt) | unit (`test_oisv_equivalence`) |
| `bvar.samplers.eq_svar_oi` | SVARSV_MH.m lines 76-87 (inline "sample alpha" block; caller computes tmpdV and keeps alpha = A(:)) | NOTHING else - the forecast fragment REWRITES this step (never-merge, below) | unit (same test) |
| `bvar.samplers.eq_tri_cs` | CS_MH.m lines 54-72 (inline "sample B" block; caller computes tmpdV and keeps beta = reshape(B',k_beta,1); the dead `zi` assignment kept verbatim) | forecast_CS_MH.m lines 45-63 (identical modulo Yt/Xt/Tt) | unit (same test) |
| `bvar.samplers.alp_tri_cs` | CS_MH.m lines 77-87 (inline "sample alp" count_alp loop; caller computes E = Y-XB and keeps A(A_id) = alp) | forecast_CS_MH.m lines 68-78 (identical modulo Tt); plus, from step 9, chan2023_joe_mlvarsv `VAR_ARSV_redu.m` 64-73 and `VAR_ARSVO_redu.m` 71-80 (the latter through the new optional `o` argument) | unit (same test; ml_varsv in `test_mlvarsv_equivalence`) |
| `bvar.samplers.horseshoe_kappa_psi` | SVARSV_MH.m lines 102-120 (psi -> z_psi -> kappa(1:2) -> z_kappa block, theta = alpha) | CS_MH.m lines 102-120 (theta = beta), forecast_SVARSV_MH.m lines 96-114, forecast_CS_MH.m lines 93-111 - all four textually identical modulo the coefficient vector's name. NEVER merge with `bvar.samplers.gig_shrinkage` (MAHP normal-gamma GIG block - different prior family). | unit (same test) |
| `bvar.priors.vtheta` (Vbeta output; row added, no new function) | - | chan_koop_yu2024_jbes_oisv `utility/getVbeta.m`: exactly vtheta's three Vbeta assignment lines on the same inputs; OISV callers use `[~,Vbeta] = bvar.priors.vtheta(...)` and discard Valp (NaN under the OI kappa(3) = NaN, never read) | diff + unit (same test) |

Reused as-is (extracted in steps 3-4, headers already list the OISV copies): `bvar.priors.resid_var_ar4`
(legacy get_resid_var), `bvar.priors.minnesota_C` (get_C), `bvar.sv.ksc_ar1_mean` (sample_SV),
`bvar.sv.sv0_params` (sample_SV0para, phi bound default .99), `bvar.sv.sv_params` (sample_SVpara,
phi bound default .999), `bvar.util.vec`, `bvar.util.build_lags` (the func_main inline lag loop),
`third_party/EvalFore.m`, `third_party/heatmap_fx.m`. `bvar.priors.impact_B0` is NOT used by OISV
(its inline Hyper.B0 = eye(n)/VB0 = ones(n) prior is a different object - see impact_B0's header).
The dead OISV utilities SURform2.m / getARh_approx1N.m stay covered by their step-3 rows
(`bvar.util.surform2`, `bvar.sv.init_approx1N`); genSV.m / gendata_SVARSV.m (simulation-study
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
draws; Psi is not read inside the block); `alp_tri_cs` returns alp preallocated
zeros(1,k_alp) instead of the legacy dynamic growth into the same 1 x k_alp row (fully
overwritten every sweep before any read); unqualified anormrnd/vec calls now
bvar.util.anormrnd/bvar.util.vec (code-identical). Everything else byte-verbatim, including
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

## Canonicalized in step 8 (Kronecker family pass, part 1: estimation + marginal likelihood, 2026-09-02)

Full-sample estimation + marginal-likelihood pipeline of chan2020_jbes_kronecker (Chan 2020,
JBES 38(1): 68-79): main_BVAR.m -> BVAR*.m (8 models) -> ml_BVAR_*.m (7 scripts; model 1's ML
is the inline cp_ml block of BVAR.m). The realtime forecasting pipeline (main_forecasting.m +
realtime_forecasts/) is part 2 - NOT canonicalized here. Equivalence test:
`tests/unit/test_kron_equivalence.m` runs the legacy cp_ml = 1 pipeline for ALL EIGHT models
from tempdir copies at small nsim - the seven MCMC estimation scripts with their ACTIVE
clock-seed lines removed as the sole patch (asserted exactly-one-occurrence each, not
commented); BVAR.m, all seven ml_*, all four intlike_* and every helper BYTE-VERBATIM
(asserted seed-line-free) - and asserts isequal on all stores, counters, script-tail
summaries, every ML piece (ML/llike/lpri/lpost/final store_lpost), and the terminal rng
state. Models 4 and 8 are compared under bugcompat (below); their corrected defaults are
additionally asserted to differ exactly where each bug lives and match everywhere else.

| Core function | Canonical source (legacy) | Also canonicalizes | Verified |
|---|---|---|---|
| `bvar.ml.lniwpdf` | chan2020_jbes_kronecker `lniwpdf.m` (single copy) | all prior/posterior NIW ordinates in the 8 ML computations | unit (`test_kron_ml_densities` bitwise + end-to-end) |
| `bvar.ml.linvgammpdf` | `linvgammpdf.m` (single copy) | the sigh2 ordinates (models 3/5/7/8) | unit (same tests) |
| `bvar.ml.llike_ma` | `llike_MA.m` (root; body verbatim incl. its `chol(Sig)'` upper-transposed Cholesky) | BVAR_MA.m + ml_BVAR_MA.m psi targets. realtime_forecasts/llike_MA.m is NOT canonicalized (its function line is named llike_MA1; part 2). | unit (same tests) |
| `bvar.ml.llike_csv_ma` | `llike_CSV_MA.m` (package ROOT copy WITH the -n/2*sum(h) term) | the psi targets of BVAR_t_MA/BVAR_CSV_MA/BVAR_CSV_t_MA and ml_BVAR_t_MA/ml_BVAR_CSV_MA/ml_BVAR_CSV_t_MA (with h := log(lam) / U pre-scaled by sqrt(lam) in the t models, exactly as the legacy calls do). The realtime/springer reduced copies stay never-merge (below). | unit (`test_kron_ml_densities`: bitwise vs root AND asserted to differ from the realtime copy by n/2*sum(h)) + end-to-end |
| `bvar.ml.intlike_csv` | `intlike_BVAR_CSV.m` (renamed; body verbatim) | (single copy) | unit (`test_kron_intlike` bitwise seeded, real data, + end-to-end model 3) |
| `bvar.ml.intlike_t_csv` | `intlike_BVAR_t_CSV.m` | (single copy) | unit (same, + end-to-end model 5) |
| `bvar.ml.intlike_csv_ma` | `intlike_BVAR_CSV_MA.m` | (single copy) | unit (same, + end-to-end model 7) |
| `bvar.ml.intlike_csv_t_ma` | `intlike_BVAR_CSV_t_MA.m` | (single copy; carries the first-observation scale quirk - see audit notes below) | unit (same, + end-to-end model 8) |
| `bvar.ml.kron_bvar` | BVAR.m lines 36-47 (inline cp_ml block) | (model 1; analytic, no rng) | unit (`test_kron_equivalence` model 1) |
| `bvar.ml.kron_bvar_t` | ml_BVAR_t.m | (model 2; deterministic given stores) - CLEAN BILL | unit (same, model 2) |
| `bvar.ml.kron_bvar_csv` | ml_BVAR_CSV.m | (model 3) - CLEAN BILL; chain-continuation leftovers made explicit (h/rho from last stored draws; countrho continues the ESTIMATION counter via est.state.countrho); its inline reduced-run h step = `bvar.sv.csv_armh(s2,rho,sigh2,h,n,isim==1,h_mean)` (NR start promoted, see below) | unit (same, model 3) |
| `bvar.ml.kron_bvar_ma` | ml_BVAR_MA.m | (model 4) - AFFECTED: line-17 leftover-psi llike term; `'bugcompat',true` reproduces it bitwise from est.state.psi, default corrects to psi_mean | unit (same, model 4 bugcompat bitwise + corrected-differs-only-in-llike teeth) |
| `bvar.ml.kron_bvar_t_csv` | ml_BVAR_t_CSV.m | (model 5) - CLEAN BILL; reduced-run rho bound .999 vs estimation .9999 kept verbatim | unit (same, model 5) |
| `bvar.ml.kron_bvar_t_ma` | ml_BVAR_t_MA.m | (model 6) - CLEAN BILL; optimizer warm start est.state.psihat made explicit; tmpden(psiidx)-centered den_psi normalization kept verbatim | unit (same, model 6) |
| `bvar.ml.kron_bvar_csv_ma` | ml_BVAR_CSV_MA.m | (model 7) - CLEAN BILL; its per-draw Hpsi rebuild is the pattern model 8 violates; dead `ht = h_mean` line kept as dead assignment | unit (same, model 7) |
| `bvar.ml.kron_bvar_csv_t_ma` | ml_BVAR_CSV_t_MA.m | (model 8) - AFFECTED twice: frozen leftover Hpsi/psi ordinate loop (lines 42-44) and leftover last-draw Sig in the reduced-run psi target (line 108); `'bugcompat',true` reproduces both bitwise from est.state.psi/est.state.Sig, default corrects (per-draw Hpsi; Sig_mean) | unit (same, model 8 bugcompat bitwise + corrected teeth: llike unchanged, lpost(1) moves, lpost(2:3) unchanged, ML moves) |

Core reuse inside the functionized estimation (`replications/chan2020_jbes_kronecker/run_all.m`):
`bvar.priors.niw('kron_script')` (construct_prior_A + the callers' S0/nu0), `bvar.util.build_lags`
(the inline X loop), `bvar.sv.csv_armh` (root sample_h; all estimation h steps and the ml reduced
runs of models 5/7/8), `bvar.sv.nu_studentt` (root sample_nu; normpdf-form MH ratio - decisions
verified identical, step-4 note), `bvar.ml.llike_ma`/`llike_csv_ma` (the estimation psi-MH
targets). The (Sig,A) joint draw, lam, sigh2, rho-MH and psi-MH blocks have no core counterpart
and live verbatim in run_all's per-model subfunctions (as in the step-6 note, functionizing
them further belongs to the springer family pass, which shares their structure).

Edits made during extraction, in full: provenance headers; scripts wrapped as functions with
data/priors/stores passed explicitly and sizes recomputed from arguments (identical integers);
`bvar.sv.csv_armh` gained an optional 7th argument `ht_start` (default h = previous behavior
bit-for-bit; ml_BVAR_CSV's inline reduced-run h step is its body with ht_start = h_mean and
is_ForcedAccept = (isim==1) - the only executable differences); leftover-workspace reads of
the legacy ml scripts became explicit `est.state.*` arguments (final chain draws
Sig/A/h/lam/rho/sigh2/psi/nu, final psi-MH mode psihat, counters - captured by run_all after
the last sweep; the legacy leftovers equal the last stored draws wherever both exist, and
Hpsi/Hrho are rebuilt from the corresponding scalar, bitwise identical since the legacy only
ever updates matrix and scalar together); `options = optimset('Display','off','LargeScale','off')`
reconstructed inside the ml functions (identical struct); the legacy clock-seed lines dropped
from run_all (mechanics as steps 5/7); figures and wall-clock displays not reproduced. Dead
Hrho rebuilds: `run_all` KEEPS all four estimation-side rebuilds verbatim (models 3/5/7/8,
annotated `%#ok<NASGU>`) - the fidelity-safe choice, and necessary for model 3, whose rebuild
IS consumed downstream when the coupled `cp_ml=1` pipeline hands its leftover Hrho to
ml_BVAR_CSV.m line 51; the ml-side rebuilds in models 5/7/8 (nothing reads Hrho there -
`bvar.sv.csv_armh` and `sample_h` rebuild it from rho internally) are not reproduced. The dead
`ht = h_mean` line of ml_BVAR_CSV_MA kept. R (IS
draws) and nsims2 (reduced-run length) exposed as options DEFAULTING to the legacy hard-coded
values (cited in preset.m). Everything else byte-verbatim, including the m6-specific
tmpden(psiidx)-centered psi-density normalization, model 8's +/-.999 psigrid, and the
rho-truncation zoo recorded in preset.m `pr.rho_mh_bnd_est` / `pr.rho_mh_bnd_ml`.

### Step-8 bug audit (2026-09-02) - all findings, including clean bills

Verified line-by-line from source; the two AFFECTED scripts reproduce bitwise under bugcompat
and are corrected by default (corrected = consistent evaluation point across ordinate pieces;
method quirks shared by both modes are listed as quirks, not corrected):

- **ml_BVAR_MA.m line 17 (defect, model 4)**: `-.5*(s2(1)/(1+psi^2) + sum(s2(2:end)))` uses
  the loop variable `psi` left over from BVAR_MA.m (the final chain draw = store_psi(nsims))
  where `psi_mean` is intended (lines 11/16 use psi_mean). The published BVAR-MA log-ML
  depends on the random final chain state through this one llike term.
- **ml_BVAR_CSV_t_MA.m lines 42-44 (defect, model 8)**: the (A,Sig) ordinate loop reads
  `Hpsi\X`, `Hpsi\shortY` and `1/(1+psi^2)` with Hpsi/psi left over from the estimation run
  and never rebuilds them from the stored psi draws - every conditional NIW term conditions
  on the single random final psi (ml_BVAR_CSV_MA.m lines 26-30 rebuild per draw).
- **ml_BVAR_CSV_t_MA.m line 108 (defect, model 8)**: the reduced-run psi target is
  `llike_CSV_MA(x,U_psi,Sig,h)` with `Sig` the final estimation DRAW (leftover) where
  Sig_mean is intended - the same reduced run's h/lam steps condition on Sig_mean (line 78),
  and ml_BVAR_CSV_MA.m line 83 uses Sig_mean.
- **ml_BVAR_t.m (model 2)**: CLEAN - all ordinates at (A_mean, Sig_mean, nu_mean); no
  leftover reads beyond stores/priors; no rng.
- **ml_BVAR_CSV.m (model 3)**: CLEAN - consistent starred point; leftover reads are chain
  continuation (h/rho/Hrho from the final state, which equals the last stored draws) plus the
  countrho counter continuation; reduced-run h step NR-starts at h_mean with first-sweep
  forced accept (unlike the later scripts, which call sample_h = NR start at current h).
- **ml_BVAR_t_CSV.m (model 5)**: CLEAN on evaluation points. Quirk: reduced-run rho MH bound
  .999 vs estimation .9999.
- **ml_BVAR_t_MA.m (model 6)**: CLEAN - its llike line 19 applies s2(1)/(1+psi_mean^2)
  correctly (the exact term model 4 gets wrong). Leftover psihat/options are optimizer
  continuation. Quirk: den_psi normalization centered at tmpden(psiidx) instead of max.
- **ml_BVAR_CSV_MA.m (model 7)**: CLEAN - per-draw Hpsi rebuild present; Sig_mean used in the
  psi target. Quirk: dead `ht = h_mean;` line 61 (vestige of model 3's inline h step).
- **intlike_BVAR_CSV_t_MA.m (QUIRK, both modes)**: deny_h receives the transformed Utld and
  applies NO (1+psi^2) scale and NO -n/2*log(1+psi^2) constant to the first observation,
  though the Gibbs sampler gives it variance (1+psi^2)exp(h_1)lam_1*Sig and the Gaussian
  intlike_BVAR_CSV_MA applies both. One-observation model/ordinate mismatch inside the
  published model-8 log-ML; a likelihood-formula property shared by both modes (documented,
  not silently changed).
- **BVAR_CSV_t_MA.m lam step (QUIRK, estimation + its ml reduced run)**: no (1+psi^2)
  first-observation correction in the lam_1 draw (model 6 corrects its lam step; model 8's h
  step corrects). Estimation and reduced run are mutually consistent; kept verbatim.
- **rho truncation zoo (VERIFIED variance of the audit's ".999 vs .99" note)**: estimation
  .9999 (models 3/5/7) vs .99 (model 8, BVAR_CSV_t_MA.m line 92); ml reduced runs .9999
  (models 3/7) vs .999 (models 5/8). All rho ordinate grids span (-.999,.999) regardless.
- **model-8 psigrid (QUIRK)**: (-.999,.999) where every other MA model uses (-.99,.99); the
  prior truncates at +/-.99, so the extra grid points carry the -1e10 penalty (~zero mass).
- **grid normalization (QUIRK, family-wide)**: inserting the starred value into a sorted
  uniform grid and normalizing by nugrid(2)-nugrid(1) treats the grid as uniform though one
  interval is split - inherent to the published method, shared by both modes everywhere.

New replication drivers (not core): `replications/chan2020_jbes_kronecker/run_all.m`
(functionized estimation, models 1-8), `run_ml.m` (estimation + ML on one stream, exposing
'bugcompat'), and `preset.m` (every hard-coded constant with per-line citations, including
the truncation-bound tables and the llike_CSV_MA path-hazard record).

## Canonicalized in step 9 (ml_varsv family pass, estimation only, 2026-09-03)

Estimation pipeline of chan2023_joe_mlvarsv (Chan 2023, JoE 235(2): 1419-1446):
main_varsv.m -> VAR_NCP.m / VAR_CSV.m / VAR_ARSV_redu.m / VAR_FSV.m / VAR_ARSVO_redu.m.
The marginal-likelihood routines `utility/ml_var_*.m` are a SEPARATE phase and are not
touched here; model 1 is the exception, its log-ML being inline and analytic in VAR_NCP.m.
Equivalence test: `tests/unit/test_mlvarsv_equivalence.m` runs all five legacy scripts from
tempdir copies at nsim/burnin = 60/20 over 15 configurations - every model at every switch
setting it reads, plus models 3/4/5 repeated at the paper's active n = 15 selection - and
asserts isequal on all stores, counters, script-tail summaries and the terminal rng state.
The four MCMC scripts' clock-seed lines are removed as the sole patch (asserted
exactly-one-occurrence and active, per file); VAR_NCP.m is byte-verbatim (asserted
seed-line-free). cp_ml = 0 for models 2-5 so no ml_var_* routine is entered.

| Core function | Canonical source (legacy) | Also canonicalizes | Verified |
|---|---|---|---|
| `bvar.samplers.eq_var_redu_tri` | chan2023_joe_mlvarsv `VAR_ARSV_redu.m` lines 44-57 (the reduced-form column-at-a-time coefficient draw) | `VAR_ARSVO_redu.m` 51-64 through the optional `o` argument - the two blocks' only textual difference is the Lambda line's `.*repmat(o,1,n-ii+1)`, mechanically diffed. NOT merged into `bvar.samplers.eq_svar_oi` (below) | unit (`test_mlvarsv_equivalence` models 3 and 5) |
| `bvar.samplers.factor_fsv` | `VAR_FSV.m` lines 38-44 (joint precision-sampler draw of the whole T x r factor path) | (single copy in the repo) | unit (same test, model 4) |
| `bvar.samplers.eq_fsv_load` | `VAR_FSV.m` lines 48-71 (per-equation joint draw of A and the free loadings L) | (single copy in the repo) | unit (same test, model 4) |
| `bvar.sv.svo_outlier` | `VAR_ARSVO_redu.m` lines 112-124 (grid draw of the outlier scales o_t, then the beta draw of po) | (single copy in the repo) | unit (same test, model 5; the test asserts at least one o_t > 1 is drawn) |
| `bvar.samplers.alp_tri_cs` (row added; optional `o` argument added) | - | `VAR_ARSV_redu.m` 64-73 and `VAR_ARSVO_redu.m` 71-80: the step-7 OISV body verbatim modulo naming (beta/Hyper.Vbeta for alp/Hyper.Valp), the SVO copy differing only by the iD line's `./o.^2` | unit (same test, models 3 and 5) |
| `bvar.sv.sv_params` (row amended; no code change) | - | chan2023_joe_mlvarsv `utility/sample_SVpara.m` at `phi_bnd = .998` - see the step-4 correction above | unit (`test_sv_params_mlvarsv`, all three ml_varsv shapes + the bound teeth; end-to-end in `test_mlvarsv_equivalence`) |

Reused as-is: `bvar.priors.niw('mlvarsv_ncp')` (legacy prior_NCP), `bvar.priors.minn` with
n0pre = 4 (prior_Minn), `bvar.priors.impact_B0` (prior_B0), `bvar.priors.minnesota_C` (get_C),
`bvar.sv.init_approx1N` (getARh_approx1N), `bvar.sv.csv_armh` (sample_CSV), `bvar.sv.ksc_ar1_mean`
(sample_SV), `bvar.util.build_lags` / `vec` / `vech` / `ldet` / `mgammaln`,
`third_party/gigrnd.m`. Not reached by the estimation scripts: `utility/getISden_ARSS.m`
(read only by the ml_var_* routines), `SURform2.m` and `tnormrnd.m` (no caller in the five
scripts; both already covered by their step-3 rows).

New replication drivers (not core): `replications/chan2023_joe_mlvarsv/run_all.m`
(`run_all(model, is_kappafixed, is_kappasym, nsim, burnin, seed, varid)` over the five models,
honouring both main_varsv switches) and `preset.m` (every hard-coded constant with per-line
citations). The four clock-seed lines are deliberately not reproduced in run_all (same
rationale and mechanics as steps 5/7/8); VAR_NCP needs no patch.

Edits made during extraction, in full: provenance headers; blocks wrapped as functions with
sizes recomputed from arguments (`[T,n] = size(Y)`, `k = size(X,2)`, `r = size(L,2)`/`size(F,2)`,
`ngrid = numel(o_grid)-1` - identical integers); `Hyper.*` passed explicitly; unqualified `vec`
now `bvar.util.vec`; `o` preallocated in `svo_outlier` instead of updated in place (every element
is overwritten before any read). The `o` arguments of `eq_var_redu_tri` and `alp_tri_cs` default
to `ones(T,1)`, which leaves the non-SVO path bit-for-bit unchanged - multiplication and
division by exactly 1.0 are identity in IEEE-754, and `test_oisv_equivalence` re-verifies the
`alp_tri_cs` default path end-to-end. Everything else byte-verbatim, including the dead
`U = zeros(T,n)` and `Hphi` assignments (kept in run_all with `%#ok<NASGU>`), the dead
`E = Y-X*A` of VAR_FSV.m 21, VAR_FSV's `store_A` running-sum convention, and the legacy
`if isim < 20` forced-accept window of the VAR-CSV h step. The wall-clock timing displays,
VAR_CSV.m's exp(h/2) figure, and main_varsv.m's dead `y = reshape(Y',n*T,1)` are not reproduced.

Deferred (NOT functionized in step 9): the VAR-CSV natural-conjugate (Sig,A) joint draw
(VAR_CSV.m 37-46) - the same block chan2020_jbes_kronecker/run_all.m keeps inline, held for
the springer family pass, which shares its structure - and the four `ml_var_*` routines with
their importance-sampling machinery (`getISden_ARSS.m`, `M`, `flag_marg`). The `ml_var_*`
routines were extracted in step 10 (below); the (Sig,A) draw is still deferred.

### Step-9 verification notes (2026-09-03)

- Teeth checks (scratch mirror of `core/`, prepended to the path; the real tree was never
  perturbed, and the mirror's reverted run passes): a 1e-7 relative perturbation of
  `eq_var_redu_tri` fails on model 3's store_alp, of `factor_fsv` and of `eq_fsv_load` on
  model 4's store_l, of `alp_tri_cs`'s `./o.^2` on model 5 only (model 3 unaffected, as the
  o == 1 identity predicts), and of `alp_tri_cs`'s shared `X_alpi'*iD*X_alpi` on BOTH
  `test_mlvarsv_equivalence` (model 3) and `test_oisv_equivalence` (CS) - which is what proves
  the two packages share that block rather than merely resembling each other. Hard-coding .999
  inside a mirrored `sv_params` fails on model 3.
- **Sensitivity floor of the outlier block.** `svo_outlier`'s o_t draw is discrete (a 32-point
  grid), so small perturbations are absorbed: 1e-7, 1e-6 and 1e-4 relative changes to its
  log-likelihood constant are NOT detected, while 2e-2 is (and so is a 2% change to the
  `-n*log(o_grid)` term). This is a property of the block, not a gap in the test - the
  surrounding continuous draws are detected at 1e-7 - but a future refactor of `svo_outlier`
  must not rely on the equivalence test to catch sub-1e-3 arithmetic drift there.
- Path nuance (same shape as steps 5-8): inside `test_mlvarsv_equivalence` the legacy side
  resolves the tempdir copies of prior_Minn/prior_NCP/prior_B0/sample_SV/sample_SVpara/
  sample_CSV/get_C/getARh_approx1N/vec/vech/ldet/mgammaln, and run_all's ONE unqualified name
  (`gigrnd`) resolves to the tempdir legacy copy too. The test now ASSERTS that the legacy
  gigrnd.m and `third_party/gigrnd.m` are code-identical (comment-stripped), converting the
  step-5 silent dependency into a checked one. Four packages define `run_all`; resolution is
  pinned with a `which` assertion, and `preset` by run_all's own cd guard.
- `test_sv_params_mlvarsv` asserts `sample_SVpara` resolves from the ml_varsv legacy copy, not
  the same-name OISV one.
- Data note: the shipped `macrodata_Q_2019Q4.csv` is 242 x 248, so T = 234 - the
  `linspace(1961,2019.75,T)` axis in VAR_CSV.m's figure implies 236 and is cosmetically off by
  two quarters. Display only; no estimate depends on it.
- There is NO constant named `nuub` (or `nu_ub`/`nub`) anywhere in this legacy package
  (grep-verified). A lead carried into this step that does not correspond to these files.

## Canonicalized in step 10 (ml_varsv family pass, part 2: marginal likelihoods, 2026-09-03)

The marginal-likelihood phase of chan2023_joe_mlvarsv: `utility/ml_var_csv.m`,
`ml_var_arsv_redu.m`, `ml_var_fsv.m`, `ml_var_arsvo_redu.m`, dispatched by each estimation
script's `if cp_ml` tail. **The method is adaptive importance sampling, not the Chib
decomposition of the Kronecker family**: for each model a small set of blocks is integrated
out analytically (always the VAR coefficients; also Sig for VAR-CSV, the latent factors for
VAR-FSV, and the log-volatility variances under `flag_marg = 2`), everything else is drawn
from importance densities fitted to the posterior draws - the log-volatility paths from
`getISden_ARSS` (a Gaussian state space fit, rho by `fminbnd` on the concentrated likelihood),
phi truncated normal, mu normal, the loadings/impact elements multivariate normal, kappa gamma
(`gamfit`), po beta (`betafit`) - and the M log weights are averaged in 50 batches, which is
also where the reported numerical standard error comes from. There is no starred point and no
posterior-ordinate telescoping, so the step-8 audit categories (evaluation point, leftover
workspace) mostly do not apply: these are functions, not workspace scripts, and everything
they read is an argument.

Equivalence test: `tests/unit/test_mlvarsv_ml.m` runs each legacy estimation script with
`cp_ml = 1` from tempdir copies (so the legacy ML routine executes inside it, on the
estimation's own stream) at nsim/burnin = 40/10, M = 100, n = 4, over 11 configurations -
every model at every switch setting its cprior/gIS branch reads - and asserts `isequal` on the
stored draws, on `lml` and `lmlstd`, on `store_w` where the legacy script keeps it (VAR-FSV),
and on the terminal rng state. The four ml routines and every density utility run BYTE-VERBATIM
(asserted seed-line-free); the estimation scripts carry the same sole clock-seed patch as step
9. VAR-SVO is compared under `'bugcompat', true`. T = 234 (the full sample) is deliberate: the
`o_hat` linear-index defect only takes its published form when T >= 32.

| Core function | Canonical source (legacy) | Also canonicalizes | Verified |
|---|---|---|---|
| `bvar.ml.lgampdf` | chan2023_joe_mlvarsv `utility/lgampdf.m` (single copy) | the kappa prior/IS ordinates of ml_var_csv (its only caller) | unit (`test_mlvarsv_ml_densities` bitwise + end-to-end) |
| `bvar.ml.ltnormpdf` | `utility/ltnormpdf.m` (single copy) | the phi prior/IS ordinates of ml_var_csv (its only caller) | unit (same) |
| `bvar.ml.lmvnpdf_pcn` | `utility/lmvnpdf_pcn.m` (single copy) | the h IS ordinate of ml_var_csv (its only caller) | unit (same) |
| `bvar.ml.isden_arss` | `utility/getISden_ARSS.m` (renamed; body and its `concen_like_h` subfunction verbatim; the only copy in the repo) | all four ml routines, once per series | unit (all five outputs bitwise) + end-to-end |
| `bvar.ml.mlvarsv_csv` | `utility/ml_var_csv.m` (incl. its `like_VAR_CSV` subfunction) | (model 2) - CLEAN BILL | unit (`test_mlvarsv_ml`, both kappa settings) |
| `bvar.ml.mlvarsv_arsv_redu` | `utility/ml_var_arsv_redu.m` | (model 3) - CLEAN BILL | unit (same, all three switch settings) |
| `bvar.ml.mlvarsv_fsv` | `utility/ml_var_fsv.m` (incl. its `deny_fsv` subfunction) | (model 4) - CLEAN BILL; the only routine implementing `flag_marg = 1` | unit (same, all three settings; `store_w` bitwise) |
| `bvar.ml.mlvarsv_arsvo_redu` | `utility/ml_var_arsvo_redu.m` | (model 5) - AFFECTED three times, all in the outlier block; `'bugcompat',true` reproduces the legacy bitwise, the default corrects | unit (same, all three settings, bugcompat bitwise + corrected teeth) |

Reused rather than re-extracted: `bvar.priors.niw('mlvarsv_ncp')` (legacy prior_NCP, called per
IS draw to refresh VA from the drawn kappa), `bvar.priors.minn` with n0pre = 4 (prior_Minn),
`bvar.priors.impact_B0` (prior_B0), `bvar.util.tnormrnd`, `bvar.util.vec`, `bvar.util.ldet`,
`bvar.util.mgammaln`, `bvar.util.surform2` (SURform2 - `deny_fsv` is the first CORE function to
call it; the springer/kronecker forecast scripts call the legacy copies). Not extracted:
`utility/ligampdf.m`, which is CODE-IDENTICAL to the step-8
`bvar.ml.linvgammpdf` (asserted in `test_mlvarsv_ml_densities`) and has no caller in this
package, and `utility/lmvnpdf.m`, which has no caller anywhere in the package. Both stay in
legacy/.

New replication driver (not core): `replications/chan2023_joe_mlvarsv/run_ml.m` -
`run_ml(model, is_kappafixed, is_kappasym, nsim, burnin, seed, varid, 'M',.., 'flag_marg',..,
'bugcompat',..)`, the legacy `cp_ml = 1` pipeline: run_all followed by the ML routine on one
continuous rng stream, as the legacy script tails do. `run_all` gained three output fields
(`out.Y`, `out.X`, `out.Y0`) so run_ml need not rebuild the design.

Edits made during extraction, in full: provenance headers; the four routines renamed
`ml_var_*` -> `mlvarsv_*` with the legacy positional outputs `[lml,lmlstd]` kept and the third
output promoted from the bare `store_w` to a detail struct (`out.store_w`, `out.bigml`, the
fitted IS parameters, and for VAR-SVO `out.store_lr_o` / `out.store_lJ_o`); helper calls
redirected to the core names above; an explicit `assert` on `flag_marg` (the legacy switch
leaves `prior`/`gIS` undefined for unimplemented values, so this converts a confusing
undefined-handle error into a clear one and changes no arithmetic); the VAR-SVO corrections
below behind `'bugcompat'`. Everything else byte-verbatim, including the dead `gamfit` block
over `1./sig2` and its commented-out consumers, the commented-out `big_sig2` blocks of the two
Cholesky-SV routines, `M = 50*ceil(M/50)`, the 50-batch reshape, and `kappa3 = 100`
re-hard-coded inside each routine.

### Step-10 bug audit (2026-09-03) - all findings, including clean bills

Verified line-by-line from source. The three defects are confined to VAR-SVO's outlier block
and are exactly the o-specific lines the routine gained when it was copied from
`ml_var_arsv_redu.m` (a comment-stripped diff of the two files returns only the o/po
additions - and `c1`, the one line the o patch should have touched and did not).

- **`ml_var_arsvo_redu.m` line 12 with line 180 (defect 1, prior mass):** line 12 sets
  `ngrid = size(o_grid,1)` = **32**, the number of ATOMS in `o_grid = [1;linspace(2,20,31)']`;
  line 180 then builds `o_lpri = log([1-po; repmat(po/ngrid,ngrid,1)])`. The sampler's own
  `ngrid` (VAR_ARSVO_redu.m line 8) is **31**, the number of OUTLIER atoms, and its line 112
  builds the same expression with it. So the ML charges every outlier draw `log(po/32)` where
  the model gives it `log(po/31)` - 0.0317487 too little per outlier period - and produces a
  33-entry vector whose last entry is never indexed (`o_idx <= 32`). The ML's implied prior
  sums to `1 - po/32 = 0.998047` over the atoms the sampler can produce. Correct value: po/31.
- **`ml_var_arsvo_redu.m` line 181 (defect 2, indexing):**
  `lr_o = sum(o_lpri(o_idx) - log(o_hat(o_idx)));`. `o_hat` is T x 32 (line 68) and `o_idx` is
  T x 1 with values in 1..32, so `o_hat(o_idx)` is a LINEAR index into a T-row matrix: for
  T >= 32 (T = 234 published) every index lands in column 1. It reads the smoothed probability
  of the atom o = 1 at PERIOD `o_idx(t)` instead of the probability of the drawn atom at period
  t - the importance density's own ordinate, evaluated at the wrong entry. Correct:
  `o_hat(sub2ind(size(o_hat),(1:T)',o_idx))`. (Below T = 32 the linear indices spill into
  column 2 and the behaviour changes; the test therefore runs at the full T.)
- **`ml_var_arsvo_redu.m` line 153 (defect 3, missing Jacobian):**
  `c1 = -n*T/2*log(2*pi) -.5*sum(sum(h)) -.5*sum(log(Hyper.Valp));` is byte-identical to
  ml_var_arsv_redu.m line 135, but line 155 scales the residuals by the outlier size -
  `diag_sqrt_D = vec(exp(h/2).*repmat(o,1,n))`. The o scaling enters the quadratic form and
  its log-determinant, `-n*sum(log(o))`, never enters the normalizing constant. The model is
  unambiguous: the sampler's own o step carries the term
  (`lliket = -n*log(o_grid) -.5*U(tt,:)*U(tt,:)'./o_grid.^2;`, VAR_ARSVO_redu.m 115). Since
  `log(o_t) >= 0`, the omission inflates the weight of every draw that contains an outlier
  (and leaves an all-ones draw alone). It is the largest of the three.
- **`ml_var_csv.m` (model 2): CLEAN.** `like_VAR_CSV` is the correct matric-t ordinate
  (`pi^{-nT/2}`, `|Omega|^{-n/2} = -n/2*sum(h)`, both `|V_A|` and `|K_A|` at power -n/2, the
  IW constants via `mgammaln`); the sigh2 integration on lines 50-52 carries its full
  normalizer including `-T/2*log(2*pi)` and `.5*log(1-phi^2)`; `Hyper.VA` is refreshed from
  each kappa draw (line 45) before every read; mu is fixed at 0 throughout, matching
  VAR_CSV.m's `sample_SVpara(h,0,phi,Hyper)`.
- **`ml_var_arsv_redu.m` (model 3): CLEAN.** Every parameter is accounted for exactly once
  (alp and sigh2 analytically, beta/h/mu/phi/kappa by importance sampling); `c1` carries the
  full Jacobian `-.5*sum(sum(h))` (|det B0| = 1); `lh_pri` and `lh_g` both carry their
  `-T*n/2*log(2*pi)`; `lr_beta` is prior-minus-IS with the two `(2*pi)^{-k_beta/2}` constants
  correctly cancelled; Valp/Vbeta are refreshed from the drawn kappas (lines 127-128) before
  every read.
- **`ml_var_fsv.m` (model 4): CLEAN.** `deny_fsv` marginalizes the factors through
  `Sy = (I kron L) Omega (I kron L') + Sig` with row-major stacking consistent with
  `SURform2(X,n)` and `reshape(Y',T*n,1)`; the free-loading index set matches the sampler's.
  **The one thing that looks like a defect and is not:** `c_hi = .5*ldet(Kh_hat)` (line 125)
  omits the `-T*(n+r)/2*log(2*pi)` of the h importance density, and `lh_prior` (lines 156-161)
  omits exactly the same constant, once per series, in BOTH `flag_marg` branches - the two
  omissions cancel in `llike + lh_prior - lh_g`, so the level is right. (The VAR-SV routine
  carries the constant on both sides instead; the two routines agree.)
- **Family-wide quirk, phi truncation:** the samplers accept a phi candidate only when
  `abs(phic) < .998` (sample_SVpara.m line 22), so the posterior they target is truncated at
  .998, while every ML routine scores the phi prior AND draws its importance density on
  (-1,1). The IS estimator remains a valid estimate of the marginal likelihood under the
  (-1,1)-truncated prior (the IS support contains the integrand's), but that prior is not
  quite the one the MCMC targets. At the published hyperparameters the two prior normalizers
  differ by log(0.6554/0.6406) = 0.0229 per log-volatility series - about 0.02 for VAR-CSV,
  0.34 for VAR-SV/VAR-SVO, 0.39 for VAR-FSV. Shared by all five models and left unchanged: it
  is a prior-definition choice, not an implementation slip.
- **Quirk, `flag_marg`:** only `ml_var_fsv` implements `flag_marg = 1`. In the other three the
  `switch flag_marg case 2` block is the only place `prior`/`gIS` are defined, so any other
  value fails later with an undefined-function error. The extractions assert on it instead.
- **Quirk, dead work:** `ml_var_arsv_redu`/`ml_var_arsvo_redu` fit a gamma to `1./sig2` for
  every series (lines 47-50) whose only consumers are commented out (line 51); `ml_var_fsv`
  draws `big_sig2` (M x (n+r) `gamrnd`) even under `flag_marg = 2`, where nothing reads it -
  rng-consuming but numerically inert. Both kept verbatim.
- **Quirk, is_kappasym prior rows:** ml_var_arsv_redu/arsvo_redu score the symmetric-kappa
  prior with `Hyper.c0(2:3,:)` while the estimation draws kappa1 from row 1. Rows 1 and 2 of
  `c0` are both `[1, 1/.2^2]`, so this is numerically a no-op - but it would stop being one if
  the two rows were ever given different hyperparameters.
- **Quirk, gamma scale round trip:** `ml_var_csv` draws `gamrnd(ckappa_hat(1),1./ckappa_hat(2),...)`
  where `ckappa_hat(2) = 1/tmp_hat(2)`, while the other three pass `tmp_hat(2)` straight
  through. Mathematically the same scale, up to one ulp different in floating point. A future
  consolidation of the four kappa blocks must not silently pick one form.

### Step-10 verification notes (2026-09-03)

- Teeth checks, all seven FAILED as required (scratch mirror of `core/` appended AFTER the real
  core so the mirror wins - `addpath` prepends, so the order matters; the real tree was never
  perturbed and the suite is green): a 1e-7 relative perturbation of the `-T/2*log(2*pi)` in
  `mlvarsv_csv` fails on model 2's lml, of `c1` in `mlvarsv_arsv_redu` on model 3, of `c_hi` in
  `mlvarsv_fsv` on model 4, of the `isden_arss` rho bound (.99 -> .98) on model 2. Applying any
  ONE of the three corrections inside the bugcompat path fails on model 5, which is what proves
  the three are materially separate: at the test scale (n = 4, T = 234, nsim = 40, M = 100) the
  legacy lml is -2044.06 and becomes -2040.88 with defect 1 alone corrected (+3.2), -1671.77
  with defect 2 alone (+372.3) and -2926.41 with defect 3 alone (-882.3).
- Corrected-mode teeth INSIDE the test: the corrected VAR-SVO run is asserted to leave the
  terminal rng state bitwise unchanged (none of the three corrections touches the stream), to
  use 31 prior atoms where bugcompat uses 32, to apply a strictly negative o Jacobian on at
  least one draw, and to move each weight by exactly `store_lJ_o + (lr_o_corrected -
  lr_o_legacy)` to 1e-6 (the residual is the re-association of `c1`).
- Path nuance (same shape as steps 5-9): inside `test_mlvarsv_ml` the legacy side resolves the
  tempdir copies of the four ml routines and every helper (asserted by `which` for the script
  and its ml routine); the core side resolves only qualified `bvar.*` names plus the package's
  run_all/run_ml/preset, with run_all's one unqualified `gigrnd` again asserted code-identical
  to `third_party/gigrnd.m`.
- Toolbox note: the ML phase needs the Statistics and Machine Learning Toolbox beyond what the
  samplers need - `gamfit` (all four routines), `betafit`/`betarnd`/`betapdf` (VAR-SVO),
  `normcdf`/`norminv` (via tnormrnd and the truncation normalizers). `fminbnd` is base MATLAB.
- Fragility worth knowing: `o_hat` is built with `sum(store_o == o_grid(ii))`, exact
  floating-point equality against the grid. It is safe because the chain assigned those values
  FROM the same array and `run_ml` hands the estimation's own `out.o_grid` back - but a
  recomputed grid with different endpoints would silently make every count zero and leave
  `o_hat` uniform at 1/(nsim+32).
- Do not read the test-scale VAR-SVO numbers as indicative. At nsim = 40 the smoothed `o_hat`
  puts ~1.6% on each outlier atom (50% outlier probability per period), so the importance
  draws are wildly over-dispersed relative to the posterior and the corrected weights - which,
  unlike the legacy ones, actually penalise implausible o - collapse. The published
  configuration (nsim = 10000) is the one that matters; see the adjudication below.

## Ranking adjudication (step 10, 2026-09-03)

Full-length runs at the published configuration (n = 15) comparing the two modes are
complete on two seeds, with a VAR-SV control; the bugcompat path reproduces the published
values within Monte Carlo error, on an affected and an unaffected model - the full-scale
end-to-end check the small-nsim equivalence test cannot give. The quantitative record of
the corrected computation is withheld pending completion of the adjudication (runs at the
remaining model dimensions are in progress) and the author's review; it will be recorded
here when settled.

- Runtime: test_mlvarsv_ml is ~16 s for its 11 configurations (both sides), and
  test_mlvarsv_ml_densities ~1 s.
- Cost warning for anyone re-running the published configuration: at n = 15, M = 10000 the
  ML phase dominates the estimation - the golden logs record 7 s (VAR-CSV), 563 s (VAR-SV),
  1769 s (VAR-FSV) and 464 s (VAR-SVO) against 11-113 s for the corresponding chains.

## NEVER MERGE - same name, numerically different

A future deduplication must not unify any of these; doing so silently changes published results.

- **`SVRW.m`**: sp_code's variant uses a DIFFUSE initial condition h_1 ~ N(0,Vh), lower-Cholesky,
  returns `[h S]`; the large_BVAR/BVAR_code/MAHP variant takes a KNOWN h0, upper-Cholesky.
  Different model, same name.
- **`sample_SV` / `sample_SVRW`**: OISV/ml_varsv stationary AR(1)-with-mean vs HYB random-walk -
  different state equations.
- **`llike_CSV_MA.m`**: the BVAR_code ROOT copy includes the `-n/2*sum(h)` term; the large_BVAR
  and realtime_forecasts copies omit it. Interchangeable inside the psi-MH at fixed h, WRONG to
  swap for marginal-likelihood ordinates (ml_BVAR_CSV_MA depends on the root version). Step 8
  canonicalized the ROOT copy as `bvar.ml.llike_csv_ma` (test_kron_ml_densities asserts the
  n/2*sum(h) gap against the realtime copy); the legacy resolution is cwd/path-order-dependent
  after main_forecasting.m's addpath('./realtime_forecasts') - recorded in the kronecker
  preset.m (`pr.llike_csv_ma_root_has_sumh`). The realtime/springer reduced copies stay
  un-canonicalized for the part-2/springer forecast passes.
- **`prior_NCP.m`**: two incompatible signatures - AD_OptHyper's 5-kappa version with lag-decay
  exponent `l^kappa2` vs ml_varsv's `(p,c1,c2,Y0,Y)` with fixed `l^2` decay.
- **`prior_Minn.m`**: large_BVAR uses `Y0(end-p+1:end,:)`, ml_varsv uses `Y0(end-4+1:end,:)` -
  different data enter the AR(4) fits when p differs from 4; ml_varsv adds a U_hat output.
- **`get_resid_var.m` vs `get_resid_var_v2.m`**: the `_v2` (HYB) regresses on 4 lags of ALL
  variables plus a 1e-4 ridge - numerically different sig2, hence different Minnesota scalings.
- **`getVtheta.m`**: HYB copy hard-codes kappa3=.2, kappa4=1 inside the body; MAHP copy takes
  them from the kappa vector.
- **`sample_SVpara.m`** (DELISTED 2026-09-03, step 9 - kept here only so the earlier verdict is
  not re-derived): the ml_varsv and OISV copies differ in (i) phi truncation bound .998
  (ml_varsv) vs .999 (OISV); (ii) OISV takes h with n+r columns and demeans only the first n,
  ml_varsv demeans all columns; (iii) OISV's mu block indexes phi(1:n)/sig2(1:n), ml_varsv uses
  full vectors. BOTH copies carry the `if mu ~= 0` vectorized gate (2026-09-01 correction) - the
  gate is not what separates them. Differences (ii) and (iii) are no-ops at r = 0, which is what
  every ml_varsv call site passes, so `bvar.sv.sv_params(h,mu,phi,Hyper,.998)` reproduces the
  ml_varsv copy draw-for-draw (`test_sv_params_mlvarsv`, plus end-to-end in
  `test_mlvarsv_equivalence`). Only difference (i) is real in practice; it is a parameter, not a
  fork. OISV additionally splits the zero-mean case into `sample_SV0para.m` with bound .99
  (`bvar.sv.sv0_params`) - that one is a genuinely separate function. Do NOT merge the .998 and
  .999 DEFAULTS: the bound changes draws (teeth-verified).
- **`macrodata_Q_2018Q4.csv`**: byte-identical between MAHP and HYB but a DIFFERENT file in
  BVAR_ACP (md5-verified). Never key a shared data folder by this filename.
- **MAHP `forecast_BVAR_NG.m` kappa/psi block**: NOT reproduced by
  `bvar.samplers.gig_shrinkage('ng',...)` at any psi_floor - its conditionals carry an extra
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
- **OISV OI alpha step, estimation vs forecast**: `bvar.samplers.eq_svar_oi` reproduces ONLY
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
- **`eq_var_redu_tri` vs `eq_svar_oi`**: siblings, not duplicates. `eq_svar_oi` (OISV
  SVARSV_MH.m 76-87) stacks ALL n rows of the rotated system - `yi = vec((Y-X*A)*B0')`,
  `Wi = kron(B0(:,ii),X)` - and has no prior-mean term. `eq_var_redu_tri` (ml_varsv
  VAR_ARSV_redu.m 44-57) stacks only rows `ii:n` (valid because B0 is unit lower triangular,
  so rows 1..ii-1 carry B0(j,ii) = 0 and contribute nothing) and adds `iValpi*alpi0`. Same
  conditional posterior in the triangular case, different dimensions and floating-point path;
  and the OISV B0 is NOT triangular, so its truncation would be wrong there. Keep separate.
- **`anormrnd.m` vs `tnormrnd.m`**: anormrnd is the OISV bimodal two-component draw for the
  first B0 rotation coordinate (one rand + one randn); tnormrnd is an inverse-cdf truncated
  normal. Same "restricted normal draw" vibe, entirely different densities and rng
  sequences - never unify.
- **`horseshoe_kappa_psi` vs `gig_shrinkage`**: the OISV horseshoe block (inverse-gamma /
  auxiliary z draws, kappa(1:2) global scales) and the MAHP normal-gamma (GIG) block are
  different prior families with different draw sequences - never unify.
- **`bvar.ml.kron_*` vs `bvar.ml.mlvarsv_*`**: two different ESTIMATORS sharing one package
  folder. `kron_*` is Chib's method (a starred point, prior and posterior ordinates, an
  integrated likelihood by importance sampling, no standard error); `mlvarsv_*` is adaptive
  importance sampling over the whole parameter vector with densities fitted to the posterior
  draws, reporting a 50-batch standard error. They share only scalar densities. Do not unify
  them behind one interface or one options struct - `'bugcompat'` denotes different defects in
  each (kron: leftover-workspace evaluation points; mlvarsv: the VAR-SVO outlier block).
- **`lmvnpdf_pcn` vs `lmvnpdf`**: precision vs covariance parameterization of the same
  Gaussian. Only the precision form is extracted (`bvar.ml.lmvnpdf_pcn`, which is what
  `getISden_ARSS` hands it); the legacy `utility/lmvnpdf.m` has no caller anywhere. If it is
  ever extracted, keep both names - handing a covariance to the precision form silently
  returns a different number.

## Verification notes (step 8 self-check, 2026-09-02)

- Teeth check passed: a +1e-7 perturbation of the flat-nu lpri constant in a scratch-mirror
  copy of `bvar.ml.kron_bvar_t` makes `test_kron_equivalence` FAIL with "model 2: ml lpri
  differs"; the real tree was never perturbed and its suite is green.
- All EIGHT models are covered end-to-end bitwise (estimation + ML on one stream, terminal
  rng state included) - the fallback subset (affected paths + models 1/3/8) was not needed:
  measured runtimes m8 ~34 s/side and m7 ~26 s/side put the whole test at ~137 s, inside the
  ~4-minute budget (test_kron_intlike + test_kron_ml_densities add ~5 s).
- Corrected-mode teeth INSIDE the test: model 4 corrected moves ONLY llike (lpri/lpost/
  den_psi bitwise unchanged; llike and ML move); model 8 corrected leaves llike (the
  intlike, drawn first at the same stream position) and lpost(2:3) bitwise unchanged while
  lpost(1) and ML move. At the test size the model-8 ML moves -8492.1 -> -8495.8; the
  golden-manifest published-run reference -8468.4 is the bugcompat path (manifest row:
  "Do NOT fix in the golden run").
- rng-decision nuance carried from step 4: the root sample_nu uses the normpdf-form MH
  ratio while `bvar.sv.nu_studentt` uses the log-form - mathematically identical, bitwise
  different in the acceptance scalar; decisions (and hence all draws) verified identical at
  the test seed here and over test_nu_studentt's 500 sweeps, but a knife-edge flip at some
  other seed is not provably impossible. Any future failure of test_kron_equivalence on
  nu-bearing stores (models 2/5/6/8) should suspect this first.
- Path nuance (same shape as steps 5-7): inside test_kron_equivalence the legacy side
  resolves tempdir copies (construct_prior_A, sample_h, sample_nu, llike_MA, the ROOT
  llike_CSV_MA - asserted - lniwpdf, linvgammpdf, intlike_*); the functionized side
  resolves only qualified bvar.* names plus the package's run_all/run_ml/preset (resolution
  pinned by a `which` assertion against the two other packages that define run_all).

## Verification notes (step 5 adversarial review, 2026-09-01)

- The perturbation "teeth" check passed: altering a single preset constant (sv_offset) in a
  scratch mirror makes test_mahp_equivalence FAIL on the stored draws - the equivalence test
  detects one-constant deviations.
- Path nuance in test_mahp_equivalence: inside the test, the unqualified gigrnd call in
  bvar.samplers.gig_shrinkage resolves to the tempdir copy of the LEGACY gigrnd.m, not
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

- Teeth check passed: a 1e-7 relative perturbation of one hierarchical-shrinkage constant in
  `bvar.samplers.horseshoe_kappa_psi` makes `test_oisv_equivalence` FAIL on store_kappa;
  suite green on revert.
- The six "textually identical modulo renaming" claims of the step-7 table (B0 block,
  CS B block, CS alp block, and the hierarchical shrinkage block across all four scripts) were verified
  mechanically: comment-stripped, whitespace-normalized diffs with Yt/Xt/Tt -> Y/X/T and
  alpha/beta -> theta come back empty; the OI alpha estimation-vs-forecast diff is
  NON-empty, confirming the never-merge entry's direction.
- Record corrected in `bvar.sv.sv_params`'s header: the step-4 note "in the OISV CS_MH run
  mu is initialized at zero, so mu is never updated there" was WRONG - CS_MH.m lines 34-37
  initialize mu(ii) = mean(log(s2i)) from the data (almost surely all-nonzero), so the
  `if mu~=0` gate passes and mu IS updated every sweep. The equivalence test's store_hpara
  comparison (columns 1:n are the mu draws) covers it draw-for-draw.
- Path nuance in test_oisv_equivalence (same shape as steps 5-6): the legacy side resolves
  the tempdir copies of anormrnd/vec/sample_SV*/get_*/getVbeta/construct_Sigt, the core
  side resolves only qualified bvar.* names - no unqualified name in any step-7 core
  function can be shadowed by the tempdir. Both packages define a `run_all`; the test
  pins resolution to the OISV package with a `which` assertion.
- CS shape quirks reproduced by construction: legacy `alp` is a 1 x k_alp ROW (dynamic
  growth) and `alp_tri_cs` returns the same row shape; `z_kappa` enters sweep 1 as a 2x1
  column and leaves every call to the hierarchical shrinkage block as a 1x2 row - kept verbatim. The equivalence
  test's isequal locks the STORES and terminal rng state (values); these internal
  orientations are draw-irrelevant (orientation-agnostic indexing) and are not
  independently asserted. Also: the test pairs OI-with-default and CS-with-flipped
  ordering; the cross pairings are unexercised end-to-end (flip only reverses var_id
  upstream of the model switch - negligible risk, swap the pairings once if paranoid).

## Verification notes (step 8 adversarial review, 2026-09-02)

- Verdict EQUIVALENT. The verifier independently re-derived the bug audit from the legacy
  sources (all 7 ml_* scripts, 7 estimation scripts, 4 intlike evaluators, both
  llike_CSV_MA copies): all three defects confirmed at the cited lines, all six clean
  bills confirmed, no additional defects found. Its teeth check perturbed a DIFFERENT
  constant than the builder's (the bugcompat leftover-Sig consumption) and the equivalence
  test failed on exactly `model 8: ml lpost differs` - proving defect 3 is materially
  separate from defect 1.
- **What "corrected" does and does not mean.** The legacy ordinate decomposition mixes
  full-run marginal ordinates (e.g. p(sigh2*|Y), p(nu*|Y)) with reduced-run conditional
  ordinates, and models 7/8 take den_rho and den_psi from ONE joint reduced run rather
  than a conditional telescoping - so it is not an exact Chib decomposition. That
  structural approximation is shared bitwise by legacy, bugcompat AND corrected modes in
  all six multi-block models; correcting it was out of scope. "Corrected" here means the
  three evaluation-point/leftover-workspace defects are fixed, not that the estimator is
  the exact Chib construction.
- The four intlike evaluators are executably verbatim, not byte-verbatim: differences are
  trailing whitespace, one stray semicolon after `while errh> 10^(-3);`, and `[T n]` ->
  `[T, n]`. Bitwise outputs and terminal rng state are asserted by `test_kron_intlike`.
- Toolbox note: models 4/6/7/8 need the Optimization Toolbox for `fminunc` (on every
  optimizer path); `fminbnd` is base MATLAB.
- Grid sizes (preset `pr.ml.ngrid`): 700-point psi grids for models 2/3/4/5, 300 for
  models 6/7, 299 for model 8.

## Ranking adjudication (step 8, 2026-09-02)

Full-length runs (nsims = 30000, burnin = 5000, full sample), both marginal likelihoods
computed from the SAME chain (rng state saved after estimation and restored before each ML
call), two seeds per affected model; log
`tests/golden/chan2020_jbes_kronecker/ml_bugcompat_comparison_20260902/`.

| Model | seed | bugcompat | corrected | delta |
|---|---|---|---|---|
| BVAR-MA (4) | 20260902 | -8703.3074 | -8703.4840 | -0.1766 |
| BVAR-MA (4) | 8177 | -8703.5934 | -8703.4922 | +0.1013 |
| BVAR-CSV-t-MA (8) | 20260902 | -8469.1532 | -8471.6073 | -2.4542 |
| BVAR-CSV-t-MA (8) | 8177 | -8473.3639 | -8473.2296 | +0.1344 |
| BVAR-CSV-MA (7), control | 20260902 | -8484.8318 | -8484.8318 | 0.0000 (bitwise) |

**The published ranking is unaffected.** The corrections are at most 2.45 log points and
straddle zero across seeds, while BVAR-CSV-t-MA leads the second-best model (BVAR-CSV-MA)
by ~17.5 points and BVAR-MA sits ~24 points above BVAR and ~161 below BVAR-t. Both
corrections are also smaller than the estimator's own seed-to-seed Monte Carlo spread
(model 8's bugcompat ML varies by 4.2 points across the two seeds), i.e. the defects move
these marginal likelihoods by less than the noise already inherent in reporting them. The
model-7 control confirms the flag is a bitwise no-op where no defect exists.
