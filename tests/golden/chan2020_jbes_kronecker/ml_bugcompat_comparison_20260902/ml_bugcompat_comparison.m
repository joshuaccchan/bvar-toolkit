% ml_bugcompat_comparison - ranking adjudication for chan2020_jbes_kronecker
% (Chan 2020, JBES 38(1): 68-79), step-8 follow-up, 2026-09-02.
%
% Question: does correcting the two known legacy ML defects (ml_BVAR_MA.m
% line 17 leftover-psi llike term; ml_BVAR_CSV_t_MA.m lines 42-44 frozen
% Hpsi/psi ordinate loop + line 108 leftover-Sig psi target) change the
% paper's model ranking?
%
% Design: for each AFFECTED model (4 = BVAR-MA, 8 = BVAR-CSV-t-MA), ONE
% full-length estimation run (nsims = 30000, burnin = 5000, full data_Q.csv
% sample - the golden-capture settings), then BOTH marginal-likelihood
% computations from the SAME chain: the rng state is saved at the end of
% estimation and restored before each ML call, so bugcompat and corrected
% ordinates consume identical random streams and the delta is purely the
% evaluation-point fix. Two seeds per affected model bound the MC error of
% the delta. Control: model 7 (BVAR-CSV-MA, clean bill) run twice through
% run_ml with bugcompat true/false from the same seed - the two MLs must be
% bitwise identical.
%
% Golden cross-check targets (clock-seeded legacy full runs, 2026-09-01;
% expect MC-error agreement, not bitwise): BVAR -8727.1, t -8542.2,
% CSV -8514.1, MA -8703.3, t-CSV -8500.5, t-MA -8522.2, CSV-MA -8485.9,
% CSV-t-MA -8468.4.

repdir = 'C:\Users\joshu\Dropbox\website\Github\bvar-toolkit\replications\chan2020_jbes_kronecker';
addpath(repdir);

nsim   = 30000;
burnin = 5000;
seeds  = [20260902 8177];

fprintf('=== ml_bugcompat_comparison: nsim=%d burnin=%d seeds=[%d %d] ===\n\n', ...
    nsim, burnin, seeds(1), seeds(2));

results = struct('model', {}, 'seed', {}, 'ML_bug', {}, 'ML_corr', {}, ...
    't_est', {}, 't_mlb', {}, 't_mlc', {});

for seed = seeds
    for model = [4 8]
        t0 = tic;
        est = run_all(model, nsim, burnin, seed);
        t_est = toc(t0);
        s = rng;    % stream state at end of estimation: both ML modes start here

        rng(s); t1 = tic;
        switch model
            case 4
                [ML_bug, db] = bvt.ml.kron_bvar_ma(est.shortY, est.X, est.pri, est, ...
                    'bugcompat', true);
            case 8
                [ML_bug, db] = bvt.ml.kron_bvar_csv_t_ma(est.shortY, est.X, est.pri, est, ...
                    'bugcompat', true);
        end
        t_mlb = toc(t1);

        rng(s); t2 = tic;
        switch model
            case 4
                [ML_corr, dc] = bvt.ml.kron_bvar_ma(est.shortY, est.X, est.pri, est);
            case 8
                [ML_corr, dc] = bvt.ml.kron_bvar_csv_t_ma(est.shortY, est.X, est.pri, est);
        end
        t_mlc = toc(t2);

        fprintf('model %d, seed %d:\n', model, seed);
        fprintf('  bugcompat  ML = %.4f   (llike %.4f, lpri %.4f, sum lpost %.4f)\n', ...
            ML_bug, db.llike, sum(db.lpri), sum(db.lpost));
        fprintf('  corrected  ML = %.4f   (llike %.4f, lpri %.4f, sum lpost %.4f)\n', ...
            ML_corr, dc.llike, sum(dc.lpri), sum(dc.lpost));
        fprintf('  delta (corrected - bugcompat) = %.4f\n', ML_corr - ML_bug);
        fprintf('  timings: est %.1fs, ml bugcompat %.1fs, ml corrected %.1fs\n\n', ...
            t_est, t_mlb, t_mlc);

        results(end+1) = struct('model', model, 'seed', seed, ...
            'ML_bug', ML_bug, 'ML_corr', ML_corr, ...
            't_est', t_est, 't_mlb', t_mlb, 't_mlc', t_mlc); %#ok<SAGROW>
    end
end

    % control: model 7 (clean bill) - the bugcompat flag must be a no-op
fprintf('=== control: model 7 (BVAR-CSV-MA), run_ml with bugcompat true vs false, seed %d ===\n', seeds(1));
t0 = tic;
o1 = run_ml(7, nsim, burnin, seeds(1), 'bugcompat', true);
t7a = toc(t0);
t0 = tic;
o2 = run_ml(7, nsim, burnin, seeds(1), 'bugcompat', false);
t7b = toc(t0);
fprintf('control model 7: bugcompat ML = %.4f, corrected ML = %.4f, bitwise equal: %d\n', ...
    o1.ML, o2.ML, isequal(o1.ML, o2.ML));
fprintf('control model 7 full-output bitwise equal (ml detail struct): %d\n', ...
    isequal(o1.ml, o2.ml));
fprintf('control timings: %.1fs, %.1fs\n\n', t7a, t7b);
assert(isequal(o1.ML, o2.ML), 'control model 7: bugcompat flag changed the ML');

    % summary table
fprintf('=== summary (full length, nsims=30000, burnin=5000) ===\n');
fprintf('%-14s %8s %14s %14s %10s\n', 'model', 'seed', 'bugcompat', 'corrected', 'delta');
mnames = containers.Map({4, 7, 8}, {'BVAR-MA', 'BVAR-CSV-MA', 'BVAR-CSV-t-MA'});
for ir = 1:numel(results)
    r = results(ir);
    fprintf('%-14s %8d %14.4f %14.4f %10.4f\n', mnames(r.model), r.seed, ...
        r.ML_bug, r.ML_corr, r.ML_corr - r.ML_bug);
end
fprintf('%-14s %8d %14.4f %14.4f %10.4f\n', mnames(7), seeds(1), o1.ML, o2.ML, o2.ML - o1.ML);
fprintf('\nGolden (clock-seeded) references: MA -8703.3, CSV-t-MA -8468.4, CSV-MA -8485.9\n');
fprintf('Done.\n');
