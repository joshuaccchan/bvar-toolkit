function test_vtheta
% deterministic function: exact equivalence with BOTH legacy getVtheta copies
% under the documented kappa settings (MAHP: 4-vector unchanged; HYB: [k1 k2 .2 1]).
rng(44, 'twister');
n = 5; p = 3;
sig2 = exp(randn(n, 1));
[C, idx_kappa1, idx_kappa2] = bvar.priors.minnesota_C(n, p, sig2);
kappa = exp(randn(4, 1)); % generic positive 4-vector (MAHP semantics)

root = getappdata(0, 'bvar_repo_root');

% --- MAHP copy: kappa 4-vector passed unchanged ---
[Va0, Vb0] = bvar.priors.vtheta(idx_kappa1, idx_kappa2, kappa, C, sig2);
leg = fullfile(root, 'replications', 'chan2021_ijf_mahp', 'legacy');
addpath(leg); c = onCleanup(@() rmpath(leg));
[Va1, Vb1] = getVtheta(idx_kappa1, idx_kappa2, kappa, C, sig2);
assert(isequal(Va1, Va0) && isequal(Vb1, Vb0), ...
    'vtheta: differs from legacy MAHP getVtheta');
clear c

% --- HYB copy: hard-codes kappa3=.2, kappa4=1; reproduce with [k1 k2 .2 1] ---
kappa_hyb = kappa(1:2); % HYB callers pass only the first two hyperparameters
[Va0, Vb0] = bvar.priors.vtheta(idx_kappa1, idx_kappa2, [kappa(1); kappa(2); .2; 1], C, sig2);
leg = fullfile(root, 'replications', 'chan2023_jbes_hybtvp', 'legacy', 'utility');
addpath(leg); c = onCleanup(@() rmpath(leg));
[Va1, Vb1] = getVtheta(idx_kappa1, idx_kappa2, kappa_hyb, C, sig2);
assert(isequal(Va1, Va0) && isequal(Vb1, Vb0), ...
    'vtheta: kappa=[k1;k2;.2;1] does not reproduce legacy HYB getVtheta');
clear c
end
