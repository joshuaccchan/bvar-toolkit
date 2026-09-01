function test_ksc_ar1_mean
% bvt.sv.ksc_ar1_mean must reproduce, draw-for-draw under one seed, BOTH legacy
% sample_SV.m copies it canonicalizes (ml_varsv and OISV) - path AND indicators.
root = getappdata(0, 'bvt_repo_root');

rng(43, 'twister');                          % fixed test data
T = 173;
e = randn(T,1).*exp(0.3*randn(T,1));
ystar = log(e.^2 + 1e-4);
h_in = -1 + 0.2*randn(T,1);
mu = -1.2; rho = 0.95; sig2 = 0.15;

rng(9, 'twister');
[h_core, S_core] = bvt.sv.ksc_ar1_mean(ystar, h_in, mu, rho, sig2);

copies = { ...
    {'chan2023_joe_mlvarsv', 'legacy', 'utility'}; ...
    {'chan_koop_yu2024_jbes_oisv', 'legacy', 'utility'}};

for ii = 1:numel(copies)
    leg = fullfile(root, 'replications', copies{ii}{:});
    addpath(leg); c = onCleanup(@() rmpath(leg));  % prepended -> legacy copy shadows
    rng(9, 'twister');
    [h_leg, S_leg] = sample_SV(ystar, h_in, mu, rho, sig2);
    clear c                                   % rmpath before the next copy
    assert(isequal(h_leg, h_core), ...
        'ksc_ar1_mean: h differs from legacy sample_SV in %s', copies{ii}{1});
    assert(isequal(S_leg, S_core), ...
        'ksc_ar1_mean: indicators differ from legacy sample_SV in %s', copies{ii}{1});
end
assert(all(ismember(S_core, 1:7)), 'ksc_ar1_mean: indicators outside 1..7');
end
