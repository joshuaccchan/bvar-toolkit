function test_ksc_rw_h0
% bvar.sv.ksc_rw_h0 must reproduce, draw-for-draw under one seed, EVERY legacy
% copy it canonicalizes: springer/kronecker/mahp SVRW.m and HYB sample_SVRW.m.
root = getappdata(0, 'bvar_repo_root');

rng(41, 'twister');                          % fixed test data
T = 137;
e = randn(T,1).*exp(0.3*randn(T,1));
ystar = log(e.^2 + 1e-4);
h_in = 0.2*randn(T,1);
sig2 = 0.12;                                 % state innovation VARIANCE
h0 = 0.4;

rng(7, 'twister');
h_core = bvar.sv.ksc_rw_h0(ystar, h_in, sig2, h0);

copies = { ...
    {'chan2020_springer_largebvar', 'legacy'},                     'SVRW'; ...
    {'chan2020_jbes_kronecker', 'legacy', 'realtime_forecasts'},   'SVRW'; ...
    {'chan2021_ijf_mahp', 'legacy'},                               'SVRW'; ...
    {'chan2023_jbes_hybtvp', 'legacy', 'utility'},                 'sample_SVRW'};

for ii = 1:size(copies, 1)
    leg = fullfile(root, 'replications', copies{ii,1}{:});
    addpath(leg); c = onCleanup(@() rmpath(leg));  % prepended -> legacy copy shadows
    rng(7, 'twister');
    h_leg = feval(copies{ii,2}, ystar, h_in, sig2, h0);
    clear c                                   % rmpath before the next copy
    assert(isequal(h_leg, h_core), ...
        'ksc_rw_h0: differs from legacy %s in %s', copies{ii,2}, copies{ii,1}{1});
end
assert(all(isfinite(h_core)), 'ksc_rw_h0: non-finite log-volatility draw');
end
