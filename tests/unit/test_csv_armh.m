function test_csv_armh
% seeded draw-for-draw equivalence of bvt.sv.csv_armh with all four legacy copies:
% sample_CSV (ml_varsv, canonical) and the three identical sample_h copies.
% Realistic harness: simulate a CSV path, form s2, run repeated MCMC sweeps.
root = getappdata(0, 'bvt_repo_root');

% --- simulate data once ---
rng(1, 'twister');
T = 120; n = 15; rho = 0.95; sigh2 = 0.05;
htrue = zeros(T,1); htrue(1) = sqrt(sigh2/(1-rho^2))*randn;
for t = 2:T, htrue(t) = rho*htrue(t-1) + sqrt(sigh2)*randn; end
s2 = exp(htrue).*chi2rnd(n,T,1);   % sum of n squared standardized errors
h0 = log(s2/n);
nrep = 50;

% --- core: default flag (sample_h behavior), and explicit false ---
rng(42, 'twister'); h = h0; Hc = zeros(T,nrep); ac = zeros(nrep,1);
for i = 1:nrep, [h,a] = bvt.sv.csv_armh(s2,rho,sigh2,h,n); Hc(:,i) = h; ac(i) = a; end
sc = rng;
rng(42, 'twister'); h = h0; Hf = zeros(T,nrep); af = zeros(nrep,1);
for i = 1:nrep, [h,a] = bvt.sv.csv_armh(s2,rho,sigh2,h,n,false); Hf(:,i) = h; af(i) = a; end
assert(isequal(Hc,Hf) && isequal(ac,af), 'csv_armh: omitted flag differs from explicit false');
assert(any(ac == 1) && numel(unique(Hc(1,:))) > 1, 'csv_armh: harness produced no accepted moves');

% --- legacy comparisons, one folder on the path at a time ---
legs = { fullfile(root,'replications','chan2023_joe_mlvarsv','legacy','utility'),  'sample_CSV'; ...
         fullfile(root,'replications','chan2020_jbes_kronecker','legacy'),          'sample_h'; ...
         fullfile(root,'replications','chan2020_jbes_kronecker','legacy','realtime_forecasts'), 'sample_h'; ...
         fullfile(root,'replications','chan2020_springer_largebvar','legacy'),      'sample_h' };
for k = 1:size(legs,1)
    addpath(legs{k,1}); c = onCleanup(@() rmpath(legs{k,1}));
    fn = str2func(legs{k,2});
    rng(42, 'twister'); h = h0; Hl = zeros(T,nrep); al = zeros(nrep,1);
    for i = 1:nrep, [h,a] = fn(s2,rho,sigh2,h,n); Hl(:,i) = h; al(i) = a; end
    sl = rng;
    assert(isequal(Hl,Hc) && isequal(al,ac), ...
        'csv_armh: differs from legacy %s in %s', legs{k,2}, legs{k,1});
    assert(isequal(sl.State, sc.State), ...
        'csv_armh: rng call sequence differs from legacy %s in %s', legs{k,2}, legs{k,1});
    clear c   % rmpath now, before the next folder shadows the same name
end
end
