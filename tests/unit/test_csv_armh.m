function test_csv_armh
% seeded draw-for-draw equivalence of bvar.sv.csv_armh with all four legacy copies:
% sample_CSV (ml_varsv, canonical) and the three identical sample_h copies.
% Realistic harness: simulate a CSV path, form s2, run repeated MCMC sweeps.
root = getappdata(0, 'bvar_repo_root');

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
for i = 1:nrep, [h,a] = bvar.sv.csv_armh(s2,rho,sigh2,h,n); Hc(:,i) = h; ac(i) = a; end
sc = rng;
rng(42, 'twister'); h = h0; Hf = zeros(T,nrep); af = zeros(nrep,1);
for i = 1:nrep, [h,a] = bvar.sv.csv_armh(s2,rho,sigh2,h,n,false); Hf(:,i) = h; af(i) = a; end
assert(isequal(Hc,Hf) && isequal(ac,af), 'csv_armh: omitted flag differs from explicit false');
assert(any(ac == 1) && numel(unique(Hc(1,:))) > 1, 'csv_armh: harness produced no accepted moves');

% --- step-8 ht_start argument: explicit h must equal the default bitwise;
%     a different NR start changes the proposal path (kron ml_BVAR_CSV's
%     inline reduced-run h step = ht_start h_mean + forced first accept,
%     verified end-to-end by test_kron_equivalence model 3) ---
rng(42, 'twister'); h = h0; Hs = zeros(T,nrep); as = zeros(nrep,1);
for i = 1:nrep, [h,a] = bvar.sv.csv_armh(s2,rho,sigh2,h,n,false,h); Hs(:,i) = h; as(i) = a; end
assert(isequal(Hc,Hs) && isequal(ac,as), 'csv_armh: ht_start = h differs from the default');
rng(42, 'twister');
h_alt = bvar.sv.csv_armh(s2,rho,sigh2,h0,n,true,zeros(T,1));
rng(42, 'twister');
h_def = bvar.sv.csv_armh(s2,rho,sigh2,h0,n,true);
assert(~isequal(h_alt,h_def), 'csv_armh: a different ht_start should change the accepted path');

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
