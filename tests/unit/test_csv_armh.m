function test_csv_armh
% Two standards, because this function has been improved beyond its legacy copies.
%
% (1) EQUIVALENCE, for the default path: seeded draw-for-draw agreement with all
%     four legacy copies - sample_CSV (ml_varsv, canonical) and the three
%     identical sample_h copies. Realistic harness: simulate a CSV path, form s2,
%     run repeated MCMC sweeps.
% (2) CORRECTNESS, for the options the legacy does not have: a Geweke
%     joint-distribution test. The block updates h given s2 in a model where both
%     conditionals are known exactly - h is a zero-mean AR(1) with
%     h_1 ~ N(0, sigh2/(1-rho^2)), and s2_t | h_t = exp(h_t)*chi2(n), which is the
%     likelihood the code implements - so the joint (h,s2) can be simulated two
%     ways. Drawing h from the prior and then s2 given h must give the same joint
%     as alternating s2 given h with this kernel. If it does, the kernel has the
%     right invariant distribution.
%
%     WHAT THIS DETECTS, measured rather than assumed. Giving the kernel a rho 2%
%     away from the one the data came from scores max|z| about 6, and 4% about 9,
%     against a threshold of 4; a wrong n scores 115. So it catches an error of a
%     few percent in a conditional, and it will NOT catch an arbitrarily small
%     bias - an earlier version of this function exposed the mode-search tolerance
%     as an option, and the bias that introduced needed a 3-variable target with
%     quadrature ground truth and 60,000 draws to see. A pass here is evidence at
%     that sensitivity, not a proof.
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

% ---------------------------------------------------------------------------
%  Correctness of the added options
% ---------------------------------------------------------------------------
Tg = 40; ng = 5; rg = .95; sg = .04; M = 4000; burn = 200; nb = 40;

    % c_reject must not change the distribution the kernel targets, only how
    % often the AR step accepts. 0.2 is far BELOW any valid envelope constant,
    % so it is the case where the MH correction has to do the work.
for c_reject = [0.2 1 3 20]
    z = geweke_z(Tg,ng,ng,rg,sg,M,burn,nb,{'c_reject',c_reject},false);
    assert(max(abs(z)) < 4, ...
        'csv_armh: joint distribution moved at c_reject = %g (max|z| = %.2f)', ...
        c_reject, max(abs(z)));
end

    % the test must be able to fail. Feed the kernel a DIFFERENT n from the one
    % the data were generated with, so its conditional is genuinely wrong;
    % changing both would just be a different, self-consistent model.
z = geweke_z(Tg,ng,ng+1,rg,sg,M,burn,nb,{},false);
assert(max(abs(z)) > 20, ...
    'csv_armh: the Geweke check did not detect a deliberately wrong n (max|z| = %.2f)', ...
    max(abs(z)));

    % is_ForcedAccept skips the MH correction. At the default c_reject the
    % envelope holds, so the AR output is already an exact draw and forcing is
    % harmless - which is why the legacy uses it at initialization. Push
    % c_reject below the bound and the MH step is what is holding the kernel
    % together, so forcing must break it.
z = geweke_z(Tg,ng,ng,rg,sg,M,burn,nb,{}, true);
assert(max(abs(z)) < 4, ...
    'csv_armh: forced accept at a valid envelope should still be exact (max|z| = %.2f)', ...
    max(abs(z)));
z = geweke_z(Tg,ng,ng,rg,sg,M,burn,nb,{'c_reject',0.2}, true);
assert(max(abs(z)) > 8, ...
    'csv_armh: forced accept below the envelope bound should NOT be exact (max|z| = %.2f)', ...
    max(abs(z)));

    % both loop caps raise rather than returning a draw that is not from the target
try
    bvar.sv.csv_armh(s2,rho,sigh2,h0,n,'MaxIterMode',1);
    error('csv_armh:testFailed', 'MaxIterMode did not raise');
catch err
    assert(strcmp(err.identifier,'bvar:sv:csv_armh:modeNotConverged'), ...
        'csv_armh: wrong identifier for the mode cap: %s', err.identifier);
end
try
    bvar.sv.csv_armh(s2,rho,sigh2,h0,n,'c_reject',1e12,'MaxIterAR',5);
    error('csv_armh:testFailed', 'MaxIterAR did not raise');
catch err
    assert(strcmp(err.identifier,'bvar:sv:csv_armh:arNotAccepted'), ...
        'csv_armh: wrong identifier for the AR cap: %s', err.identifier);
end
try
    bvar.sv.csv_armh(s2,rho,sigh2,h0,n,'ARbound',3);
    error('csv_armh:testFailed', 'an unknown option did not raise');
catch err
    assert(strcmp(err.identifier,'bvar:sv:csv_armh:badOption'), ...
        'csv_armh: wrong identifier for an unknown option: %s', err.identifier);
end

    % Values that used to fail SILENTLY, which is worse than failing loudly:
    % c_reject = 0 froze the chain (log(0) makes alpMH NaN, so nothing is ever
    % accepted and h comes back unchanged with is_accept = 0, indistinguishable
    % from a legitimately sticky chain), and a negative c_reject was treated as
    % its magnitude because log() of it is complex and MATLAB's relational
    % operators compare real parts only. Both are now rejected.
bad = {{'c_reject',0}, {'c_reject',-3}, {'c_reject',Inf}, {'c_reject',NaN}, ...
       {'c_reject',[1 10]}, {'MaxIterAR',0}, {'MaxIterAR',Inf}, {'MaxIterMode',2.5}};
for ib = 1:numel(bad)
    try
        bvar.sv.csv_armh(s2,rho,sigh2,h0,n,bad{ib}{:});
        error('csv_armh:testFailed', 'csv_armh accepted %s = %s', ...
            bad{ib}{1}, mat2str(bad{ib}{2}));
    catch err
        assert(strcmp(err.identifier,'bvar:sv:csv_armh:badArgument'), ...
            'csv_armh: %s = %s gave %s, expected badArgument', ...
            bad{ib}{1}, mat2str(bad{ib}{2}), err.identifier);
    end
end

    % A vector in position 6 is ht_start supplied without the flag. It used to
    % be read as is_ForcedAccept and then crash intermittently, on whichever
    % sweep the MH ratio first rejected; it now fails immediately and says so.
try
    bvar.sv.csv_armh(s2,rho,sigh2,h0,n,h0);
    error('csv_armh:testFailed', 'a vector is_ForcedAccept did not raise');
catch err
    assert(strcmp(err.identifier,'bvar:sv:csv_armh:badArgument'), ...
        'csv_armh: vector in position 6 gave %s, expected badArgument', err.identifier);
end

    % the parser must not mistake a positional optional for an option name
rng(42,'twister'); p1 = bvar.sv.csv_armh(s2,rho,sigh2,h0,n,false,h0,'c_reject',3);
rng(42,'twister'); p2 = bvar.sv.csv_armh(s2,rho,sigh2,h0,n,false,h0);
assert(isequal(p1,p2), 'csv_armh: name-value after the positional args changed the default path');
rng(42,'twister'); p3 = bvar.sv.csv_armh(s2,rho,sigh2,h0,n,'c_reject',3);
rng(42,'twister'); p4 = bvar.sv.csv_armh(s2,rho,sigh2,h0,n);
assert(isequal(p3,p4), 'csv_armh: options in position 6 were misread as is_ForcedAccept');
end

% ---------------------------------------------------------------------------
function z = geweke_z(T,n_gen,n_ker,rho,sigh2,M,burn,nb,opts,forced)
% z-statistics comparing four functionals of h between the marginal-conditional
% simulator (h from the prior, s2 given h) and the successive-conditional one
% (s2 given h, then h through the kernel). Batch means give the chain's variance.
rng(20260910,'twister');
A = zeros(M,4);
for i = 1:M, A(i,:) = geweke_funcs(ar1_prior(T,rho,sigh2)); end

rng(20260910,'twister');
h = ar1_prior(T,rho,sigh2);
B = zeros(M,4);
for i = 1:(M+burn)
    s2 = exp(h).*chi2rnd(n_gen,T,1);
    h = bvar.sv.csv_armh(s2,rho,sigh2,h,n_ker,forced,h,opts{:});
    if i > burn, B(i-burn,:) = geweke_funcs(h); end
end

mA = mean(A); vA = var(A)/M;
K = floor(M/nb);
bm = squeeze(mean(reshape(B(1:K*nb,:),nb,K,[]),1));
z = (mA - mean(bm))./sqrt(vA + var(bm)/K);
end

function h = ar1_prior(T,rho,sigh2)
h = zeros(T,1);
h(1) = sqrt(sigh2/(1-rho^2))*randn;
for t = 2:T, h(t) = rho*h(t-1) + sqrt(sigh2)*randn; end
end

function f = geweke_funcs(h)
f = [mean(h), h(1), h(end), mean(h.^2)];
end
