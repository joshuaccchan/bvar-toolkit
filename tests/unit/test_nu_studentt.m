function test_nu_studentt
% seeded draw-for-draw equivalence of bvt.sv.nu_studentt with all three legacy
% sample_nu copies (springer canonical; kronecker root uses the normpdf-ratio
% acceptance arithmetic - decisions must still match draw-for-draw).
root = getappdata(0, 'bvt_repo_root');

% latent scales lam_t ~ IG(nu/2, nu/2): heavier-tail regime so the MH step
% both accepts and rejects (acceptance ~0.9 over these seeds)
rng(3, 'twister');
lam = 1./gamrnd(2.5/2, 2/2.5, 60, 1);
nu_ub = 50; nrep = 500;

rng(11, 'twister'); nu = 8; Nc = zeros(nrep,1); Fc = zeros(nrep,1);
for i = 1:nrep, [nu,f,f_nu] = bvt.sv.nu_studentt(lam,nu,nu_ub); Nc(i) = nu; Fc(i) = f; end
sc = rng;
assert(any(Fc == 1) && any(Fc == 0), 'nu_studentt: harness did not exercise both MH branches');

legs = { fullfile(root,'replications','chan2020_springer_largebvar','legacy'); ...
         fullfile(root,'replications','chan2020_jbes_kronecker','legacy'); ...
         fullfile(root,'replications','chan2020_jbes_kronecker','legacy','realtime_forecasts') };
for k = 1:numel(legs)
    addpath(legs{k}); c = onCleanup(@() rmpath(legs{k}));
    rng(11, 'twister'); nu = 8; Nl = zeros(nrep,1); Fl = zeros(nrep,1);
    for i = 1:nrep, [nu,f,f_nu_l] = sample_nu(lam,nu,nu_ub); Nl(i) = nu; Fl(i) = f; end
    sl = rng;
    assert(isequal(Nl,Nc) && isequal(Fl,Fc), 'nu_studentt: differs from legacy in %s', legs{k});
    assert(isequal(sl.State, sc.State), 'nu_studentt: rng call sequence differs from legacy in %s', legs{k});
    % third output: the log target kernel handle must evaluate identically
    xs = [2.5 5 12 30];
    assert(isequal(f_nu(xs), f_nu_l(xs)), 'nu_studentt: f_nu handle evaluates differently vs %s', legs{k});
    clear c   % rmpath now, before the next folder shadows the same name
end
end
