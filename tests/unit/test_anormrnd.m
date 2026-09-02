function test_anormrnd
% exact draw-for-draw equivalence of bvt.util.anormrnd with the legacy OISV
% copy under one seed, across (mu,rho) pairs hitting both mixture branches
% (rho = 1/T with T = 706 is the value the B0 row sampler passes), plus a
% bimodality sanity check at mu = 0 (w = .5: both signs must appear).
root = getappdata(0, 'bvt_repo_root');
leg = fullfile(root, 'replications', 'chan_koop_yu2024_jbes_oisv', 'legacy', 'utility');
addpath(leg); c = onCleanup(@() rmpath(leg));

mus = [-3.2, -0.4, 0, 0.7, 5.1];
rhos = [1/706, .01, 1];
ndraw = 200;

d1 = zeros(ndraw, numel(mus)*numel(rhos));
rng(20260902, 'twister');
col = 0;
for mu = mus
    for rho = rhos
        col = col + 1;
        for i = 1:ndraw
            d1(i, col) = bvt.util.anormrnd(mu, rho);
        end
    end
end

d0 = zeros(ndraw, numel(mus)*numel(rhos));
rng(20260902, 'twister');
col = 0;
for mu = mus
    for rho = rhos
        col = col + 1;
        for i = 1:ndraw
            d0(i, col) = anormrnd(mu, rho);     % legacy (path-shadowed name)
        end
    end
end
assert(isequal(d0, d1), 'anormrnd: differs from legacy under the same seed');

zero_cols = d1(:, (find(mus == 0, 1)-1)*numel(rhos) + (1:numel(rhos)));
assert(any(zero_cols(:) > 0) && any(zero_cols(:) < 0), ...
    'anormrnd: mu = 0 draws must visit both modes');
end
