function test_minnesota_C
% deterministic function: exact equivalence with all four legacy get_C copies
rng(43, 'twister');
n = 5; p = 3;
sig2 = exp(randn(n, 1));
[C0, i10, i20] = bvar.priors.minnesota_C(n, p, sig2);

root = getappdata(0, 'bvar_repo_root');
folders = { ...
    fullfile(root, 'replications', 'chan2021_ijf_mahp', 'legacy'); ...
    fullfile(root, 'replications', 'chan2023_jbes_hybtvp', 'legacy', 'utility'); ...
    fullfile(root, 'replications', 'chan2023_joe_mlvarsv', 'legacy', 'utility'); ...
    fullfile(root, 'replications', 'chan_koop_yu2024_jbes_oisv', 'legacy', 'utility')};
for jj = 1:numel(folders)
    addpath(folders{jj});
    c = onCleanup(@() rmpath(folders{jj}));
    [C1, i11, i21] = get_C(n, p, sig2);
    assert(isequal(C1, C0) && isequal(i11, i10) && isequal(i21, i20), ...
        'minnesota_C: differs from legacy get_C in %s', folders{jj});
    clear c
end
end
