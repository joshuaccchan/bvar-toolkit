function test_resid_var_ar4
% deterministic function: exact equivalence with all four legacy get_resid_var copies
rng(41, 'twister');
n = 5;
Y0 = randn(9, n);
Y  = randn(60, n);
s_core = bvar.priors.resid_var_ar4(Y0, Y);

root = getappdata(0, 'bvar_repo_root');
folders = { ...
    fullfile(root, 'replications', 'chan2019wp_acp', 'legacy'); ...
    fullfile(root, 'replications', 'chan2021_ijf_mahp', 'legacy'); ...
    fullfile(root, 'replications', 'chan2022_qe_acp', 'legacy', 'utility'); ...
    fullfile(root, 'replications', 'chan_koop_yu2024_jbes_oisv', 'legacy', 'utility')};
for jj = 1:numel(folders)
    addpath(folders{jj});
    c = onCleanup(@() rmpath(folders{jj}));
    s_leg = get_resid_var(Y0, Y);
    assert(isequal(s_leg, s_core), ...
        'resid_var_ar4: differs from legacy get_resid_var in %s', folders{jj});
    clear c
end
end
