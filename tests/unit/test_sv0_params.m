function test_sv0_params
% seeded draw-for-draw equivalence with the legacy OISV sample_SV0para under defaults
root = getappdata(0, 'bvt_repo_root');
leg = fullfile(root, 'replications', 'chan_koop_yu2024_jbes_oisv', 'legacy', 'utility');
addpath(leg); c = onCleanup(@() rmpath(leg));

T = 70; n = 4;
rng(51, 'twister');
h = 0.5*randn(T, n);
phi_in = 0.9*ones(n, 1);
Hyper = struct('nuh', 3*ones(n,1), 'Sh', 0.1*ones(n,1), ...
    'phi0', 0.95*ones(n,1), 'Vphi', 0.05^2*ones(n,1));
rng(52, 'twister');
[phi1, sig21, f1] = bvt.sv.sv0_params(h, phi_in, Hyper);
rng(52, 'twister');
[phi2, sig22, f2] = sample_SV0para(h, phi_in, Hyper);
assert(isequal(phi2,phi1) && isequal(sig22,sig21) && isequal(f2,f1), ...
    'sv0_params: differs from legacy sample_SV0para');
assert(any(f1 == 1), 'sv0_params: degenerate test, no phi candidate ever accepted');

% explicit default bound gives the identical draw
rng(52, 'twister');
[phi3, sig23, f3] = bvt.sv.sv0_params(h, phi_in, Hyper, .99);
assert(isequal(phi3,phi1) && isequal(sig23,sig21) && isequal(f3,f1), ...
    'sv0_params: explicit phi_bnd=.99 differs from default');
end
