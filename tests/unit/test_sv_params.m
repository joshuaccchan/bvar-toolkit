function test_sv_params
% seeded draw-for-draw equivalence with the legacy OISV sample_SVpara under defaults
root = getappdata(0, 'bvt_repo_root');
leg = fullfile(root, 'replications', 'chan_koop_yu2024_jbes_oisv', 'legacy', 'utility');
addpath(leg); c = onCleanup(@() rmpath(leg));

% --- case 1: r = 0, all-nonzero mu so the mu block runs ---
T = 60; n = 3;
rng(41, 'twister');
h = -1 + 0.3*randn(T, n);
mu_in = [-1; 0.5; 2];
phi_in = 0.95*ones(n, 1);
Hyper = struct('nuh', 3*ones(n,1), 'Sh', 0.1*ones(n,1), ...
    'phi0', 0.95*ones(n,1), 'Vphi', 0.05^2*ones(n,1), ...
    'mu0', zeros(n,1), 'Vmu', 100*ones(n,1));
rng(42, 'twister');
[mu1, phi1, sig21, f1] = bvt.sv.sv_params(h, mu_in, phi_in, Hyper);
rng(42, 'twister');
[mu2, phi2, sig22, f2] = sample_SVpara(h, mu_in, phi_in, Hyper);
assert(isequal(mu2,mu1) && isequal(phi2,phi1) && isequal(sig22,sig21) && isequal(f2,f1), ...
    'sv_params: differs from legacy sample_SVpara (r = 0)');
assert(any(f1 == 1), 'sv_params: degenerate test, no phi candidate ever accepted (r = 0)');

% explicit default bound gives the identical draw
rng(42, 'twister');
[mu3, phi3, sig23, f3] = bvt.sv.sv_params(h, mu_in, phi_in, Hyper, .999);
assert(isequal(mu3,mu1) && isequal(phi3,phi1) && isequal(sig23,sig21) && isequal(f3,f1), ...
    'sv_params: explicit phi_bnd=.999 differs from default');

% --- case 2: r = 2 extra zero-mean columns (h is T x (n+r)) ---
r = 2;
rng(43, 'twister');
h = [-1 + 0.3*randn(T, n), 0.4*randn(T, r)];
phi_in = 0.9*ones(n+r, 1);
Hyper = struct('nuh', 3*ones(n+r,1), 'Sh', 0.1*ones(n+r,1), ...
    'phi0', 0.95*ones(n+r,1), 'Vphi', 0.05^2*ones(n+r,1), ...
    'mu0', zeros(n,1), 'Vmu', 100*ones(n,1));
rng(44, 'twister');
[mu1, phi1, sig21, f1] = bvt.sv.sv_params(h, mu_in, phi_in, Hyper);
rng(44, 'twister');
[mu2, phi2, sig22, f2] = sample_SVpara(h, mu_in, phi_in, Hyper);
assert(isequal(mu2,mu1) && isequal(phi2,phi1) && isequal(sig22,sig21) && isequal(f2,f1), ...
    'sv_params: differs from legacy sample_SVpara (r = 2)');

% --- case 3: mu contains an exact zero -> the mu block is gated off, verbatim ---
mu_in = zeros(n, 1);
rng(45, 'twister');
[mu1, phi1, sig21, f1] = bvt.sv.sv_params(h, mu_in, phi_in, Hyper);
rng(45, 'twister');
[mu2, phi2, sig22, f2] = sample_SVpara(h, mu_in, phi_in, Hyper);
assert(isequal(mu2,mu1) && isequal(phi2,phi1) && isequal(sig22,sig21) && isequal(f2,f1), ...
    'sv_params: differs from legacy sample_SVpara (mu gated at zero)');
assert(isequal(mu1, zeros(n,1)), 'sv_params: mu block ran despite the mu~=0 gate');
end
