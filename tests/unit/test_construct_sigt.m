function test_construct_sigt
% exact equivalence of bvt.structural.construct_Sigt with the legacy OISV
% utility copy on a full impact matrix (OI-style B0) and a unit-lower-
% triangular one (CS-style A). The second legacy copy - the private
% subfunction inside func_main_SVAR_v2.m lines 67-73 - cannot be called
% directly; it is comment-stripped IDENTICAL to the utility copy (diff,
% 2026-09-02) and its numbers are covered end-to-end by the Sig_mean
% assertion in test_oisv_equivalence.
root = getappdata(0, 'bvt_repo_root');
leg = fullfile(root, 'replications', 'chan_koop_yu2024_jbes_oisv', 'legacy', 'utility');
addpath(leg); c = onCleanup(@() rmpath(leg));

rng(11, 'twister');
T = 17; n = 5;
h = .3*randn(T, n);
B0 = eye(n) + .2*randn(n);                  % full impact matrix (OI)
S1 = bvt.structural.construct_Sigt(h, B0);
S0 = construct_Sigt(h, B0);                 % legacy (path-shadowed name)
assert(isequal(S0, S1), 'construct_Sigt: OI-style B0 differs from legacy');
assert(isequal(size(S1), [T, n, n]), 'construct_Sigt: wrong output shape');

A = eye(n); A(3, 1) = .5; A(5, 2) = -.7;    % unit-lower-triangular (CS)
assert(isequal(construct_Sigt(h, A), bvt.structural.construct_Sigt(h, A)), ...
    'construct_Sigt: CS-style A differs from legacy');

    % spot-check the algebra at one t against the dense formula
t = 4;
ref = (B0\diag(exp(h(t, :))))/(B0');
assert(max(abs(squeeze(S1(t, :, :)) - ref), [], 'all') < 1e-14, ...
    'construct_Sigt: disagrees with the dense reference formula');
end
