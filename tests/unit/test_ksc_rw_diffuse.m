function test_ksc_rw_diffuse
% bvt.sv.ksc_rw_diffuse must reproduce the sp_code SVRW.m (diffuse h_1 ~ N(0,Vh),
% returns indicators) draw-for-draw under one seed - path AND indicators.
root = getappdata(0, 'bvt_repo_root');

rng(42, 'twister');                          % fixed test data
T = 211;
e = randn(T,1).*exp(0.3*randn(T,1));
ystar = log(e.^2 + 1e-4);
h_in = 0.2*randn(T,1);
omega2h = 0.09;
Vh = 9;

rng(8, 'twister');
[h_core, S_core] = bvt.sv.ksc_rw_diffuse(ystar, h_in, omega2h, Vh);

leg = fullfile(root, 'replications', 'chan_jeliazkov2009_statespace', 'legacy', 'sp_code');
addpath(leg); c = onCleanup(@() rmpath(leg));  % prepended -> legacy copy shadows
rng(8, 'twister');
[h_leg, S_leg] = SVRW(ystar, h_in, omega2h, Vh);

assert(isequal(h_leg, h_core), 'ksc_rw_diffuse: h differs from legacy sp_code SVRW');
assert(isequal(S_leg, S_core), 'ksc_rw_diffuse: indicators differ from legacy sp_code SVRW');
assert(all(ismember(S_core, 1:7)), 'ksc_rw_diffuse: indicators outside 1..7');
end
