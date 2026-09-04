function test_igrnd
% bvar.util.igrnd must draw identically to the inline 1/gamrnd(...) expression
% the legacy samplers use, under the same seed.

nu = 5; S = 0.04;

rng(11, 'twister');
x_core = bvar.util.igrnd(nu, S);
rng(11, 'twister');
x_inline = 1/gamrnd(nu, 1/S);
assert(isequal(x_core, x_inline), 'igrnd: differs from the inline expression');

% the conjugate variance step as written in the SV samplers
rng(12, 'twister');
T = 200; e = 0.3*randn(T,1); nuh0 = 5; Sh0 = 0.01*(nuh0-1);
rng(13, 'twister');
s_core = bvar.util.igrnd(nuh0 + T/2, Sh0 + sum(e.^2)/2);
rng(13, 'twister');
s_inline = 1/gamrnd(nuh0 + T/2, 1/(Sh0 + sum(e.^2)/2));
assert(isequal(s_core, s_inline), 'igrnd: conjugate step differs from inline');

% vectorized: one draw per element, same stream as the vector gamrnd call
nuv = [3; 4; 5]; Sv = [0.1; 0.2; 0.3];
rng(14, 'twister');
v_core = bvar.util.igrnd(nuv, Sv);
rng(14, 'twister');
v_inline = 1./gamrnd(nuv, 1./Sv);
assert(isequal(v_core, v_inline), 'igrnd: vector form differs from inline');
assert(isequal(size(v_core), [3 1]), 'igrnd: wrong output size');
assert(all(v_core > 0), 'igrnd: draws must be positive');

% distributional sanity: IG(nu,S) has mean S/(nu-1) for nu > 1
rng(15, 'twister');
nu_m = 8; S_m = 3;
draws = bvar.util.igrnd(nu_m*ones(2e5,1), S_m*ones(2e5,1));
assert(abs(mean(draws) - S_m/(nu_m-1)) < 0.02, 'igrnd: sample mean off target');

% bad input rejected
try
    bvar.util.igrnd(-1, 1);
    error('test:noThrow', 'igrnd should reject nu <= 0');
catch err
    assert(strcmp(err.identifier, 'bvar:util:igrnd:badParam'), 'igrnd: wrong error id');
end
end
