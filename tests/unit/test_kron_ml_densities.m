function test_kron_ml_densities
% deterministic bitwise equivalence of the step-8 density/likelihood
% utilities with their chan2020_jbes_kronecker legacy copies:
%   bvar.ml.lniwpdf      = legacy lniwpdf.m
%   bvar.ml.linvgammpdf  = legacy linvgammpdf.m
%   bvar.ml.llike_ma     = legacy llike_MA.m (root)
%   bvar.ml.llike_csv_ma = legacy llike_CSV_MA.m (package ROOT copy)
% plus the never-merge direction check: the realtime_forecasts copy of
% llike_CSV_MA omits the -n/2*sum(h) term, so the two must DIFFER by
% n/2*sum(h) (up to one-rounding tolerance - the term is folded into the
% root copy's constant before the common subtraction).
root = getappdata(0, 'bvar_repo_root');
leg = fullfile(root, 'replications', 'chan2020_jbes_kronecker', 'legacy');
legrt = fullfile(leg, 'realtime_forecasts');

% --- synthetic evaluation points (deterministic given the seed) ---
rng(8, 'twister');
T = 25; n = 4; k = 9;
A = randn(k, n); A0 = randn(k, n)/5;
Q = randn(n, n+3); Sig = Q*Q'/(n+3); Sig = (Sig+Sig')/2;
Q = randn(n, n+5); S0 = Q*Q'/(n+5); S0 = (S0+S0')/2;
VA0 = 1./gamrnd(3, 1, k, 1);
iVA0 = sparse(1:k, 1:k, 1./VA0);
nu0 = n+3;
U = randn(T, n)*chol(Sig, 'lower')';
h = .3*randn(T, 1);
psi = .15;

% --- core values ---
c_lniw = bvar.ml.lniwpdf(A, Sig, A0, iVA0, nu0, S0);
yv = [.5 1.2 3]; av = [2 3.5 4]; bv = [.5 .8 2];
c_ligp = bvar.ml.linvgammpdf(yv, av, bv);
c_lma = bvar.ml.llike_ma(psi, U, Sig);
c_lcsvma = bvar.ml.llike_csv_ma(psi, U, Sig, h);

% --- legacy root copies, bitwise ---
addpath(leg); c1 = onCleanup(@() rmpath(leg));
assert(strncmpi(which('lniwpdf'), leg, numel(leg)), 'lniwpdf must resolve from legacy');
assert(isequal(lniwpdf(A, Sig, A0, iVA0, nu0, S0), c_lniw), ...
    'bvar.ml.lniwpdf differs from legacy lniwpdf');
assert(isequal(linvgammpdf(yv, av, bv), c_ligp), ...
    'bvar.ml.linvgammpdf differs from legacy linvgammpdf');
assert(isequal(llike_MA(psi, U, Sig), c_lma), ...
    'bvar.ml.llike_ma differs from legacy llike_MA');
l_root = llike_CSV_MA(psi, U, Sig, h);
assert(isequal(l_root, c_lcsvma), ...
    'bvar.ml.llike_csv_ma differs from the legacy ROOT llike_CSV_MA');
clear c1

% --- realtime copy: must differ by exactly the -n/2*sum(h) term ---
addpath(legrt); c2 = onCleanup(@() rmpath(legrt));
assert(strncmpi(which('llike_CSV_MA'), legrt, numel(legrt)), ...
    'llike_CSV_MA must now resolve from realtime_forecasts');
l_rt = llike_CSV_MA(psi, U, Sig, h);
assert(~isequal(l_rt, l_root), ...
    'root and realtime llike_CSV_MA should differ (never-merge)');
gap = (l_rt - l_root) - n/2*sum(h);
assert(abs(gap) < 1e-9*max(1, abs(l_root)), ...
    'root vs realtime llike_CSV_MA difference is not the n/2*sum(h) term (gap %.3g)', gap);
clear c2
end
