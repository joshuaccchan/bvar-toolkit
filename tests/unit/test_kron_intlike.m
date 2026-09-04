function test_kron_intlike
% seeded draw-for-draw bitwise equivalence of the four extracted integrated-
% likelihood evaluators with their chan2020_jbes_kronecker legacy copies:
%   bvar.ml.intlike_csv      = legacy intlike_BVAR_CSV.m
%   bvar.ml.intlike_t_csv    = legacy intlike_BVAR_t_CSV.m
%   bvar.ml.intlike_csv_ma   = legacy intlike_BVAR_CSV_MA.m
%   bvar.ml.intlike_csv_t_ma = legacy intlike_BVAR_CSV_t_MA.m
% at a small importance-sample size R, on the REAL package data with an OLS
% (A, Sig) evaluation point (so the Newton-Raphson / EM mode searches run
% exactly as in production). Asserts isequal on [intlike, store_llike] and
% on the terminal rng state (identical randn call sequence).
root = getappdata(0, 'bvar_repo_root');
leg = fullfile(root, 'replications', 'chan2020_jbes_kronecker', 'legacy');

    % design exactly as main_BVAR.m lines 26-31 (read-only)
data_Q = load(fullfile(leg, 'data_Q.csv'));
data = data_Q(:, [1:3 6:15 17 19:24]);
Y0 = data(1:4, :); shortY = data(5:end, :);
[T, n] = size(shortY); p = 4;
tmpY = [Y0; shortY];
X = zeros(T, n*p);
for i = 1:p
    X(:, (i-1)*n+1:i*n) = tmpY(p-i+1:end-i, :);
end
X = [ones(T, 1) X];

    % OLS evaluation point + working-model parameter values
A = (X'*X)\(X'*shortY);
U = shortY - X*A;
Sig = U'*U/T; Sig = (Sig+Sig')/2;
rho = .95; sigh2 = .05; psi = .1; nu = 10;
R = 25; seed = 20260902;

addpath(leg); c = onCleanup(@() rmpath(leg));
assert(strncmpi(which('intlike_BVAR_CSV'), leg, numel(leg)), ...
    'intlike_BVAR_CSV must resolve from the legacy folder');

cases = { ...
    @() bvar.ml.intlike_csv(shortY,X,A,Sig,rho,sigh2,R),          @() intlike_BVAR_CSV(shortY,X,A,Sig,rho,sigh2,R),          'intlike_csv'; ...
    @() bvar.ml.intlike_t_csv(shortY,X,A,Sig,rho,sigh2,nu,R),     @() intlike_BVAR_t_CSV(shortY,X,A,Sig,rho,sigh2,nu,R),     'intlike_t_csv'; ...
    @() bvar.ml.intlike_csv_ma(shortY,X,A,Sig,psi,rho,sigh2,R),   @() intlike_BVAR_CSV_MA(shortY,X,A,Sig,psi,rho,sigh2,R),   'intlike_csv_ma'; ...
    @() bvar.ml.intlike_csv_t_ma(shortY,X,A,Sig,psi,rho,sigh2,nu,R), @() intlike_BVAR_CSV_t_MA(shortY,X,A,Sig,psi,rho,sigh2,nu,R), 'intlike_csv_t_ma'};

for kc = 1:size(cases, 1)
    rng(seed, 'twister');
    [ic, sc] = cases{kc, 1}();
    rc = rng;
    rng(seed, 'twister');
    [il, sl] = cases{kc, 2}();
    rl = rng;
    assert(isequal(ic, il), '%s: intlike differs from legacy', cases{kc, 3});
    assert(isequal(sc, sl), '%s: store_llike differs from legacy', cases{kc, 3});
    assert(isequal(rc.State, rl.State), '%s: rng call sequence differs', cases{kc, 3});
    assert(isfinite(ic), '%s: intlike not finite - harness broken', cases{kc, 3});
end
end
