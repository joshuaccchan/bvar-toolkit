function test_var_ncp
% bvar.models.var_ncp against the frozen draw loop of guides/var_forecast, draw for draw under one
% seed; its prior mean, its choice of c1 and its log marginal likelihood against bvar.ml.ncp; the
% rows that 'omit' leaves out of the likelihood; and the errors
rng(20261122, 'twister');
n = 3; R = 40; T = 80;
A1 = [.5 .1 0; .1 .4 .1; .2 0 .6];
Sigtrue = .5*eye(n) + .2*ones(n);

%% 1. bitwise against the frozen loop of guides/var_forecast, which has four lags
p = 4;
Atrue = [.1*ones(1, n); A1'; zeros(n*(p-1), n)];
data = bvar.util.simulate_var(Atrue, Sigtrue, T + 4);
rng(3, 'twister');
[A_d, Sig_d, c1] = forecast_draws(data, p, R);
rng(3, 'twister');
res = bvar.models.var_ncp(data(1:4, :), data(5:end, :), p, 'nsim', R);
assert(isequal(res.c1, c1), 'var_ncp: c1 differs');
assert(isequal(res.A, reshape(permute(A_d, [3 1 2]), R, [])), 'var_ncp: A differs');
assert(isequal(res.Sig, permute(Sig_d, [3 1 2])), 'var_ncp: Sig differs');

%% 2. a prior mean of one on the own first lag of the variables lev, against bvar.ml.ncp
lev = [2 3]; c3 = 50;
for p = [1 3 5]
    k = 1 + n*p; q0 = max(p, 4);
    Atrue = [.1*ones(1, n); A1'; zeros(n*(p-1), n)];
    Yfull = bvar.util.simulate_var(Atrue, Sigtrue, T + q0);
    res = bvar.models.var_ncp(Yfull(1:q0, :), Yfull(q0+1:end, :), p, 'lev', lev, 'c3', c3, ...
        'nsim', 5);
    [Ys, Zs] = bvar.util.build_lags(Yfull(q0-p+1:end, :), p);
    s2 = bvar.priors.resid_var_ar4(Yfull(1:q0, :), Ys);
    B0 = zeros(k, n); B0(1 + lev(1), lev(1)) = 1; B0(1 + lev(2), lev(2)) = 1;
    VA = @(c1) [c3; kron(1./(1:p)'.^2, c1./s2)];
    assert(isequal(res.s2, s2) && isequal(res.B0, B0) && isequal(res.VA, VA(res.c1)), ...
        'var_ncp: the prior differs, p = %d', p);
    lml = @(c1) bvar.ml.ncp(Ys, Zs, VA(c1), n + 2, diag(s2), B0);
    assert(isequal(res.lml, lml(res.c1)), ...
        'var_ncp: the log marginal likelihood differs, p = %d', p);
    assert(res.lml >= lml(1.01*res.c1) && res.lml >= lml(res.c1/1.01), ...
        'var_ncp: c1 is not a maximum of the marginal likelihood, p = %d', p);
    assert(isequal(size(res.A), [5, k*n]) && isequal(size(res.Sig), [5, n, n]), ...
        'var_ncp: the draws have the wrong size, p = %d', p);
    st = rng;
    r04 = bvar.models.var_ncp(Yfull(1:q0, :), Yfull(q0+1:end, :), p, 'lev', lev, 'c3', c3, ...
        'c1', .04, 'nsim', 0);
    assert(isequal(rng, st), 'var_ncp: nsim = 0 uses the random number generator');
    assert(r04.c1 == .04 && isequal(r04.lml, lml(.04)), ...
        'var_ncp: the log marginal likelihood at a given c1 differs, p = %d', p);
    assert(isequal(size(r04.A), [0, k*n]) && isequal(size(r04.Sig), [0, n, n]), ...
        'var_ncp: nsim = 0 returns draws');
end

%% 3. rows left out of the likelihood, which remain lags of the rows after them
p = 2; k = 1 + n*p;
Atrue = [.1*ones(1, n); A1'; zeros(n*(p-1), n)];
Yfull = bvar.util.simulate_var(Atrue, Sigtrue, T + 4);
Y0 = Yfull(1:4, :); Y = Yfull(5:end, :);
rng(5, 'twister'); r1 = bvar.models.var_ncp(Y0, Y, p, 'lev', lev, 'nsim', 3);
rng(5, 'twister'); r2 = bvar.models.var_ncp(Y0, Y, p, 'lev', lev, 'omit', [], 'nsim', 3);
assert(isequal(r1, r2), 'var_ncp: an empty omit changes the result');
omit = [10 11 12 40];
st = rng;
res = bvar.models.var_ncp(Y0, Y, p, 'lev', lev, 'omit', omit, 'nsim', 7);
[Ys, Zs] = bvar.util.build_lags([Y0(end-p+1:end, :); Y], p);
keep = setdiff(1:T, omit);
s2 = bvar.priors.resid_var_ar4(Y0, Y);
B0 = zeros(k, n); B0(1 + lev(1), lev(1)) = 1; B0(1 + lev(2), lev(2)) = 1;
VA = @(c1) [100; kron(1./(1:p)'.^2, c1./s2)];
lml = @(c1) bvar.ml.ncp(Ys(keep, :), Zs(keep, :), VA(c1), n + 2, diag(s2), B0);
assert(isequal(res.lml, lml(res.c1)) && res.lml >= lml(1.01*res.c1) ...
    && res.lml >= lml(res.c1/1.01), 'var_ncp: the marginal likelihood with omit differs');
rng(st);
[~, Bhat, K, Shat] = bvar.ml.ncp(Ys(keep, :), Zs(keep, :), VA(res.c1), n + 2, diag(s2), B0);
C = chol(K, 'lower');
for r = 1:7
    Sr = iwishrnd(Shat, n + 2 + numel(keep));
    Ar = Bhat + (C'\randn(k, n))*chol(Sr, 'lower')';
    assert(isequal(res.A(r, :), Ar(:)') && isequal(squeeze(res.Sig(r, :, :)), Sr), ...
        'var_ncp: the draws with omit differ, draw %d', r);
end

%% 4. the errors
expect_error('bvar:models:var_ncp:badData', { ...
    @() bvar.models.var_ncp(Y0(:, 1:2), Y, 2), 'Y0 with too few columns';
    @() bvar.models.var_ncp(Y0(1:3, :), Y, 2), 'Y0 with too few rows';
    @() bvar.models.var_ncp(Y0, [Y; NaN(1, n)], 2), 'NaN in Y'});
expect_error('bvar:models:var_ncp:badOption', { ...
    @() bvar.models.var_ncp(Y0, Y, 2, 'shrink', 1), 'an unknown option';
    @() bvar.models.var_ncp(Y0, Y, 2, 'nsim'), 'an option without a value';
    @() bvar.models.var_ncp(Y0, Y, 2, 'lev', 4), 'lev outside the columns';
    @() bvar.models.var_ncp(Y0, Y, 2, 'c1', 0), 'c1 = 0';
    @() bvar.models.var_ncp(Y0, Y, 2, 'c3', -1), 'c3 < 0';
    @() bvar.models.var_ncp(Y0, Y, 2, 'nsim', 2.5), 'a fractional nsim';
    @() bvar.models.var_ncp(Y0, Y, 2, 'omit', T + 1), 'omit outside the rows';
    @() bvar.models.var_ncp(Y0, Y, 2, 'omit', 1:T), 'omit of every row'});
fprintf('PASS test_var_ncp\n');
end

function [Adraw, Sdraw, c1] = forecast_draws(data, p, R)
% frozen from section 2 of guides/var_forecast/guide.m at 479b50d, with R draws and the seed set by
% the caller
[Y, Z] = bvar.util.build_lags(data, p);           % the first four quarters are initial values
[T, n] = size(Y); k = size(Z, 2);
s2 = bvar.priors.resid_var_ar4(data(1:4, :), Y);  % AR(4) residual variances
nu0 = n + 2; S0 = diag(s2);                       % Sig ~ IW(nu0, S0)
c3 = 10^2;                                        % the intercepts
VA = @(c1) [c3; kron(1./(1:p)'.^2, c1./s2)];      % the diagonal of V_A
lc = fminbnd(@(lc) -bvar.ml.ncp(Y, Z, VA(exp(lc)), nu0, S0), log(1e-3), log(10));
c1 = exp(lc);                                     % the shrinkage chosen by the marginal likelihood
[~, Ahat, KA, Shat] = bvar.ml.ncp(Y, Z, VA(c1), nu0, S0);
nuhat = nu0 + T;
LK = chol(KA, 'lower');                           % KA = LK*LK'
Adraw = zeros(k, n, R); Sdraw = zeros(n, n, R);
for r = 1:R
    Sdraw(:, :, r) = iwishrnd(Shat, nuhat);       % Sig ~ IW(nuhat, Shat)
    Adraw(:, :, r) = Ahat + (LK'\randn(k, n))*chol(Sdraw(:, :, r), 'lower')';   % A given Sig
end
end

function expect_error(id, cases)
% each call in the first column of cases must raise the error id
for ib = 1:size(cases, 1)
    try
        cases{ib, 1}();
        error('test_var_ncp:noError', '%s should have errored', cases{ib, 2});
    catch err
        assert(strcmp(err.identifier, id), 'wrong identifier for %s: %s', cases{ib, 2}, ...
            err.identifier);
    end
end
end
