function test_var_soe
% bvar.models.var_soe, bvar.util.soe_reduced_form and bvar.ml.ncp with a prior mean against the
% frozen local functions of the guides they replaced, draw for draw under one seed;
% soe_reduced_form against the projection that defines the domestic block; and the errors
rng(20261121, 'twister');
n = 4; nf = 2; nd = n - nf; us = 1:nf; au = nf+1:n;
lev = [1 3]; c3 = 100; R = 40; T = 60;
A1 = [.5 .1 0 0; .1 .4 0 0; .2 .1 .5 .1; .1 .2 .1 .3];   % no domestic lags in the foreign equations
Sigtrue = .5*eye(n) + .2*ones(n);

%% 1. bitwise against the frozen functions, with four and with five presample rows
for p = [2 5]
    k = 1 + n*p; q0 = max(p, 4);
    Atrue = [.1*ones(1, n); A1'; zeros(n*(p-1), n)];
    Yfull = bvar.util.simulate_var(Atrue, Sigtrue, T + q0);
    rng(7, 'twister');
    [A_d, Sig_d, c1] = block_draws(Yfull, p, us, au, lev, c3, R);
    rng(7, 'twister');
    res = bvar.models.var_soe(Yfull(1:q0, :), Yfull(q0+1:end, :), p, nf, 'lev', lev, 'c3', c3, ...
        'nsim', R);
    assert(isequal(res.c1, c1), 'var_soe: c1 differs, p = %d', p);
    assert(isequal(res.A, reshape(permute(A_d, [3 1 2]), R, [])), 'var_soe: A differs, p = %d', p);
    assert(isequal(res.Sig, permute(Sig_d, [3 1 2])), 'var_soe: Sig differs, p = %d', p);
    [Ys, Zs] = bvar.util.build_lags(Yfull(q0-p+1:end, :), p);
    s2 = bvar.priors.resid_var_ar4(Yfull(1:q0, :), Ys);
    [VA, B0] = soe_prior(p, n, lev, c1, c3, s2);
    assert(isequal(res.s2, s2) && isequal(res.VA, VA) && isequal(res.B0, B0), ...
        'var_soe: the prior differs, p = %d', p);
    assert(isequal(res.lml, soe_lml(Ys, Zs, us, au, lev, c1, c3, s2)), ...
        'var_soe: the log marginal likelihood differs, p = %d', p);
    st = rng;
    r04 = bvar.models.var_soe(Yfull(1:q0, :), Yfull(q0+1:end, :), p, nf, 'lev', lev, 'c3', c3, ...
        'c1', .04, 'nsim', 0);
    assert(isequal(rng, st), 'var_soe: nsim = 0 uses the random number generator');
    assert(isequal(r04.lml, soe_lml(Ys, Zs, us, au, lev, .04, c3, s2)) && r04.c1 == .04, ...
        'var_soe: the log marginal likelihood at a given c1 differs, p = %d', p);
    assert(isequal(size(r04.A), [0, k*n]) && isequal(size(r04.Sig), [0, n, n]), ...
        'var_soe: nsim = 0 returns draws');
    % bvar.ml.ncp with a prior mean, on the domestic block
    Xd = [Zs, Ys(:, us)]; VAd = [VA; 1./s2(us)]; B0d = [B0(:, au); zeros(nf, nd)];
    [l1, B1, K1, S1] = bvar.ml.ncp(Ys(:, au), Xd, VAd, n + 2, diag(s2(au)), B0d);
    [B2, K2, S2] = ncp_post(Ys(:, au), Xd, B0d, VAd, diag(s2(au)));
    l2 = ncp_lml(Ys(:, au), Xd, B0d, VAd, n + 2, diag(s2(au)));
    assert(isequal(l1, l2) && isequal(B1, B2) && isequal(K1, K2) && isequal(S1, S2), ...
        'ncp: differs from the frozen ncp_post and ncp_lml, p = %d', p);
end

%% 2. soe_reduced_form, bitwise against the frozen reduced_form and against the projection
p = 3; k = 1 + n*p;
ius = [1, reshape(1 + us' + (0:p-1)*n, 1, [])];
iau = setdiff(1:k, ius);
Bf = randn(numel(ius), nf); Bd = randn(k + nf, nd);
Sf = [1 .3; .3 .8]; Om = [.6 -.2; -.2 .9];
[A, Sig] = bvar.util.soe_reduced_form(Bf, Bd, Sf, Om);
[A2, Sig2] = reduced_form(Bf, Bd, Sf, Om, ius, k);
assert(isequal(A, A2) && isequal(Sig, Sig2), 'soe_reduced_form: differs from reduced_form');
G0 = Sig(au, us)/Sig(us, us);             % the projection of the domestic errors on the foreign
assert(isequal(A(ius, us), Bf) && all(A(iau, us) == 0, 'all') && isequal(Sig(us, us), Sf), ...
    'soe_reduced_form: the foreign block is not the foreign VAR');
assert(max(abs([G0 - Bd(k+1:end, :)', Sig(au, au) - G0*Sig(us, au) - Om]), [], 'all') < 1e-12, ...
    'soe_reduced_form: the domestic errors do not have the projection of the domestic block');
assert(max(abs(A(:, au) - A(:, us)*G0' - Bd(1:k, :)), [], 'all') < 1e-12, ...
    'soe_reduced_form: the domestic coefficients do not match the domestic block');

%% 3. the errors
Y0 = Yfull(1:4, :); Y = Yfull(5:end, :);
expect_error('bvar:models:var_soe:badData', {
    @() bvar.models.var_soe(Y0, Y, 2, 0), 'nf = 0'
    @() bvar.models.var_soe(Y0, Y, 2, n), 'nf = n'
    @() bvar.models.var_soe(Y0, Y, 2, 1.5), 'nf = 1.5'
    @() bvar.models.var_soe(Y0(2:end, :), Y, 2, nf), 'three presample rows'
    @() bvar.models.var_soe(Y0(:, 1:3), Y, 2, nf), 'a narrow Y0'
    @() bvar.models.var_soe(Y0, [Y; nan(1, n)], 2, nf), 'a NaN'});
expect_error('bvar:models:var_soe:badOption', {
    @() bvar.models.var_soe(Y0, Y, 2, nf, 'zzz', 1), 'an unknown option'
    @() bvar.models.var_soe(Y0, Y, 2, nf, 'nsim'), 'an option with no value'
    @() bvar.models.var_soe(Y0, Y, 2, nf, 'lev', 5), 'lev = 5'
    @() bvar.models.var_soe(Y0, Y, 2, nf, 'c1', 0), 'c1 = 0'
    @() bvar.models.var_soe(Y0, Y, 2, nf, 'c3', -1), 'c3 = -1'
    @() bvar.models.var_soe(Y0, Y, 2, nf, 'nsim', 1.5), 'nsim = 1.5'});
expect_error('bvar:util:soe_reduced_form:badDims', {
    @() bvar.util.soe_reduced_form(Bf, Bd(2:end, :), Sf, Om), 'a short Bd'
    @() bvar.util.soe_reduced_form(Bf, Bd, Om(1, 1), Om), 'a scalar Sf'});
end

function expect_error(id, cases)
% each call in the first column of cases must raise the error id
for ib = 1:size(cases, 1)
    try
        cases{ib, 1}();
        error('test_var_soe:noError', '%s should have errored', cases{ib, 2});
    catch err
        assert(strcmp(err.identifier, id), 'wrong identifier for %s: %s', cases{ib, 2}, ...
            err.identifier);
    end
end
end

% -------------------------------------------------------------------------
% FROZEN copies of the local functions of guides/soe_var, soe_exogeneity_test, soe_foreign_shock
% and soe_conditional_forecast, as they stood before bvar.models.var_soe, bvar.util.soe_reduced_form
% and bvar.ml.ncp existed. Do not edit to follow the library.
function [A_d, Sig_d, c1] = block_draws(Yfull, p, us, au, lev, c3, R)
% R draws of the reduced form from the posterior of the two blocks, as in guides/soe_var, from a
% sample whose first max(p, 4) rows are the presample; c1 maximizes the marginal likelihood
nus = numel(us); nau = numel(au); n = nus + nau;
q0 = max(p, 4);
[Ys, Zs] = bvar.util.build_lags(Yfull(q0-p+1:end, :), p);
[T, k] = size(Zs);
ius = [1, reshape(1 + us' + (0:p-1)*n, 1, [])];
s2 = bvar.priors.resid_var_ar4(Yfull(1:q0, :), Ys);
lc = fminbnd(@(lc) -soe_lml(Ys, Zs, us, au, lev, exp(lc), c3, s2), log(1e-3), log(10));
c1 = exp(lc);
[VA, B0] = soe_prior(p, n, lev, c1, c3, s2);
[Bh_us, K_us, Sh_us] = ncp_post(Ys(:, us), Zs(:, ius), B0(ius, us), VA(ius), diag(s2(us)));
[Bh_au, K_au, Sh_au] = ncp_post(Ys(:, au), [Zs, Ys(:, us)], [B0(:, au); zeros(nus, nau)], ...
    [VA; 1./s2(us)], diag(s2(au)));
Cus = chol(K_us, 'lower'); Cau = chol(K_au, 'lower');
A_d = zeros(k, n, R); Sig_d = zeros(n, n, R);
for r = 1:R
    Sus = iwishrnd(Sh_us, nus + 2 + T);                 % the US block
    Bus = Bh_us + (Cus'\randn(numel(ius), nus))*chol(Sus, 'lower')';
    Om = iwishrnd(Sh_au, n + 2 + T);                    % Australia, independent of the US block
    Bau = Bh_au + (Cau'\randn(k + nus, nau))*chol(Om, 'lower')';
    [A_d(:, :, r), Sig_d(:, :, r)] = reduced_form(Bus, Bau, Sus, Om, ius, k);
end
end

function [VA, B0] = soe_prior(p, n, lev, c1, c3, s2)
% the prior of the coefficients on x_t = (1, y_{t-1}', ..., y_{t-p}')': their variances VA, and
% their means B0, one on the own first lag of the variables lev, in levels, and zero elsewhere
VA = [c3; kron(1./(1:p)'.^2, c1./s2)];
B0 = zeros(1 + n*p, n);
B0(sub2ind(size(B0), 1 + lev, lev)) = 1;
end

function lml = soe_lml(Y, Z, us, au, lev, c1, c3, s2)
% the log marginal likelihood of the VAR with block exogeneity: the sum over the two blocks
n = size(Y, 2); p = (size(Z, 2) - 1)/n; nus = numel(us); nau = numel(au);
ius = [1, reshape(1 + us' + (0:p-1)*n, 1, [])];
[VA, B0] = soe_prior(p, n, lev, c1, c3, s2);
lml = ncp_lml(Y(:, us), Z(:, ius), B0(ius, us), VA(ius), nus + 2, diag(s2(us))) ...
    + ncp_lml(Y(:, au), [Z, Y(:, us)], [B0(:, au); zeros(nus, nau)], [VA; 1./s2(us)], ...
    n + 2, diag(s2(au)));
end

function [B, K, Shat] = ncp_post(Y, X, B0, VA, S0)
% the posterior of the multivariate regression Y = X*B + E under vec(B) | Sig ~ N(vec(B0), Sig
% kron diag(VA)) and Sig ~ IW(nu0, S0): vec(B) | Sig, Y ~ N(vec(B), Sig kron inv(K)) and
% Sig | Y ~ IW(nu0 + T, Shat)
K = diag(1./VA) + X'*X;
B = K\(B0./VA + X'*Y);
E = Y - X*B; D = B - B0;
Shat = S0 + E'*E + D'*(D./VA);           % S0 + Y'*Y + B0'*(B0./VA) - B'*K*B, computed stably
Shat = (Shat + Shat')/2;
end

function lml = ncp_lml(Y, X, B0, VA, nu0, S0)
% the log marginal likelihood of the multivariate regression Y = X*B + E under the prior of
% ncp_post, in closed form, as in guides/natural_conjugate_ml
[T, m] = size(Y);
[~, K, Shat] = ncp_post(Y, X, B0, VA, S0);
lml = -m*T/2*log(pi) - m/2*(sum(log(VA)) + bvar.util.ldet(K)) ...
    + nu0/2*bvar.util.ldet(S0) - (nu0 + T)/2*bvar.util.ldet(Shat) ...
    + bvar.util.mgammaln(m, (nu0 + T)/2) - bvar.util.mgammaln(m, nu0/2);
end

function [A, Sig] = reduced_form(Bus, Bau, Sus, Om, ius, k)
% the reduced form of the two blocks: the US block yf_t' = x_t(ius)'*Bus + uf_t', and Australia
% yd_t' = (x_t', yf_t')*Bau + e_t', with e_t independent of uf_t
nus = size(Bus, 2);
G0 = Bau(k+1:end, :)';                   % the current US coefficients
A = zeros(k, nus + size(Bau, 2));
A(ius, 1:nus) = Bus;                     % zeros for the Australian lags
A(:, nus+1:end) = Bau(1:k, :) + A(:, 1:nus)*G0';
Sig = [Sus, Sus*G0'; G0*Sus, G0*Sus*G0' + Om];
end
