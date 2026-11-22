function test_ncp
% bvar.ml.ncp against the frozen ncp_lml of the guides it replaced, against Bayes' theorem at two
% points of the parameter space and bvar.ml.kron_bvar, with a prior mean, and its errors
rng(20261121, 'twister');
T = 40; m = 3; k = 7;
X = [ones(T, 1), randn(T, k - 1)];
Y = X*randn(k, m)*.5 + randn(T, m);
VA = [100; .05 + .3*rand(k - 1, 1)];
nu0 = m + 2; S0 = diag(.5 + rand(m, 1));
B0 = .3*randn(k, m);

%% 1. bitwise against the frozen function of guides/natural_conjugate_ml, which has no B0
[l1, B1, K1, S1] = ncp_lml(Y, X, VA, nu0, S0);
[l2, B2, K2, S2] = bvar.ml.ncp(Y, X, VA, nu0, S0);
assert(isequal(l1, l2) && isequal(B1, B2) && isequal(K1, K2) && isequal(S1, S2), ...
    'ncp: differs from the frozen ncp_lml of guides/natural_conjugate_ml');
[l3, B3, K3, S3] = bvar.ml.ncp(Y, X, VA', nu0, S0, zeros(k, m));   % VA as a row, B0 = 0
assert(isequal(l2, l3) && isequal(B2, B3) && isequal(K2, K3) && isequal(S2, S3), ...
    'ncp: a row VA or an explicit zero B0 changes the result');

%% 2. Bayes' theorem at two points: log p(Y) = log p(Y | B, Sig) + log p(B, Sig) - log p(B, Sig | Y)
[lml, Bhat, K, Shat] = bvar.ml.ncp(Y, X, VA, nu0, S0, B0);
pts = {Bhat, Shat/(nu0 + T + m + 1); B0 + .2, S0 + .1*ones(m)};
for ip = 1:2
    B = pts{ip, 1}; Sig = pts{ip, 2};
    E = Y - X*B;
    ll = -T*m/2*log(2*pi) - T/2*bvar.util.ldet(Sig) - .5*trace(Sig\(E'*E));
    id = ll + bvar.ml.lniwpdf(B, Sig, B0, sparse(diag(1./VA)), nu0, S0) ...
        - bvar.ml.lniwpdf(B, Sig, Bhat, K, nu0 + T, Shat);
    assert(abs(id - lml) < 1e-8, 'ncp: Bayes'' theorem fails at point %d by %.1e', ip, id - lml);
end
pri = struct('A0', B0, 'VA0', VA, 'nu0', nu0, 'S0', S0);
est = struct('Ahat', Bhat, 'Shat', Shat, 'KA', K);
assert(abs(bvar.ml.kron_bvar(Y, X, pri, est) - lml) < 1e-8, 'ncp: differs from kron_bvar');

%% 3. the errors
bad = { @() bvar.ml.ncp(Y(2:end, :), X, VA, nu0, S0), 'bvar:ml:ncp:badDims', 'rows of X and Y'
        @() bvar.ml.ncp(Y, X, VA(2:end), nu0, S0), 'bvar:ml:ncp:badPrior', 'a short VA'
        @() bvar.ml.ncp(Y, X, -VA, nu0, S0), 'bvar:ml:ncp:badPrior', 'a negative VA'
        @() bvar.ml.ncp(Y, X, VA, m - 1, S0), 'bvar:ml:ncp:badPrior', 'nu0 = m - 1'
        @() bvar.ml.ncp(Y, X, VA, nu0, eye(m + 1)), 'bvar:ml:ncp:badPrior', 'a large S0'
        @() bvar.ml.ncp(Y, X, VA, nu0, S0, B0'), 'bvar:ml:ncp:badPrior', 'a transposed B0'};
for ib = 1:size(bad, 1)
    try
        bad{ib, 1}();
        error('test_ncp:noError', '%s should have errored', bad{ib, 3});
    catch err
        assert(strcmp(err.identifier, bad{ib, 2}), 'wrong identifier for %s: %s', bad{ib, 3}, ...
            err.identifier);
    end
end
end

% -------------------------------------------------------------------------
% FROZEN copy of the local function of guides/natural_conjugate_ml/guide.m, var_forecast,
% var_forecast_evaluation and var_lag_shrinkage, as it stood before bvar.ml.ncp existed. Do not
% edit to follow the library.
function [lml, Ahat, KA, Shat] = ncp_lml(Y, Z, VA, nu0, S0)
% the log marginal likelihood under vec(A) | Sig ~ N(0, Sig kron diag(VA)) and Sig ~ IW(nu0, S0),
% with the posterior quantities Ahat, KA and Shat
[T, n] = size(Y);
KA = diag(1./VA) + Z'*Z;
Ahat = KA\(Z'*Y);
E = Y - Z*Ahat;                        % the residuals at the posterior mean
Shat = S0 + E'*E + Ahat'*(Ahat./VA);   % S0 + Y'*Y - Ahat'*KA*Ahat, without the cancellation
Shat = (Shat + Shat')/2;
lml = -n*T/2*log(pi) - n/2*(sum(log(VA)) + bvar.util.ldet(KA)) ...
    + nu0/2*bvar.util.ldet(S0) - (nu0 + T)/2*bvar.util.ldet(Shat) ...
    + bvar.util.mgammaln(n, (nu0 + T)/2) - bvar.util.mgammaln(n, nu0/2);
end
