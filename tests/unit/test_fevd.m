function test_fevd
% fevd against the forecast error variances of the reduced form for two impact matrices, the
% stationary variances at long horizons, a recursive impact matrix on impact, posterior draws in
% a fourth dimension, a single horizon, and its error
n = 3; p = 2; H = 12;
B = [.5 .1 0; .2 .4 .1; 0 .1 .3; -.1 0 .05; .05 -.1 0; 0 .05 .1];   % lag 1, lag 2 ((n p) x n)
Sig = [1 .3 .2; .3 .8 .1; .2 .1 .5];
F = [B'; eye(n*(p-1)), zeros(n*(p-1), n)];
assert(max(abs(eig(F))) < 1, 'test_fevd: the test VAR must be stable');
L = chol(Sig, 'lower');
[Q, ~] = qr([1 2 0; -1 1 3; 2 0 1]);         % an orthogonal rotation
Phi = bvar.structural.irf_redu(B, eye(n), H, n);   % reduced-form moving-average coefficients

% against the h-step forecast error covariance sum_s Phi_s Sig Phi_s', for L and for L*Q
for M = {L, L*Q}
    Lm = M{1};
    share = bvar.structural.fevd(bvar.structural.irf_redu(B, Lm, H, n));
    want = zeros(n, n, H);
    V = zeros(n, n);                          % h-step forecast error covariance
    C = zeros(n, n);                          % squared responses summed over horizons, by shock
    for h = 1:H
        V = V + Phi(:, :, h)*Sig*Phi(:, :, h)';
        C = C + (Phi(:, :, h)*Lm).^2;
        want(:, :, h) = C./diag(V);
    end
    assert(max(abs(share - want), [], 'all') < 1e-12, 'fevd: differs from the reduced form');
    assert(max(abs(sum(share, 2) - 1), [], 'all') < 1e-12, 'fevd: shares do not sum to one');
end

% at a long horizon, the shares of the stationary variance: the process driven by shock j alone
% has the covariance V_j with vec(V_j) = (I - F kron F)^{-1} vec(Q_j)
share = bvar.structural.fevd(bvar.structural.irf_redu(B, L, 400, n));
Vj = zeros(n, n);
for j = 1:n
    Qj = blkdiag(L(:, j)*L(:, j)', zeros(n*(p-1)));
    Vfull = reshape((eye((n*p)^2) - kron(F, F))\Qj(:), n*p, n*p);
    Vj(:, j) = diag(Vfull(1:n, 1:n));
end
assert(max(abs(share(:, :, end) - Vj./sum(Vj, 2)), [], 'all') < 1e-10, ...
    'fevd: long-horizon shares differ from the stationary variances');

% recursive impact matrix: on impact the first variable is all shock 1, and shock j does not
% move the variables ordered before it
share = bvar.structural.fevd(bvar.structural.irf_redu(B, L, H, n));
assert(share(1, 1, 1) == 1 && all(triu(share(:, :, 1), 1) == 0, 'all'), ...
    'fevd: wrong shares on impact for a recursive impact matrix');

% posterior draws in a fourth dimension, and a single horizon passed as an n x n matrix
R = cat(4, bvar.structural.irf_redu(B, L, H, n), bvar.structural.irf_redu(B, L*Q, H, n));
S = bvar.structural.fevd(R);
assert(isequal(S(:, :, :, 1), share) && ...
    isequal(S(:, :, :, 2), bvar.structural.fevd(bvar.structural.irf_redu(B, L*Q, H, n))), ...
    'fevd: draws in a fourth dimension');
assert(max(abs(bvar.structural.fevd(L) - share(:, :, 1)), [], 'all') < 1e-15, ...
    'fevd: a single horizon');

% error: responses to fewer than n shocks
assert_error(@() bvar.structural.fevd(bvar.structural.irf_redu(B, L, H, 2)), 'fevd:notAllShocks');
end

function assert_error(f, id)
try
    f();
catch err
    assert(strcmp(err.identifier, id), 'expected %s, got %s', id, err.identifier);
    return
end
error('test_fevd:noError', 'expected the error %s', id);
end
