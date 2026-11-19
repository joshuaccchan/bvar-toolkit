% bvar.structural.hist_decomp - historical decomposition: each variable's sample path split into
% the contribution of each structural shock and a base path from the intercepts and the initial
% values. Written for this toolkit.
%
%   [contrib, base, shocks] = bvar.structural.hist_decomp(Yfull, A, L)
%
%   Yfull   : T0 x n data; the first p rows are the initial values
%   A       : k x n coefficients, intercepts first, then the lag-1 block, ..., the lag-p block,
%             with k = 1 + n*p (the layout of bvar.util.build_lags)
%   L       : n x n impact matrix with L*L' the error covariance, or n x n x T with one impact
%             matrix per period, as in a VAR with stochastic volatility
%   contrib : T x n x n with T = T0 - p; contrib(t,i,j) is the contribution of shock j to
%             variable i in period t
%   base    : T x n path of the VAR from the initial values with every shock set to zero
%   shocks  : T x n structural shocks, inv(L)*u_t for the residuals u_t
%
% The pieces add up to the data: Yfull(p+1:end,:) = base + sum(contrib, 3), up to rounding. The
% contribution of shock j runs the VAR without intercepts from zero initial values, driven by
% shock j alone, so it equals the sum over s of the response at horizon s times the shock in
% period t - s.

function [contrib, base, shocks] = hist_decomp(Yfull, A, L)
[k, n] = size(A);
p = (k - 1)/n;
if p < 1 || p ~= round(p)
    error('hist_decomp:badA', 'A must have 1 + n*p rows for its n = %d columns', n);
end
[Y, Z] = bvar.util.build_lags(Yfull, p);
T = size(Y, 1);
if size(L, 1) ~= n || size(L, 2) ~= n || ~any(size(L, 3) == [1, T])
    error('hist_decomp:badL', 'L must be n x n or n x n x T, with n = %d and T = %d', n, T);
end
U = Y - Z*A;                                % reduced-form residuals
shocks = zeros(T, n);
X = zeros(n, n, p + T);                     % X(:,j,p+t): contribution of shock j in period t
for t = 1:T
    Lt = L(:, :, min(t, size(L, 3)));
    shocks(t, :) = (Lt\U(t, :)')';
    x = Lt.*shocks(t, :);                   % column j: L(:,j) times shock j
    for l = 1:p
        x = x + A(1 + (l-1)*n + (1:n), :)'*X(:, :, p + t - l);
    end
    X(:, :, p + t) = x;
end
contrib = permute(X(:, :, p+1:end), [3 1 2]);
Yb = [Yfull(1:p, :); zeros(T, n)];          % the base path: every shock at zero
for t = p+1:p+T
    Yb(t, :) = [1, reshape(Yb(t-1:-1:t-p, :)', 1, [])]*A;
end
base = Yb(p+1:end, :);
end
