% bvar.util.build_lags - lagged design matrix with the toolkit-wide convention:
% intercept first, then the lag-1 block, ..., then the lag-p block.
%
% [Y, Z] = bvar.util.build_lags(Yfull, p)
%
%   Yfull : T0 x n data matrix (the first p rows serve as initial conditions)
%   p     : lag length, a positive integer strictly less than T0
%   Y     : (T0-p) x n left-hand-side observations, Yfull(p+1:end,:)
%   Z     : (T0-p) x (1+n*p) regressor matrix [1, y_{t-1}, ..., y_{t-p}]
%
% If a driver holds its initial conditions separately as Y0, pass them in the
% same matrix: build_lags([Y0(end-p+1:end,:); Y], p).
%
% Provenance and the legacy copies this stands in for: tests/variant_map.md.
% Equivalence: tests/unit/test_build_lags.m.

function [Y, Z] = build_lags(Yfull, p)
    [T0, n] = size(Yfull);
    if p < 1 || p ~= round(p)
        error('build_lags:badLag', 'p must be a positive integer');
    end
    if T0 <= p
        error('build_lags:tooShort', 'need more than p = %d rows of data', p);
    end
    T = T0 - p;
    Y = Yfull(p+1:end, :);
    Z = ones(T, 1 + n*p);
    for ii = 1:p
        Z(:, 1 + (ii-1)*n + (1:n)) = Yfull(p+1-ii:T0-ii, :);
    end
end
