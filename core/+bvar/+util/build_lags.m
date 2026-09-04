% bvar.util.build_lags - lagged design matrix with the toolkit-wide convention:
% intercept first, then the lag-1 block, ..., then the lag-p block.
%
% [Y, Z] = bvar.util.build_lags(Yfull, p)
%
%   Yfull : T0 x n data matrix (the first p rows serve as initial conditions)
%   Y     : (T0-p) x n left-hand-side observations, Yfull(p+1:end,:)
%   Z     : (T0-p) x (1+n*p) regressor matrix [1, y_{t-1}, ..., y_{t-p}]
%
% New in the consolidated toolkit (2026-09-01, step 3): codifies the inline
% construction repeated in every legacy package (verified identical convention
% across all 12 packages in the 2026-09-01 audit). Where a legacy driver keeps
% separate initial conditions Y0, call build_lags([Y0(end-p+1:end,:); Y], p).

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
