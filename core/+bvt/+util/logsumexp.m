% bvt.util.logsumexp - numerically stable log(sum(exp(x))) along a dimension.
%
% y = bvt.util.logsumexp(x)        operates along dim 1
% y = bvt.util.logsumexp(x, dim)
%
% For averaging log predictive-likelihood draws, use
%   logsumexp(logdraws) - log(M)
% which replaces the ad-hoc max-shift blocks inlined in the legacy forecasting
% scripts. New in the consolidated toolkit (2026-09-01, step 3).

function y = logsumexp(x, dim)
    if nargin < 2
        dim = 1;
    end
    m = max(x, [], dim);
    ms = m;
    ms(~isfinite(ms)) = 0;                    % avoid Inf - Inf = NaN below
    y = ms + log(sum(exp(x - ms), dim));
    y(m == -Inf) = -Inf;                      % all-(-Inf) slices
end
