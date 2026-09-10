% bvar.util.logsumexp - numerically stable log(sum(exp(x))) along a dimension.
%
% y = bvar.util.logsumexp(x)        operates along dim 1
% y = bvar.util.logsumexp(x, dim)
%
%   x   : array of log-scale values; a slice that is entirely -Inf returns
%         -Inf, with no NaN
%   dim : dimension to reduce (default 1)
%
% For averaging M log predictive-likelihood draws, use
%   logsumexp(logdraws) - log(M)
%
% Tests: tests/unit/test_logsumexp.m.

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
