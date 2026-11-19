% bvar.structural.fevd - forecast error variance decomposition: the share of each variable's
% h-step-ahead forecast error variance due to each structural shock. Written for this toolkit.
%
%   share = bvar.structural.fevd(response)
%
%   response : n x n x H impulse responses to all n shocks, response(:,:,h) at horizon h-1, as
%              bvar.structural.irf_redu returns them for an impact matrix L with L*L' the error
%              covariance. A fourth dimension may hold posterior draws
%   share    : the same size as response; share(i,j,h) is the fraction of the h-step-ahead
%              forecast error variance of variable i due to shock j, and share(i,:,h) sums to one
%
% The h-step-ahead forecast error variance of variable i is the sum of its squared responses to
% all shocks at horizons 0 to h-1. This sum equals the variance of the VAR only when response
% holds all n shocks of an impact matrix with L*L' equal to the error covariance, so a response
% with fewer columns is refused. In a VAR with stochastic volatility, pass the responses for the
% impact matrix of the period of interest.

function share = fevd(response)
n = size(response, 1);
if size(response, 2) ~= n
    error('fevd:notAllShocks', 'response must hold all n = %d shocks, not %d', n, ...
        size(response, 2));
end
mse = cumsum(response.^2, 3);   % mse(i,j,h): squared responses summed over horizons 0 to h-1
share = mse./sum(mse, 2);
end
