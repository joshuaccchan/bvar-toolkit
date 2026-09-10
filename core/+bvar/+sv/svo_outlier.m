% bvar.sv.svo_outlier - conditional draw of the outlier scales o_t (from a
% discrete grid, one independent draw per period) and of the outlier
% probability po (beta conjugate) in the VAR-SVO model. Period t is an outlier
% when o_t > 1; the grid's first point is 1, the remaining ngrid points share
% the prior mass po.
%
%   [o,po] = bvar.sv.svo_outlier(Y, X, A, B0, h, o_grid, po, p0a, p0b)
%
%   Y, X   : T x n data and T x k regressors
%   A      : k x n VAR coefficient matrix
%   B0     : n x n impact matrix; (Y-X*A)*B0' are the orthogonalized errors
%   h      : T x n log-volatility paths
%   o_grid : (ngrid+1) x 1, [1; grid of outlier sizes]
%   po     : current outlier probability (in), new draw (out)
%   p0a,p0b: beta prior parameters
%   o      : T x 1 draw, each element one of the o_grid points
%
% rng consumption: rand once per period t = 1:T, then one betarnd.
%
% Provenance and the legacy copies this stands in for: tests/variant_map.md.
% Equivalence: tests/unit/test_mlvarsv_equivalence.m.
%
% See:
% Chan, J.C.C. (2023). Comparing stochastic volatility specifications for large
% Bayesian VARs, Journal of Econometrics, 235(2), 1419-1446.

function [o,po] = svo_outlier(Y,X,A,B0,h,o_grid,po,p0a,p0b)
[T,n] = size(Y);
ngrid = numel(o_grid) - 1;
o = zeros(T,1);
o_lpri = log([1-po; repmat(po/ngrid,ngrid,1)]);
U = ((Y-X*A)*B0')./exp(h/2);
for tt=1:T
    lliket = -n*log(o_grid) -.5*U(tt,:)*U(tt,:)'./o_grid.^2;
    o_post = exp(lliket + o_lpri - max(lliket));
    o_post = o_post/sum(o_post);
    idx = find(rand<cumsum(o_post),1);
    o(tt) = o_grid(idx);
end

    % sample po
tmp = sum(o>1);
po = betarnd(p0a + tmp, p0b + T-tmp);
end
