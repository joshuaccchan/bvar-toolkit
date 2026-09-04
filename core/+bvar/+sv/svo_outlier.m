% bvar.sv.svo_outlier - conditional draw of the outlier scales o_t (from a
% discrete grid, one independent draw per period) and of the outlier
% probability po (beta conjugate) in the VAR-SVO model. Period t is an outlier
% when o_t > 1; the grid's first point is 1, the remaining ngrid points share
% the prior mass po.
%
%   [o,po] = bvar.sv.svo_outlier(Y, X, A, B0, h, o_grid, po, p0a, p0b)
%
%   o_grid : (ngrid+1) x 1, [1; grid of outlier sizes]
%   po     : current outlier probability (in), new draw (out)
%   p0a,p0b: beta prior parameters
%
% rng: rand once per period t = 1:T, then one betarnd.
%
% Extracted 2026-09-03 (step 9). Canonical source, body verbatim:
% chan2023_joe_mlvarsv/legacy/VAR_ARSVO_redu.m lines 112-124 (single copy in
% the repo). Wrapped as a function with T, n, ngrid taken from the arguments
% and Hyper.p0a/p0b passed explicitly; o is preallocated instead of updated in
% place (every element is overwritten before any read, so the values match).
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
