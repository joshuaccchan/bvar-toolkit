% This script returns the log density of the inverse-gamma distribuiton
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

function lden = linvgammpdf(y, a, b)

lden = a.*log(b) - gammaln(a) - (a+1) .* log(y) - b./y;

end