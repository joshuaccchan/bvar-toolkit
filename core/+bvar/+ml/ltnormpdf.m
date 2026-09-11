% bvar.ml.ltnormpdf - log density of the normal N(mu,sig2) truncated to (lb,ub);
% the density counterpart of bvar.util.tnormrnd.
%
% See:
% Chan, J.C.C. (2023). Comparing stochastic volatility specifications for large
% Bayesian VARs, Journal of Econometrics, 235(2), 1419-1446.

function lden = ltnormpdf(x, mu, sig2, lb, ub)
c = -.5*log(2*pi*sig2) - log(normcdf((ub-mu)/sqrt(sig2))-normcdf((lb-mu)/sqrt(sig2)));
lden = c -.5/sig2*(x-mu).^2 ;
end
