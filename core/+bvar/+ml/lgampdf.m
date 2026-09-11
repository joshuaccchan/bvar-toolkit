% bvar.ml.lgampdf - log density of the gamma distribution at x (shape a, RATE
% b; elementwise over array inputs).
%
% b is a rate. gamfit returns a scale, so a caller fitting with gamfit must
% invert it first (ckappa_hat = [a; 1/scale]).
%
% See:
% Chan, J.C.C. (2023). Comparing stochastic volatility specifications for large
% Bayesian VARs, Journal of Econometrics, 235(2), 1419-1446.

function lden = lgampdf(x,a,b)
    lden = a.*log(b) -gammaln(a) +(a-1).*log(x) -b.*x;
end
