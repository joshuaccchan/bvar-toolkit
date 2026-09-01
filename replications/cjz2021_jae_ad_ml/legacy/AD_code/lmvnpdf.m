% % =======================================================================
% % support function: evaluate the log-density of a multivariate normal
% % distribution
% %
% % See Chan, J.C.C. and Eisenstat, E. (2015). "Marginal Likelihood Estimation
% % with the Cross-Entropy Method," Econometric Reviews, 34(3), 256-285.
% %
% % (c) 2013, Joshua Chan. Email: joshuacc.chan@gmail.com
% % =======================================================================
function den = lmvnpdf(X,mu,sigma2)
N = length(X);
err = X-mu;
kernel = -err'*(sparse(sigma2)\err)/2;
den = (-N/2)*log(2*pi) - 1/2*log(det(sigma2))+kernel;
end
