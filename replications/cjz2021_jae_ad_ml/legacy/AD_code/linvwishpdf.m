% % =======================================================================
% % support function: evaluate the log-density of a inverse-Wishart
% % distribution
% %
% % See Chan, J.C.C. and Eisenstat, E. (2015). "Marginal Likelihood Estimation
% % with the Cross-Entropy Method," Econometric Reviews, 34(3), 256-285.
% %
% % (c) 2013, Joshua Chan. Email: joshuacc.chan@gmail.com
% % =======================================================================
function den = linvwishpdf(W,nu,S)
K = size(S,1);
sum1 = sum(gammaln((nu+1-(1:K))/2));
const = nu*K/2*log(2) + K*(K-1)/4*log(pi) + sum1;
den = -const + nu/2 * log(det(S)) - (nu + K + 1)/2 * log(det(W)) - trace(S/W)/2;
end

