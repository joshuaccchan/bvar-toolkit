% % =======================================================================
% % support function: evaluate the log-density of an inverse-gamma
% % distribution
% %
% % See Chan, J.C.C. and Eisenstat, E. (2015). "Marginal Likelihood Estimation
% % with the Cross-Entropy Method," Econometric Reviews, 34(3), 256-285.
% %
% % (c) 2013, Joshua Chan. Email: joshuacc.chan@gmail.com
% % =======================================================================

function lden = linvgammpdf(y,a,b)
lden = a.*log(b) - gammaln(a) - (a+1) .* log(y) - b./y;
end
