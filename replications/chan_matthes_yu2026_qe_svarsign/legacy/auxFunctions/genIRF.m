function [irf,gradIrf] = genIRF(q,vma,neg)
% Generate IRF of the ith variable to the first structural shock at the hth
% horizon. Also generates analytical gradient of the IRF.
% Inputs:
% q: first column of Q
% vma: row vector of coefficients in orthogonal reduced-form VMA 
% representation for ith variable at hth horizon
% neg: if neg = -1, return negative of IRF; if neg = 1, return positive

irf = neg*vma*q; % IRF of ith variable to first shock
gradIrf = neg*vma'; % Gradient 

end