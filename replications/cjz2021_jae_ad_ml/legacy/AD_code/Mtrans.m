% This is written by Dan Zhu(dan.zhu@monash.edu) on the 6th of June 2019
% This function computes matrix transporse, 
% 
function [XT] = Mtrans(X, K_nq)
XT.v=X.v';
if isempty(K_nq)
[r,c]=size(X.v);
K_nq=commutation_matrix(r,c);
end
XT.d = K_nq * X.d;