% The function performs the QR decomposition such that the diagonals of R 
% are normalized to be positive
%
% See:
% Chan, J.C.C. (2021). Asymmetric conjugate priors for large Bayesian VARs,
% Quantitative Economics, forthcoming.

function [Q,R] = QR(A)
m = size(A,1);
[Q,R] = qr(A);
Q = Q*sparse(1:m,1:m,sign(diag(R)));
R = sparse(1:m,1:m,sign(diag(R)))*R;
end
    