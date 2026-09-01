% The function performs the QR decomposition such that the diagonals of R 
% are normalized to be positive
%
% See:
% Chan, J.C.C. (2022). Asymmetric Conjugate Priors for Large Bayesian VARs,
% Quantitative Economics, 13(3): 1145-1169

function [Q,R] = QR(A)
m = size(A,1);
[Q,R] = qr(A);
Q = Q*sparse(1:m,1:m,sign(diag(R)));
R = sparse(1:m,1:m,sign(diag(R)))*R;
end
    