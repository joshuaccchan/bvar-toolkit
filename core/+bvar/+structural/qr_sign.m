% bvar.structural.qr_sign - QR decomposition with the diagonal of R normalized
% to be positive, which makes the decomposition unique and gives the orthogonal
% factor Q used to draw rotations in sign-restricted SVARs.
%
%   [Q,R] = bvar.structural.qr_sign(A)
%
% MATLAB's qr returns a decomposition unique only up to the signs of the columns
% of Q and the rows of R. Multiplying both by diag(sign(diag(R))) fixes those
% signs, so a random A gives a rotation drawn from the uniform (Haar) measure on
% the orthogonal group rather than an arbitrary member of a sign-equivalence
% class. Callers draw A = randn(n,n) and use Q to rotate the Cholesky factor of
% the reduced-form covariance.
%
% See:
% Chan, J.C.C. (2022). Asymmetric Conjugate Priors for Large Bayesian VARs,
% Quantitative Economics, 13(3): 1145-1169.

function [Q,R] = qr_sign(A)
m = size(A,1);
[Q,R] = qr(A);
Q = Q*sparse(1:m,1:m,sign(diag(R)));
R = sparse(1:m,1:m,sign(diag(R)))*R;
end
