% bvar.util.surform - T x Tk BLOCK-DIAGONAL sparse SUR expansion (TVP state
% stacking).
%
%   Xout = bvar.util.surform(X)
%
%   X    : T x k matrix of regressors
%   Xout : T x T*k sparse matrix whose row t holds X(t,:) in columns
%          (t-1)*k + (1:k), so that Xout*[b_1; b_2; ...; b_T] returns the
%          fitted series of the TVP regression y_t = X(t,:)*b_t
%
% NOT the same operator as bvar.util.surform2, which is the n-row Kronecker
% expansion kron(speye(n), X(t,:)) for stacked-vector VARs.
%
% Provenance and the legacy copies this stands in for: tests/variant_map.md.
% Equivalence: tests/unit/test_surform.m.
%
% See:
% Chan, J.C.C. (2023). Large Hybrid Time-Varying Parameter VARs, Journal of
% Business and Economic Statistics, 41(3): 890-905

function Xout = surform( X )
[r,c] = size( X );
idi = kron((1:r)',ones(c,1));
idj = (1:r*c)';
Xout = sparse(idi,idj,reshape(X',r*c,1));
end