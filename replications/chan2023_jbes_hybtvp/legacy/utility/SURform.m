% Support function for estimating the hybrid TVP-VAR in Chan (2022)
%
% See:
% Chan, J.C.C. (2022). Large Hybrid Time-Varying Parameter VARs, Journal of
% Business and Economic Statistics, forthcoming

function Xout = SURform( X )
[r,c] = size( X );
idi = kron((1:r)',ones(c,1));
idj = (1:r*c)';
Xout = sparse(idi,idj,reshape(X',r*c,1));
end