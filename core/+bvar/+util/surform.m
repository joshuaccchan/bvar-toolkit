% bvar.util.surform - T x Tk BLOCK-DIAGONAL sparse SUR expansion (TVP state stacking).
% NOT the same operator as bvar.util.surform2 (n-row Kronecker expansion) - see headers.
% Extracted 2026-09-01 (step 3, zero-risk core): verified identical (modulo comments/whitespace) across chan2023_jbes_hybtvp/legacy/utility/SURform.m (canonical, this copy)
% and chan_jeliazkov2009_statespace/legacy/sp_code/SURform.m. Function renamed SURform -> surform.
% Support function for estimating the hybrid TVP-VAR in Chan (2022)
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