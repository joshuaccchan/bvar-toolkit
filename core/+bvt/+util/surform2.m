% bvt.util.surform2 - sparse SUR expansion with rows kron(speye(n), X(i,:)) for stacked-vector
% VAR sampling. NOT the same operator as bvt.util.surform (block-diagonal T x Tk).
% Extracted 2026-09-01 (step 3, zero-risk core): verified identical (modulo comments/whitespace) across chan2023_joe_mlvarsv/legacy/utility/SURform2.m (canonical, this copy),
% chan2020_springer_largebvar/legacy/SURform2.m, chan2020_jbes_kronecker/legacy/
% realtime_forecasts/SURform2.m, and chan_koop_yu2024_jbes_oisv/legacy/utility/SURform2.m
% (differences: output variable spelling, [r c] vs [r,c], indentation).
% Function renamed SURform2 -> surform2.
function Xout = surform2( X, n )
repX = kron(X,ones(n,1));
[r,c] = size( X );
idi = kron((1:r*n)',ones(c,1));
idj = repmat((1:n*c)',r,1);
Xout = sparse(idi,idj,reshape(repX',n*r*c,1));
end