% bvar.util.surform2 - sparse SUR expansion with rows kron(speye(n), X(i,:)) for stacked-vector
% VAR sampling. NOT the same operator as bvar.util.surform (block-diagonal T x Tk).
% Body from chan2023_joe_mlvarsv/legacy/utility/SURform2.m, renamed; stands in
% for the springer_largebvar, jbes_kronecker and oisv copies, which differ only
% in variable spelling and whitespace.
% Equivalence: tests/unit/test_surform2.m. Record: tests/variant_map.md.
function Xout = surform2( X, n )
repX = kron(X,ones(n,1));
[r,c] = size( X );
idi = kron((1:r*n)',ones(c,1));
idj = repmat((1:n*c)',r,1);
Xout = sparse(idi,idj,reshape(repX',n*r*c,1));
end