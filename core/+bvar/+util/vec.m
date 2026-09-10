% bvar.util.vec - column-stacking operator: vec(Y) = Y(:).
%
% Provenance and the legacy copies this stands in for: tests/variant_map.md.
% Equivalence: tests/unit/test_util_small.m.
function y=vec(Y)
y=Y(:);
end
