% bvar.util.vec - column-stacking operator: vec(Y) = Y(:).
%
% Body from chan2023_joe_mlvarsv/legacy/utility/vec.m.
% Equivalence: tests/unit/test_util_small.m. Record: tests/variant_map.md.
function y=vec(Y)
y=Y(:);
end
