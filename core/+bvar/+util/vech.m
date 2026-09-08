% Body from chan2023_joe_mlvarsv/legacy/utility/vech.m.
% Equivalence: tests/unit/test_util_small.m. Record: tests/variant_map.md.
function y=vech(Y)
y = nonzeros(tril(Y));
end
