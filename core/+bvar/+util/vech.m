% bvar.util.vech - half-vectorization: the lower triangle of Y, stacked by column.
%
% NOTE the body is `nonzeros(tril(Y))`, so an exact zero in the lower triangle is
% DROPPED rather than stacked and the result is shorter than n(n+1)/2. Kept verbatim;
% every caller passes a matrix whose lower triangle has no exact zeros.
%
% Body from chan2023_joe_mlvarsv/legacy/utility/vech.m.
% Equivalence: tests/unit/test_util_small.m. Record: tests/variant_map.md.
function y=vech(Y)
y = nonzeros(tril(Y));
end
