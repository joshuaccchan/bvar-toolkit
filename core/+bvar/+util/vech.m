% bvar.util.vech - half-vectorization: the lower triangle of Y, stacked by column.
%
% NOTE the body is `nonzeros(tril(Y))`, so an exact zero in the lower triangle is
% DROPPED instead of stacked and the result is shorter than n(n+1)/2. Every caller
% in the toolkit passes a matrix whose lower triangle has no exact zeros.
function y=vech(Y)
y = nonzeros(tril(Y));
end
