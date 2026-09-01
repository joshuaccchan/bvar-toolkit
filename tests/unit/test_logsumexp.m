function test_logsumexp
x = [0 1; 2 3];
assert(abs(bvt.util.logsumexp(x(:,1)) - log(exp(0)+exp(2))) < 1e-12, 'logsumexp: wrong value');
y = bvt.util.logsumexp(x, 2);
assert(abs(y(1) - log(exp(0)+exp(1))) < 1e-12, 'logsumexp: wrong along dim 2');

big = [1000; 1001];   % naive overflow case
assert(abs(bvt.util.logsumexp(big) - (1001 + log(1+exp(-1)))) < 1e-12, 'logsumexp: overflow handling');
assert(bvt.util.logsumexp([-Inf; -Inf]) == -Inf, 'logsumexp: all -Inf slice');
end
