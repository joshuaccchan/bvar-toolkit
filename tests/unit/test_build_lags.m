function test_build_lags
% hand-checked small case reproducing the legacy inline construction
Yfull = [1 10; 2 20; 3 30; 4 40; 5 50];
[Y, Z] = bvar.util.build_lags(Yfull, 2);
assert(isequal(Y, [3 30; 4 40; 5 50]), 'build_lags: wrong Y');
Zexp = [1 2 20 1 10;
        1 3 30 2 20;
        1 4 40 3 30];
assert(isequal(Z, Zexp), 'build_lags: wrong Z (intercept-first, lag-1-block-first)');

% agrees with the inline pattern used across the legacy packages
rng(7, 'twister');
Yr = randn(30, 4); p = 3; [T0, n] = size(Yr);
tmpY = Yr; T = T0 - p; Z2 = zeros(T, n*p);
for ii = 1:p
    Z2(:, (ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii, :);
end
Z2 = [ones(T,1) Z2];
[~, Zc] = bvar.util.build_lags(Yr, p);
assert(isequal(Zc, Z2), 'build_lags: differs from the legacy inline construction');
end
