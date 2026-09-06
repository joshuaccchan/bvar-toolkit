function test_eq_var_oi
% bvar.samplers.eq_var_oi draws the same conditional as eq_svar_oi with the
% same rng consumption: under a common seed the two draws must agree to
% floating-point precision, across several (n, T, p) shapes, with
% heterogeneous log-variances, a non-diagonal B0 and heterogeneous prior
% variances. Then: the prior-mean reparameterization its header documents, the
% input guards, and a loose floor on the speedup that is the function's reason
% for existing.

shapes = [3 60 1; 6 120 2; 10 90 4; 24 241 4];     % n, T, p
for s = 1:size(shapes, 1)
    n = shapes(s, 1); T = shapes(s, 2); p = shapes(s, 3); k = 1 + n*p;
    rng(20 + s, 'twister');
    [Y, X] = bvar.util.build_lags(randn(T + p, n), p);
    B0 = eye(n) + 0.3*randn(n);
    h = 0.6*randn(T, n);
    A = 0.1*randn(k, n);
    V = 0.05*rand(k, n) + 0.01; V(1, :) = 100;

    rng(100 + s, 'twister'); A1 = bvar.samplers.eq_svar_oi(Y, X, B0, h, A, V(:)); r1 = rng;
    rng(100 + s, 'twister'); A2 = bvar.samplers.eq_var_oi(Y, X, B0, h, A, V(:));  r2 = rng;
    d = max(abs(A1(:) - A2(:)))/max(abs(A1(:)));
    % measured 4e-16 to 2e-15 on these shapes; 1e-12 leaves three orders of
    % headroom and still catches a defect far smaller than any real one
    assert(d < 1e-12, 'eq_var_oi: differs from eq_svar_oi (rel %.1e) at n = %d, T = %d, p = %d', d, n, T, p);
    assert(isequal(r1.State, r2.State), 'eq_var_oi: rng consumption differs from eq_svar_oi');
    assert(isequal(size(A2), [k n]), 'eq_var_oi: wrong output size');
end

% every column is updated (a stale column would still match A)
assert(~any(all(A2 == A, 1)), 'eq_var_oi: some column of A was not redrawn');

% --- the prior-mean reparameterization the header documents, and that the
%     out-of-tree consumer relies on: drawing A - A0 on Y - X*A0 and adding A0
%     back must equal an explicit nonzero-prior-mean sweep under the same seed
n = 6; T = 120; p = 2; k = 1 + n*p;
rng(7, 'twister');
[Y, X] = bvar.util.build_lags(randn(T + p, n), p);
B0 = eye(n) + 0.3*randn(n); h = 0.6*randn(T, n);
A = 0.1*randn(k, n); V = 0.05*rand(k, n) + 0.01; V(1, :) = 100;
A0 = zeros(k, n); A0(2:n+1, :) = 0.8*eye(n);        % prior mean on the first own lag
rng(11, 'twister');
Ashift = A0 + bvar.samplers.eq_var_oi(Y - X*A0, X, B0, h, A - A0, V(:));
rng(11, 'twister');                                  % explicit level-parameterization sweep
Aref = A; eh_inv = exp(-h);
for ii = 1:n
    Aref(:, ii) = 0;
    Etil = (Y - X*Aref)*B0'; b = B0(:, ii);
    w = eh_inv*(b.^2); c = (eh_inv.*Etil)*b; iVi = 1./V(:, ii);
    Kai = X'*(w.*X); Kai(1:k+1:end) = Kai(1:k+1:end) + iVi';
    CK = chol(Kai, 'lower');
    Aref(:, ii) = CK'\(CK\(X'*c + iVi.*A0(:, ii))) + CK'\randn(k, 1);
end
d = max(abs(Ashift(:) - Aref(:)))/max(abs(Aref(:)));
assert(d < 1e-12, 'eq_var_oi: the documented prior-mean shift is not exact (rel %.1e)', d);

% --- input guards ---
for bad = {{'badH', @() bvar.samplers.eq_var_oi(Y, X, B0, h', A, V(:))}, ...
           {'badB0', @() bvar.samplers.eq_var_oi(Y, X, B0(:, 1:end-1), h, A, V(:))}, ...
           {'badV', @() bvar.samplers.eq_var_oi(Y, X, B0, h, A, V(1:end-1))}}
    try
        bad{1}{2}();
        error('test:noThrow', 'eq_var_oi: %s was accepted', bad{1}{1});
    catch err
        assert(strcmp(err.identifier, ['bvar:samplers:eq_var_oi:' bad{1}{1}]), ...
            'eq_var_oi: wrong error id for %s (got %s)', bad{1}{1}, err.identifier);
    end
end

% --- speed: a loose floor, so that delegating back to the stacked form, or
%     otherwise reintroducing the O(T n k^2) arithmetic, fails the suite.
%     Measured ratios at this shape run 9x to 20x depending on machine load;
%     3x cannot be reached by the stacked form and cannot be flaky.
n = 24; T = 241; p = 4; k = 1 + n*p;
rng(31, 'twister');
[Y, X] = bvar.util.build_lags(randn(T + p, n), p);
B0 = eye(n) + 0.3*randn(n); h = 0.6*randn(T, n);
A = 0.1*randn(k, n); V = 0.05*rand(k, n) + 0.01; V(1, :) = 100;
nrep = 7;
for r = 1:3, bvar.samplers.eq_svar_oi(Y,X,B0,h,A,V(:)); bvar.samplers.eq_var_oi(Y,X,B0,h,A,V(:)); end
t1 = zeros(nrep,1); t2 = zeros(nrep,1);
for r = 1:nrep
    tic; bvar.samplers.eq_svar_oi(Y, X, B0, h, A, V(:)); t1(r) = toc;
    tic; bvar.samplers.eq_var_oi(Y, X, B0, h, A, V(:));  t2(r) = toc;
end
ratio = median(t1)/median(t2);
fprintf('      eq_var_oi at n = 24, T = 241, p = 4: eq_svar_oi %.3f s, eq_var_oi %.3f s per sweep (x%.1f)\n', ...
        median(t1), median(t2), ratio);
assert(ratio > 3, 'eq_var_oi: speedup collapsed to %.1fx - has the fast path been lost?', ratio);
end
