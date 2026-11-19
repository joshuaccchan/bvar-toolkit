function test_hist_decomp
% hist_decomp: the pieces add up to the data, the shocks are the simulated ones, the
% contributions are sums of responses times shocks, the base path, one impact matrix per period,
% a scalar AR(1), and its errors
n = 3; p = 2; T = 60;
A = [.2 -.1 .3; .5 .1 0; .2 .4 .1; 0 .1 .3; -.1 0 .05; .05 -.1 0; 0 .05 .1];   % k x n
Sig = [1 .3 .2; .3 .8 .1; .2 .1 .5];
L = chol(Sig, 'lower');
Y0 = [.5 -.2 .1; .3 .1 -.4];                % y_{-1}, y_0

% simulate_var turns randn(n,1) draws into errors L*z_t, so the shocks are those draws
rng(21, 'twister');
Y = bvar.util.simulate_var(A, Sig, T, Y0);
rng(21, 'twister');
E = randn(n, T)';
[contrib, base, shocks] = bvar.structural.hist_decomp([Y0; Y], A, L);
assert(isequal(size(contrib), [T n n]) && isequal(size(base), [T n]), 'hist_decomp: sizes');
assert(max(abs(shocks - E), [], 'all') < 1e-12, 'hist_decomp: shocks differ from the draws');
assert(max(abs(base + sum(contrib, 3) - Y), [], 'all') < 1e-10, ...
    'hist_decomp: the pieces do not add up to the data');

% each contribution is the sum over s of the response at horizon s times the shock at t - s
resp = bvar.structural.irf_redu(A(2:end, :), L, T, n);
want = zeros(T, n, n);
for t = 1:T
    for s = 0:t-1
        want(t, :, :) = want(t, :, :) + reshape(resp(:, :, s+1).*E(t-s, :), 1, n, n);
    end
end
assert(max(abs(contrib - want), [], 'all') < 1e-10, 'hist_decomp: differs from the responses');

% the base path runs the VAR from the initial values with no errors
Yb = [Y0; zeros(T, n)];
for t = p+1:p+T
    Yb(t, :) = A(1, :) + Yb(t-1, :)*A(2:4, :) + Yb(t-2, :)*A(5:7, :);
end
assert(max(abs(base - Yb(p+1:end, :)), [], 'all') < 1e-12, 'hist_decomp: wrong base path');

% one impact matrix per period: errors L_t*z_t with L_t = s_t*L
s = exp(.3*sin((1:T)'/5));
Lt = L.*reshape(s, 1, 1, T);
Ys = [Y0; zeros(T, n)];
for t = p+1:p+T
    Ys(t, :) = A(1, :) + Ys(t-1, :)*A(2:4, :) + Ys(t-2, :)*A(5:7, :) + (Lt(:, :, t-p)*E(t-p, :)')';
end
[contrib, base, shocks] = bvar.structural.hist_decomp(Ys, A, Lt);
Phi = bvar.structural.irf_redu(A(2:end, :), eye(n), T, n);
want = zeros(T, n, n);
for t = 1:T
    for r = 0:t-1
        want(t, :, :) = want(t, :, :) + ...
            reshape(Phi(:, :, r+1)*Lt(:, :, t-r).*E(t-r, :), 1, n, n);
    end
end
assert(max(abs(shocks - E), [], 'all') < 1e-12, 'hist_decomp: shocks with L_t');
assert(max(abs(contrib - want), [], 'all') < 1e-10, 'hist_decomp: contributions with L_t');
assert(max(abs(base + sum(contrib, 3) - Ys(p+1:end, :)), [], 'all') < 1e-10, ...
    'hist_decomp: the pieces do not add up with L_t');

% a scalar AR(1): contribution sig*sum_s phi^s z_{t-s}, base c(1 - phi^t)/(1 - phi) + phi^t y_0
c = .5; phi = .8; sig = .7; y0 = 2;
z = E(:, 1);
y = filter(1, [1 -phi], c + sig*z, phi*y0);
[contrib, base] = bvar.structural.hist_decomp([y0; y], [c; phi], sig);
t = (1:T)';
assert(max(abs(contrib - filter(1, [1 -phi], sig*z))) < 1e-12, 'hist_decomp: AR(1) contribution');
assert(max(abs(base - (c*(1 - phi.^t)/(1 - phi) + phi.^t*y0))) < 1e-12, 'hist_decomp: AR(1) base');

% errors: A of the wrong size, L with the wrong number of periods
assert_error(@() bvar.structural.hist_decomp([Y0; Y], A(1:6, :), L), 'hist_decomp:badA');
assert_error(@() bvar.structural.hist_decomp([Y0; Y], A, Lt(:, :, 1:T-1)), 'hist_decomp:badL');
end

function assert_error(f, id)
try
    f();
catch err
    assert(strcmp(err.identifier, id), 'expected %s, got %s', id, err.identifier);
    return
end
error('test_hist_decomp:noError', 'expected the error %s', id);
end
