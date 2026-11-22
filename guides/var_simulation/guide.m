%% guide - the code of guides/var_simulation/README.md: a VAR(4) fitted to US GDP growth, PCE
% inflation and the unemployment rate, the true model of the series Forecasting with a VAR, and
% samples simulated from it. Sections 1 to 3 are the code to adapt. Saves fig_simulated.png next
% to this file.

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(fileparts(here)), 'setup.m'));
names = {'GDP growth', 'PCE inflation', 'unemployment'};

%% 1. the VAR(4) fitted to the data: the code to adapt
D = readmatrix(fullfile(fileparts(here), 'data', 'macro3_Q.csv'));
data = D(:, 2:4);                         % GDP growth, PCE inflation, unemployment, 1960Q1-2019Q4
p = 4;
[Y, Z] = bvar.util.build_lags(data, p);   % y_t from 1961Q1 on; x_t = (1, y_{t-1}', ..., y_{t-4}')'
[T, n] = size(Y); k = size(Z, 2);
A = Z\Y;                                  % least squares, A = (Z'Z)^{-1} Z'Y
U = Y - Z*A;
Sig = U'*U/(T - k);                       % Sig = U'U/(T - k)

%% 2. stability and the stationary distribution: the code to adapt
F = [A(2:end, :)'; eye(n*(p-1)), zeros(n*(p-1), n)];   % s_t = c + F s_{t-1} + e_t
rho = max(abs(eig(F)));                                % the largest modulus; stable if below one
c = [A(1, :)'; zeros(n*(p-1), 1)];
mu = (eye(n*p) - F)\c;                                 % mu = (I - F)^{-1} c
Q = blkdiag(Sig, zeros(n*(p-1)));
V = reshape((eye((n*p)^2) - kron(F, F))\Q(:), n*p, n*p);   % vec(V) = (I - F kron F)^{-1} vec(Q)

%% 3. a simulated sample: the code to adapt
rng(2, 'twister');
[Ysim, Y0sim] = bvar.util.simulate_var(A, Sig, 240);   % 240 quarters, as many as in the data

%% the results
pop = [mu(1:n)'; sqrt(diag(V(1:n, 1:n)))'; diag(F(1:n, :)*V(:, 1:n))'./diag(V(1:n, 1:n))'];
fprintf('\nVAR(%d) of %d variables by least squares, T = %d\n', p, n, T);
fprintf('largest modulus of the companion eigenvalues: %.3f\n', rho);
fprintf(['residual sd %s; correlations: GDP-inflation %.2f, GDP-unemployment %.2f, ' ...
    'inflation-unemployment %.2f\n'], ...
    mat2str(sqrt(diag(Sig))', 3), corr_of(Sig, 1, 2), corr_of(Sig, 1, 3), corr_of(Sig, 2, 3));
fprintf('own lag-1 coefficients: %s\n', mat2str(diag(A(2:n+1, :))', 3));
fprintf('\n%-14s %22s %22s %22s\n', '', 'mean', 'sd', 'autocorrelation');
fprintf('%-14s %22s %22s %22s\n', '', 'data/sim/stationary', 'data/sim/stationary', ...
    'data/sim/stationary');
md = moments(data); ms = moments(Ysim);
for i = 1:n
    fprintf('%-14s %6.2f %6.2f %6.2f    %6.2f %6.2f %6.2f    %6.2f %6.2f %6.2f\n', ...
        names{i}, ...
        md(1, i), ms(1, i), pop(1, i), md(2, i), ms(2, i), pop(2, i), md(3, i), ms(3, i), pop(3, i));
end

fprintf('lowest simulated unemployment rate: %.2f, in quarter %d\n', min(Ysim(:, 3)), ...
    find(Ysim(:, 3) == min(Ysim(:, 3))));
qname = @(i) sprintf('%dQ%d', 1960 + floor((i - 1)/4), mod(i - 1, 4) + 1);
iu = find(data(:, 3) == min(data(:, 3)));
fprintf('lowest unemployment rate in the data: %.2f, in %s\n', min(data(:, 3)), ...
    strjoin(arrayfun(qname, iu', 'UniformOutput', false), ' and '));
irf = zeros(13, 1); S = zeros(n*p, 1); S(3) = 1;        % unemployment's response to its own shock
for h = 0:12
    irf(h+1) = S(3); S = F*S;
end
[~, hpk] = max(irf);
fprintf('response of unemployment to its own shock: peak %.2f after %d quarters, %.2f after 12\n', ...
    max(irf), hpk - 1, irf(end));

%% the lag length: information criteria on a common sample
pmax = 8; crit = zeros(pmax, 2);
for q = 1:pmax
    [Yq, Zq] = bvar.util.build_lags(data(pmax-q+1:end, :), q);   % the same quarters for every q
    Uq = Yq - Zq*(Zq\Yq); Tq = size(Yq, 1); kq = size(Zq, 2);
    ld = log(det(Uq'*Uq/Tq));                                    % log|Sig_hat|, Sig_hat = U'U/T
    crit(q, :) = [ld + 2*n*kq/Tq, ld + log(Tq)*n*kq/Tq];         % AIC, BIC
end
[~, best] = min(crit);
fprintf('\nlag length chosen on 1962Q1-2019Q4: AIC %d, BIC %d\n', best);

%% 4. 1,000 additional samples of 240 quarters
rng(1, 'twister');
nrep = 1000; mrep = zeros(3, n, nrep);
for r = 1:nrep
    mrep(:, :, r) = moments(bvar.util.simulate_var(A, Sig, 240));
end
lo = quantile(mrep, .05, 3); hi = quantile(mrep, .95, 3);
fprintf('\n90%% ranges over %d simulated samples of 240 quarters:\n', nrep);
for i = 1:n
    fprintf(['%-14s mean %5.2f to %5.2f (data %5.2f), sd %4.2f to %4.2f (data %4.2f), ' ...
        'autocorrelation %4.2f to %4.2f (data %4.2f)\n'], names{i}, lo(1, i), hi(1, i), md(1, i), ...
        lo(2, i), hi(2, i), md(2, i), lo(3, i), hi(3, i), md(3, i));
end
fprintf('every moment of the data inside its 90%% range: %d\n', all(md >= lo & md <= hi, 'all'));

%% 5. the checks
gapV = max(abs(V - F*V*F' - Q), [], 'all')/max(abs(V(:)));
rng(3, 'twister');
Ylong = bvar.util.simulate_var(A, Sig, 200000);
ml = moments(Ylong);
fprintf('\nV against F*V*F'' + Q: %.1e (relative)\n', gapV);
fprintf(['200,000 simulated quarters against the stationary values: means %.3f, sds %.3f, ' ...
    'autocorrelations %.4f (largest gaps)\n'], ...
    max(abs(ml(1, :) - pop(1, :))), max(abs(ml(2, :) - pop(2, :))), max(abs(ml(3, :) - pop(3, :))));
assert(gapV < 1e-10, 'the stationary covariance does not solve V = F*V*F'' + Q');
assert(max(abs(ml(1, :) - pop(1, :))) < .1 && max(abs(ml(2, :)./pop(2, :) - 1)) < .03, ...
    'the long simulation does not match the stationary distribution');

%% 6. 20 quarters from the data of 2019: the adaptation
H = 20; R = 10000;
paths = zeros(H, n, R);
rng(4, 'twister');
for r = 1:R
    paths(:, :, r) = bvar.util.simulate_var(A, Sig, H, data(end-p+1:end, :));   % from 2019Q1-2019Q4
end
fprintf('\n%d paths from the data of 2019 (2019Q4: %s):\n', R, mat2str(data(end, :), 3));
for h = [4 20]
    qh = quantile(squeeze(paths(h, :, :))', [.05 .5 .95]);
    fprintf('  h = %2d: medians %s, 90%% intervals %s to %s\n', h, mat2str(qh(2, :), 3), ...
        mat2str(qh(1, :), 3), mat2str(qh(3, :), 3));
end
fprintf('  mean at h = 20 %s against the stationary means %s\n', ...
    mat2str(mean(squeeze(paths(H, :, :)), 2)', 3), mat2str(mu(1:n)', 3));

%% the figure
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 18 13]);
q = 1960 + (0:239)/4;
for i = 1:n
    yl = [min([data(:, i); Ysim(:, i)]), max([data(:, i); Ysim(:, i)])];
    subplot(n, 2, 2*i - 1);
    plot(q, data(:, i), 'k'); box off; xlim(q([1 end])); ylim(yl);
    title([names{i} ', 1960-2019'], 'FontWeight', 'normal');
    subplot(n, 2, 2*i);
    plot(1:240, Ysim(:, i), 'k'); box off; xlim([1 240]); ylim(yl);
    title([names{i} ', simulated'], 'FontWeight', 'normal');
end
exportgraphics(fig, fullfile(here, 'fig_simulated.png'), 'Resolution', 150); close(fig);

function m = moments(X)
% the mean, standard deviation and first-order autocorrelation of each column
Xc = X - mean(X);
m = [mean(X); std(X); sum(Xc(2:end, :).*Xc(1:end-1, :))./sum(Xc.^2)];
end

function r = corr_of(S, i, j)
r = S(i, j)/sqrt(S(i, i)*S(j, j));
end
