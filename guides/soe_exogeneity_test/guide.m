%% guide - the code of guides/soe_exogeneity_test/README.md: the log Bayes factor for block
% exogeneity in the VAR of guides/soe_var, in closed form. Section 3 is the code to adapt.

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(fileparts(here)), 'setup.m'));
rng(1, 'twister');
names = {'AU GDP growth', 'AU CPI inflation', 'AU interest rate', 'AUD depreciation'};

%% 1. the true model of guides/soe_var: the VAR with block exogeneity fitted to the data
D = readmatrix(fullfile(fileparts(here), 'data', 'soe_Q.csv'));
data = D(:, 2:9);                         % 1984Q1-2019Q4: the US block, then Australia
nus = 4; nau = 4; n = nus + nau;
us = 1:nus; au = nus+1:n;
lev = [3 4 7];                            % the variables in levels: the two rates and the premium
c3 = 10^2; p = 3;                         % the intercepts, and the lag length of guides/soe_var
[Y, Z] = bvar.util.build_lags(data, p);   % x_t = (1, y_{t-1}', ..., y_{t-p}')'
[T, k] = size(Z);
ius = [1, reshape(1 + us' + (0:p-1)*n, 1, [])];   % the intercept and the US lags in x_t
iau = setdiff(2:k, ius);                  % the Australian lags in x_t
Xus = Z(:, ius); Xau = [Z, Y(:, us)];
Bus = Xus\Y(:, us);                       % least squares, block by block
Bau = Xau\Y(:, au);
Uus = Y(:, us) - Xus*Bus;  Sus = Uus'*Uus/(T - size(Xus, 2));
Uau = Y(:, au) - Xau*Bau;  Om = Uau'*Uau/(T - size(Xau, 2));
[Atrue, Sigtrue] = bvar.util.soe_reduced_form(Bus, Bau, Sus, Om);

%% 2. the sample of guides/soe_var
q0 = max(p, 4);                           % presample quarters: p lags, four for the AR(4) fits
Yfull = bvar.util.simulate_var(Atrue, Sigtrue, T + q0);     % the presample, then y_1, ..., y_T
[Ys, Zs] = bvar.util.build_lags(Yfull(q0-p+1:end, :), p);   % y_1, ..., y_T and their regressors

%% 3. the two marginal likelihoods of the US block: the code to adapt
res = bvar.models.var_soe(Yfull(1:q0, :), Ys, p, nus, 'lev', lev, 'c3', c3, 'nsim', 0);
c1 = res.c1;                              % the shrinkage of guides/soe_var
s2 = res.s2;                              % AR(4) variances of the sample
VA = res.VA; B0 = res.B0;                 % the prior variances and means on x_t
nu0 = nus + 2; S0 = diag(s2(us));         % Sig_us ~ IW(nu0, S0) in both models
lml_r = bvar.ml.ncp(Ys(:, us), Zs(:, ius), VA(ius), nu0, S0, B0(ius, us));   % block exogeneity
lml_u = bvar.ml.ncp(Ys(:, us), Zs, VA, nu0, S0, B0(:, us));     % with the Australian lags
lbf = lml_r - lml_u;                      % the log Bayes factor for block exogeneity
assert(abs(lbf - log_bf(Yfull, p, nus, lev, c3, 1)) < 1e-10);

%% 4. the prior variances of the Australian lags
fs = 10.^(-2:3);                          % multiples of the Australian lags' prior variances
lbf_f = log_bf(Yfull, p, nus, lev, c3, fs);             % the sample
lbf_data_f = log_bf(data, p, nus, lev, c3, fs);         % the data
lbf_data = lbf_data_f(fs == 1);

%% 5. one Australian variable at a time: the adaptation
lbf_j = zeros(1, nau);
for j = 1:nau
    keep = setdiff(1:k, 1 + au(j) + (0:p-1)*n);         % x_t without the lags of variable j
    lbf_j(j) = bvar.ml.ncp(Ys(:, us), Zs(:, keep), VA(keep), nu0, S0, B0(keep, us)) - lml_u;
end

%% the results
fprintf(['\nthe sample: %d quarters; each US equation has %d coefficients, %d with the ' ...
    'Australian lags; c1 = %.3f\n'], T, numel(ius), k, c1);
fprintf('log marginal likelihood of the US block: %.2f with block exogeneity, %.2f without\n', ...
    lml_r, lml_u);
fprintf('log Bayes factor for block exogeneity %.2f; posterior probability 1 - %.1e\n', lbf, ...
    1/(1 + exp(lbf)));
fprintf('on the data: log Bayes factor %.2f\n', lbf_data);
fprintf('\nmultiples of the prior variances of the Australian lags: %s\n', mat2str(fs));
fprintf('  log Bayes factor, the sample: %s\n', mat2str(round(lbf_f, 2)));
fprintf('  log Bayes factor, the data:   %s\n', mat2str(round(lbf_data_f, 2)));
fprintf('  rise per tenfold increase, the sample: %s\n', mat2str(round(diff(lbf_f), 2)));
fprintf('  rise per tenfold increase, the data:   %s\n', mat2str(round(diff(lbf_data_f), 2)));
fprintf('  its limit, half of the %d coefficients times log(10): %.2f\n', nus*numel(iau), ...
    nus*numel(iau)/2*log(10));
fprintf('\none Australian variable at a time, log Bayes factor for excluding its lags:\n');
for j = 1:nau
    fprintf('  %-17s %6.2f\n', names{j}, lbf_j(j));
end

%% the checks
% Bayes' theorem at the true values:
% log p(Y) = log p(Y | B, Sig) + log p(B, Sig) - log p(B, Sig | Y)
[~, Bh_r, K_r, Sh_r] = bvar.ml.ncp(Ys(:, us), Zs(:, ius), VA(ius), nu0, S0, B0(ius, us));
[~, Bh_u, K_u, Sh_u] = bvar.ml.ncp(Ys(:, us), Zs, VA, nu0, S0, B0(:, us));
B_r = Atrue(ius, us); B_u = Atrue(:, us); Sig0 = Sigtrue(us, us);
id_r = llike(Ys(:, us), Zs(:, ius), B_r, Sig0) ...
    + bvar.ml.lniwpdf(B_r, Sig0, B0(ius, us), sparse(diag(1./VA(ius))), nu0, S0) ...
    - bvar.ml.lniwpdf(B_r, Sig0, Bh_r, K_r, nu0 + T, Sh_r);
id_u = llike(Ys(:, us), Zs, B_u, Sig0) ...
    + bvar.ml.lniwpdf(B_u, Sig0, B0(:, us), sparse(diag(1./VA)), nu0, S0) ...
    - bvar.ml.lniwpdf(B_u, Sig0, Bh_u, K_u, nu0 + T, Sh_u);
fprintf('\nBayes'' theorem at the true values against the closed form: gaps %.1e and %.1e\n', ...
    id_r - lml_r, id_u - lml_u);
% the unrestricted VAR of guides/soe_var, under one natural conjugate prior for the
% eight equations with Sig ~ IW(n + 2, diag(s2)): the posterior of its US equations has the
% mean and the scale matrix of the US block with the Australian lags
[~, Bh_var, ~, Sh_var] = bvar.ml.ncp(Ys, Zs, VA, n + 2, diag(s2), B0);
gap_var = [max(abs(Bh_var(:, us) - Bh_u), [], 'all'), max(abs(Sh_var(us, us) - Sh_u), [], 'all')];
fprintf(['the unrestricted VAR of guides/soe_var, US equations: posterior mean ' ...
    'gap %.1e, scale matrix gap %.1e\n'], gap_var);
assert(max(abs([id_r - lml_r, id_u - lml_u, gap_var])) < 1e-8, ...
    'the closed form disagrees with the identities of the marginal likelihood');

%% local functions
function lbf = log_bf(Yfull, p, nus, lev, c3, f)
% the log Bayes factor for block exogeneity, one for each entry of f, from a sample whose first
% max(p, 4) rows are the presample: the log marginal likelihood of the US block without the
% Australian lags, minus that with them when f multiplies their prior variances; c1 maximizes
% the marginal likelihood of the VAR with block exogeneity
q0 = max(p, 4);
res = bvar.models.var_soe(Yfull(1:q0, :), Yfull(q0+1:end, :), p, nus, 'lev', lev, 'c3', c3, ...
    'nsim', 0);
[Ys, Zs] = bvar.util.build_lags(Yfull(q0-p+1:end, :), p);
n = size(Ys, 2); k = size(Zs, 2); us = 1:nus;
VA = res.VA; B0 = res.B0;
nu0 = nus + 2; S0 = diag(res.s2(us));
ius = [1, reshape(1 + us' + (0:p-1)*n, 1, [])];
iau = setdiff(2:k, ius);
lml_r = bvar.ml.ncp(Ys(:, us), Zs(:, ius), VA(ius), nu0, S0, B0(ius, us));
lbf = zeros(size(f));
for j = 1:numel(f)
    VAu = VA; VAu(iau) = f(j)*VA(iau);
    lbf(j) = lml_r - bvar.ml.ncp(Ys(:, us), Zs, VAu, nu0, S0, B0(:, us));
end
end

function l = llike(Y, X, B, Sig)
% the log likelihood of the multivariate regression Y = X*B + E with rows of E ~ N(0, Sig)
[T, m] = size(Y);
C = chol(Sig, 'lower');
E = (Y - X*B)/C';
l = -T*m/2*log(2*pi) - T*sum(log(diag(C))) - .5*sum(E(:).^2);
end
