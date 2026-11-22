%% guide - the code of guides/var_stacking/README.md: a VAR(2) of three variables on generated
% data, written as Y = Z*A + U and as y = X*beta + eps. The sections called "the code to adapt"
% contain everything a sampler uses. Saves fig_forms.png next to this file.

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(fileparts(here)), 'setup.m'));
rng(1, 'twister');
rel = @(a, b) max(abs(a(:) - b(:)))/max(abs(b(:)));

%% 1. the data: a VAR(2) of three variables
n = 3; p = 2; T = 200;
a = [.5; .2; .1];                    % intercepts
A1 = [.5 .1 0; .2 .4 .1; 0 .2 .5];   % lag 1
A2 = [.2 0 0; 0 .1 0; .1 0 .2];      % lag 2
Sig = [1 .3 .2; .3 .8 .1; .2 .1 .6];
A = [a'; A1'; A2'];   % k x n, A = (a, A_1, A_2)'
C = chol(Sig, 'lower'); burn = 100;
Yfull = zeros(burn + p + T, n);
for t = p+1:burn+p+T
    Yfull(t,:) = a' + Yfull(t-1,:)*A1' + Yfull(t-2,:)*A2' + (C*randn(n,1))';
end
Yfull = Yfull(burn+1:end, :);   % y_{-1}, y_0, y_1, ..., y_T

%% 2. the two forms: the code to adapt
[Y, Z] = bvar.util.build_lags(Yfull, p);   % T x n and T x k, k = 1 + n*p
k = size(Z, 2);
Ahat = Z\Y;   % least squares, one column per equation
U = Y - Z*Ahat;
Sighat = U'*U/T;

y = reshape(Y', T*n, 1);        % y = vec(Y'): y_1, ..., y_T stacked
X = bvar.util.surform2(Z, n);   % rows (t-1)*n+1:t*n hold kron(I_n, x_t')
beta = A(:);                    % beta = vec(A)

%% 2. the checks
loglik_mat = @(A, Sig) -T*n/2*log(2*pi) - T/2*bvar.util.ldet(Sig) ...
    - .5*trace(Sig\((Y - Z*A)'*(Y - Z*A)));   % the matrix form
loglik_vec = @(b, Sig) -T*n/2*log(2*pi) - T/2*bvar.util.ldet(Sig) ...
    - .5*(y - X*b)'*kron(speye(T), inv(Sig))*(y - X*b);   % the stacked form
Zbyhand = [ones(T,1), Yfull(p:end-1,:), Yfull(p-1:end-2,:)];
lam = abs(eig([A1 A2; eye(n) zeros(n)]));   % the companion matrix
fprintf('\na VAR(%d) of %d variables, T = %d: k = %d regressors per equation, %d coefficients\n', ...
    p, n, T, k, n*k);
fprintf('largest modulus of the companion eigenvalues: %.3f\n', max(lam));
dYZ = [max(abs(Y - Yfull(p+1:end,:)), [], 'all'), max(abs(Z - Zbyhand), [], 'all')];
fprintf('Y and Z against direct indexing: %g, %g (largest difference)\n', dYZ);
fprintf('X: %d x %d with %d nonzero entries, %d per period\n', size(X,1), size(X,2), nnz(X), nnz(X)/T);
gapX = rel(X*beta, reshape((Z*A)', [], 1));
fprintf('X*vec(A) against vec((Z*A)''): %.1e (relative)\n', gapX);
ll = [loglik_mat(A, Sig), loglik_vec(beta, Sig); ...         % matrix form, stacked form
      loglik_mat(Ahat, Sighat), loglik_vec(Ahat(:), Sighat)];
fprintf('log likelihood at the true parameters:     %.6f (matrix form), %.6f (stacked form)\n', ...
    ll(1, :));
fprintf('log likelihood at the least squares values: %.6f (matrix form), %.6f (stacked form)\n', ...
    ll(2, :));
assert(all(dYZ == 0) && gapX < 1e-12 && max(abs(ll(:, 1) - ll(:, 2))) < 1e-8, ...
    'the matrix form and the stacked form disagree');
fprintf('least squares estimates: largest error %.3f, mean absolute error %.3f\n', ...
    max(abs(Ahat - A), [], 'all'), mean(abs(Ahat - A), 'all'));
fprintf('the intercepts and lag-1 coefficients of equation 1 (true, then estimated):\n');
disp([A(1:4,1), Ahat(1:4,1)]');

%% 3. generalized least squares on the system: the code to adapt
W = kron(speye(T), inv(Sig));   % the precision of eps, I_T kron inv(Sig)
XWX = X'*W*X;                   % nk x nk
XWy = X'*W*y;
bgls = XWX\XWy;

%% 3. the checks
fprintf('\ngeneralized least squares on the stacked form\n');
gapK = [rel(XWX, kron(inv(Sig), Z'*Z)), rel(XWy, reshape(Z'*Y/Sig, [], 1))];
fprintf('X''*W*X against kron(inv(Sig), Z''*Z): %.1e; X''*W*y against vec(Z''*Y*inv(Sig)): %.1e\n', ...
    gapK);
SigR = cov(randn(10, n));
Slist = {Sig, eye(n), SigR}; gapG = zeros(1, 3);
for j = 1:3
    WS = kron(speye(T), inv(Slist{j}));
    gapG(j) = rel((X'*WS*X)\(X'*WS*y), Ahat(:));
    fprintf('GLS with %s against least squares by equation: %.1e (relative)\n', ...
        ifelse_name(Slist{j}, Sig), gapG(j));
end
Xd = full(X); Wd = full(W);
gapD = rel(full(XWX), Xd'*Wd*Xd);
fprintf('the shortcut against dense algebra: %.1e\n', gapD);
assert(all(gapK < 1e-12) && all(gapG < 1e-10) && gapD < 1e-12, 'the Kronecker shortcut fails');

%% 4. one equation with fewer regressors: the adaptation
keep = setdiff(1:n*k, [3 4 6 7]);   % equation 1 drops the lags of variables 2, 3
Xr = X(:, keep);
br_gls = (Xr'*W*Xr)\(Xr'*W*y);
br_ols = [Z(:, [1 2 5])\Y(:,1); reshape(Ahat(:, 2:n), [], 1)];   % equation by equation
Wdiag = kron(speye(T), diag(1./diag(Sig)));   % the same variances, no correlation
br_gls_diag = (Xr'*Wdiag*Xr)\(Xr'*Wdiag*y);
fprintf('\nequation 1 without the lags of variables 2 and 3: %d coefficients\n', numel(keep));
gapR = [max(abs(br_gls(1:3) - br_ols(1:3))), max(abs(br_gls(4:end) - br_ols(4:end))), ...
    max(abs(br_gls_diag - br_ols))];
fprintf('GLS against least squares by equation, largest difference: %.1e in equation 1, %.3f in equations 2 and 3\n', ...
    gapR(1:2));
fprintf('the same with a diagonal Sig: %.1e\n', gapR(3));
assert(gapR(1) < 1e-10 && gapR(3) < 1e-10, 'GLS differs from least squares where the two agree');

%% the figure
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 16 6.5]);
subplot(1,2,1); spy(Z(1:4,:), 'k', 8); box off
title('Z, t = 1, ..., 4', 'FontWeight', 'normal'); xlabel('');
subplot(1,2,2); spy(X(1:4*n,:), 'k', 6); box off
title('X, t = 1, ..., 4', 'FontWeight', 'normal'); xlabel('');
exportgraphics(fig, fullfile(here, 'fig_forms.png'), 'Resolution', 150); close(fig);

function s = ifelse_name(S, Sig)
% the label of a covariance matrix in the printed checks
if isequal(S, Sig), s = 'the true Sig'; elseif isequal(S, eye(size(S))), s = 'Sig = I'; else, s = 'a random Sig'; end
end
