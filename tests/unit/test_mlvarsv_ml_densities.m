function test_mlvarsv_ml_densities
% bitwise checks of the ml_varsv marginal-likelihood helper functions against
% the legacy copies (run from tempdir), plus two executable statements of the
% VAR-SVO outlier-block defects that bvar.ml.mlvarsv_arsvo_redu corrects.
%
% Covered: bvar.ml.lgampdf, bvar.ml.ltnormpdf, bvar.ml.lmvnpdf_pcn (the three
% density utilities the ML phase actually calls - only from ml_var_csv.m) and
% bvar.ml.isden_arss (getISden_ARSS, called by all four routines). Also asserts
% the cross-package identity bvar.ml.linvgammpdf == legacy ligampdf, which is why
% ligampdf is not re-extracted. utility/lmvnpdf.m has no caller anywhere in the
% package and is deliberately left in legacy/.
root = getappdata(0, 'bvar_repo_root');
leg = fullfile(root, 'replications', 'chan2023_joe_mlvarsv', 'legacy', 'utility');

tmp = tempname; mkdir(tmp);
ctmp = onCleanup(@() cleanup_tmp(tmp));
copies = {'lgampdf.m', 'ligampdf.m', 'ltnormpdf.m', 'lmvnpdf_pcn.m', 'getISden_ARSS.m'};
for k = 1:numel(copies)
    copyfile(fullfile(leg, copies{k}), fullfile(tmp, copies{k}));
end
addpath(tmp);

rng(20260903, 'twister');

% --- lgampdf: shape a, RATE b ---
x = gamrnd(2, 3, 7, 1); a = 1.7; b = 1/.2^2;
assert(isequal(lgampdf(x, a, b), bvar.ml.lgampdf(x, a, b)), 'lgampdf differs');
assert(isequal(lgampdf(x, [1;2;3;4;5;6;7], b), bvar.ml.lgampdf(x, [1;2;3;4;5;6;7], b)), ...
    'lgampdf differs on a vector shape');

% --- ligampdf == bvar.ml.linvgammpdf (the step-8 Kronecker extraction) ---
assert(isequal(ligampdf(x, a, b), bvar.ml.linvgammpdf(x, a, b)), ...
    'legacy ligampdf and bvar.ml.linvgammpdf differ');

% --- ltnormpdf on (-1,1), the phi prior/IS ordinate of ml_var_csv ---
xr = -.9 + 1.8*rand(5,1);
assert(isequal(ltnormpdf(xr, .98, .05^2, -1, 1), bvar.ml.ltnormpdf(xr, .98, .05^2, -1, 1)), ...
    'ltnormpdf differs');

% --- lmvnpdf_pcn: precision parameterization ---
T = 12;
Hr = speye(T) - .9*sparse(2:T, 1:(T-1), ones(1,T-1), T, T);
K = Hr'*sparse(1:T, 1:T, 1 + rand(T,1))*Hr;
mu = randn(T,1); xh = mu + randn(T,1);
assert(isequal(lmvnpdf_pcn(xh, mu, K), bvar.ml.lmvnpdf_pcn(xh, mu, K)), 'lmvnpdf_pcn differs');

% --- isden_arss (getISden_ARSS): all five outputs, on draws with real serial
%     dependence so the fminbnd rho search has something to find ---
R = 60; Tb = 40;
store_h = zeros(R, Tb);
for ii = 1:R
    e = .3*randn(Tb,1); hh = zeros(Tb,1); hh(1) = e(1)/sqrt(1-.95^2);
    for t = 2:Tb, hh(t) = .95*hh(t-1) + e(t); end
    store_h(ii,:) = hh' - 1.2;
end
[h1,K1,r1,m1,v1] = getISden_ARSS(store_h);
[h2,K2,r2,m2,v2] = bvar.ml.isden_arss(store_h);
assert(isequal(h1,h2) && isequal(K1,K2) && isequal(r1,r2) && isequal(m1,m2) && isequal(v1,v2), ...
    'isden_arss differs from getISden_ARSS');
assert(abs(r1) < .99 && all(v1 > 0), 'isden_arss returned a degenerate fit');

% --- defect 1: the ML spreads the outlier prior mass over numel(o_grid) = 32
%     atoms where the sampler uses 31, so the ML prior does not sum to one over
%     the atoms the sampler can produce, and its o_lpri has an unread 33rd entry
ngrid_sampler = 31;                                  % VAR_ARSVO_redu.m line 8
o_grid = [1; linspace(2,20,ngrid_sampler)'];         % line 9
ngrid_ml = size(o_grid,1);                           % ml_var_arsvo_redu.m line 12
po = 1/16;
lpri_sampler = log([1-po; repmat(po/ngrid_sampler,ngrid_sampler,1)]);   % VAR_ARSVO_redu.m 112
lpri_ml      = log([1-po; repmat(po/ngrid_ml,ngrid_ml,1)]);             % ml_var_arsvo_redu.m 180
assert(numel(lpri_sampler) == 32 && numel(lpri_ml) == 33, 'o_lpri lengths changed');
assert(abs(sum(exp(lpri_sampler)) - 1) < 1e-12, 'the sampler o prior must sum to 1');
assert(abs(sum(exp(lpri_ml(1:32))) - (1 - po/32)) < 1e-12, ...
    'the ML o prior must leak po/32 onto the unreachable 33rd atom');
assert(abs((lpri_ml(2) - lpri_sampler(2)) - log(31/32)) < 1e-14, ...
    'the per-outlier-atom prior gap must be log(31/32)');

% --- defect 2: o_hat(o_idx) linear-indexes a T x 32 matrix, so for T >= 32
%     every index falls in column 1 ---
Tp = 234;                                            % the published sample length
o_hat = rand(Tp, ngrid_ml); o_hat = o_hat./sum(o_hat,2);
o_idx = randi(ngrid_ml, Tp, 1);
assert(isequal(o_hat(o_idx), o_hat(o_idx,1)), ...
    'the legacy linear index is expected to read column 1 of o_hat');
correct = o_hat(sub2ind(size(o_hat), (1:Tp)', o_idx));
assert(~isequal(o_hat(o_idx), correct), 'the corrected index must differ from the legacy one');
end

% -------------------------------------------------------------------------
function cleanup_tmp(tmp)
if any(strcmpi(strsplit(path, pathsep), tmp))
    rmpath(tmp);
end
if exist(tmp, 'dir')
    rmdir(tmp, 's');
end
end
