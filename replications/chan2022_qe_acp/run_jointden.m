% chan2022_qe_acp/run_jointden - main_ACP_jointden.m, functionized: evaluates
% the marginal likelihood of the 15-variable VAR over a grid of the two
% shrinkage hyperparameters, which is the joint posterior density of
% kappa1_tilde and kappa2_tilde up to a constant.
%
%   out = run_jointden()
%   out = run_jointden(kappa1_grid, kappa2_grid)
%
% No simulation is involved. bvar.ml.acp gives
% the marginal likelihood in closed form, so the surface is evaluated directly,
% once per grid point, and the optimum found by bvar.priors.acp_opt_kappa can be
% placed on it. The default grid is the paper's: kappa1 over 0.01:0.001:0.2 and
% kappa2 over 0.001:0.0002:0.012, which is 191 x 56 = 10696 evaluations.
%
% Output: store_lml (the log marginal likelihood at each grid point), store_ml
% (exponentiated and normalized by its maximum, as the legacy plots it), the
% meshgrid arrays Kappa1 and Kappa2, and kappa_Sym / ml_Sym, the symmetric-prior
% optimum marked on the paper's figure. The subjective prior's point,
% [.04, .0016], is in preset as pr.jd.subjective.
%
% NO PLOTTING, as in run_all: the contour, the axis limits and the two
% annotations are left to the caller.
%
% Core used: bvar.priors.resid_var_ar4 (get_resid_var), bvar.priors.acp_redu
% (prior_ACP_redu), bvar.ml.acp (ml_VAR_ACP), and bvar.priors.acp_opt_kappa with
% 'symmetric', true (get_OptSymKappa).
%
% Functionized 2026-09-07 (step 12). The pieces are covered by
% tests/unit/test_acp_equivalence.m; this driver is a loop over them.
%
% See:
% Chan, J.C.C. (2022). Asymmetric Conjugate Priors for Large Bayesian VARs,
% Quantitative Economics, 13(3): 1145-1169.

function out = run_jointden(kappa1_grid, kappa2_grid)
thisdir = fileparts(mfilename('fullpath'));
if isempty(which('bvar.ml.acp'))
    root = fileparts(fileparts(thisdir));
    addpath(fullfile(root, 'core'));
end
od = cd(thisdir);
guard = onCleanup(@() cd(od));
pr = preset();
clear guard

if nargin < 1 || isempty(kappa1_grid), kappa1_grid = pr.jd.kappa1_grid; end
if nargin < 2 || isempty(kappa2_grid), kappa2_grid = pr.jd.kappa2_grid; end

p = pr.p;
idx_ns = pr.jd.idx_ns;
data = xlsread(fullfile(thisdir, 'legacy', pr.data_file)); %#ok<XLSRD>
Y0 = data(1:pr.n_init, pr.jd.var_id);
Y  = data(pr.n_init+1:end, pr.jd.var_id);
[T,n] = size(Y);
tmpY = [Y0(end-p+1:end,:); Y];
Z = zeros(T,n*p);
for ii = 1:p
    Z(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
end
Z = [ones(T,1) Z];
sig2 = bvar.priors.resid_var_ar4(Y0, Y);
kappa = pr.jd.kappa;

[Kappa1,Kappa2] = meshgrid(kappa1_grid, kappa2_grid);
store_lml = zeros(size(Kappa1,1), size(Kappa2,2));
fprintf('Evaluating the marginal likelihood over a %d x %d grid...\n', ...
    size(store_lml,1), size(store_lml,2));
for ii = 1:size(Kappa1,1)
    for ij = 1:size(Kappa2,2)
        store_lml(ii,ij) = bvar.ml.acp(p,Y,Z, ...
            bvar.priors.acp_redu(n,p,[Kappa1(ii,ij),Kappa2(ii,ij),kappa(3),kappa(4)],sig2,idx_ns));
    end
end
store_ml = exp(store_lml-max(max(store_lml)));
[ml_Sym,kappa_Sym] = bvar.priors.acp_opt_kappa(Y0,Y,Z,p,[],'redu',idx_ns,'symmetric',true);

out = struct('Kappa1',Kappa1, 'Kappa2',Kappa2, 'store_lml',store_lml, ...
    'store_ml',store_ml, 'ml_Sym',ml_Sym, 'kappa_Sym',kappa_Sym, ...
    'Y',Y, 'Y0',Y0, 'Z',Z, 'sig2',sig2, 'p',p, 'n',n, 'T',T, ...
    'idx_ns',idx_ns, 'kappa',kappa, 'preset',pr);
end
