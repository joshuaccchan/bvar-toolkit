% chan2023_joe_mlvarsv/run_ml - the legacy cp_ml = 1 pipeline of main_varsv.m,
% functionized: run the estimation (run_all in this folder) and then the
% marginal-likelihood computation (the extracted bvar.ml.mlvarsv_* functions) on
% its output, in the legacy order and on one continuous rng stream - what the
% legacy scripts do when their tail dispatches the matching utility/ml_var_*
% routine in the same workspace.
%
%   out = run_ml(model, is_kappafixed, is_kappasym, nsim, burnin, seed, varid)
%   out = run_ml(..., 'M', 10000, 'flag_marg', 2, 'bugcompat', true)
%
%   model .. varid - exactly as run_all (seed seeds rng ONCE, before the
%       estimation; the ML computation continues the same stream)
%   'M'         - importance-sampling draws (default 10000, main_varsv.m 28);
%       every routine rounds it up to a multiple of 50, the batch count behind
%       the reported numerical standard error
%   'flag_marg' - default 2 (main_varsv.m 30). Only VAR-FSV implements 1
%   'bugcompat' - default false; forwarded to bvar.ml.mlvarsv_arsvo_redu, the
%       only affected routine (three defects in the VAR-SVO outlier block: the
%       prior mass split over 32 atoms instead of 31, the linear indexing of
%       o_hat, and the missing -n*sum(log(o)) Jacobian). true reproduces
%       ml_var_arsvo_redu.m bitwise; false runs the corrected computation. For
%       the other models the flag is accepted and ignored - clean bills.
%       See tests/variant_map.md for the audit.
%
% VAR-NCP (model 1) computes its log marginal likelihood inline and
% analytically, so run_all already returns it and run_ml only reports it.
%
% Functionized 2026-09-03 (step 10). Output: the run_all output plus out.lml,
% out.lmlstd, out.ml (the routine's detail struct: store_w, the 50 batch
% values, the fitted importance-sampling parameters) and out.ml_bugcompat. The
% legacy 'Computing the marginal likelihood of ...' and 'log marginal
% likelihood of ...' displays are reproduced; the ML functions print nothing.
%
% See:
% Chan, J.C.C. (2023). Comparing stochastic volatility specifications for large
% Bayesian VARs, Journal of Econometrics, 235(2), 1419-1446.

function out = run_ml(model, is_kappafixed, is_kappasym, nsim, burnin, seed, varid, varargin)
if nargin < 1, model = []; end
if nargin < 2, is_kappafixed = []; end
if nargin < 3, is_kappasym = []; end
if nargin < 4, nsim = []; end
if nargin < 5, burnin = []; end
if nargin < 6, seed = []; end
if nargin < 7, varid = []; end
M = [];                     % filled from preset below (main_varsv.m 28)
flag_marg = [];             % (main_varsv.m 30)
bugcompat = false;
for iv = 1:2:numel(varargin)
    switch lower(varargin{iv})
        case 'm', M = varargin{iv+1};
        case 'flag_marg', flag_marg = varargin{iv+1};
        case 'bugcompat', bugcompat = varargin{iv+1};
        otherwise, error('run_ml:badOption', 'unknown option ''%s''', varargin{iv});
    end
end

    % estimation (seeds the stream; the ML computation continues it)
out = run_all(model, is_kappafixed, is_kappasym, nsim, burnin, seed, varid);
if isempty(M),         M = out.preset.ml.M;                 end
if isempty(flag_marg), flag_marg = out.preset.ml.flag_marg; end

Y = out.Y; X = out.X; Y0 = out.Y0; Hyper = out.Hyper;
kfix = out.is_kappafixed; ksym = out.is_kappasym;
detail = struct();
if out.model_num == 5
    if bugcompat
        fprintf(['note: bugcompat is on - reproducing the published computation, including\n' ...
                 'its three documented defects. The corrected version (the default) differs.\n']);
    else
        fprintf(['note: the published ml_var_arsvo_redu.m has three defects in the density\n' ...
                 'evaluation; this run uses the corrected version, so the value will differ\n' ...
                 'from Table 6 of the paper. Pass ''bugcompat'', true to reproduce the\n' ...
                 'published computation. Audit: core/+bvar/+ml/mlvarsv_arsvo_redu.m header.\n']);
    end
end
if out.model_num ~= 1
    disp(['Computing the marginal likelihood of ' out.model_name '... ']);
end
switch out.model_num
    case 1
        lml = out.lml; lmlstd = [];       % inline and analytic in VAR_NCP.m
    case 2
        [lml,lmlstd,detail] = bvar.ml.mlvarsv_csv(X,Y,Y0,M,Hyper, ...
            out.store_h,out.store_hpara,out.store_kappa,kfix);
    case 3
        [lml,lmlstd,detail] = bvar.ml.mlvarsv_arsv_redu(X,Y,Y0,M,Hyper,flag_marg, ...
            out.store_h,out.store_beta,out.store_hpara,out.store_kappa,kfix,ksym);
    case 4
        [lml,lmlstd,detail] = bvar.ml.mlvarsv_fsv(X,Y,Y0,M,Hyper,flag_marg, ...
            out.store_h,out.store_hpara,out.store_l,out.store_kappa,kfix,ksym);
    case 5
        [lml,lmlstd,detail] = bvar.ml.mlvarsv_arsvo_redu(X,Y,Y0,M,Hyper,flag_marg, ...
            out.store_h,out.store_beta,out.store_hpara,out.store_kappa, ...
            out.store_o,out.store_po,out.o_grid,kfix,ksym,'bugcompat',bugcompat);
end

out.lml = lml;
out.lmlstd = lmlstd;
out.ml = detail;
out.ml_M = M;
out.ml_flag_marg = flag_marg;
out.ml_bugcompat = bugcompat;

fprintf(['log marginal likelihood of ' out.model_name ': %.1f \n'], lml);  % main_varsv.m 140
end
