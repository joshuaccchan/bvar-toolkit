% chan2020_jbes_kronecker/run_ml - the legacy cp_ml = 1 pipeline of
% main_BVAR.m, functionized: run the estimation (run_all in this folder) and
% then the marginal-likelihood computation (the extracted bvt.ml.kron_bvar_*
% functions) on its output, in the legacy order and on one continuous rng
% stream - exactly what the legacy scripts do when the estimation script's
% tail dispatches its ml_BVAR_* companion in the same workspace.
%
%   out = run_ml(model, nsim, burnin, seed)
%   out = run_ml(model, nsim, burnin, seed, 'bugcompat', true)
%
%   model/nsim/burnin/seed - exactly as run_all (seed seeds rng ONCE, before
%       estimation; the ML computation continues the same stream).
%   'bugcompat' (default false) - forwarded to the two affected ML functions:
%       model 4 (bvt.ml.kron_bvar_ma: the ml_BVAR_MA.m line-17 leftover-psi
%       llike term) and model 8 (bvt.ml.kron_bvar_csv_t_ma: the frozen
%       leftover Hpsi/psi ordinate loop and the line-108 leftover-Sig psi
%       target). true reproduces each legacy ml script bitwise from the same
%       seed and chain, consuming run_all's out.state exactly as the legacy
%       script consumes its workspace leftovers; false (default) runs the
%       corrected computation (consistent evaluation point across all
%       ordinate pieces). For every other model the flag is accepted and
%       ignored - their ml scripts have no divergent modes (clean bills,
%       step-8 audit; see tests/variant_map.md).
%
% Functionized 2026-09-02 (step 8, Kronecker family pass, part 1). Output:
% the run_all output plus out.ML, out.ml (the ML function's detail struct:
% llike, lpri, lpost, ordinate stores), and out.ml_bugcompat. The legacy
% fprintf('log marginal likelihood of ...') line is reproduced; the ML
% functions themselves print nothing.
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics,
% 38(1), 68-79.

function out = run_ml(model, nsim, burnin, seed, varargin)
if nargin < 1, model = []; end
if nargin < 2, nsim = []; end
if nargin < 3, burnin = []; end
if nargin < 4, seed = []; end
bugcompat = false;
for iv = 1:2:numel(varargin)
    switch lower(varargin{iv})
        case 'bugcompat', bugcompat = varargin{iv+1};
        otherwise, error('run_ml:badOption', 'unknown option ''%s''', varargin{iv});
    end
end

    % estimation (seeds the stream; the ML computation continues it)
out = run_all(model, nsim, burnin, seed);

switch out.model
    case 1
        [ML, detail] = bvt.ml.kron_bvar(out.shortY, out.X, out.pri, out);
        name = 'BVAR';
    case 2
        disp('Computing the marginal likelihood of BVAR-t... ');
        [ML, detail] = bvt.ml.kron_bvar_t(out.shortY, out.X, out.pri, out);
        name = 'BVAR-t';
    case 3
        disp('Computing the marginal likelihood of BVAR-CSV... ');
        [ML, detail] = bvt.ml.kron_bvar_csv(out.shortY, out.X, out.pri, out);
        name = 'BVAR-CSV';
    case 4
        disp('Computing the marginal likelihood of BVAR-MA... ');
        [ML, detail] = bvt.ml.kron_bvar_ma(out.shortY, out.X, out.pri, out, ...
            'bugcompat', bugcompat);
        name = 'BVAR-MA';
    case 5
        disp('Computing the marginal likelihood of BVAR-t-CSV... ');
        [ML, detail] = bvt.ml.kron_bvar_t_csv(out.shortY, out.X, out.pri, out);
        name = 'BVAR-t-CSV';
    case 6
        disp('Computing the marginal likelihood of BVAR-t-MA... ');
        [ML, detail] = bvt.ml.kron_bvar_t_ma(out.shortY, out.X, out.pri, out);
        name = 'BVAR-t-MA';
    case 7
        disp('Computing the marginal likelihood of BVAR-CSV-MA... ');
        [ML, detail] = bvt.ml.kron_bvar_csv_ma(out.shortY, out.X, out.pri, out);
        name = 'BVAR-CSV-MA';
    case 8
        disp('Computing the marginal likelihood of BVAR-CSV-t-MA... ');
        [ML, detail] = bvt.ml.kron_bvar_csv_t_ma(out.shortY, out.X, out.pri, out, ...
            'bugcompat', bugcompat);
        name = 'BVAR-CSV-t-MA';
end

fprintf('log marginal likelihood of %s: %.1f \n', name, ML);
out.ML = ML;
out.ml = detail;
out.ml_bugcompat = bugcompat;
end
