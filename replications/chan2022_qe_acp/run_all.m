% chan2022_qe_acp/run_all - main_ACP_apps.m, functionized: estimates the VAR of
% Chan (2022) under the asymmetric conjugate prior and computes impulse
% responses identified by sign restrictions and row inequalities.
%
%   out = run_all()                     % dataset 1 (n = 6), the legacy default
%   out = run_all(dataset, nsim, seed)
%   out = run_all(dataset, nsim, seed, 'nbatch', 50000, 'horizon', 36)
%
%   dataset - 1 for the 6-variable VAR (kappa fixed at [1 1 1 100]), 2 for the
%             15-variable one (kappa chosen by maximizing the closed-form
%             marginal likelihood). Default 1.
%   nsim    - accepted draws to collect; default 5000 for dataset 1 and 1000 for
%             dataset 2, as in the legacy file. THIS IS A REJECTION SAMPLER: the
%             cost is nsim divided by the acceptance rate, and the paper's
%             15-variable run takes days. Start small.
%   seed    - seeds the stream once, before the first batch. The legacy script
%             does not seed at all, so a run of it is irreproducible; passing a
%             seed here is the only way to get a repeatable answer, and the
%             equivalence test relies on it.
%
% Output: store_response (nsim x n x m x horizon), the accepted draws' impulse
% responses; response_median and response_CI (the 16th and 84th percentiles) as
% the legacy computes them; count_total, the number of posterior draws examined,
% and the implied acceptance rate; plus the data, prior and settings used.
%
% NO PLOTTING. The legacy script ends by drawing figures. run_all returns the
% arrays those figures are made from and leaves the drawing to the caller, so
% that it can be used from a script, a test, or a headless session. preset
% carries the axis limits and titles if you want to reproduce the figures.
%
% Core used: bvar.priors.resid_var_ar4 (get_resid_var), bvar.priors.acp_redu and
% acp_stru (prior_ACP_redu / prior_ACP_stru), bvar.priors.acp_opt_kappa
% (get_OptKappa), bvar.samplers.acp_theta_sig (sample_ThetaSig),
% bvar.structural.reduced_form (getReducedForm), bvar.structural.qr_sign (QR),
% bvar.structural.sign_restrict (the inline check at main_ACP_apps.m 103-131)
% and bvar.structural.irf_redu (IRredu).
%
% Functionized 2026-09-07 (step 12). Draw-for-draw equivalence against the
% unmodified legacy script: tests/unit/test_acp_equivalence.m.
%
% See:
% Chan, J.C.C. (2022). Asymmetric Conjugate Priors for Large Bayesian VARs,
% Quantitative Economics, 13(3): 1145-1169.

function out = run_all(dataset, nsim, seed, varargin)
thisdir = fileparts(mfilename('fullpath'));

if isempty(which('bvar.samplers.acp_theta_sig'))
    root = fileparts(fileparts(thisdir));
    addpath(fullfile(root, 'core'));
end

od = cd(thisdir);
guard = onCleanup(@() cd(od));
pr = preset();
clear guard

if nargin < 1 || isempty(dataset), dataset = 1; end
if nargin < 3, seed = []; end
assert(any(dataset == [1 2]), 'run_all:badDataset', 'dataset must be 1 or 2');
cfg = pr.(sprintf('d%d', dataset));
if nargin < 2 || isempty(nsim), nsim = cfg.nsim; end

nbatch = pr.nbatch;
horizon = pr.horizon;
for iv = 1:2:numel(varargin)
    switch lower(varargin{iv})
        case 'nbatch',  nbatch = varargin{iv+1};
        case 'horizon', horizon = varargin{iv+1};
        otherwise, error('run_all:badOption', 'unknown option ''%s''', varargin{iv});
    end
end
if ~isempty(seed)
    rng(seed, 'twister');
end

p = pr.p;

    % ---- data (main_ACP_apps.m 28-48) ----
    % xlsread, not readmatrix: on this file readmatrix returns [141 16] where
    % xlsread returns [140 15], keeping a header row and an index column that
    % the legacy numeric read drops.
data = xlsread(fullfile(thisdir, 'legacy', pr.data_file)); %#ok<XLSRD>
var_id = cfg.var_id;
idx_ns = cfg.idx_ns;
Y0 = data(1:pr.n_init, var_id);
Y  = data(pr.n_init+1:end, var_id);
[T,n] = size(Y);
tmpY = [Y0(end-p+1:end,:); Y];
Z = zeros(T, n*p);
for ii = 1:p
    Z(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
end
Z = [ones(T,1) Z];

    % ---- prior (49-57) ----
if dataset == 1
    kappa = cfg.kappa;
    ml_opt = [];
else
        % find the optimal kappa values
    [ml_opt, kappa] = bvar.priors.acp_opt_kappa(Y0, Y, Z, p, cfg.kappa_init, 'redu', idx_ns);
end
sig2 = bvar.priors.resid_var_ar4(Y0, Y);
prior_stru = bvar.priors.acp_stru(n, p, kappa, sig2, idx_ns);   %#ok<NASGU> % legacy 56: computed, never used
prior_redu = bvar.priors.acp_redu(n, p, kappa, sig2, idx_ns);

    % ---- identification (59-85) ----
S = cfg.S;
m = size(S,2);
Rineq = cfg.Rineq;
Ridx = cfg.Ridx;

    % ---- rejection sampling (87-144) ----
start_time = clock; %#ok<CLOCK>
count_sat = 0;
count_total = 0;
fprintf('Computing impulse responses from a %d-variable VAR\n', n);
disp('    to an one-standard-deviation financial shock...');
store_response = zeros(nsim,n,m,horizon);
while count_sat < nsim
        %  sample nbatch draws from the posterior
    [store_alp,store_beta,store_Sig] = bvar.samplers.acp_theta_sig(Y0,Y,p,prior_redu,nbatch);
    count_total = count_total + nbatch;

        % obtain the reduced-form parameters
    [store_Btilde,store_Sigtilde] = bvar.structural.reduced_form(store_alp,store_beta,store_Sig);

    for isim = 1:nbatch
        Sigtilde = squeeze(store_Sigtilde(isim,:,:));
        [Q,R] = bvar.structural.qr_sign(randn(n,n)); %#ok<ASGLU>
        L0 = chol(Sigtilde,'lower');
        L = L0*Q;
        [ok,L] = bvar.structural.sign_restrict(L,S,Rineq,Ridx);
            % No count_sat < nsim guard, deliberately: the legacy loop finishes
            % the batch it is in, so the final batch overshoots and
            % store_response ends with MORE than nsim rows. The median and the
            % percentiles below are taken over all of them, so truncating here
            % would change the reported responses.
        if ok
            count_sat = count_sat + 1;
            Btilde = reshape(store_Btilde(isim,:),n*p+1,n);
            store_response(count_sat,:,:,:) = ...
                bvar.structural.irf_redu(Btilde(2:end,:),L,horizon,m);
        end
    end

    if (mod(count_total, pr.progress_every) == 0)
        fprintf('Out of %d million posterior draws, %d satisfy all the restrictions.\n', ...
            count_total/1e6, count_sat);
    end
end
fprintf('%d posterior draws that satisfy all the restrictions are obtained.\n', nsim);
fprintf('The simulation took %g minutes.\n', etime(clock,start_time)/60); %#ok<CLOCK,DETIM>

response_median = squeeze(median(store_response));
response_CI = squeeze(quantile(store_response,[.16,.84]));

out = struct();
out.dataset = dataset; out.nsim = nsim; out.seed = seed;
out.nbatch = nbatch; out.horizon = horizon;
out.Y = Y; out.Y0 = Y0; out.Z = Z; out.var_id = var_id; out.idx_ns = idx_ns;
out.T = T; out.n = n; out.p = p; out.m = m;
out.kappa = kappa; out.ml_opt = ml_opt; out.sig2 = sig2;
out.prior_redu = prior_redu;
out.S = S; out.Rineq = Rineq; out.Ridx = Ridx;
out.store_response = store_response;
out.response_median = response_median;
out.response_CI = response_CI;
out.count_total = count_total;
out.acceptance_rate = nsim / count_total;
out.preset = pr;
end
