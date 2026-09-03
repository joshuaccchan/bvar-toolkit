%% ex05 - Marginal likelihoods and model comparison
%
% THE QUANTITY. The marginal likelihood of a model M is
%
%       p(y | M) = INT p(y | theta, M) p(theta | M) d theta,
%
% the prior-predictive density of the data actually observed. Ratios of
% marginal likelihoods are Bayes factors, so this one number is what ranks the
% competing error specifications in Chan (2020, JBES). It is not a by-product
% of the MCMC: computing it takes a second pass over the stored draws. The
% method here is Chib's - evaluate the likelihood and the prior at one
% "starred" parameter point (the posterior means), then divide by the
% posterior ordinate at that same point, which the chain lets you estimate:
%
%       log p(y) = log p(y | theta*) + log p(theta*) - log p(theta* | y).
%
% Each piece is a density evaluated at theta*, so the accounting is exact and
% every term is worth inspecting separately - this example prints them.
%
% WHAT IS IN THE TOOLKIT. core/+bvt/+ml/ holds one function per model of that
% paper (kron_bvar, kron_bvar_t, ..., kron_bvar_csv_t_ma), the shared density
% pieces (lniwpdf, linvgammpdf, llike_ma, llike_csv_ma) and four
% importance-sampling evaluators (intlike_*) for the models whose likelihood
% has no closed form because the volatility path must be integrated out.
% replications/chan2020_jbes_kronecker/run_ml.m runs the estimation and then
% the matching ML computation on ONE continuous random-number stream, exactly
% as the legacy scripts do.
%
% THREE MODELS. We compare the error specification, holding the VAR and the
% prior fixed:
%   1  BVAR      iid Gaussian errors
%   2  BVAR-t    Student-t errors (fat tails, no time variation)
%   3  BVAR-CSV  a common stochastic volatility factor (time variation)
% Model 1's marginal likelihood is available in closed form; models 2 and 3
% need the chain.
%
% Sizes here are deliberately tiny (a few hundred draws). The published run
% uses 30000 - do not read the numbers below as results. Their ORDER is
% already informative, but a real comparison needs the full settings.
%
% See: Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics, 38(1), 68-79.

run(fullfile(fileparts(fileparts(mfilename('fullpath'))),'setup.m'))

repo = fileparts(fileparts(mfilename('fullpath')));
fprintf('\n=== ex05: marginal likelihoods and model comparison ===\n');

nsim = 200; burnin = 50; seed = 20260903;
fprintf('\nsettings: nsim = %d, burnin = %d, seed = %d\n', nsim, burnin, seed);
fprintf('(the published run uses nsim = 30000, burnin = 5000)\n');

%% ------------------------------------------------------------------
%  Run the driver from its own folder, so run_all, preset and the legacy
%  data file resolve unambiguously (three replication packages define a
%  function called run_all). The working directory is restored below
%  whatever happens.
%  ------------------------------------------------------------------
kdir  = fullfile(repo, 'replications', 'chan2020_jbes_kronecker');
oldwd = cd(kdir);
try
    fprintf('\n--- model 1: BVAR (iid Gaussian) ---\n');
    o1 = run_ml(1, nsim, burnin, seed);
    fprintf('\n--- model 2: BVAR-t (Student-t errors) ---\n');
    o2 = run_ml(2, nsim, burnin, seed);
    fprintf('\n--- model 3: BVAR-CSV (common stochastic volatility) ---\n');
    o3 = run_ml(3, nsim, burnin, seed);
catch err
    cd(oldwd);
    rethrow(err);
end
cd(oldwd);

%% ------------------------------------------------------------------
%  1. Where the number comes from: the three pieces of Chib's identity
%  ------------------------------------------------------------------
fprintf('\nlog p(y) = log p(y|theta*) + log p(theta*) - log p(theta*|y), for BVAR-t:\n');
fprintf('  log likelihood at theta*                      : %12.4f\n', o2.ml.llike);
fprintf('  log prior at theta*                           : %12.4f\n', o2.ml.lpri);
fprintf('  log posterior ordinate, (A,Sig) block         : %12.4f\n', o2.ml.lpost(1));
fprintf('  log posterior ordinate, nu block              : %12.4f\n', o2.ml.lpost(2));
fprintf('  ------------------------------------------------------------\n');
fprintf('  log marginal likelihood                       : %12.4f\n', o2.ML);
fprintf('  (posterior mean of the t degrees of freedom nu : %.2f - well below\n', o2.ml.nu_mean);
fprintf('   the ~30 at which a t is indistinguishable from a normal, so the\n');
fprintf('   data are asking for fat tails.)\n');

fprintf('\nNote the two ordinate blocks are SUBTRACTED. A model can raise its\n');
fprintf('likelihood by fitting the sample more closely and still lose, because\n');
fprintf('a sharper posterior means a larger ordinate at theta*. That subtraction\n');
fprintf('is the Ockham factor: the marginal likelihood prices complexity\n');
fprintf('automatically, with no penalty term put in by hand.\n');

%% ------------------------------------------------------------------
%  2. The comparison
%  ------------------------------------------------------------------
names = {'BVAR      (iid Gaussian)', 'BVAR-t    (Student-t)', 'BVAR-CSV  (common SV)'};
MLs   = [o1.ML, o2.ML, o3.ML];
[~, best] = max(MLs);

fprintf('\n%s\n', repmat('-', 1, 66));
fprintf('%-26s %14s %12s %10s\n', 'model', 'log ML', 'log BF vs 1', 'rank');
fprintf('%s\n', repmat('-', 1, 66));
[~, ord] = sort(MLs, 'descend');
rank = zeros(1,3); rank(ord) = 1:3;
for ii = 1:3
    fprintf('%-26s %14.2f %12.2f %10d%s\n', names{ii}, MLs(ii), MLs(ii) - MLs(1), ...
        rank(ii), repmat('  <- best', 1, ii == best));
end
fprintf('%s\n', repmat('-', 1, 66));

fprintf('\nHow to read the middle column: it is a log Bayes factor against the\n');
fprintf('plain Gaussian VAR. On Kass and Raftery''s scale a log BF above 5 is\n');
fprintf('"very strong" evidence, and these are far larger - allowing for fat\n');
fprintf('tails or time-varying volatility buys an enormous amount, which is\n');
fprintf('the paper''s point. The full-length ranking puts the model that does\n');
fprintf('BOTH (plus an MA term) on top; see tests/golden/ for those tables.\n');

%% ------------------------------------------------------------------
%  3. Using the ml functions on your own model
%  ------------------------------------------------------------------
fprintf('\nTo compute an ML for a different model, the pattern is:\n');
fprintf('  out = run_ml(model, nsim, burnin, seed)   %% estimation + ML, one stream\n');
fprintf('  out.ML          the log marginal likelihood\n');
fprintf('  out.ml.llike / .lpri / .lpost    the three pieces above\n');
fprintf('  out.est         the draws, if you want them\n');
fprintf('The bvt.ml.* functions can also be called directly on stored draws;\n');
fprintf('their headers document the pri/est structs they expect.\n');

fprintf('\nOne footnote for replication: two of the eight legacy ML scripts\n');
fprintf('(models 4 and 8) evaluate one term at a leftover chain draw rather\n');
fprintf('than the posterior mean. The core functions use the posterior mean by\n');
fprintf('default and reproduce the published computation under\n');
fprintf('''bugcompat'', true. It does not affect the ranking; the audit and the\n');
fprintf('comparison are in tests/variant_map.md.\n');

fprintf('\nex05 done. For a published number use the full settings:\n');
fprintf('  cd replications/chan2020_jbes_kronecker\n');
fprintf('  out = run_ml(3, 30000, 5000);\n');
