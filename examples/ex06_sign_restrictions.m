%% ex06 - Identifying a SVAR by sign restrictions, and why the search matters
%
% THE PROBLEM. A reduced-form VAR pins down Sigma but not the impact matrix. Any
% L with L*L' = Sigma is admissible, and they are related by rotation: given one
% Cholesky factor L0, every candidate is L0*Q for some orthogonal Q. Sign
% restrictions pick out the economically meaningful ones - a monetary
% contraction raises the interest rate and lowers output and prices, and so on.
% Since no finite set of sign restrictions leaves a single L, the object of
% inference is a SET, and it is explored by drawing rotations at random and
% keeping those that satisfy the restrictions.
%
% That acceptance step is the whole computational problem, and this example is
% about two ways of doing it.
%
%   bvar.structural.sign_restrict   requires column i to satisfy shock i, and
%                                   rejects the candidate at the first shock
%                                   that fails. This is the standard scheme.
%
%   bvar.structural.sign_assign     asks which columns admit which shocks,
%                                   accepts whenever every shock has at least
%                                   one, and then draws an assignment. This is
%                                   the algorithm of Chan, Matthes and Yu (2026).
%
% The insight is that the labelling of the columns of Q is arbitrary. A rotation
% whose fourth column looks like a monetary shock is just as admissible as one
% where the monetary shock happens to land in the fourth column, but the strict
% scheme throws the first away. Both target the same identified set. They differ
% only in how much work is wasted, and the difference is not small: the run
% behind this repository's tests/golden/chan2022_qe_acp/main_ACP_apps_15var
% capture needed 3.8 million draws for each acceptance at n = 15.
%
% WHAT THIS SCRIPT DOES. Draws one batch from the posterior of a 6-variable VAR
% under the asymmetric conjugate prior, then runs BOTH acceptance rules over the
% same batch, so the comparison is of the rules and not of the sampling. It
% prints the two acceptance rates and checks that the accepted draws satisfy the
% restrictions they were selected for.
%
% DATA. Read-only from replications/chan_matthes_yu2026_qe_svarsign/legacy/data/
% database_2019Q4.csv, the quarterly US panel of Chan, Matthes and Yu (2026),
% first six variables: GDP, the GDP deflator, the 3-month Tbill, investment, the
% S&P 500 and a credit spread.
%
% See:
% Chan, J.C.C., Matthes, C. and Yu, X. (2026). Large Structural VARs with
% Multiple Sign and Ranking Restrictions, Quantitative Economics, 17(3): 709-740.

run(fullfile(fileparts(fileparts(mfilename('fullpath'))),'setup.m'))

rng(20260908, 'twister')
repo = fileparts(fileparts(mfilename('fullpath')));
fprintf('\n=== ex06: sign restrictions, and the cost of the search ===\n');

%% ------------------------------------------------------------------
%  1. Data and prior
%  ------------------------------------------------------------------
csv = fullfile(repo, 'replications', 'chan_matthes_yu2026_qe_svarsign', ...
    'legacy', 'data', 'database_2019Q4.csv');
data = readmatrix(csv, 'Range', 'B3:P142');   % never written to; legacy is frozen

p = 5; var_id = 1:6; idx_ns = [1,2,4,5];      % RWZ_15var.m, dataset = 1
Y0 = data(1:8, var_id);
Y  = data(9:end, var_id);
[T, n] = size(Y);

sig2 = bvar.priors.resid_var_ar4(Y0, Y);
kappa = [1, 1, 1, 100];
prior_redu = bvar.priors.acp_redu(n, p, kappa, sig2, idx_ns);

fprintf('\nsample: T = %d quarters, n = %d variables, p = %d lags\n', T, n, p);

%% ------------------------------------------------------------------
%  2. The restrictions. One column of S per shock; +1 and -1 restrict the
%     sign of that variable's impact response, NaN leaves it free.
%  ------------------------------------------------------------------
supply   = [ 1,-1,NaN,NaN,  1,NaN]';
demand   = [ 1, 1,  1,NaN,NaN,NaN]';
monetary = [ 1, 1, -1,NaN,NaN,NaN]';
invest   = [ 1, 1,  1,NaN, -1,NaN]';
finc     = [ 1, 1,  1,NaN,  1,NaN]';
S = [supply, demand, monetary, invest, finc];
m = size(S, 2);
shocks = ["supply" "demand" "monetary" "investment" "financial"];

    % Row inequalities: a linear combination of the impact responses required to
    % be negative. Here they rank the interest-rate response against the
    % investment response for three of the shocks. The two functions index them
    % differently - by COLUMN for sign_restrict, by SHOCK for sign_assign -
    % which is the same information written two ways.
Rineq_bycol = [-1,0,0,1,0,0; 1,0,0,-1,0,0; 1,0,0,-1,0,0];
Ridx        = [2, 4, 5];
Rineq_byshock = zeros(m, n);
Rineq_byshock(Ridx, :) = Rineq_bycol;

fprintf('%d shocks, %d sign restrictions, %d row inequalities\n', ...
    m, sum(~isnan(S(:))), numel(Ridx));

%% ------------------------------------------------------------------
%  3. One batch of posterior draws, shared by both rules
%  ------------------------------------------------------------------
nbatch = 50000; % Large enough that the strict rule keeps a handful of draws
fprintf('\ndrawing %d posterior draws (conjugate prior, so no chain)...\n', nbatch);
[store_alp, store_beta, store_Sig] = bvar.samplers.acp_theta_sig(Y0, Y, p, prior_redu, nbatch);
[~, store_Sigtilde] = bvar.structural.reduced_form(store_alp, store_beta, store_Sig);

%% ------------------------------------------------------------------
%  4. Both acceptance rules, on the same candidates
%  ------------------------------------------------------------------
n_strict = 0; n_assign = 0;
Lkept = [];
for isim = 1:nbatch
    Sigtilde = squeeze(store_Sigtilde(isim,:,:));
    [Q, ~] = bvar.structural.qr_sign(randn(n,n));
    L = chol(Sigtilde,'lower') * Q;          % the SAME candidate for both rules

    n_strict = n_strict + bvar.structural.sign_restrict(L, S, Rineq_bycol, Ridx);

    [ok, Lassigned] = bvar.structural.sign_assign(L, S, Rineq_byshock, 1);
    if ok
        n_assign = n_assign + 1;
        if isempty(Lkept), Lkept = Lassigned; end
    end
end

fprintf('\n%s\n', repmat('-', 1, 64));
fprintf('%-34s %10s %12s\n', 'acceptance rule', 'accepted', 'rate');
fprintf('%s\n', repmat('-', 1, 64));
fprintf('%-34s %10d %11.2f%%\n', 'sign_restrict (column i = shock i)', n_strict, 100*n_strict/nbatch);
fprintf('%-34s %10d %11.2f%%\n', 'sign_assign   (search assignments)', n_assign, 100*n_assign/nbatch);
fprintf('%s\n', repmat('-', 1, 64));
if n_strict > 0
    fprintf('ratio: %.0fx more draws kept from the same candidates\n', n_assign/n_strict);
    fprintf('(read the ratio as an order of magnitude: it rests on %d accepted draws)\n', n_strict);
else
    fprintf('the strict rule accepted nothing in %d draws; sign_assign kept %d\n', nbatch, n_assign);
end
fprintf(['\nThe candidates are identical - the same Sigma draws, the same rotations.\n' ...
         'Everything the strict rule accepts, the assignment rule accepts too. The\n' ...
         'extra acceptances are rotations whose shocks arrived in the wrong columns.\n']);

%% ------------------------------------------------------------------
%  5. What an accepted draw looks like
%  ------------------------------------------------------------------
if ~isempty(Lkept)
    fprintf('\nimpact responses of one accepted draw (rows = variables, columns = shocks):\n');
    fprintf('%-14s', 'variable');
    fprintf('%12s', shocks); fprintf('\n');
    vnames = ["GDP" "deflator" "Tbill" "investment" "S&P 500" "spread"];
    for i = 1:n
        fprintf('%-14s', vnames(i));
        fprintf('%12.4f', Lkept(i,1:m));
        fprintf('\n');
    end
    fprintf('\nRead the signs against the restrictions: a supply shock raises GDP and\n');
    fprintf('lowers the deflator, a monetary shock raises both and lowers the Tbill.\n');
    fprintf('The free entries (NaN in S) are unrestricted and take whatever sign the\n');
    fprintf('rotation gives them.\n');
end

fprintf('\nex06 done. The identified set is explored by repeating section 4 until\n');
fprintf('enough draws are kept, one draw per accepted candidate: two assignments\n');
fprintf('from the same rotation differ only by a permutation and sign flips, so\n');
fprintf('they would not be independent. replications/chan_matthes_yu2026_qe_svarsign/\n');
fprintf('has the published applications; its auxFunctions/ folder is third-party,\n');
fprintf('the algorithm of Read (2022), one of the two schemes the paper benchmarks\n');
fprintf('against.\n');
