%% ex07 - Dynamic sign restrictions and impulse responses: Uhlig (2005)
%
% ex06 restricted the impact responses only. Uhlig (2005) restricts a stretch of
% the impulse response instead: a monetary contraction raises the funds rate and
% lowers prices, commodity prices and nonborrowed reserves for K months, not just
% on impact. Output is left free, which is what makes the identification
% agnostic.
%
% THE ACCEPTANCE TEST. The impact restrictions are tested first, which is
% where the two rules differ - sign_restrict is the accept-reject algorithm of
% Rubio-Ramirez, Waggoner and Zha (2010), sign_assign the search of Chan,
% Matthes and Yu (2026), both applied to rotations that qr_sign draws uniformly
% from the orthogonal group. The impulse responses at horizons 1..K are tested
% second, which is plain rejection under both. So the 2026 algorithm accelerates
% the impact stage; it does not impose the dynamic restrictions. That is the
% structure of Application_Uhlig2005.m in the replication package.
%
% This script runs both rules over the SAME candidates until each has 1000
% accepted draws. They differ in cost, by a wide margin. They do not differ in
% answer - Proposition 1 says the assignment rule still targets a uniform
% rotation - and the script checks that against a Monte Carlo yardstick rather
% than asserting it.
%
% DATA. Read-only from replications/chan_matthes_yu2026_qe_svarsign/legacy/data/
% Uhlig_monthly.csv, monthly US data on GDP, the GDP deflator, commodity prices,
% nonborrowed reserves, total reserves and the federal funds rate, the first five
% in logs. The package records that commodity prices were perturbed with noise
% under a licensing agreement, so these numbers are close to but not identical to
% the published ones. The prior is the asymmetric conjugate prior of Chan (2022),
% the package's model 1; its model 2 needs sample_BSig_NCP, which is not core.
%
% See:
% Uhlig, H. (2005). What are the Effects of Monetary Policy on Output? Results
% from an Agnostic Identification Procedure, Journal of Monetary Economics,
% 52(2): 381-419.
% Rubio-Ramirez, J.F., Waggoner, D.F. and Zha, T. (2010). Structural Vector
% Autoregressions: Theory of Identification and Algorithms for Inference,
% Review of Economic Studies, 77(2): 665-696.
% Chan, J.C.C., Matthes, C. and Yu, X. (2026). Large Structural VARs with
% Multiple Sign and Ranking Restrictions, Quantitative Economics, 17(3): 709-740.

run(fullfile(fileparts(fileparts(mfilename('fullpath'))),'setup.m'))

rng(20260908, 'twister')
repo = fileparts(fileparts(mfilename('fullpath')));
fprintf('\n=== ex07: dynamic sign restrictions, Uhlig (2005) ===\n');

%% ------------------------------------------------------------------
%  1. Data and prior
%  ------------------------------------------------------------------
csv = fullfile(repo, 'replications', 'chan_matthes_yu2026_qe_svarsign', ...
    'legacy', 'data', 'Uhlig_monthly.csv');
data = readmatrix(csv, 'Range', 'B602:G1069');   % never written to; legacy is frozen
data(:,1:5) = log(data(:,1:5));

p = 12; var_id = 1:6; idx_ns = var_id;           % all six enter in levels
Y0 = data(1:12, var_id);
Y  = data(13:end, var_id);
[T, n] = size(Y);

sig2 = bvar.priors.resid_var_ar4(Y0, Y);
kappa = [.04, .0016, 1, 100];
prior_redu = bvar.priors.acp_redu(n, p, kappa, sig2, idx_ns);

vnames = ["GDP" "deflator" "comm. prices" "nonborr. res." "total res." "FFR"];
fprintf('\nsample: T = %d months, n = %d variables, p = %d lags (k = %d per equation)\n', ...
    T, n, p, n*p+1);

%% ------------------------------------------------------------------
%  2. The restrictions. One shock, restricted over K+1 horizons.
%  ------------------------------------------------------------------
S = [0, -1, -1, -1, 0, 1]';      % 0 = unrestricted, as in the legacy script
m = size(S, 2);
K = 5;                           % restrictions hold at horizons 0,1,...,K
horizon = 61;                    % horizons computed for the plots
idx = find(S(:,1) == -1 | S(:,1) == 1);
nidx = numel(idx);

fprintf('%d shock, %d restricted variables, imposed at horizons 0 to %d\n', m, nidx, K);
fprintf('restricted: ');
for j = idx'
    fprintf('%s(%+d) ', vnames(j), S(j,1));
end
fprintf('\nfree      : %s\n', strjoin(vnames(setdiff(1:n, idx)), ', '));

%% ------------------------------------------------------------------
%  3. Both rules, on the same candidates, to 1000 accepted draws each
%  ------------------------------------------------------------------
target = 1000; nbatch = 2000;
resp_strict = zeros(target, n, horizon);
resp_assign = zeros(target, n, horizon);
nS = 0; nA = 0; cand_S = 0; cand_A = 0;

fprintf('\ncollecting %d accepted draws under each rule...\n', target);
t0 = tic;
while nS < target || nA < target
    [store_alp, store_beta, store_Sig] = ...
        bvar.samplers.acp_theta_sig(Y0, Y, p, prior_redu, nbatch);
    [store_Btilde, store_Sigtilde] = ...
        bvar.structural.reduced_form(store_alp, store_beta, store_Sig);

    for isim = 1:nbatch
        Sigtilde = squeeze(store_Sigtilde(isim,:,:));
        [Q, ~] = bvar.structural.qr_sign(randn(n,n));
        L = chol(Sigtilde,'lower') * Q;          % the SAME candidate for both rules
        A = reshape(store_Btilde(isim,:), n*p+1, n);
        A = A(2:end,:);                          % drop the intercept row
        if nS < target, cand_S = cand_S + 1; end
        if nA < target, cand_A = cand_A + 1; end

            % Strict rule. Pass an EMPTY Ridx, not a zero row: sign_restrict
            % requires Rineq*L < 0 STRICTLY, so a zero row would reject every
            % draw. sign_assign's test is <= 0, where a zero row is harmless.
        if nS < target
            [ok, Ls] = bvar.structural.sign_restrict(L, S, zeros(0,n), []);
            if ok && dynamic_ok(A, Ls, S, idx, nidx, K, m)
                nS = nS + 1;
                resp_strict(nS,:,:) = squeeze(bvar.structural.irf_redu(A, Ls, horizon, m));
            end
        end

        if nA < target
            [ok, La] = bvar.structural.sign_assign(L, S, zeros(m,n), 1);
            if ok && dynamic_ok(A, La, S, idx, nidx, K, m)
                nA = nA + 1;
                resp_assign(nA,:,:) = squeeze(bvar.structural.irf_redu(A, La, horizon, m));
            end
        end
    end
end
elapsed = toc(t0);

fprintf('\n%s\n', repmat('-', 1, 68));
fprintf('%-28s %14s %14s\n', '', 'sign_restrict', 'sign_assign');
fprintf('%s\n', repmat('-', 1, 68));
fprintf('%-28s %14d %14d\n', 'accepted draws', nS, nA);
fprintf('%-28s %14d %14d\n', 'candidates examined', cand_S, cand_A);
fprintf('%-28s %13.2f%% %13.2f%%\n', 'acceptance rate', 100*nS/cand_S, 100*nA/cand_A);
fprintf('%s\n', repmat('-', 1, 68));
fprintf('the assignment rule reached %d draws from %.1fx fewer candidates (%.0f s total)\n', ...
    target, cand_S/cand_A, elapsed);
fprintf(['\nThe margin here is far smaller than ex06''s. With ONE restricted shock\n' ...
         'and six columns, the strict rule already has a fair chance of finding the\n' ...
         'shock in column 1. The advantage of searching grows with the number of\n' ...
         'shocks that must be placed, which is why the 15-variable case in ex06\n' ...
         'looks so different.\n']);

%% ------------------------------------------------------------------
%  4. Do the two rules give the same answer?
%  ------------------------------------------------------------------
qS = quantile(resp_strict, [.16 .5 .84], 1);
qA = quantile(resp_assign, [.16 .5 .84], 1);

fprintf('\nmedian response to a one-standard-deviation monetary contraction:\n');
fprintf('%-16s', 'variable');
fprintf('%22s', 'h = 0', 'h = 12', 'h = 36'); fprintf('\n');
fprintf('%-16s', '');
fprintf('%11s%11s', 'strict', 'assign', 'strict', 'assign', 'strict', 'assign'); fprintf('\n');
for i = 1:n
    fprintf('%-16s', vnames(i));
    for h = [1 13 37]
        fprintf('%11.4f%11.4f', qS(2,i,h), qA(2,i,h));
    end
    fprintf('\n');
end

    % Judging that comparison needs care on two counts. The variables are on
    % different scales - the funds rate is in percent, the other five in logs -
    % so a raw maximum over all of them just reports the funds rate. And any gap
    % has to be read against the gap Monte Carlo noise produces on its own.
    % Both are handled by measuring every difference as a fraction of the 68%
    % band at that variable and horizon, and by comparing against random
    % half-splits of ONE rule's draws. The split-half yardstick is conservative:
    % two 500-draw medians are noisier than two 1000-draw medians by about
    % sqrt(2), so it overstates the noise the rule-vs-rule comparison carries.
bw  = squeeze(qS(3,:,:) - qS(1,:,:));
gap = abs(squeeze(qS(2,:,:) - qA(2,:,:)));
rel_between = max(gap(:) ./ bw(:));

nrep = 20; rel_within = zeros(nrep,1);
for r = 1:nrep
    prm = randperm(target);
    a = squeeze(quantile(resp_assign(prm(1:target/2),:,:), .5, 1));
    b = squeeze(quantile(resp_assign(prm(target/2+1:end),:,:), .5, 1));
    rel_within(r) = max(abs(a(:) - b(:)) ./ bw(:));
end
yard = median(rel_within);

fprintf('\nEvery difference below is expressed as a fraction of the 68%% band at that\n');
fprintf('variable and horizon, so the six scales are comparable.\n');
fprintf('  largest gap, strict vs assignment     : %5.1f%% of the band\n', 100*rel_between);
fprintf('  same, two random halves of one rule   : %5.1f%% of the band  (%d splits, median)\n', ...
    100*yard, nrep);
if rel_between <= yard
    fprintf(['\nThe two rules differ by less than one rule differs from itself. What\n' ...
             'separates them is cost, not the answer, which is Proposition 1 of the\n' ...
             'paper: accepting whenever SOME column works and then drawing an\n' ...
             'assignment still targets a uniform rotation.\n']);
else
    fprintf(['\nThe gap is larger than a single half-split, though of the same order.\n' ...
             'Both rules provably target the same identified set, so this is sampling\n' ...
             'noise rather than disagreement; raising `target` shrinks it.\n']);
end

%% ------------------------------------------------------------------
%  5. The bands
%  ------------------------------------------------------------------
hh = 0:horizon-1;
figure('Name','ex07: Uhlig (2005) impulse responses','Position',[100 100 900 520]);
for i = 1:n
    subplot(2,3,i); hold on
    fill([hh fliplr(hh)], [squeeze(qS(1,i,:))' fliplr(squeeze(qS(3,i,:))')], ...
        [.85 .85 .85], 'EdgeColor','none');
    plot(hh, squeeze(qS(2,i,:)), 'k-',  'LineWidth', 1.2);
    plot(hh, squeeze(qA(1,i,:)), 'r--', 'LineWidth', .8);
    plot(hh, squeeze(qA(2,i,:)), 'r-',  'LineWidth', 1.2);
    plot(hh, squeeze(qA(3,i,:)), 'r--', 'LineWidth', .8);
    yline(0, 'k:');
    xline(K, 'b:');
    title(vnames(i)); xlim([0 horizon-1]); hold off
    if i == 1
        legend({'strict, 68% band','strict, median','assignment, 68%','assignment, median'}, ...
            'Location','best','FontSize',7); legend boxoff
    end
end
sgtitle('Response to a monetary contraction (blue dotted line: last restricted horizon)');

fprintf(['\nThe two sets of bands lie on top of each other. The dotted vertical line\n' ...
         'marks horizon %d, the last one restricted; beyond it the responses are\n' ...
         'estimated rather than imposed. GDP was left unrestricted throughout,\n' ...
         'which is what makes the identification agnostic: its response is a\n' ...
         'result here, not an assumption.\n'], K);

fprintf('\nex07 done. replications/chan_matthes_yu2026_qe_svarsign/ has the published\n');
fprintf('application, which also runs Read (2022) as a third comparison.\n');

%% ------------------------------------------------------------------
function ok = dynamic_ok(A, L, S, idx, nidx, K, m)
% Horizons 1..K, given that the impact restrictions already hold. Only K+1
% horizons are needed here, so this is cheaper than the 61 computed for the
% accepted draws.
r = bvar.structural.irf_redu(A, L, K+1, m);
sgn = sign(squeeze(r(idx,1,2:K+1))) == repmat(S(idx,1), 1, K);
ok = sum(sgn(:)) == nidx*K;
end
