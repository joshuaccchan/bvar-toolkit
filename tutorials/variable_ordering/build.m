%% build - regenerate the figures and numbers of tutorials/variable_ordering/README.md
%
% Part 1 estimates both models on four FRED-MD series, in the published and the
% reverse order and with a second seed, at the published chain length of 30,000
% draws after 5,000 burn-in. Part 2 reads the 20-variable posterior means
% that the Chan, Koop and Yu (2024) package ships in its results_mat folder. Part 3
% runs the package's Table3_forecasting.m on a temporary copy of the package, which
% recomputes the forecast comparison from the shipped forecasts, and checks the
% result against the capture in tests/golden. Everything printed goes to a log in
% tempdir and the figures are written next to this file.
%
% Usage, from anywhere:  run tutorials/variable_ordering/build.m

tut_dir = fileparts(mfilename('fullpath'));
tut_repo = fileparts(fileparts(tut_dir));
tut_log = fullfile(tempdir, 'bvar_variable_ordering_build_log.txt');
if exist(tut_log, 'file'), delete(tut_log); end
diary(tut_log);
fprintf('tutorials/variable_ordering/build.m, %s, MATLAB %s\n', ...
    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), version);
tut_t0 = tic;

tut_style = {{'-',  [0 0.447 0.741], 1.0}, ...    % Cholesky SV, published order
             {'--', [0.85 0.325 0.098], 1.0}, ... % Cholesky SV, reversed order
             {'-',  [0.62 0.62 0.62], 2.6}, ...   % order-invariant SV, published order
             {':',  [0 0 0], 1.4}};               % order-invariant SV, reversed order
tut_legend = {'Cholesky SV, published order', 'Cholesky SV, reversed order', ...
              'order-invariant SV, published order', 'order-invariant SV, reversed order'};

%% ------------------------------------------------------------------
%  Part 1. Four variables: both models, both orders and a second seed
%  ------------------------------------------------------------------
tut_cky = fullfile(tut_repo, 'replications', 'chan_koop_yu2024_jbes_oisv');
tut_od = cd(tut_cky); tut_guard = onCleanup(@() cd(tut_od));    % preset.m resolves in its own folder
tut_pr = preset();
clear tut_guard
tut_data = load(fullfile(tut_cky, 'legacy', tut_pr.data_file));  % never written to; legacy is frozen
tut4.vnames = ["IP" "unemployment" "PCE inflation" "fed funds"];  % columns tut_pr.var_id1 = [4 6 12 13]
tut4.cfg = struct('model', {'OI','OI','OI','CS','CS','CS'}, ...
                  'rev',   {false, true, false, false, true, false}, ...
                  'seed',  {20260915, 20260915, 31415926, 20260915, 20260915, 31415926});
tut_chain = struct('nsim', 30000, 'burnin', 5000);                % the package's full-sample chains
fprintf('\n=== Part 1: four variables, %d draws after %d burn-in ===\n', tut_chain.nsim, tut_chain.burnin);
tut4.res = cell(1, numel(tut4.cfg));
for q = 1:numel(tut4.cfg)
    tut_ord = tut_pr.var_id1;
    if tut4.cfg(q).rev, tut_ord = tut_ord(end:-1:1); end
    tut_t1 = tic;
    tut4.res{q} = bvar.models.var_sv(tut_data(1:tut_pr.n0, tut_ord), tut_data(tut_pr.n0+1:end, tut_ord), ...
        tut_pr.p, 'model', tut4.cfg(q).model, 'nsim', tut_chain.nsim, ...
        'burnin', tut_chain.burnin, 'seed', tut4.cfg(q).seed);
    if tut4.cfg(q).rev                            % back to the published order for comparison
        tut4.res{q}.Sig_mean = tut4.res{q}.Sig_mean(:, end:-1:1, end:-1:1);
    end
    fprintf('%s, order %d, seed %d: %.0f s\n', tut4.cfg(q).model, tut4.cfg(q).rev + 1, tut4.cfg(q).seed, toc(tut_t1));
end

relgap = @(a,b) 100*mean(abs(a-b)) / mean(abs(a)+abs(b))*2;        % change as a share of the level
cpath = @(S,i,j) S(:,i,j)./sqrt(S(:,i,i).*S(:,j,j));
tut_pairs = [4 5; 4 6; 1 2; 1 3];                 % CS orderings, CS seeds, OI orderings, OI seeds
prs = nchoosek(1:numel(tut4.vnames), 2);
gap = zeros(size(prs,1), 4);
for q = 1:size(prs,1)
    for c = 1:4
        gap(q,c) = mean(abs(cpath(tut4.res{tut_pairs(c,1)}.Sig_mean, prs(q,1), prs(q,2)) ...
            - cpath(tut4.res{tut_pairs(c,2)}.Sig_mean, prs(q,1), prs(q,2))));
    end
end
[~, tut_ip] = max(gap(:,1));                      % the pair the figure shows: largest Cholesky change
tut4.i = prs(tut_ip,1);  tut4.j = prs(tut_ip,2);
tut4.T = size(tut4.res{1}.Sig_mean, 1);
tut4.dates = datetime(1961,3,1) + calmonths(0:tut4.T-1);
tut4.order = [4 5 1 2];                           % CS 1, CS 2, OI 1, OI 2 in tut4.cfg
fprintf('\n');
fprintf('| | Cholesky SV, two orderings | Cholesky SV, two seeds | order-invariant SV, two orderings | order-invariant SV, two seeds |\n');
fprintf('|---|---|---|---|---|\n');
for q = 1:4
    fprintf('| var(%s) |', tut4.vnames(q));
    for c = 1:4
        fprintf(' %.1f%% |', relgap(tut4.res{tut_pairs(c,1)}.Sig_mean(:,q,q), tut4.res{tut_pairs(c,2)}.Sig_mean(:,q,q)));
    end
    fprintf('\n');
end
for q = 1:size(prs,1)
    fprintf('| corr(%s, %s) |', tut4.vnames(prs(q,1)), tut4.vnames(prs(q,2)));
    fprintf(' %.3f |', gap(q,:));
    fprintf('\n');
end

tut_v = zeros(tut4.T, 4);  tut_c = zeros(tut4.T, 4);
for q = 1:4
    S = tut4.res{tut4.order(q)}.Sig_mean;
    tut_v(:,q) = S(:,4,4);
    tut_c(:,q) = cpath(S, tut4.i, tut4.j);
end
fprintf('range over time of the %s x %s correlation:\n', tut4.vnames(tut4.i), tut4.vnames(tut4.j));
for q = 1:4, fprintf('  %-36s %.2f to %.2f\n', tut_legend{q}, min(tut_c(:,q)), max(tut_c(:,q))); end
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 18 13]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; draw_four(tut4.dates, tut_v, tut_style, tut_legend);
title('Variance of the federal funds rate equation');
nexttile; draw_four(tut4.dates, tut_c, tut_style, {});
title(sprintf('Correlation between the %s and %s equations', tut4.vnames(tut4.i), tut4.vnames(tut4.j)));
exportgraphics(fig, fullfile(tut_dir, 'fig_4var.png'), 'Resolution', 150);
close(fig)

%% ------------------------------------------------------------------
%  Part 2. The 20-variable VAR: the posterior means the package ships
%  ------------------------------------------------------------------
tut_pkg = fullfile(tut_repo, 'replications', 'chan_koop_yu2024_jbes_oisv', 'legacy');
tut20.names = ["real personal income" "real PCE" "real M&T sales" "IP" "capacity utilization" ...
    "unemployment" "payrolls" "hours" "hourly earnings" "PPI finished goods" ...
    "PPI commodities" "PCE inflation" "fed funds" "housing starts" "S&P 500" "USD/GBP" ...
    "1-year spread" "10-year spread" "BAA spread" "ISM new orders"];   % the paper's data appendix, in file order
tut20.var_id = [4,6,12,13, 1:3,5,7:11,14:20];                       % func_main_SVAR_v2.m lines 15-17
tut20.lab = tut20.names(tut20.var_id);                               % labels in the published order
tut20.models = ["CS1" "CS2" "OI1" "OI2"];
for q = 1:4
    L = load(fullfile(tut_pkg, 'results_mat', "fullsample" + tut20.models(q) + ".mat"), 'Sig_mean');
    if endsWith(tut20.models(q), "2")                               % reversed order: back to the published one
        L.Sig_mean = L.Sig_mean(:, end:-1:1, end:-1:1);
    end
    tut20.S{q} = L.Sig_mean;
end
[tut20.T, tut20.n, ~] = size(tut20.S{1});
tut20.dates = datetime(1961,3,1) + calmonths(0:tut20.T-1);

tut20.vg = zeros(tut20.n, 2);
for q = 1:tut20.n
    tut20.vg(q,1) = relgap(tut20.S{1}(:,q,q), tut20.S{2}(:,q,q));
    tut20.vg(q,2) = relgap(tut20.S{3}(:,q,q), tut20.S{4}(:,q,q));
end
tut20.prs = nchoosek(1:tut20.n, 2);
tut20.cg = zeros(size(tut20.prs,1), 2);
for q = 1:size(tut20.prs,1)
    a = tut20.prs(q,1); b = tut20.prs(q,2);
    tut20.cg(q,1) = mean(abs(cpath(tut20.S{1},a,b) - cpath(tut20.S{2},a,b)));
    tut20.cg(q,2) = mean(abs(cpath(tut20.S{3},a,b) - cpath(tut20.S{4},a,b)));
end

fprintf('\n=== Part 2: 20 variables, shipped posterior means (T = %d, n = %d) ===\n', tut20.T, tut20.n);
fprintf('| | Cholesky SV | order-invariant SV |\n|---|---|---|\n');
fprintf('| Variances, median change | %.1f%% | %.1f%% |\n', median(tut20.vg));
fprintf('| Variances, largest change | %.1f%% | %.1f%% |\n', max(tut20.vg));
fprintf('| Correlations, median change | %.3f | %.3f |\n', median(tut20.cg));
fprintf('| Correlations, 90th percentile | %.3f | %.3f |\n', prctile(tut20.cg, 90));
fprintf('| Correlations, largest change | %.3f | %.3f |\n', max(tut20.cg));
fprintf('variances where the Cholesky change exceeds the largest order-invariant change: %d of %d\n', ...
    sum(tut20.vg(:,1) > max(tut20.vg(:,2))), tut20.n);
fprintf('correlations where the Cholesky change exceeds the largest order-invariant change: %d of %d\n', ...
    sum(tut20.cg(:,1) > max(tut20.cg(:,2))), size(tut20.prs,1));
[~, o] = sort(tut20.vg(:,1), 'descend');
fprintf('largest variance changes (Cholesky / order-invariant):\n');
for q = 1:5, fprintf('  %-22s %6.1f%%  %5.1f%%\n', tut20.lab(o(q)), tut20.vg(o(q),:)); end
[~, o] = sort(tut20.cg(:,1), 'descend');
fprintf('largest correlation changes (Cholesky / order-invariant):\n');
for q = 1:5
    fprintf('  %-22s x %-22s %.3f  %.3f\n', tut20.lab(tut20.prs(o(q),1)), tut20.lab(tut20.prs(o(q),2)), tut20.cg(o(q),:));
end
tut20.top = tut20.prs(o(1),:);
fprintf('the four core variables inside the 20-variable VAR (Cholesky / order-invariant):\n');
for q = find(tut20.prs(:,2) <= 4)'
    fprintf('  corr(%s, %s): %.3f  %.3f\n', tut20.lab(tut20.prs(q,1)), tut20.lab(tut20.prs(q,2)), tut20.cg(q,:));
end
for q = 1:4, fprintf('  var(%s): %.1f%%  %.1f%%\n', tut20.lab(q), tut20.vg(q,:)); end

tut_v = zeros(tut20.T, 4);  tut_c = zeros(tut20.T, 4);
for q = 1:4
    tut_v(:,q) = tut20.S{q}(:,4,4);
    tut_c(:,q) = cpath(tut20.S{q}, tut20.top(1), tut20.top(2));
end
fprintf('time average of the %s x %s correlation (Cholesky 1, 2, order-invariant 1, 2): %.2f %.2f %.2f %.2f\n', ...
    tut20.lab(tut20.top(1)), tut20.lab(tut20.top(2)), mean(tut_c));
fprintf('its range over time under order-invariant SV, published order: %.2f to %.2f\n', min(tut_c(:,3)), max(tut_c(:,3)));
fprintf('peak of the fed funds variance (Cholesky 1, 2, order-invariant 1, 2): %.2f %.2f %.2f %.2f\n', max(tut_v));
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 18 13]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; draw_four(tut20.dates, tut_v, tut_style, tut_legend);
title('Variance of the federal funds rate equation, 20-variable VAR');
nexttile; draw_four(tut20.dates, tut_c, tut_style, {});
title(sprintf('Correlation between the %s and %s equations', tut20.lab(tut20.top(1)), tut20.lab(tut20.top(2))));
exportgraphics(fig, fullfile(tut_dir, 'fig_20var_paths.png'), 'Resolution', 150);
close(fig)

fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 18 8]);
[~, o] = sort(tut20.cg(:,1), 'descend');
hold on; box off
plot(1:numel(o), tut20.cg(o,1), 'o', 'Color', tut_style{1}{2}, 'MarkerSize', 3.5, 'MarkerFaceColor', tut_style{1}{2});
plot(1:numel(o), tut20.cg(o,2), 'o', 'Color', [0 0 0], 'MarkerSize', 3.5, 'MarkerFaceColor', [0 0 0]);
xlabel('the 190 correlations, ranked by the change under Cholesky SV');
ylabel('mean absolute change');
legend({'Cholesky SV', 'order-invariant SV'}, 'Location', 'northeast');
xlim([0 numel(o)+1]);
exportgraphics(fig, fullfile(tut_dir, 'fig_20var_gaps.png'), 'Resolution', 150);
close(fig)

%% ------------------------------------------------------------------
%  Part 3. The forecast comparison, recomputed from the shipped forecasts
%  ------------------------------------------------------------------
tut_tmp = tempname;
copyfile(tut_pkg, tut_tmp);                        % the legacy folder is never written to
tut_here = pwd;
run(fullfile(tut_tmp, 'Table3_forecasting.m'));    % leaves Table_re and Table_DM behind
cd(tut_here);
rmpath(fullfile(tut_tmp, 'results_mat'));
rmdir(tut_tmp, 's');
tut_fc.re = Table_re;  tut_fc.dm = Table_DM;

    % the recomputed table must match the capture in tests/golden, number for number
tut_gold = fileread(fullfile(tut_repo, 'tests', 'golden', 'chan_koop_yu2024_jbes_oisv', ...
    'Table3_forecasting_20260901_1300', 'golden_log_Table3_forecasting.txt'));
tut_gold = extractBetween(tut_gold, 'Table_re_withDM =', 'Elapsed');     % the table alone
tut_gnum = str2double(regexp(tut_gold{1}, '-?\d+\.\d{3}(?!\d)', 'match'));
tut_new = str2double(compose('%.3f', reshape(tut_fc.re', 1, [])));
assert(isequal(tut_gnum, tut_new), 'recomputed forecast table differs from the golden capture');
fprintf('\nforecast table recomputed from the shipped forecasts: matches the golden capture (%d numbers)\n', numel(tut_new));

tut_fc.vars = ["IP" "unemployment" "PCE inflation" "fed funds"];
tut_fc.rows = ["Cholesky SV, published order" "Cholesky SV, reversed order" ...
               "order-invariant SV, published order" "order-invariant SV, reversed order"];
tut_star = ["" "*" "**" "***"];
fprintf('\n=== Part 3: forecasts, 1970:03 onward ===\n');
fprintf('| Variable | Model | RMSFE h=1 | RMSFE h=6 | RMSFE h=12 | ALPL h=1 | ALPL h=6 | ALPL h=12 |\n');
fprintf('|---|---|---|---|---|---|---|---|\n');
for v = 1:4
    for m = 1:4
        r = 4*(v-1) + m;
        cells = strings(1, 6);
        for c = 1:6
            s = tut_fc.dm(r,c);
            if isnan(s), s = 0; end
            cells(c) = sprintf('%.4g', tut_fc.re(r,c));
            if c > 3, cells(c) = sprintf('%.3f', tut_fc.re(r,c)); end
            cells(c) = cells(c) + tut_star(s + 1);
        end
        fprintf('| %s | %s | %s |\n', tut_fc.vars(v), tut_fc.rows(m), strjoin(cells, ' | '));
    end
end
fprintf(['RMSFE: lower is better. ALPL: average log predictive likelihood, higher is better. ' ...
    'Stars: two-sided Diebold-Mariano test against Cholesky SV in the published order, 10/5/1%%.\n']);
tut_rng = zeros(4, 3);
for v = 1:4
    blk = tut_fc.re(4*(v-1)+(1:4), 1:3);
    tut_rng(v,:) = 100*(max(blk) - min(blk))./min(blk);
end
fprintf('largest spread of RMSFE across the four models, by variable and horizon: %.1f%%\n', max(tut_rng(:)));
tut_pair = @(a,b) 100*abs(tut_fc.re(a,1:3) - tut_fc.re(b,1:3))./min(tut_fc.re(a,1:3), tut_fc.re(b,1:3));
tut_cs = []; tut_oi = [];
for v = 1:4
    r = 4*(v-1);
    tut_cs = [tut_cs, tut_pair(r+1, r+2)]; %#ok<AGROW>
    tut_oi = [tut_oi, tut_pair(r+3, r+4)]; %#ok<AGROW>
end
fprintf('largest RMSFE difference between the two orderings: Cholesky SV %.1f%%, order-invariant SV %.1f%%\n', ...
    max(tut_cs), max(tut_oi));

fprintf('\nbuild finished in %.1f minutes\n', toc(tut_t0)/60);
diary off

function draw_four(x, Y, style, labels)
% The four paths in the columns of Y (Cholesky 1 and 2, order-invariant 1 and 2),
% order-invariant first so the thinner Cholesky lines stay visible on top.
hold on; box off
h = gobjects(1, 4);
for q = [3 4 2 1]
    h(q) = plot(x, Y(:,q), style{q}{1}, 'Color', style{q}{2}, 'LineWidth', style{q}{3});
end
if ~isempty(labels), legend(h, labels, 'Location', 'northeast'); end
end
