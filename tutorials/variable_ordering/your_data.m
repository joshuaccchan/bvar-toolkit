%% your_data - the Cholesky and order-invariant VARs with stochastic volatility on your data
%
% Set the file, the columns, their names, the lag length and the chain length
% below, then run the script. It estimates both models with the variables in the
% order given and in the reverse order, prints how much each model's correlations
% change when the order is reversed, and plots the correlation that changes most
% under the Cholesky model. The defaults use the four FRED-MD series of the
% tutorial with short chains; Chan, Koop and Yu (2024) use 30,000 draws after
% 5,000 burn-in.
%
% The file is read with readmatrix, so it may have a header row, and the column
% numbers count every column of the file. The selected columns must be stationary,
% transformed as needed, and have no missing values inside the sample; rows missing
% at either end are dropped. Setting datecol checks that the rows kept are evenly
% spaced. The first n0 rows serve as initial conditions: at least max(p,4) of them.

repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(repo, 'setup.m'))

%% ---- settings ----
file   = fullfile(repo, 'replications', 'chan_koop_yu2024_jbes_oisv', 'legacy', 'FRED_MD_20vars.csv');
cols   = [4 6 12 13];                          % columns of the file, in the order to use
names  = ["IP" "unemployment" "PCE inflation" "fed funds"];
datecol = [];                                  % the column holding the dates; [] to skip the check
p      = 13;                                   % lags
n0     = 24;                                   % rows used as initial conditions
nsim   = 1000;                                 % draws kept
burnin = 200;                                  % draws discarded first
seed   = 1;

%% ---- data ----
raw = readmatrix(file);
sel = raw(:, cols);
assert(numel(names) == numel(cols), 'names and cols must have the same length');
keep = all(isfinite(sel), 2);
assert(any(keep), 'every row of the selected columns has a missing value');
lo = find(keep, 1);  hi = find(keep, 1, 'last');
assert(all(keep(lo:hi)), 'the selected columns have missing values inside the sample');
data = sel(lo:hi, :);
flat = std(data) == 0;
assert(~any(flat), 'these series are constant: %s', strjoin(names(flat), ', '));
if ~isempty(datecol)
    step = diff(raw(lo:hi, datecol));
    assert(all(step > 0) && max(abs(step - median(step))) <= 1e-6*max(1, abs(median(step))), ...
        'the dates of the rows kept are not evenly spaced');
end

%% ---- estimation: both models, the order given and its reverse ----
Y0 = data(1:n0, :);
Y  = data(n0+1:end, :);
n = numel(cols);
rev = n:-1:1;
runs = {'CS', false; 'CS', true; 'OI', false; 'OI', true};
orders = ["given" "reverse"];
res = cell(4, 1);
for r = 1:4
    ord = 1:n;
    if runs{r,2}, ord = rev; end
    t0 = tic;
    res{r} = bvar.models.var_sv(Y0(:,ord), Y(:,ord), p, 'model', runs{r,1}, ...
        'nsim', nsim, 'burnin', burnin, 'seed', seed);
    if runs{r,2}                               % back to the order given
        res{r}.Sig_mean = res{r}.Sig_mean(:, rev, rev);
    end
    fprintf('%s, %s order: %.0f s\n', runs{r,1}, orders(runs{r,2} + 1), toc(t0));
end

%% ---- how much reversing the order changes each correlation ----
corrpath = @(S,i,j) S(:,i,j)./sqrt(S(:,i,i).*S(:,j,j));
prs = nchoosek(1:n, 2);
gap = zeros(size(prs,1), 2);
fprintf('\naverage absolute change in each correlation when the order is reversed\n');
fprintf('%-40s %10s %16s\n', '', 'Cholesky', 'order-invariant');
for q = 1:size(prs,1)
    i = prs(q,1); j = prs(q,2);
    gap(q,1) = mean(abs(corrpath(res{1}.Sig_mean,i,j) - corrpath(res{2}.Sig_mean,i,j)));
    gap(q,2) = mean(abs(corrpath(res{3}.Sig_mean,i,j) - corrpath(res{4}.Sig_mean,i,j)));
    fprintf('%-40s %10.3f %16.3f\n', names(i) + ", " + names(j), gap(q,:));
end
fprintf(['\nThe order-invariant model is the same model in both orders, so its column is\n' ...
    'Monte Carlo error. Correlations are computed from the posterior mean of Sigma_t.\n']);

%% ---- the correlation that changes most under the Cholesky model ----
[~, q] = max(gap(:,1));
i = prs(q,1); j = prs(q,2);
lbl = {'Cholesky, order given', 'Cholesky, reverse order', ...
       'order-invariant, order given', 'order-invariant, reverse order'};
sty = {'-', '--', '-', ':'};
col = {[0 0.447 0.741], [0.85 0.325 0.098], [0.62 0.62 0.62], [0 0 0]};
wid = [1 1 2.6 1.4];
figure; hold on; box off
h = gobjects(1, 4);
for r = [3 4 2 1]
    h(r) = plot(corrpath(res{r}.Sig_mean, i, j), sty{r}, 'Color', col{r}, 'LineWidth', wid(r));
end
legend(h, lbl, 'Location', 'best');
title(sprintf('Correlation between the %s and %s equations', names(i), names(j)));
xlabel('observation'); hold off
