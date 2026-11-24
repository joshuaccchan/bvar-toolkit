%% your_data - a mixed-frequency VAR with common stochastic volatility on your data
%
% Set the file, the columns, their names, which of them are quarterly, how each is
% transformed, the sample and the chain length below, then run the script. It estimates
% the VAR of the tutorial, with the priors of Chan, Poon and Zhu (2023), prints the
% monthly estimates of each quarterly series over the last twelve months, plots the
% first quarterly series with its 68 percent credible band, and writes the monthly
% estimates to a csv and a mat file. The defaults use the tutorial's data from 1990 with
% short chains and take about a minute; the tutorial uses 20,000 draws after a burn-in
% period of 2,000.
%
% The file is read with readtable and has one row per month, with the dates in datecol
% in the format datefmt, such as 2026-07 for yyyy-MM. A quarterly series holds its level
% in the last month of each quarter and nothing in the other months. It may start late
% or end early. The script stops if a quarter is missing between its first and last
% values in the file. 'dlog' turns a series into 100 times its log change, quarter on
% quarter for a quarterly series, so that its monthly values are monthly growth rates
% tied to it by the aggregation of Mariano and Murasawa (2003); 'level' uses a monthly
% series as it is, after multiplying by scale. A monthly series may have missing values
% anywhere, and when one in levels is missing inside the sample, the growth rates around
% it keep their known sum. The first p months serve as initial conditions.

repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(repo, 'setup.m'))
lastwarn('');       % the report records the last warning raised from here on

%% ---- settings ----
file      = fullfile(repo, 'tutorials', 'mixed_frequency', 'mf_data.csv');
datecol   = "date";
datefmt   = "yyyy-MM";                         % format of the dates, first and last
cols      = ["INDPRO" "CPIAUCSL" "UNRATE" "PAYEMS" "AWHMAN" "GDPC1" "FPIx"];
names     = ["IP" "CPI" "Unemployment" "Payrolls" "Hours" "GDP" "Investment"];
quarterly = [false false false false false true true];
transform = ["dlog" "dlog" "level" "dlog" "level" "dlog" "dlog"];
scale     = [1 1 1 1 .1 1 1];                  % multiplies a series used in levels
first     = "1990-01";                         % first month of the sample; "" for all
last      = "";                                % last month; "" for all
p         = 12;                                % lags
nsim      = 1000;                              % draws kept
burnin    = 500;                               % draws discarded first
seed      = 1;
outdir    = tempdir;                           % where the report goes; '' for none

%% ---- data ----
n = numel(cols);
assert(numel(names) == n && numel(quarterly) == n && numel(transform) == n && numel(scale) == n, ...
    'cols, names, quarterly, transform and scale must have the same length');
assert(all(ismember(transform, ["dlog" "level"])), 'transform must be "dlog" or "level"');
assert(~any(quarterly & transform ~= "dlog"), 'a quarterly series must be in levels, with transform "dlog"');
tbl = readtable(file, 'TextType', 'string', 'VariableNamingRule', 'preserve');
dates = datetime(tbl.(datecol), 'InputFormat', datefmt);
mi = 12*year(dates) + month(dates);
assert(all(diff(mi) == 1), 'the file must have one row per month, in order');
Lv = tbl{:, cols};
assert(all(Lv(:, transform == "dlog") > 0 | isnan(Lv(:, transform == "dlog")), 'all'), ...
    'a series in log changes must be positive');

% the window, with one month before it for the first growth rate
span = sprintf('the file runs from %s to %s', string(dates(1), 'yyyy-MM'), string(dates(end), 'yyyy-MM'));
lo = 2;  hi = numel(dates);
if strlength(first) > 0
    lo = find(dates >= datetime(first, 'InputFormat', datefmt), 1);
    assert(~isempty(lo), 'first = %s is after the last month in the file; %s', first, span);
end
if strlength(last) > 0
    hi = find(dates <= datetime(last, 'InputFormat', datefmt), 1, 'last');
    assert(~isempty(hi), 'last = %s is before the first month in the file; %s', last, span);
end
assert(lo >= 2, 'first must be %s or later, since the first growth rate needs the month before it; %s', ...
    string(dates(2), 'yyyy-MM'), span);
assert(hi - lo > p + 12, 'the sample from %s to %s is too short for %d lags; %s', ...
    string(dates(lo), 'yyyy-MM'), string(dates(hi), 'yyyy-MM'), p, span);
w = lo-1:hi;
d = dates(lo:hi);
T = numel(d);

X = nan(T, n);
% quarterly series: 100 x log change between consecutive quarterly levels in the file
for j = find(quarterly)
    iobs = find(~isnan(Lv(:,j)));
    assert(all(mod(month(dates(iobs)), 3) == 0), ...
        'series %s must hold its values in the last month of each quarter', names(j));
    k = find(diff(iobs) ~= 3, 1);
    assert(isempty(k), 'series %s has no value for %s, between its first and last values in the file', ...
        names(j), string(dates(iobs(k)) + calmonths(3), 'yyyy-MM'));
    g = nan(size(Lv,1), 1);
    g(iobs(2:end)) = 100*diff(log(Lv(iobs,j)));
    X(:,j) = g(lo:hi);
end
% monthly series in log changes, with a restriction across each gap
jd = find(~quarterly & transform == "dlog");
[G, Mg, zg] = bvar.util.dlog_gaps(Lv(w, jd));
X(:, jd) = G;
jl = find(~quarterly & transform == "level");
X(:, jl) = Lv(lo:hi, jl).*scale(jl);
[Mq, zq, Y] = bvar.util.mm_constraint(X, quarterly);
[r, c] = find(Mg);
r = r(:);  c = c(:);
tg = ceil(c/numel(jd));  jg = c - (tg-1)*numel(jd);
M = [Mq; sparse(r, (tg-1)*n + reshape(jd(jg), [], 1), 1, size(Mg,1), T*n)];
z = [zq; zg];

% the Minnesota scale of a quarterly series: 9/19 of the AR(4) residual variance of its
% quarterly growth rates, the monthly variance that gives it under the aggregation
sig2 = nan(n,1);
for j = find(quarterly)
    g = X(~isnan(X(:,j)), j);
    assert(numel(g) >= 12, 'series %s has fewer than 12 quarterly values in the sample', names(j));
    Zq = [ones(numel(g)-4,1) g(4:end-1) g(3:end-2) g(2:end-3) g(1:end-4)];
    e = g(5:end) - Zq*(Zq\g(5:end));
    sig2(j) = 9/19*mean(e.^2);
end
fprintf('\n%s to %s, %d months, the first %d as initial conditions; %d restrictions\n', ...
    string(d(1), 'yyyy-MM'), string(d(end), 'yyyy-MM'), T, p, numel(z));
fprintf('missing values: %s\n', strjoin(compose('%s %d', names', sum(isnan(Y))'), ', '));

%% ---- estimation ----
t0 = tic;
res = bvar.models.mfvar_csv(Y, p, 'M', M, 'z', z, 'sig2', sig2, 'nsim', nsim, ...
    'burnin', burnin, 'seed', seed);
fprintf('%d draws after %d burn-in: %.0f s; share of h blocks accepted %.2f\n', ...
    nsim, burnin, toc(t0), res.accept_rate);

%% ---- the monthly estimates of the quarterly series ----
jq = find(quarterly);
for j = jq
    fprintf('\n%s, monthly growth in percent: posterior mean (90%% band)\n', names(j));
    for t = T-11:T
        fprintf('  %s %7.2f (%6.2f, %6.2f)\n', string(d(t), 'yyyy-MM'), res.Y_mean(t,j), ...
            res.Y_q(t,j,1), res.Y_q(t,j,end));
    end
end

j = jq(1);
k = d >= d(end) - calyears(10);
figure; hold on; box off; grid on
bvar.util.shaded_band(d(k), res.Y_q(k,j,2), res.Y_q(k,j,4));
plot(d(k), res.Y_mean(k,j), 'k', 'LineWidth', 1.2);
title(sprintf('%s, monthly growth: posterior mean and 68%% credible band', names(j)));
ylabel('percent'); hold off

%% ---- the report ----
% One row per month and quarterly series, and the settings in the mat file.
if ~isempty(outdir)
    rows = numel(jq)*T;
    series = reshape(repmat(names(jq), T, 1), rows, 1);
    mon = repmat(string(d, 'yyyy-MM'), numel(jq), 1);
    est = reshape(res.Y_mean(:, jq), rows, 1);
    qs = reshape(res.Y_q(:, jq, :), rows, []);
    monthly = [table(series, mon, est, 'VariableNames', {'series', 'month', 'mean'}), ...
        array2table(qs, 'VariableNames', compose('q%02d', round(100*res.quantiles)))];
    meta = struct('file', file, 'columns', cols, 'names', names, 'quarterly', quarterly, ...
        'transform', transform, 'scale', scale, 'first', string(d(1), 'yyyy-MM'), ...
        'last', string(d(end), 'yyyy-MM'), 'p', p, 'nsim', nsim, 'burnin', burnin, ...
        'seed', seed, 'sig2', res.sig2', 'accept_rate', res.accept_rate);
    bvar.util.report('mixed_frequency_report', monthly, meta, outdir);
end
