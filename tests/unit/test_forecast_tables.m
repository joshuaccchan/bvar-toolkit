function test_forecast_tables
% deterministic equivalence of bvar.forecast.tables with the legacy
% accumulation + RMSFE/ALPL table code, on synthetic accumulated arrays.
% The legacy side is NOT retyped: the exact line ranges are sliced out of the
% two frozen legacy main_forecasting.m files at test time (content-asserted)
% into tempdir scripts and dispatched in a primed workspace:
%   chan2021_ijf_mahp        lines 98-105 (accumulation), 111-116 (tables)
%   chan2020_springer_largebvar lines 147-154 (accumulation), 160-172 (tables,
%     both the model==1 all-variable form and the model~=1 var_core form)
% One vintage's synthetic tmpyhat is complex-typed with zero imaginary part,
% locking in the legacy magnitude-based max() behavior documented in
% bvar.forecast.iterate/tables.
root = getappdata(0, 'bvar_repo_root');
mahp_main = fullfile(root, 'replications', 'chan2021_ijf_mahp', 'legacy', 'main_forecasting.m');
spr_main = fullfile(root, 'replications', 'chan2020_springer_largebvar', 'legacy', 'main_forecasting.m');

tmpdir = tempname; mkdir(tmpdir);
ctmp = onCleanup(@() cleanup_tmp(tmpdir));   %#ok<NASGU>
slice(mahp_main, 98, 105, 'tmpmax = max(tmpyhat1(:,n+1:end));', fullfile(tmpdir, 'legacy_accum_mahp.m'));
slice(mahp_main, 111, 116, 'RMSFE_1 = sqrt(mean((yhat1(4:end,1:n)', fullfile(tmpdir, 'legacy_tables_mahp.m'));
slice(spr_main, 147, 154, 'tmpmax = max(tmpyhat0(:,n+1:end));', fullfile(tmpdir, 'legacy_accum_springer.m'));
slice(spr_main, 160, 172, 'if model==1', fullfile(tmpdir, 'legacy_tables_springer.m'));
addpath(tmpdir);

check_mahp(tmpdir);
check_springer(tmpdir);
end

% ---------------------------------------------------------------------------
function cleanup_tmp(tmpdir)
if any(strcmpi(strsplit(path, pathsep), tmpdir))
    rmpath(tmpdir);
end
if exist(tmpdir, 'dir')
    rmdir(tmpdir, 's');
end
end

function slice(srcfile, a, b, first_line_marker, outfile)
% write lines a..b of the frozen legacy file to outfile, guarding the range
% against drift by asserting the marker appears in line a
txt = fileread(srcfile);
lines = regexp(txt, '\r?\n', 'split');
assert(numel(lines) >= b, 'slice: %s has fewer than %d lines', srcfile, b);
assert(~isempty(strfind(lines{a}, first_line_marker)), ...
    'slice: line %d of %s does not contain "%s"', a, srcfile, first_line_marker); %#ok<STREMP>
fid = fopen(outfile, 'w');
fprintf(fid, '%s\n', lines{a:b});
fclose(fid);
end

% ---------------------------------------------------------------------------
function check_mahp(tmpdir) %#ok<INUSD>
% MAHP variant: yhat1 (h=1) + yhat4 (h=4, guarded t<=T-4), then the
% RMSFE_1/RMSFE_4 + aveprelike tables over all n variables and the joint
n = 5; T0 = 6; T = 16; nsim = 9;
rng(20260901, 'twister');
Y = randn(T+1, n);      % outturns; rows t+1..t+4 are read
C1 = cell(T, 1); C4 = cell(T, 1);
for t = T0:T-1
    C1{t} = randn(nsim, 2*n+1) - 1;
    C4{t} = randn(nsim, 2*n+1) - 2;
end
C1{T0+2} = complex(C1{T0+2}, 0);    % zero-imag complex vintage (magnitude max)
C4{T0+2} = complex(C4{T0+2}, 0);

    % legacy side: primed workspace + sliced scripts
yhat1 = zeros(T-T0, 3*n+1); yhat4 = zeros(T-4-T0+1, 3*n+1);
RMSFE_1 = []; RMSFE_4 = []; RMSFE = []; aveprelike_1 = []; aveprelike_4 = []; aveprelike = [];
for t = T0:T-1
    tmpyhat1 = C1{t}; tmpyhat4 = C4{t};                     %#ok<NASGU>
    legacy_accum_mahp;
end
legacy_tables_mahp;
Lyhat1 = yhat1; Lyhat4 = yhat4;

    % functionized side
yhat1 = zeros(T-T0, 3*n+1); yhat4 = zeros(T-4-T0+1, 3*n+1);
for t = T0:T-1
    yhat1(t-T0+1, :) = bvar.forecast.tables('accum_row', C1{t}, Y(t+1, :));
    if t <= T-4
        yhat4(t-T0+1, :) = bvar.forecast.tables('accum_row', C4{t}, Y(t+4, :));
    end
end
S = bvar.forecast.tables('mahp', yhat1, yhat4);

assert(isequal(Lyhat1, yhat1), 'mahp: accumulated yhat1 differs');
assert(isequal(Lyhat4, yhat4), 'mahp: accumulated yhat4 differs');
assert(isequal(RMSFE_1, S.RMSFE_1) && isequal(RMSFE_4, S.RMSFE_4) ...
    && isequal(RMSFE, S.RMSFE), 'mahp: RMSFE tables differ');
assert(isequal(aveprelike_1, S.aveprelike_1) && isequal(aveprelike_4, S.aveprelike_4) ...
    && isequal(aveprelike, S.aveprelike), 'mahp: ALPL tables differ');
end

% ---------------------------------------------------------------------------
function check_springer(tmpdir) %#ok<INUSD>
% springer variant: yhat0 (h=0) + yhat1 (h=1, guarded t<=T-1), then both
% table forms: model==1 (all variables + joint ALPL) and model~=1 (var_core
% RMSFE, per-variable-only ALPL)
n = 6; T0 = 5; T = 14; nsim = 8;
var_core = [1 3 6]';    % column, like the legacy [1 7 8 12]'
rng(20260902, 'twister');
O0 = cell(T, 1); O1 = cell(T, 1); C0 = cell(T, 1); C1 = cell(T, 1);
for t = T0:T
    O0{t} = randn(1, n); O1{t} = randn(1, n);
    C0{t} = randn(nsim, 2*n+1) - 1;
    C1{t} = randn(nsim, 2*n+1) - 2;
end
C0{T0+1} = complex(C0{T0+1}, 0);    % zero-imag complex vintage
C1{T0+1} = complex(C1{T0+1}, 0);

    % legacy side
yhat0 = zeros(T-T0+1, 3*n+1); yhat1 = zeros(T-1-T0+1, 3*n+1);
for t = T0:T
    data_tp0 = O0{t}; data_tp1 = O1{t};                     %#ok<NASGU>
    tmpyhat0 = C0{t}; tmpyhat1 = C1{t};                     %#ok<NASGU>
    legacy_accum_springer;
end
Lyhat0 = yhat0; Lyhat1 = yhat1;
Lt = cell(2, 1);
for model = [1 2]
    RMSFE_0 = []; RMSFE_1 = []; RMSFE = []; aveprelike_0 = []; aveprelike_1 = []; aveprelike = [];
    legacy_tables_springer;
    Lt{1+(model~=1)} = struct('RMSFE_0', RMSFE_0, 'RMSFE_1', RMSFE_1, 'RMSFE', RMSFE, ...
        'aveprelike_0', aveprelike_0, 'aveprelike_1', aveprelike_1, 'aveprelike', aveprelike);
end

    % functionized side
yhat0 = zeros(T-T0+1, 3*n+1); yhat1 = zeros(T-1-T0+1, 3*n+1);
for t = T0:T
    yhat0(t-T0+1, :) = bvar.forecast.tables('accum_row', C0{t}, O0{t});
    if t <= T-1
        yhat1(t-T0+1, :) = bvar.forecast.tables('accum_row', C1{t}, O1{t});
    end
end
S1 = bvar.forecast.tables('springer', yhat0, yhat1, []);         % model==1 form
S2 = bvar.forecast.tables('springer', yhat0, yhat1, var_core);   % model~=1 form

assert(isequal(Lyhat0, yhat0), 'springer: accumulated yhat0 differs');
assert(isequal(Lyhat1, yhat1), 'springer: accumulated yhat1 differs');
for kv = 1:2
    if kv == 1, S = S1; else, S = S2; end
    assert(isequal(Lt{kv}.RMSFE_0, S.RMSFE_0) && isequal(Lt{kv}.RMSFE_1, S.RMSFE_1) ...
        && isequal(Lt{kv}.RMSFE, S.RMSFE), 'springer form %d: RMSFE tables differ', kv);
    assert(isequal(Lt{kv}.aveprelike_0, S.aveprelike_0) && isequal(Lt{kv}.aveprelike_1, S.aveprelike_1) ...
        && isequal(Lt{kv}.aveprelike, S.aveprelike), 'springer form %d: ALPL tables differ', kv);
end
end
