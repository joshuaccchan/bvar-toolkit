% run_unit_tests - run every unit test in this folder; error on any failure.
% Usage (from repo root or anywhere):  matlab -batch "run('tests/unit/run_unit_tests.m')"
%
% Set the environment variable BVAR_SKIP_TESTS to a comma-separated list of test
% names to leave out, e.g.
%
%   BVAR_SKIP_TESTS=test_forecast_iterate_springer
%
% Skipped tests are named individually in the output and counted in the summary,
% so a run that covers less than the whole suite says so rather than reporting a
% clean sweep. The CI workflow uses this for the one test that cannot run on a
% hosted runner (see .github/workflows/unit-tests.yml). Skipping is for tests
% that CANNOT run in an environment, never for tests that fail in it.

thisdir = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(thisdir));
addpath(fullfile(root, 'core'));
addpath(fullfile(root, 'third_party'));
setappdata(0, 'bvar_repo_root', root);

skip = strtrim(strsplit(getenv('BVAR_SKIP_TESTS'), ','));
skip = skip(~cellfun(@isempty, skip));

tests = dir(fullfile(thisdir, 'test_*.m'));
nfail = 0; nrun = 0; nskip = 0;
for ii = 1:numel(tests)
    [~, name] = fileparts(tests(ii).name);
    if any(strcmp(name, skip))
        nskip = nskip + 1;
        fprintf('SKIP  %s (listed in BVAR_SKIP_TESTS)\n', name);
        continue
    end
    try
        feval(name);
        nrun = nrun + 1;
        fprintf('PASS  %s\n', name);
    catch err
        nfail = nfail + 1;
        fprintf('FAIL  %s: %s\n', name, err.message);
    end
end

% A name in BVAR_SKIP_TESTS that matches nothing is a typo, and would silently
% widen coverage claims rather than narrow them - fail loudly instead.
names = cellfun(@(f) erase(f, '.m'), {tests.name}, 'UniformOutput', false);
unknown = setdiff(skip, names);
if ~isempty(unknown)
    error('BVAR_SKIP_TESTS names no such test: %s', strjoin(unknown, ', '));
end

if nfail > 0
    error('%d unit test(s) failed', nfail);
end
if nskip > 0
    fprintf('All %d unit tests passed (%d skipped).\n', nrun, nskip);
else
    fprintf('All %d unit tests passed.\n', nrun);
end
