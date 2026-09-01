% run_unit_tests - run every unit test in this folder; error on any failure.
% Usage (from repo root or anywhere):  matlab -batch "run('tests/unit/run_unit_tests.m')"

thisdir = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(thisdir));
addpath(fullfile(root, 'core'));
addpath(fullfile(root, 'third_party'));
setappdata(0, 'bvt_repo_root', root);

tests = dir(fullfile(thisdir, 'test_*.m'));
nfail = 0;
for ii = 1:numel(tests)
    [~, name] = fileparts(tests(ii).name);
    try
        feval(name);
        fprintf('PASS  %s\n', name);
    catch err
        nfail = nfail + 1;
        fprintf('FAIL  %s: %s\n', name, err.message);
    end
end
if nfail > 0
    error('%d unit test(s) failed', nfail);
end
fprintf('All %d unit tests passed.\n', numel(tests));
