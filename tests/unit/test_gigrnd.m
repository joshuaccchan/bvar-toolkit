function test_gigrnd
% third_party/gigrnd must reproduce the legacy copy draw-for-draw under one seed
root = getappdata(0, 'bvt_repo_root');

rng(5, 'twister');
y_core = gigrnd(0.5, 2, 3, 1);               % third_party copy (on path via runner)

leg = fullfile(root, 'replications', 'chan2021_ijf_mahp', 'legacy');
addpath(leg); c = onCleanup(@() rmpath(leg));  % prepended -> legacy copy shadows ours
rng(5, 'twister');
y_leg = gigrnd(0.5, 2, 3, 1);
assert(isequal(y_leg, y_core), 'gigrnd: third_party copy differs from legacy');
assert(y_core > 0, 'gigrnd: GIG draw must be positive');
end
