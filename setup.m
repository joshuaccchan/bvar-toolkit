% setup.m - put bvar-toolkit on the MATLAB path for this session.
%
%   run setup.m          (from anywhere; the script locates the repo itself)
%
% Adds the repo's core/ and third_party/ folders to the path. core/ contains
% the MATLAB package folder +bvar, so functions are called with the package
% prefix once core/ is on the path:
%
%   [Y, Z] = bvar.util.build_lags(Yraw, 4);
%   h      = bvar.sv.ksc_rw_h0(ystar, h, sigh2, h0);
%
% Nothing is added permanently - the path change lasts for the session. Add
% this line to your own startup.m if you want it every time:
%
%   run('<path-to-repo>/setup.m')
%
% Replication drivers are NOT added to the path: each lives beside the paper
% it reproduces (replications/<paper>/run_all.m) and expects to be called
% with that folder on the path or as the working directory, so that its
% legacy data files resolve. See the README.

bvar_root = fileparts(mfilename('fullpath'));
addpath(fullfile(bvar_root, 'core'));
addpath(fullfile(bvar_root, 'third_party'));

if verLessThan('matlab', '9.1')      % R2016b
    warning('bvar:setup:oldMATLAB', ...
        ['This toolkit is tested on R2025b. Releases before R2016b lack ' ...
         'features used here (implicit expansion, string handling).']);
end
if isempty(ver('stats'))
    warning('bvar:setup:noStats', ...
        ['The Statistics and Machine Learning Toolbox is not available. ' ...
         'Most samplers need gamrnd/iwishrnd/normpdf and will fail.']);
end

fprintf('bvar-toolkit: core/ and third_party/ added to the path (%s).\n', bvar_root);
fprintf('  quick start:  cd examples; ex01_precision_sampler\n');
fprintf('  run tests:    run(fullfile(''%s'',''tests'',''unit'',''run_unit_tests.m''))\n', bvar_root);

clear bvar_root
