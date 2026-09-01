function test_util_small
% vec, vech, ldet, mgammaln: known values + equivalence with legacy copies
A = [1 2; 3 4];
assert(isequal(bvt.util.vec(A), [1; 3; 2; 4]), 'vec: wrong value');
assert(isequal(bvt.util.vech(A), [1; 3; 4]), 'vech: wrong value');

rng(3, 'twister');
B = randn(5); S = B*B' + 5*eye(5);
assert(abs(bvt.util.ldet(S) - log(det(S))) < 1e-10, 'ldet: disagrees with log(det)');
assert(abs(bvt.util.mgammaln(1, 3.7) - gammaln(3.7)) < 1e-12, 'mgammaln: p=1 should equal gammaln');

root = getappdata(0, 'bvt_repo_root');
leg = fullfile(root, 'replications', 'chan2023_joe_mlvarsv', 'legacy', 'utility');
addpath(leg); c = onCleanup(@() rmpath(leg));
assert(isequal(vec(A), bvt.util.vec(A)), 'vec: differs from legacy');
assert(isequal(vech(A), bvt.util.vech(A)), 'vech: differs from legacy');
assert(isequal(ldet(S), bvt.util.ldet(S)), 'ldet: differs from legacy');
assert(isequal(mgammaln(3, 4.2), bvt.util.mgammaln(3, 4.2)), 'mgammaln: differs from legacy');
end
