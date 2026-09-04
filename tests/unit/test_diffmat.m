function test_diffmat
% bvar.util.diffmat must reproduce the inline state-equation matrices that the
% legacy-derived core functions build by hand, exactly.

T = 40;

% random-walk form, as built in bvar.sv.ksc_rw_h0 (spdiags spelling)
H_rw = speye(T) - spdiags(ones(T-1,1), -1, T, T);
assert(isequal(bvar.util.diffmat(T), H_rw), 'diffmat: RW form differs');
assert(isequal(bvar.util.diffmat(T, 1), H_rw), 'diffmat: a=1 must equal the RW form');

% AR(1) form, as built in bvar.sv.csv_armh (sparse spelling)
rho = 0.937;
H_ar = speye(T) - rho*sparse(2:T, 1:(T-1), ones(1,T-1), T, T);
assert(isequal(bvar.util.diffmat(T, rho), H_ar), 'diffmat: AR(1) form differs');

% MA(1) transform Hpsi = I + psi*L, as built in bvar.ml.kron_bvar_ma
psi = -0.42;
H_ma = speye(T) + psi*sparse(2:T, 1:(T-1), ones(1,T-1), T, T);
assert(isequal(bvar.util.diffmat(T, -psi), H_ma), 'diffmat: MA(1) form differs');

% structure: unit diagonal, one subdiagonal, sparse, right size
H = bvar.util.diffmat(T, rho);
assert(issparse(H), 'diffmat: must be sparse');
assert(isequal(size(H), [T T]), 'diffmat: wrong size');
assert(nnz(H) == T + (T-1), 'diffmat: wrong sparsity');
assert(all(full(diag(H)) == 1), 'diffmat: diagonal must be 1');
assert(isequal(full(diag(H,-1)), -rho*ones(T-1,1)), 'diffmat: wrong subdiagonal');

% the defining property: H*x differences a path
x = cumsum(randn(T,1));
assert(norm(H_rw*x - [x(1); diff(x)]) < 1e-12, 'diffmat: H*x must difference x');

% edge case
assert(isequal(bvar.util.diffmat(1, rho), speye(1)), 'diffmat: T=1');

% bad input rejected
try
    bvar.util.diffmat(0);
    error('test:noThrow', 'diffmat should reject T = 0');
catch err
    assert(strcmp(err.identifier, 'bvar:util:diffmat:badT'), 'diffmat: wrong error id');
end
end
