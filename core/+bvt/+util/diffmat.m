% bvt.util.diffmat - the sparse first-order difference matrix of a state
% equation, H = I_T - a*L, where L is the first subdiagonal (the lag operator
% on a stacked path).
%
%   H = bvt.util.diffmat(T)        % a = 1: the random-walk difference matrix
%   H = bvt.util.diffmat(T, a)
%
% This is the matrix that turns a state equation into a linear system in the
% whole path at once, which is what makes the precision-based samplers work:
% for x_t = a*x_{t-1} + u_t with u ~ N(0, S), H*x = u, so the path has
% precision H'*inv(S)*H - banded, hence an O(T) Cholesky. See examples/
% ex01_precision_sampler.m.
%
% Sign conventions in this toolkit, both of which appear in the papers:
%   AR(1) / random-walk state:  H_rho  = I - rho*L   ->  diffmat(T, rho)
%   MA(1) error transform:      H_psi  = I + psi*L   ->  diffmat(T, -psi)
%
% NOTE for maintainers: the functions under core/ that were extracted from the
% published packages keep their inline `speye(T) - rho*sparse(2:T,...)`
% expressions verbatim, so that they remain line-by-line diffable against
% replications/*/legacy/. This helper is for new code; it is not retrofitted
% into them. It returns exactly the same matrix (asserted in
% tests/unit/test_diffmat.m).

function H = diffmat(T, a)
if nargin < 2
    a = 1;
end
if ~isscalar(T) || T < 1 || T ~= round(T)
    error('bvt:util:diffmat:badT', 'T must be a positive integer');
end
if ~isscalar(a)
    error('bvt:util:diffmat:badA', 'a must be a scalar');
end
if T == 1
    H = speye(1);
    return
end
H = speye(T) - a*sparse(2:T, 1:(T-1), ones(1, T-1), T, T);
end
