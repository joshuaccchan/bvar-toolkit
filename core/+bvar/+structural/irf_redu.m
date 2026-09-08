% bvar.structural.irf_redu - impulse responses of a reduced-form VAR to the
% structural shocks identified by an impact matrix L, computed by powering the
% companion matrix.
%
%   response = bvar.structural.irf_redu(A, L, nstep, nshock)
%
%   A      : (n p) x n reduced-form coefficients, EXCLUDING the intercept -
%            column ii holds equation ii's coefficients, lag-1 block first
%   L      : n x n impact matrix with L*L' the reduced-form covariance; the
%            first nshock columns are the identified shocks
%   nstep  : number of horizons returned, counting the impact period
%   response(:,:,it) : responses at horizon it-1, so response(:,:,1) = impact
%
% HORIZON INDEXING. response(:,:,1) is the impact response L(:,1:nshock), and
% response(:,:,it) for it >= 2 applies the (1:n,1:n) block of the companion
% matrix raised to the (it-1)th power. Two other copies of this routine in the
% wider codebase start the loop at the wrong index and so report every response
% shifted one period; this package's copy is the corrected one and is the reason
% it, rather than another copy, is canonical here. Anyone comparing figures with
% an older run should check which copy produced them.
%
% Body from chan2022_qe_acp/legacy/utility/IRredu.m (the BVAR_ACP_R1
% zip of 2026-08-27, the copy carrying the fix), renamed.
% Equivalence: tests/unit/test_acp_equivalence.m. Record: tests/variant_map.md.
%
% See:
% Chan, J.C.C. (2022). Asymmetric Conjugate Priors for Large Bayesian VARs,
% Quantitative Economics, 13(3): 1145-1169.

function response = irf_redu(A,L,nstep,nshock)
[np, n] = size(A);
p = np/n;
response = zeros(n,nshock,nstep);
Acomp = [A'; sparse(1:n*(p-1),1:n*(p-1),ones(n*(p-1),1),n*(p-1),np)];
response(:,:,1) = L(:,1:nshock);
Apower = Acomp;
for it = 2:nstep
    response(:,:,it) = Apower(1:n,1:n)*L(:,1:nshock);
    Apower = Apower*Acomp;
end
end
