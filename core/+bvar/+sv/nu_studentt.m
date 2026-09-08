% bvar.sv.nu_studentt - independence-MH draw of the Student-t degree-of-freedom
% parameter nu given the latent mixture scales lam (lam_t ~ IG(nu/2, nu/2)),
% under a uniform prior on (2, nu_ub): Newton-Raphson mode nut, Gaussian
% N(nut, -1/H) proposal, MH accept-reject.
%
% Body from chan2020_springer_largebvar/legacy/sample_nu.m, renamed, with one
% deviation: a duplicated recomputation of T/sum1/sum2 - three lines reassigning
% identical values, no rng calls - was removed. Also stands in for the two
% sample_nu copies in chan2020_jbes_kronecker; one of those forms the MH ratio
% from normpdf ratios rather than in logs. That is mathematically identical and
% its accept/reject decisions matched on every seed tested, but the two are not
% bitwise-identical in alpha; the log-form is kept because it has no Inf*0 = NaN
% hazard when |nu - nut| is extreme.
% Equivalence: tests/unit/test_nu_studentt.m. Record: tests/variant_map.md.
%
% Inputs:  lam   - T x 1 latent scale draws
%          nu    - current df draw
%          nu_ub - upper bound of the uniform prior support
% Outputs: nu, flag (1 if accepted), f_nu (handle: log target kernel of nu)
%
% See:
% Chan, J.C.C. (2020). Large Bayesian Vector Autoregressions. In: P. Fuleky (Eds),
% Macroeconomic Forecasting in the Era of Big Data, 95-125, Springer, Cham.
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics, 38(1), 68-79.

function [nu, flag, f_nu] = nu_studentt(lam,nu,nu_ub)
flag = 0;
T = size(lam,1);
sum1 = sum(log(lam));
sum2 = sum(1./lam);
f_nu = @(x) T*(x/2.*log(x/2)-gammaln(x/2)) - (x/2+1)*sum1 - x/2*sum2;
df_nu = @(x) T/2*(log(x/2) + 1 - psi(x/2)) - .5*(sum1+sum2);
d2f_nu = @(x) T/(2*x) - T/4*psi(1,x/2);
S_nu = 1;
nut = nu;
while abs(S_nu) > 10^(-5)   % stopping criteria
    S_nu = df_nu(nut);
    H_nu = d2f_nu(nut);
    nut = nut - H_nu\S_nu;
    if nut<2
        nut = 5;
        H_nu = d2f_nu(nut);
        break;
    end
end
Dnu = -1/H_nu;
nuc = nut + sqrt(Dnu)*randn;
if nuc > 2 && nuc < nu_ub
    lalp_MH = f_nu(nuc) - f_nu(nu) - .5*(nu-nut)^2/Dnu + .5*(nuc-nut)^2/Dnu;
    if exp(lalp_MH) > rand
        nu = nuc;
        flag = 1;
    end
end
end
