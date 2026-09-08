% bvar.samplers.acp_theta_sig - draw the structural VAR parameters under the
% asymmetric conjugate prior, equation by equation and INDEPENDENTLY across
% draws: the prior is conjugate, so each equation's posterior is
% normal-inverse-gamma in closed form and nsim draws come out in one call
% without a Markov chain.
%
%   [Alp,Beta,Sig] = bvar.samplers.acp_theta_sig(Y0, Y, p, prior, nsim)
%
%   Y0, Y : initial conditions and the estimation sample; the lag matrix is
%           rebuilt here rather than passed in, so this call is self-contained
%   prior : the struct from bvar.priors.acp_redu or acp_stru (fields beta0,
%           Vbeta, alp0, Valp, nu, S)
%   Alp   : nsim x n(n-1)/2 free elements of the unit-lower-triangular A
%   Beta  : nsim x (n^2 p + n) coefficients, equation by equation
%   Sig   : nsim x n structural innovation variances
%
% Equation ii regresses y_ii on the lags and on the CONTEMPORANEOUS values of
% the preceding equations, Xi = [Z, -Y(:,1:ii-1)], which is what makes A unit
% lower triangular and the whole system a set of independent regressions. Use
% bvar.structural.reduced_form to map the draws back to the reduced form.
%
% rng consumption per equation ii, in order: gamrnd(...,nsim,1) then
% randn(nsim,ki) with ki = n p + ii. The whole block for one equation is drawn
% at once, so the stream position after the call depends on nsim as well as on
% n and p.
%
% Body from chan2022_qe_acp/legacy/utility/sample_ThetaSig.m, renamed.
% Equivalence: tests/unit/test_acp_equivalence.m. Record: tests/variant_map.md.
%
% See:
% Chan, J.C.C. (2022). Asymmetric Conjugate Priors for Large Bayesian VARs,
% Quantitative Economics, 13(3): 1145-1169.

function [Alp,Beta,Sig] = acp_theta_sig(Y0,Y,p,prior,nsim)
[T,n] = size(Y);
tmpY = [Y0(end-p+1:end,:); Y];
Z = zeros(T,n*p);
for ii=1:p
    Z(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
end
Z = [ones(T,1) Z];
k_beta = n^2*p+n;
k_alp = n*(n-1)/2;
Beta = zeros(nsim,k_beta);
Alp = zeros(nsim,k_alp);
Sig = zeros(nsim,n);
count_alp = 0;
for ii = 1:n
    yi = Y(:,ii);
    ki = n*p+ii;
    mi = [prior.beta0(:,ii);prior.alp0(count_alp+1:count_alp+ii-1)];
    Vi = sparse(1:ki,1:ki,[prior.Vbeta(:,ii);prior.Valp(count_alp+1:count_alp+ii-1)]);
    nui = prior.nu(ii);
    Si = prior.S(ii);
    Xi = [Z -Y(:,1:ii-1)];
        % compute the parameters of the posterior distribution
    iVi = Vi\speye(ki);
    Kthetai = iVi + Xi'*Xi;
    CKthetai = chol(Kthetai,'lower');
    thetai_hat = (CKthetai')\(CKthetai\(iVi*mi + Xi'*yi));
    Si_hat = Si + (yi'*yi + mi'*iVi*mi - thetai_hat'*Kthetai*thetai_hat)/2;
        % sample sig and theta
    Sigi = 1./gamrnd(nui+T/2,1./Si_hat,nsim,1);
    U = randn(nsim,ki).*repmat(sqrt(Sigi),1,ki);
    Thetai = repmat(thetai_hat',nsim,1) + U/CKthetai;

    Sig(:,ii) = Sigi;
    Beta(:,(ii-1)*(n*p+1)+1:ii*(n*p+1)) =  Thetai(:,1:n*p+1);
    Alp(:,count_alp+1:count_alp+ii-1) =  Thetai(:,n*p+2:end);
    count_alp = count_alp + ii -1;
end
end
