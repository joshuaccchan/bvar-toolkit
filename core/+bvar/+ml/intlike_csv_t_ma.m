% bvar.ml.intlike_csv_t_ma - integrated likelihood of the BVAR-CSV-t-MA model
% at (A,Sig,psi,rho,sigh2,nu): the Student-t scales lam are integrated
% analytically and the common log-volatility path h by importance sampling
% (R draws), with the EM-within-Newton mode search of the t-CSV case run on
% the MA(1)-transformed residuals Utld = Hpsi\U. Consumes R*T randn calls.
%
% Extracted 2026-09-02 (step 8, Kronecker family pass). Canonical source:
%   chan2020_jbes_kronecker/legacy/intlike_BVAR_CSV_t_MA.m  (body verbatim;
%   renamed intlike_BVAR_CSV_t_MA -> intlike_csv_t_ma). Called by the
% BVAR-CSV-t-MA marginal likelihood with R = 10000 (legacy
% ml_BVAR_CSV_t_MA.m lines 20-21). Verified by tests/unit/test_kron_intlike.m
% (bitwise, seeded) and end-to-end through the ML equivalence tests.
%
% Known legacy quirk (kept verbatim in both bugcompat and corrected ML modes -
% it is a likelihood-formula property, not an evaluation-point inconsistency):
% unlike the Gaussian intlike_csv_ma, whose deny_h scales the first
% observation by (1+psi^2)exp(h_1) and carries the -n/2*log(1+psi^2)
% constant, this function's deny_h receives the already-transformed Utld and
% applies NO (1+psi^2) first-observation correction and NO -n/2*log(1+psi^2)
% term - the first transformed observation is treated as t with scale
% exp(h_1)*Sig, whereas the Gibbs sampler it pairs with (BVAR_CSV_t_MA.m
% line 59) gives it variance (1+psi^2)*exp(h_1)*lam_1*Sig. A one-observation
% (out of T) mismatch between the estimated model and this likelihood
% ordinate; the published BVAR-CSV-t-MA marginal likelihood includes it.
% Recorded in tests/variant_map.md (step-8 audit).
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics, 38(1), 68-79.

function [intlike,store_llike] = intlike_csv_t_ma(shortY,X,A,Sig,psi,rho,sigh2,nu,R)
[T,n] = size(shortY);
Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);

% obtain the mode and negative Hessian of the conditional density of h
CSig = chol(Sig,'lower');
U = shortY - X*A;
Utld = Hpsi\U;
tmp = (Utld/CSig');
ht = log(mean(sum(tmp.^2,2)))*ones(T,1);
HiSH = Hrho'*sparse(1:T,1:T,[(1-rho^2)/sigh2; 1/sigh2*ones(T-1,1)])*Hrho;
errh_out = 1;
while errh_out> 10^(-3)
        % E-step
    tmps2 = sum(tmp.^2,2)./exp(ht);
    Eilam = (nu+n)./(tmps2+nu);
    s2 = sum(tmp.^2,2).*Eilam;
        % M-step
    htt = ht;
    errh_in = 1;
    while errh_in> 10^(-3)
        eht = exp(htt);
        sieht = s2./eht;
        fh = -n/2 + .5*sieht;
        Gh = .5*sieht;
        Kh = HiSH + sparse(1:T,1:T,Gh);
        newht = Kh\(fh+Gh.*htt);
        errh_in = max(abs(newht-htt));
        htt = newht;
    end
    errh_out = max(abs(ht-htt));
    ht = htt;
end
    % compute negative Hessian
s2 = sum(tmp.^2,2);
Gh = (nu+n)/(2*nu)*(s2.*exp(ht))./((exp(ht)+s2/nu).^2);
Kh = HiSH + sparse(1:T,1:T,Gh);
CKh = chol(Kh,'lower');

%% evaluate the importance weights
c_pri = -T/2*log(2*pi*sigh2) -.5*log(1/(1-rho^2));
c_IS = -T/2*log(2*pi) + sum(log(diag(CKh)));
pri_den = @(x) c_pri -.5*x'*HiSH*x;
IS_den = @(x) c_IS -.5*(x-ht)'*Kh*(x-ht);
store_llike = zeros(R,1);
for i=1:R
    hc = ht + CKh'\randn(T,1);
    store_llike(i) = deny_h(Utld,hc,CSig,nu) + pri_den(hc) - IS_den(hc);
end
maxllike = max(store_llike);
intlike = log(mean(exp(store_llike-maxllike))) + maxllike;
end

function llike = deny_h(U,h,CSig,nu)
[T,n] = size(U);
c = T*(gammaln((nu+n)/2) - gammaln(nu/2) - n/2*log(nu*pi))...
    - T*sum(log(diag(CSig))) - n/2*sum(h);
tmp = U/CSig';
s2 = sum(tmp.^2,2);
llike = c -(nu+n)/2*sum(log(1+(s2./exp(h))/nu));
end
