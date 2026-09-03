% bvt.ml.intlike_csv_ma - integrated likelihood of the BVAR-CSV-MA model at
% (A,Sig,psi,rho,sigh2): the common log-volatility path h is integrated out
% by importance sampling with R draws around the Newton-Raphson mode of its
% conditional density given the MA(1)-transformed residuals. The first
% observation carries the (1+psi^2) initialization variance both in the mode
% search and inside the deny_h likelihood evaluation. Consumes R*T randn calls.
%
% Extracted 2026-09-02 (step 8, Kronecker family pass). Canonical source:
%   chan2020_jbes_kronecker/legacy/intlike_BVAR_CSV_MA.m  (body verbatim;
%   renamed intlike_BVAR_CSV_MA -> intlike_csv_ma). Called by the BVAR-CSV-MA
% marginal likelihood with R = 5000 (legacy ml_BVAR_CSV_MA.m lines 10-11).
% Verified by tests/unit/test_kron_intlike.m (bitwise, seeded).
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics, 38(1), 68-79.

function [intlike,store_llike] = intlike_csv_ma(shortY,X,A,Sig,psi,rho,sigh2,R)
[T,n] = size(shortY);
Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);

% obtain the mode and negative Hessian of the conditional density of h
CSig = chol(Sig,'lower');
U = shortY - X*A;
Utld = Hpsi\U;
tmp = (Utld/CSig');
s2 = sum(tmp.^2,2);
s2(1) = s2(1)/(1+psi^2);
HiSH = Hrho'*sparse(1:T,1:T,[(1-rho^2)/sigh2; 1/sigh2*ones(T-1,1)])*Hrho;
ht = log(mean(s2))*ones(T,1);
errh = 1;
while errh> 10^(-3)
    eht = exp(ht);
    sieht = s2./eht;
    fh = -n/2 + .5*sieht;
    Gh = .5*sieht;
    Kh = HiSH + sparse(1:T,1:T,Gh);
    newht = Kh\(fh+Gh.*ht);
    errh = max(abs(newht-ht));
    ht = newht;
end
CKh = chol(Kh,'lower');

%% evaluate the importance weights
c_pri = -T/2*log(2*pi*sigh2) -.5*log(1/(1-rho^2));
c_IS = -T/2*log(2*pi) + sum(log(diag(CKh)));
pri_den = @(x) c_pri -.5*x'*HiSH*x;
IS_den = @(x) c_IS -.5*(x-ht)'*Kh*(x-ht);
store_llike = zeros(R,1);
for i=1:R
    hc = ht + CKh'\randn(T,1);
    store_llike(i) = deny_h(U,hc,CSig,psi) + pri_den(hc) - IS_den(hc);
end
maxllike = max(store_llike);
intlike = log(mean(exp(store_llike-maxllike))) + maxllike;
end

function llike = deny_h(U,h,CSig,psi)
[T, n] = size(U);
Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
Utld = Hpsi\U;

c = -T*n/2*log(2*pi) - T*sum(log(diag(CSig))) - n/2*log(1+psi^2) -n/2*sum(h);
tmp = (Utld/CSig');
s2 = sum(tmp.^2,2);
llike = c -.5*(s2(1)/((1+psi^2)*exp(h(1))) + sum(s2(2:end)./exp(h(2:end))));
end
