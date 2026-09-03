% bvt.ml.kron_bvar - log marginal likelihood of the standard (homoskedastic,
% natural-conjugate) BVAR of Chan (2020, JBES): fully analytic - the ML
% identity llike + lpri - lpost is evaluated exactly at the posterior mean
% pair A = Ahat, Sig = Shat/(T+nu0), with both NIW ordinates in closed form.
% No rng is consumed.
%
% Extracted 2026-09-02 (step 8, Kronecker family pass). Canonical source:
%   chan2020_jbes_kronecker/legacy/BVAR.m lines 36-47 (the inline cp_ml
%   block; model 1 has no separate ml_* script). Body verbatim; lniwpdf ->
%   bvt.ml.lniwpdf (code-identical). No evaluation-point defect (clean bill,
%   step-8 audit): every piece is evaluated at the same (A,Sig).
%
%   [ML, out] = bvt.ml.kron_bvar(shortY, X, pri, est)
%
%   pri: A0, VA0 (k x 1), nu0, S0        [replication preset, cited there]
%   est: Ahat, Shat, KA                  [run_all model 1 output]
%   out: llike, A, Sig                   [the evaluation point used]
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics, 38(1), 68-79.

function [ML, out] = kron_bvar(shortY, X, pri, est)
[T, n] = size(shortY);
k = size(X, 2);
A0 = pri.A0; VA0 = pri.VA0; nu0 = pri.nu0; S0 = pri.S0;
Ahat = est.Ahat; Shat = est.Shat; KA = est.KA;

    % [BVAR.m lines 37-44, verbatim]
A = Ahat;
Sig = Shat/(T+nu0);
CSig = chol(Sig,'lower');
tmp = (shortY-X*A)/CSig';
s2 = sum(tmp.^2,2);
llike =  -T*n/2*log(2*pi) - T*sum(log(diag(CSig))) -.5*sum(s2);
ML = llike + bvt.ml.lniwpdf(A,Sig,A0,sparse(1:k,1:k,1./VA0),nu0,S0) ...
    - bvt.ml.lniwpdf(A,Sig,Ahat,KA,nu0+T,Shat);

out = struct('llike', llike, 'A', A, 'Sig', Sig);
end
