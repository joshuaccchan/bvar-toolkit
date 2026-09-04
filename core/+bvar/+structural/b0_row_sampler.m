% bvar.structural.b0_row_sampler - one full row-wise pass of the posterior draw
% of the impact matrix B0 in the order-invariant SVAR-SV model. For each row
% ii = 1:n it forms the row's conditional Gaussian precision Kbi from the
% residuals U weighted by exp(-h(:,ii)), rotates onto an orthonormal basis
% aligned with the orthogonal complement of the other rows (null-space
% construction), draws the first rotation coordinate from the bimodal
% two-component approximation bvar.util.anormrnd and the remaining n-1
% coordinates as zetaj_hat + randn/sqrt(T), maps back, and fixes the sign of
% the ii-th element to be positive.
%
%   B0 = bvar.structural.b0_row_sampler(U, h, B0, B00, VB0)
%
%   U   : T x n residual matrix Y - X*A (computed by the CALLER, verbatim
%         legacy position - estimation U = Y-X*A, forecast U = Yt-Xt*A)
%   h   : T x n log-volatilities
%   B0  : current n x n impact matrix (rows updated in place, in order)
%   B00 : n x n prior mean of B0   (legacy Hyper.B0,  = eye(n) in the paper)
%   VB0 : n x n prior variances    (legacy Hyper.VB0, = ones(n) in the paper)
%
% rng consumption, per row: one rand then one randn (inside anormrnd), then
% n-1 further randn - rows in order ii = 1:n.
%
% Extracted 2026-09-02 (step 7, OISV family pass). Canonical source (body
% verbatim, including comments): chan_koop_yu2024_jbes_oisv/legacy/SVARSV_MH.m
% lines 49-72 (the inline "sammple B0" loop). Also canonicalizes
% forecast_SVARSV_MH.m lines 43-66, which are textually identical modulo the
% Y/X/T -> Yt/Xt/Tt renaming (U and T enter only through the arguments here).
% Edits made, in full: wrapped as a function with [T,n] = size(U) replacing the
% workspace T,n (identical integers); Hyper.B0/Hyper.VB0 renamed B00/VB0; the
% unqualified anormrnd call now bvar.util.anormrnd (code-identical to the legacy
% utility copy). Everything else byte-verbatim. Draw-for-draw equivalence:
% tests/unit/test_oisv_equivalence.m.
%
% See:
% Chan, J.C.C., Koop, G. and Yu, X. (2024). Large Order-Invariant Bayesian
% VARs with Stochastic Volatility, Journal of Business and Economic
% Statistics, 42(2): 825-837.

function B0 = b0_row_sampler(U,h,B0,B00,VB0)
[T,n] = size(U);
for ii=1:n
    EiOhi = U'*sparse(1:T,1:T,exp(-h(:,ii)));
    Kbi = sparse(1:n,1:n,1./VB0(ii,:)) + EiOhi*U;
    mui = Kbi\(B00(ii,:)./VB0(ii,:))';
    Ci = chol(Kbi,'lower')/sqrt(T);
    Gam_mi = B0([1:ii-1 ii+1:end],:)';
    Gam_miperp = null(Gam_mi');

    V = zeros(n,n); zeta = zeros(n,1);
    for jj=1:n
        if jj==1
            v1 = Ci\Gam_miperp; v1 = v1/norm(v1);
            V = [v1 null(v1')];
            zetaj_hat = mui'*(Ci*v1);
            zeta(1) = bvar.util.anormrnd(zetaj_hat,1/T);
        else
            zetaj_hat = mui'*(Ci*V(:,jj));
            zeta(jj) = zetaj_hat + 1/sqrt(T)*randn;
        end
    end
    phii = (Ci')\sum(V.*repmat(zeta',n,1),2);
    % B0(ii,:) = phii;
    B0(ii,:) = phii*sign(phii(ii)); % fix the sign of the i-th element to be positive
end
end
