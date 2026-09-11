% bvar.structural.reduced_form - map stored structural-form draws of the
% asymmetric-conjugate VAR to their reduced form, one draw at a time.
%
%   [store_Btilde,store_Sigtilde] = bvar.structural.reduced_form(store_alp, ...
%       store_beta, store_Sig)
%
%   store_alp  : nsim x n(n-1)/2 free elements of the unit-lower-triangular A,
%                stacked equation by equation
%   store_beta : nsim x (n^2 p + n) structural coefficients, equation by equation
%   store_Sig  : nsim x n structural innovation variances (diagonal)
%   store_Btilde   : nsim x (n^2 p + n), the reduced-form coefficients stacked
%                    the same way
%   store_Sigtilde : nsim x n x n reduced-form covariances
%
% The structural form is A y_t = B x_t + eps_t with A unit lower triangular and
% Var(eps_t) = diag(sig). The reduced form follows by inversion:
% Sigtilde = A^{-1} diag(sig) A^{-T} and Btilde = (A^{-1} B')'. Both are formed
% with backslash rather than an explicit inverse.
%
% See:
% Chan, J.C.C. (2022). Asymmetric Conjugate Priors for Large Bayesian VARs,
% Quantitative Economics, 13(3): 1145-1169.

function [store_Btilde,store_Sigtilde] = reduced_form(store_alp,store_beta,store_Sig)
[nsim,n] = size(store_Sig);
k_beta = size(store_beta,2);
p = (k_beta/n-1)/n;
A_id = nonzeros(tril(reshape(1:n^2,n,n),-1)');
A = eye(n);
store_Btilde = zeros(nsim,n^2*p+n,1);
store_Sigtilde = zeros(nsim,n,n);
    % compute reduced-form parameters
for isim = 1:nsim
    alp = store_alp(isim,:)';
    beta = store_beta(isim,:)';
    sig = store_Sig(isim,:)';

        % trasnform the parameters into reduced-form
    A(A_id) = alp;
    Sig = (A\sparse(1:n,1:n,sig))/A';
    Btilde = (A\(reshape(beta,n*p+1,n)'))';

    store_Btilde(isim,:) = Btilde(:); % stack by equations
    store_Sigtilde(isim,:,:) = Sig;
end
end
