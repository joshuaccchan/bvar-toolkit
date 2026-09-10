% bvar.sv.csv_armh - one accept-reject Metropolis-Hastings sweep for the common
% stochastic volatility log-volatility path h, marginal of the mixture indicators:
% Newton-Raphson mode-finding, a Gaussian proposal at the mode, then MH correction.
%
%   [h,is_accept] = bvar.sv.csv_armh(s2, rho, sigh2, h, n)
%   [h,is_accept] = bvar.sv.csv_armh(..., is_ForcedAccept, ht_start)
%   [h,is_accept] = bvar.sv.csv_armh(..., 'c_reject', 3)
%
%   s2              : T x 1, sum over the n series of squared (orthogonalized,
%                     lambda-scaled) errors at each t
%   rho, sigh2      : AR(1) coefficient and innovation variance of h
%   h, n            : current log-volatility path (T x 1), and number of series
%   is_ForcedAccept : take the proposal regardless of the MH ratio; default false
%   ht_start        : starting point of the mode search; default h
%   is_accept       : 1 if the MH proposal was accepted
%
%   'c_reject'      c in the envelope target <= c * proposal; default 3. Larger c
%                   costs more proposals per sweep and does not move the target
%   'MaxIterMode', 'MaxIterAR'   caps on the two loops; defaults 500 and 1000
%
% An independence sampler, so a bad starting h can leave the chain stuck there
% silently; is_accept stays 0 when it does, so check it rather than the path.
%
% rng consumption: randn(T,1) and one rand per accept-reject proposal until one is
% accepted, so the count is data-dependent, then one rand for the MH step.
%
% Provenance, the legacy copies this stands in for, and the design notes:
% tests/variant_map.md. Equivalence: tests/unit/test_csv_armh.m.
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics, 38(1), 68-79.
% Chan, J.C.C. (2023). Comparing stochastic volatility specifications for large
% Bayesian VARs, Journal of Econometrics, 235(2), 1419-1446.

function [h,is_accept] = csv_armh(s2,rho,sigh2,h,n,varargin)
    % positional optionals first, then any name-value pairs
is_ForcedAccept = false;
ht_start = [];
iv = 1;
if numel(varargin) >= 1 && ~(ischar(varargin{1}) || isstring(varargin{1}))
    is_ForcedAccept = varargin{1};
    iv = 2;
    if numel(varargin) >= 2 && ~(ischar(varargin{2}) || isstring(varargin{2}))
        ht_start = varargin{2};
        iv = 3;
    end
end
c_reject = 3; maxit_mode = 500; maxit_ar = 1000;
tol_mode = 1e-3;                        % NOT an option: see the header
while iv <= numel(varargin)
    if iv == numel(varargin)
        error('bvar:sv:csv_armh:badOption', ...
            'option ''%s'' has no value', name_of(varargin{iv}));
    end
    switch lower(name_of(varargin{iv}))
        case 'c_reject',    c_reject   = varargin{iv+1};
        case 'maxitermode', maxit_mode = varargin{iv+1};
        case 'maxiterar',   maxit_ar   = varargin{iv+1};
        otherwise, error('bvar:sv:csv_armh:badOption', ...
                'unknown option ''%s''', name_of(varargin{iv}));
    end
    iv = iv + 2;
end
if ~(isscalar(is_ForcedAccept) && (islogical(is_ForcedAccept) || isnumeric(is_ForcedAccept)))
    error('bvar:sv:csv_armh:badArgument', ...
        ['is_ForcedAccept must be a logical scalar; got %s. To supply ht_start ' ...
         'while leaving the flag alone, call csv_armh(...,false,ht_start).'], ...
        name_of(is_ForcedAccept));
end
if ~(isscalar(c_reject) && isnumeric(c_reject) && isreal(c_reject) ...
        && isfinite(c_reject) && c_reject > 0)
    error('bvar:sv:csv_armh:badArgument', ...
        ['c_reject must be a finite positive real scalar; got %s. Zero freezes ' ...
         'the chain and a negative value silently acts as its magnitude.'], ...
        name_of(c_reject));
end
maxit_mode = check_cap(maxit_mode, 'MaxIterMode');
maxit_ar   = check_cap(maxit_ar,   'MaxIterAR');
if isempty(ht_start), ht_start = h; end

is_accept = 0;
T = size(s2,1);
Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
HiSH = Hrho'*sparse(1:T,1:T,[(1-rho^2)/sigh2; 1/sigh2*ones(T-1,1)])*Hrho;
errh = Inf; ht = ht_start; it_mode = 0;
while ~(errh <= tol_mode)             % NaN fails this, so it hits the cap
    it_mode = it_mode + 1;
    if it_mode > maxit_mode
        error('bvar:sv:csv_armh:modeNotConverged', ...
            ['mode search did not converge to %g in %d iterations ' ...
             '(last max|dh| = %g)'], tol_mode, maxit_mode, errh);
    end
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
% AR-step
hstar = ht;
logc = -.5*hstar'*HiSH*hstar - n/2*sum(hstar) - .5*exp(-hstar)'*s2 + log(c_reject);
flag = 0; it_ar = 0;
while flag == 0
    it_ar = it_ar + 1;
    if it_ar > maxit_ar
        error('bvar:sv:csv_armh:arNotAccepted', ...
            ['no accept-reject proposal accepted in %d draws at c_reject = ' ...
             '%g; lower c_reject, or check s2, rho and sigh2'], ...
            maxit_ar, c_reject);
    end
    hc = ht + CKh'\randn(T,1);
    alpARc = -.5*hc'*HiSH*hc - n/2*sum(hc) - .5*exp(-hc)'*s2 ...
        + .5*(hc-ht)'*Kh*(hc-ht) - logc;
    if alpARc > log(rand)
        flag = 1;
    end
end
% MH-step
alpAR = -.5*h'*HiSH*h - n/2*sum(h) -.5*exp(-h)'*s2 + .5*(h-ht)'*Kh*(h-ht) - logc;
if alpAR < 0
    alpMH = 1;
elseif alpARc < 0
    alpMH = - alpAR;
else
    alpMH = alpARc - alpAR;
end
if alpMH > log(rand) || is_ForcedAccept
    h = hc;
    is_accept = 1;
end
end

function s = name_of(v)
if ischar(v) || isstring(v), s = char(v); else, s = mat2str(v); end
end

function v = check_cap(v, nm)
if ~(isscalar(v) && isnumeric(v) && isreal(v) && isfinite(v) && v >= 1 && v == fix(v))
    error('bvar:sv:csv_armh:badArgument', ...
        ['%s must be a finite positive integer; got %s. Inf would restore the ' ...
         'unbounded loop this argument exists to prevent.'], nm, name_of(v));
end
end
