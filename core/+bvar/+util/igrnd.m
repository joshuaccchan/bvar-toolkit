% bvar.util.igrnd - draws from the inverse-gamma distribution in the (nu, S)
% parameterization used throughout this toolkit and the papers:
%
%       x ~ IG(nu, S)   with density proportional to x^(-nu-1) exp(-S/x),
%       mean S/(nu-1) for nu > 1.
%
%   x = bvar.util.igrnd(nu, S)
%
% nu and S may be scalars or conformable arrays; the result has the size of
% the expansion, one independent draw per element. Implemented as
% 1./gamrnd(nu, 1./S), which is the expression the legacy samplers write
% inline, so a seeded call here consumes the random stream identically.
%
% The conjugate variance step of a Gaussian model - with prior IG(nu0, S0)
% and residuals e - is
%
%       sig2 = bvar.util.igrnd(nu0 + T/2, S0 + sum(e.^2)/2);
%
% which is the update repeated in the state-equation variance draws of the SV
% and TVP samplers.
%
% NOTE for maintainers: the functions under core/ extracted from the published
% packages keep their inline `1/gamrnd(...)` calls verbatim, so they remain
% line-by-line diffable against replications/*/legacy/. This helper is for new
% code and is not retrofitted into them; it draws identically under the same
% seed (asserted in tests/unit/test_igrnd.m).
%
% Requires the Statistics and Machine Learning Toolbox (gamrnd).

function x = igrnd(nu, S)
if nargin < 2
    error('bvar:util:igrnd:nargin', 'igrnd needs both nu and S');
end
if any(nu(:) <= 0) || any(S(:) <= 0)
    error('bvar:util:igrnd:badParam', 'nu and S must be positive');
end
x = 1./gamrnd(nu, 1./S);
end
