% bvar.util.igrnd - draws from the inverse-gamma distribution in the (nu, S)
% parameterization used throughout this toolkit and the papers:
%
%       x ~ IG(nu, S)   with density proportional to x^(-nu-1) exp(-S/x),
%       mean S/(nu-1) for nu > 1.
%
%   x = bvar.util.igrnd(nu, S)
%
% nu and S must be positive, and may be scalars or conformable arrays; the
% result has the size of the expansion, one independent draw per element.
%
% rng consumption: implemented as 1./gamrnd(nu, 1./S), the same expression the
% samplers write inline, so a seeded call here advances the random stream
% identically (asserted in tests/unit/test_igrnd.m).
%
% The conjugate variance step of a Gaussian model - with prior IG(nu0, S0)
% and residuals e - is
%
%       sig2 = bvar.util.igrnd(nu0 + T/2, S0 + sum(e.^2)/2);
%
% which is the update repeated in the state-equation variance draws of the SV
% and TVP samplers.
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
