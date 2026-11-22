% bvar.util.soe_reduced_form - the reduced form of a VAR with block exogeneity from its two blocks.
% Written for this toolkit.
%
%   [A, Sig] = bvar.util.soe_reduced_form(Bf, Bd, Sf, Om)
%
% The VAR has n = nf + nd variables, the nf foreign ones first, and p lags. The foreign block is
% the VAR yf_t' = xf_t'*Bf + uf_t', on the intercept and the foreign lags in xf_t. The domestic
% block is the regression yd_t' = (x_t', yf_t')*Bd + e_t', on all lags and the current foreign
% variables, with e_t independent of uf_t.
%
%   Bf  : (1 + nf*p) x nf coefficients of the foreign block, intercept first
%   Bd  : (k + nf) x nd coefficients of the domestic block: on x_t = (1, y_{t-1}', ...,
%         y_{t-p}')', k = 1 + n*p, then on the current foreign variables
%   Sf  : nf x nf error covariance matrix of the foreign block
%   Om  : nd x nd error covariance matrix of the domestic block
%   A   : k x n reduced-form coefficients, in the layout of bvar.util.build_lags, with zeros on
%         the domestic lags in the foreign equations
%   Sig : n x n reduced-form error covariance matrix
%
% See:
% Zha, T. (1999). Block Recursion and Structural Vector Autoregressions. Journal of
% Econometrics, 90(2): 291-316.

function [A, Sig] = soe_reduced_form(Bf, Bd, Sf, Om)
nf = size(Bf, 2); nd = size(Bd, 2); n = nf + nd;
k = size(Bd, 1) - nf;
p = (k - 1)/n;
if p < 1 || p ~= fix(p) || size(Bf, 1) ~= 1 + nf*p
    error('bvar:util:soe_reduced_form:badDims', ['Bd must have 1 + n*p + nf rows and Bf ' ...
        '1 + nf*p, for n = %d variables and an integer p'], n);
end
if ~isequal(size(Sf), [nf nf]) || ~isequal(size(Om), [nd nd])
    error('bvar:util:soe_reduced_form:badDims', 'Sf must be %d x %d and Om %d x %d', nf, nf, ...
        nd, nd);
end
iff = [1, reshape(1 + (1:nf)' + (0:p-1)*n, 1, [])];   % the intercept and the foreign lags in x_t
G0 = Bd(k+1:end, :)';                    % the current foreign coefficients
A = zeros(k, n);
A(iff, 1:nf) = Bf;                       % zeros for the domestic lags
A(:, nf+1:end) = Bd(1:k, :) + A(:, 1:nf)*G0';
Sig = [Sf, Sf*G0'; G0*Sf, G0*Sf*G0' + Om];
end
