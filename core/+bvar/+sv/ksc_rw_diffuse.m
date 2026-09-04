% bvar.sv.ksc_rw_diffuse - KSC auxiliary-mixture sampler for the log-volatility
% path, random-walk state equation with a DIFFUSE-style initial condition
% h_1 ~ N(0,Vh) (Vh a variance, supplied by the caller):
%   ystar_t = h_t + eps_t,  eps_t approximated by the Kim-Shephard-Chib (1998)
%             7-component normal mixture,
%   h_t = h_{t-1} + v_t,    v_t ~ N(0,omega2h),   h_1 ~ N(0,Vh).
% Returns the new path h AND the mixture indicators S.
% Consumes rand(T,1) then randn(T,1) - one of each per call.
%
% Extracted 2026-09-01 (step 4, SV/prior core). Canonical body:
% chan_jeliazkov2009_statespace/legacy/sp_code/SVRW.m, verbatim (single copy).
% Function renamed SVRW -> ksc_rw_diffuse; output list written [h,S] instead of
% the legacy space-separated [h S] (syntax only, no numeric effect); nothing
% parameterized. NEVER merge with bvar.sv.ksc_rw_h0: different initial condition
% (h_1 ~ N(0,Vh) here vs h_1 ~ N(h0,sig) there), hence different draws.
%
% See Chan, J.C.C. (2013). Moving Average Stochastic Volatility Models
%     with Application to Inflation Forecast, Journal of Econometrics, 176(2): 162-172
% (c) 2012, Joshua Chan. Email: joshuacc.chan@gmail.com

function [h,S] = ksc_rw_diffuse(ystar,h,omega2h,Vh)

T = length(h);
%% parameters for the Gaussian mixture
pi = [0.0073 .10556 .00002 .04395 .34001 .24566 .2575];
mui = [-10.12999 -3.97281 -8.56686 2.77786 .61942 1.79518 -1.08819] - 1.2704;
sigma2i = [5.79596 2.61369 5.17950 .16735 .64009 .34023 1.26261];
sigmai = sqrt(sigma2i);

%% sample S from a 7-point distrete distribution
temprand = rand(T,1);
q = repmat(pi,T,1).*normpdf(repmat(ystar,1,7),repmat(h,1,7) ...
    +repmat(mui,T,1),repmat(sigmai,T,1));
q = q./repmat(sum(q,2),1,7);
S = 7 - sum(repmat(temprand,1,7)<cumsum(q,2),2)+1;

%% sample h
H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);
invOmegah = spdiags([1/Vh; 1/omega2h*ones(T-1,1)],0,T,T);
d = mui(S)'; invSigystar = spdiags(1./sigma2i(S)',0,T,T);
Kh = H'*invOmegah*H + invSigystar;
Ch = chol(Kh,'lower');              % so that Ch*Ch' = Kh
hhat = Kh\(invSigystar*(ystar-d));
h = hhat + Ch'\randn(T,1);          % note the transpose
