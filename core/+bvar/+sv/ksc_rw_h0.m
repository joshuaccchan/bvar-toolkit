% bvar.sv.ksc_rw_h0 - KSC auxiliary-mixture sampler for the log-volatility path,
% random-walk state equation with KNOWN initial value h0 (h_1 ~ N(h0,sig)):
%   ystar_t = h_t + eps_t,  eps_t = log of a chi^2_1 draw, approximated by the
%             Kim-Shephard-Chib (1998) 7-component normal mixture,
%   h_t = h_{t-1} + v_t,    v_t ~ N(0,sig),   with sig a VARIANCE.
%
%   h = bvar.sv.ksc_rw_h0(Ystar, h, sig, h0)
%
%   Ystar is T x 1; h on input is the current path and on output the new draw;
%   h0 is the known log-volatility at time 0. Only the path is returned: the
%   mixture indicators are drawn internally and discarded.
%
% Distinct from ksc_rw_diffuse: different initial condition, hence different
% draws.
%
% When all seven mixture densities of a period underflow, its indicator probabilities are
% computed in logs.
%
% See:
% Kim, S., Shephard, N. and Chib, S. (1998). Stochastic Volatility: Likelihood
% Inference and Comparison with ARCH Models, Review of Economic Studies, 65(3):
% 361-393.
% Chan, J.C.C. and Jeliazkov, I. (2009). Efficient Simulation and Integrated
% Likelihood Estimation in State Space Models, International Journal of
% Mathematical Modelling and Numerical Optimisation, 1(1/2): 101-120.
% Chan, J.C.C. (2020). Large Bayesian Vector Autoregressions. In: P. Fuleky (Eds),
% Macroeconomic Forecasting in the Era of Big Data, 95-125, Springer, Cham

function h = ksc_rw_h0(Ystar,h,sig,h0)
T = length(h);
%% normal mixture
pi = [0.0073 .10556 .00002 .04395 .34001 .24566 .2575];
mi = [-10.12999 -3.97281 -8.56686 2.77786 .61942 1.79518 -1.08819] - 1.2704;  %% means already adjusted!! %%
sigi = [5.79596 2.61369 5.17950 .16735 .64009 .34023 1.26261];
sqrtsigi = sqrt(sigi);

%% sample S from a 7-point discrete distribution
temprand = rand(T,1);
q = repmat(pi,T,1).*normpdf(repmat(Ystar,1,7),repmat(h,1,7)+repmat(mi,T,1), repmat(sqrtsigi,T,1));
bad = ~(sum(q,2) > 0);    % all seven densities underflow: weigh them in logs
if any(bad)
    lq = log(pi) - log(sqrtsigi) - .5*((Ystar(bad) - h(bad) - mi)./sqrtsigi).^2;
    q(bad,:) = exp(lq - max(lq, [], 2));
end
q = q./repmat(sum(q,2),1,7);
S = 7 - sum(repmat(temprand,1,7)<cumsum(q,2),2)+1;

%% sample h
% Given the indicators S: Ystar = h + d + e with e ~ N(0,Omega), where d_t and
% Omega_tt are the mean and variance of mixture component S_t; and
% Hh*h = [h0;0] + v with v ~ N(0,sig*I_T), so h has prior mean alph = Hh\[h0;0]
% and precision Kh.
Hh =  speye(T) - spdiags(ones(T-1,1),-1,T,T);
invSh = spdiags(1/sig*ones(T,1),0,T,T);
dconst = mi(S)'; invOmega = spdiags(1./sigi(S)',0,T,T);
alph = Hh\[h0;sparse(T-1,1)];
Kh = Hh'*invSh*Hh;
Ph = Kh + invOmega;
Ch = chol(Ph,'lower');              % so that Ch*Ch' = Ph
hhat = (Ch')\(Ch\(Kh*alph + invOmega*(Ystar-dconst)));
h = hhat + Ch'\randn(T,1);
end
