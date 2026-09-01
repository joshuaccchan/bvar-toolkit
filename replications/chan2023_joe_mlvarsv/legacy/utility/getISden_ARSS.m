% This function obtains the importance sampling density for h
% using the parameterization
% h_t = a_t + rho h_{t-1} + u_t, u_t ~ N(0,b_t)
%    
% See:
% Chan, J.C.C. (2022). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, forthcoming.


function [h_hat,Kh_hat,r_hat,m_hat,v_hat] = getISden_ARSS(store_h)
T = size(store_h,2);
r_hat = fminbnd(@(x) -concen_like_h(store_h,x),-.99,.99);
Hr = speye(T) -r_hat*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
[~,m_hat,v_hat] = concen_like_h(store_h,r_hat);
h_hat = Hr\m_hat;
Kh_hat = Hr'*sparse(1:T,1:T,1./v_hat)*Hr;

end

function [concen_like,m_hat,v_hat] = concen_like_h(store_h,rho)
[R,T] = size(store_h);
Hrho = speye(T) -rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
Hrhoh = Hrho*store_h';
m_hat = mean(Hrhoh,2);
E = Hrhoh - repmat(m_hat,1,R);
v_hat = mean(E.^2,2);

concen_like = -T*R/2*log(2*pi) -R/2*sum(log(v_hat)) ...
    -.5*sum(sum(E.^2./repmat(v_hat,1,R)));
end

