% Support function for estimating the hybrid TVP-VAR in Chan (2022)
%
% See:
% Chan, J.C.C. (2022). Large Hybrid Time-Varying Parameter VARs, Journal of
% Business and Economic Statistics, forthcoming
    
    % initialize alp0, beta0, Sigalp, Sigbeta
beta0 = zeros(k_beta,1);
alp0 = zeros(k_alp,1);
count1 = 0;
U = zeros(T,n);
for ii=1:n         
        ki = n*p+1+ii-1;
        Xi = [X2 -Y(:,1:ii-1)];
        theta0 = (Xi'*Xi)\(Xi'*Y(:,ii));
        beta0((ii-1)*k_beta/n+1:ii*k_beta/n) = theta0(1:k_beta/n);
        alp0(count1+1:count1+ii-1) = theta0(k_beta/n+1:ki);        
        U(:,ii) = Y(:,ii)-Xi*theta0;        
        count1 = count1 + ii-1;        
end
Sigalp = Vsigalp;
Sigbeta = Vsigbeta;

    % initialize p0
p0 = .5*ones(n,2);  % fixed at 0.5 so the prior for gam is symmetric and does not affect the ML value

    % initialize h0 and Sigh    
h0 = mean(log(U.^2))';
Sigh = Sh0;

    % initialize h
h = repmat(h0',T,1);
for ii=1:n
    if size(Sigh,2) == 1
        h(:,ii) = sample_SVRW(log(U(:,ii).^2+1e-4),h(:,ii),Sigh(ii),h0(ii));   
    else
        h(:,ii) = sample_SVRW(log(U(:,ii).^2+1e-4),h(:,ii),Sigh(ii,ii),h0(ii));   
    end
end