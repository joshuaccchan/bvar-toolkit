function [Y0,Y,B0,B,h,phi,sig2,A] = gendata_SVARSV(T,n,p)
T0 = 8; % number of initial obs

    % generate parameters 
phi = .95*ones(n,1);
sig2 = .05*ones(n,1);
B0_idl = nonzeros(tril(reshape(1:n^2,n,n),-1)');
B0_idu = nonzeros(triu(reshape(1:n^2,n,n),1)');
B0_id = [B0_idl;B0_idu];
B0 = diag(.5 + 1*rand(n,1));
B0(B0_id) = randn(length(B0_id),1)/sqrt(n);
% B0 = eye(n);
% B0(B0_idu) = -.8;
% B0(B0_idl) = .8;

    % generate reduced-form VAR coefficients
A = zeros(n,n*p+1);
A(:,1) = -10 + 20*rand(n,1);
for ii=1:p
    if ii == 1
        Ai = -.2 + .4*rand(n,n);
        Ai(1:n+1:end) = .5*rand(n,1);    
    else
        Ai = .1*randn(n,n)/ii^2; 
    end       
    A(:,(ii-1)*n+2:ii*n+1) = Ai;
end
B = B0*A;  % structural-form parameters

   % generate SV
h = zeros(T+T0,n);
for ii=1:n
    h(:,ii) = genSV(0,phi(ii),sig2(ii),T+T0);
end
    
    % generate the data
Uy = (exp(h/2).*randn(T+T0,n))/(B0');
Y = zeros(T+T0,n);
for tt=1:T+T0
    if tt <= p 
        xt = [1 reshape(Y(tt-1:-1:1,:)',n*(tt-1),1)'];
        Y(tt,:) = xt*A(:,1:n*(tt-1)+1)' + Uy(tt,:);
    else
        xt = [1 reshape(Y(tt-1:-1:tt-p,:)',n*p,1)'];
        Y(tt,:) = xt*A' + Uy(tt,:);
    end    
end
Y0 = Y(1:T0,:);
Y = Y(T0+1:end,:);
h = h(T0+1:end,:);
end