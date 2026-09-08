function [Y0,Y,L,A] = Func_gendata_VARSign(T,n,p,Restr)
T0 = 8; % number of initial obs
flag = true;
while flag % return only if the series does not explode
    % generate parameters 
L_idl = nonzeros(tril(reshape(1:n^2,n,n),-1)');
L_idu = nonzeros(triu(reshape(1:n^2,n,n),1)');
L_id = [L_idl;L_idu];
L = diag(.5 + 1*rand(n,1));
L(L_id) = randn(length(L_id),1); %L(L_id) = randn(length(L_id),1)*0.1;

Restr_full = [Restr, zeros(n,n-length(Restr(1,:)))];
signL2 = sign(L).*Restr_full; signL2(find(signL2==0))=1; L = L.*signL2;

    % generate reduced-form VAR coefficients
A = zeros(n,n*p+1);
A(:,1) = -1 + 2*rand(n,1);
for ii=1:p
    if ii == 1
        Ai = -.2 + .4*rand(n,n);
        Ai(1:n+1:end) = .5*rand(n,1);    
    else
        Ai = .1*randn(n,n)/ii^2; 
    end       
    A(:,(ii-1)*n+2:ii*n+1) = Ai;
end
   
    % generate the data
Uy = randn(T+T0,n)*L';  %note: this should be L'
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
if max(abs(max(Y)-min(Y))) < 1e3
        flag = false;
end
end