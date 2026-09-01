% See:
% Chan, J.C.C., Jacobi, L. and Zhu, D. (2018). How Sensitive Are VAR 
% Forecasts to Prior Hyperparameters? An Automated Sensitivity Analysis,
% CAMA Working Paper 25/2018.

function yf=forecastVAR(ind1, ind2, ind3, N,beta,sigma)
global T n  Y
C.v=chol(sigma.v)';
C.d=d_cholasky(C.v,sigma.d);
C.v=C.v';

Z=randn(N,n);
yf.v=zeros(N,n);
yf.vm=zeros(N,n);
yf.dm=zeros(N,n,3);
yf.d=zeros(N,n,3);
%% First step
%obtain mean sensitivities
yf.vm(1,:)=beta.v(1,:)+Y(T,:)*beta.v(2:4,:)+Y(T-1,:)*beta.v(5:7,:);
yf.dm(1,:,:)=beta.d(ind1,:)+kron(eye(n),Y(T,:))*beta.d(ind2,:)+kron(eye(n),Y(T-1,:))*beta.d(ind3,:);
%forecast sensitivity
yf.v(1,:)=yf.vm(1,:)+Z(1,:)*C.v;
yf.d(1,:,:)=squeeze(yf.dm(1,:,:))+kron(Z(1,:),eye(n))*C.d;

%% Second time step
yf.vm(2,:)=beta.v(1,:)+yf.v(1,:)*beta.v(2:4,:)+Y(T,:)*beta.v(5:7,:);
yf.dm(2,:,:)=beta.d(ind1,:)+kron(eye(n),yf.v(1,:))*beta.d(ind2,:)+kron(eye(n),Y(T,:))*beta.d(ind3,:)...
             +beta.v(2:4,:)'*squeeze(yf.d(1,:,:));
yf.v(2,:)=yf.vm(2,:)+Z(2,:)*C.v;
yf.d(2,:,:)=squeeze(yf.dm(2,:,:))+kron(Z(2,:),eye(n))*C.d;


%% Further into the future
for i=3:N
   yf.vm(i,:)=beta.v(1,:)+yf.v(i-1,:)*beta.v(2:4,:)+yf.v(i-2,:)*beta.v(5:7,:);
   
    yf.dm(i,:,:)=beta.d(ind1,:)+kron(eye(n),yf.v(i-1,:))*beta.d(ind2,:)+kron(eye(n),yf.v(i-2,:))*beta.d(ind3,:)...
                 +beta.v(2:4,:)'*squeeze(yf.d(i-1,:,:))+beta.v(5:7,:)'*squeeze(yf.d(i-2,:,:));

    
    yf.v(i,:)=yf.vm(i,:)+Z(i,:)*C.v;
    yf.d(i,:,:)=squeeze(yf.dm(i,:,:))+kron(Z(i,:),eye(n))*C.d;
end  

