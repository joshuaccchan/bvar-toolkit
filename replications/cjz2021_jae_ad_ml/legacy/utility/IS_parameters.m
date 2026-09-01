

function [a_bar,Da_bar,CDa_bar,...
    nusig2_bar,Ssig2_bar,...
    nuomega2_bar,Somega2_bar]=IS_parameters(store_a,store_Sig,store_Omega, R)
n=size(store_Sig.v,2); r=size(store_Omega.v,2);
a_bar.v = mean(store_a.v)';
a_bar.d=mean(store_a.d,3);
Da_bar = d_Cov2(store_a);
CDa_bar = Cholasky(Da_bar,[],[],[],[]);
dim=size(a_bar.d,2);
nusig2_bar.v = zeros(n,1);nusig2_bar.d=zeros(n,dim);
Ssig2_bar.v = zeros(n,1);Ssig2_bar.d=zeros(n,dim);
for i=1:n
    m =gamfit(store_Sig.v(:,i));
    d=d_gamMLE(m(1),m(2),store_Sig.v(:,i)',squeeze(store_Sig.d(i,:,:))',R);
    nusig2_bar.v(i)=m(1);nusig2_bar.d(i,:)=d(1,:);
    Ssig2_bar.v(i)=m(2);Ssig2_bar.d(i,:)=d(2,:);
end

nuomega2_bar.v = zeros(r,1); Somega2_bar.v = zeros(r,1);
nuomega2_bar.d=zeros(r,dim);Somega2_bar.d=zeros(r,dim);
for jj=1:r
    m = gamfit(store_Omega.v(:,jj),"Gamma");
    d=d_gamMLE(m(1),m(2),store_Omega.v(:,jj)',squeeze(store_Omega.d(jj,:,:))',R);
    nuomega2_bar.v(jj) = m(1);
    Somega2_bar.v(jj) =m(2);  
    nuomega2_bar.d(jj,:)=d(1,:);
    Somega2_bar.d(jj,:)=d(2,:);
end
