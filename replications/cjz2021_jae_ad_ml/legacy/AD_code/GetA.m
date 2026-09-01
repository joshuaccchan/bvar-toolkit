function C=GetA(df,M)
A=zeros(2,2,M);At=zeros(2,2,M);

C1=sqrt(chi2rnd(df,1,M));
C2=sqrt(chi2rnd(df-1,1,M));
Z=randn(1,M);

A(1,1,:)=C1;
A(2,2,:)=C2;

At(1,1,:)=C1;
At(2,2,:)=C2;
A(2,1,:)=Z;
At(1,2,:)=Z;

C=tmult(A,At);
