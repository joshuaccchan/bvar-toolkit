% This script generates data from a reduced-form VAR with SV (AR state equations)

clc; clear all; 
% Check whether the restriction matrix satisfy the identifying assumptions first
for i=[1:6]
    eval(['load',' SignRestr',num2str(i),'.mat;'])
    check_result(i,1:2) = [Func_check_restrictions(Restr), sum(sum(abs(Restr)))];
end
check_result

clc; clear all; T=200; p=5; 
n=10; m=5; i=1
eval(['load', ' SignRestr', num2str(i), '.mat', ' Restr']);
rng('default'); rng(1000);
[Y0,Y,L,A] = Func_gendata_VARSign(T,n,p,Restr); data = Y; L
%always, n-th/1-th, if 1-th is bigger then set it -1 
[Rineq_newalg, Rineq_RWZ, Ridx] = Func_Gen_Rineq(i,L,m,n); save('RankRestr6.mat','Rineq_newalg','Rineq_RWZ','Ridx');
eval(['save', ' DataT200_', num2str(i), '.mat', ' data']);

clc; clear all; T=200; p=5; 
n=30; m=5; i=2
eval(['load', ' SignRestr', num2str(i), '.mat', ' Restr']);
rng('default'); rng(1000);
[Y0,Y,L,A] = Func_gendata_VARSign(T,n,p,Restr); data = Y; L
%always, n-th/1-th, if 1-th is bigger then set it -1 
[Rineq_newalg, Rineq_RWZ, Ridx] = Func_Gen_Rineq(i,L,m,n); save('RankRestr7.mat','Rineq_newalg','Rineq_RWZ','Ridx');
eval(['save', ' DataT200_', num2str(i), '.mat', ' data']);

clc; clear all; T=200; p=5; 
n=50; m=5; i=3
eval(['load', ' SignRestr', num2str(i), '.mat', ' Restr']);
rng('default'); rng(1000);
[Y0,Y,L,A] = Func_gendata_VARSign(T,n,p,Restr); data = Y; L
%always, n-th/1-th, if 1-th is bigger then set it -1 
[Rineq_newalg, Rineq_RWZ, Ridx] = Func_Gen_Rineq(i,L,m,n); save('RankRestr8.mat','Rineq_newalg','Rineq_RWZ','Ridx');
eval(['save', ' DataT200_', num2str(i), '.mat', ' data']);

clc; clear all; T=200; p=5; 
n=10; m=8; i=4
eval(['load', ' SignRestr', num2str(i), '.mat', ' Restr']);
rng('default'); rng(1000);
[Y0,Y,L,A] = Func_gendata_VARSign(T,n,p,Restr); data = Y; L
%always, n-th/1-th, if 1-th is bigger then set it -1 
[Rineq_newalg, Rineq_RWZ, Ridx] = Func_Gen_Rineq(i,L,m,n); save('RankRestr10.mat','Rineq_newalg','Rineq_RWZ','Ridx');
eval(['save', ' DataT200_', num2str(i), '.mat', ' data']);

clc; clear all; T=200; p=5; 
n=30; m=8; i=5
eval(['load', ' SignRestr', num2str(i), '.mat', ' Restr']);
rng('default'); rng(1000);
[Y0,Y,L,A] = Func_gendata_VARSign(T,n,p,Restr); data = Y; L
%always, n-th/1-th, if 1-th is bigger then set it -1 
[Rineq_newalg, Rineq_RWZ, Ridx] = Func_Gen_Rineq(i,L,m,n); save('RankRestr11.mat','Rineq_newalg','Rineq_RWZ','Ridx');
eval(['save', ' DataT200_', num2str(i), '.mat', ' data']);

clc; clear all; T=200; p=5; 
n=50; m=8; i=6
eval(['load', ' SignRestr', num2str(i), '.mat', ' Restr']);
rng('default'); rng(100);
[Y0,Y,L,A] = Func_gendata_VARSign(T,n,p,Restr); data = Y; L
%always, n-th/1-th, if 1-th is bigger then set it -1 
[Rineq_newalg, Rineq_RWZ, Ridx] = Func_Gen_Rineq(i,L,m,n); save('RankRestr12.mat','Rineq_newalg','Rineq_RWZ','Ridx');
eval(['save', ' DataT200_', num2str(i), '.mat', ' data']);
