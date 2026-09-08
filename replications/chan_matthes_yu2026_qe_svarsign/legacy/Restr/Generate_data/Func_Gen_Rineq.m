function [Rineq_newalg, Rineq_RWZ, Ridx] = Func_Gen_Rineq(i,L,m,n)
S1 = [1,2,5,6,10]; S2 =[3,4,7,8,11,12];
if ismember(i,S1)
Ridx = [1, 2, 3];
elseif ismember(i,S2)
Ridx = [1, 2, 3];
end
Rineq = zeros(m,n);
L1n = L([1,n],Ridx); 
Rineq(Ridx,1) = -sign(L1n(1,:) - L1n(2,:)); 
Rineq(Ridx,n) = sign(L1n(1,:) - L1n(2,:)); 
Rineq_newalg = Rineq; Rineq = Rineq(Ridx,:); Rineq_RWZ = Rineq;
end
