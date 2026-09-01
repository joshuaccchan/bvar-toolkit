% Support script for Chan, Jacobi and Zhu (2019)
%
% See: 
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019
function L=GettheL(S,n)
    L.v=chol(S.inv)';
    L.d=d_cholasky(L.v,S.invd);
    L.dt=d_trans(L.v,L.d);
    L.inv=(L.v\speye(n));
    L.tinv=(L.v'\speye(n));
end 