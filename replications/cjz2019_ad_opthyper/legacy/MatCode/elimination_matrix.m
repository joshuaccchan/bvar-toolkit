% Support script for Chan, Jacobi and Zhu (2019)
%
% See: 
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019
function [m0] = elimination_matrix(n)

   %src = 1;
    %tgt = 1;
    %for j = 1:n
     %   for i = 1:n           
      %      if i >= j
       %         m0(tgt, src) = 1;
        %        tgt = tgt + 1;
         %   end
          %  src = src + 1;
        %end
    %end
    ind1=repmat((1:n)',n,1);
    ind2=repmat(1:n,n,1);    
    ind2=ind2(:);
    ind3=find(ind1>=ind2);
    
    m0 = sparse((1:0.5*n*(n+1))',ind3, 1);
    
    
