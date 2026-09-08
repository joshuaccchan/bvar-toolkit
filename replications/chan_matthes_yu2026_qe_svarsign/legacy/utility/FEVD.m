
A = Btilde;
m = 8;
[k,n] = size(A);
A_comp = [A(2:end,:)';[eye((p-1)*n) zeros((p-1)*n,n)]];
J = [speye(n), sparse(n,n*(p-1))];

MSPE = zeros(n,m,3); % 3 different horizons: 1, 5, 20

for jj=1:20
    if jj == 1
        A_powerj = A_comp;
    else 
        A_powerj = A_comp*A_comp;
    end
    Phij = J*A_powerj*J';
    Thetaj = Phij*L;
end
