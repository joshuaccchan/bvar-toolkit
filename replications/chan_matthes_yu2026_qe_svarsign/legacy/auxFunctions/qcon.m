function [c,ceq,DC,DCeq] = qcon(q)
% Generate non-linear equality constraint on the first column of Q
% (q'q = 1), along with analytical gradient of the constraint.
% Inputs:
% q: first column of orthonormal matrix, Q.

c = []; % No nonlinear inequality constraint
DC = []; % No associated gradient

ceq = q'*q - 1; % Equality constraint in form C(Q) = 0
DCeq = 2*q; % Gradient of equality constraint with respect to q

end