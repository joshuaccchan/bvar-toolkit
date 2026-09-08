function [irfMin,irfMax] = numericalBounds(restr,phi,q0,opt,optimOptions)
% Function uses a constrained optimisation routine to find the upper and
% lower bounds of identified set for impulse responses.
% Inputs are:
% - restr: structure representing restrictions
% - phi: structure containing reduced-form VAR parameters
% - q0: value of q satisfying identifying restrictions
% - opt: structure containing model information and options
% - optimOptions: structure containing optimisation option

F = restr.F;
S = restr.S;
vma = phi.vma;
H = opt.H;
ivar = opt.ivar;

% Set equality constraint for optimisation (Aeq*q = beq).
Aeq = F;
beq = zeros(size(F,1),1);

% Set inequality constraints for optimisation (Aineq*q <= b).
Aineq = -S;
bineq = zeros(size(S,1),1);

irfMin = zeros(H+1,length(opt.ivar));
irfMax = irfMin;
LB = [];
UB = [];

parfor hh = 1:H+1 % For each horizon
    
    [irfMin(hh,:),irfMax(hh,:)] = genBounds(q0,ivar,vma(:,:,hh),...
        Aineq,bineq,Aeq,beq,LB,UB,optimOptions)

end

irfMax = -irfMax; % Max obtained by minimising negative of objective

% Take optima over different initial values.
irfMin = min(irfMin,[],3);
irfMax = max(irfMax,[],3);

end