function [irfMin,irfMax] = genBounds(q0,ivar,vma,Aineq,bineq,Aeq,...
    beq,LB,UB,optimOptions)

irfMin = zeros(1,length(ivar));
irfMax = irfMin;

for ii = 1:length(ivar) % For each variable of interest

    % Minimise IRF subject to sign and equality restrictions.
    [~,irfMin(ii)] = fmincon(@(q) genIRF(q,vma(ivar(ii),:),1),...
        q0,Aineq,bineq,Aeq,beq,LB,UB,@(q) qcon(q),optimOptions);
    % Maximise IRF subject to sign and equality restrictions.
    [~,irfMax(ii)] = fmincon(@(q) genIRF(q,vma(ivar(ii),:),-1),...
        q0,Aineq,bineq,Aeq,beq,LB,UB,@(q) qcon(q),optimOptions);

end