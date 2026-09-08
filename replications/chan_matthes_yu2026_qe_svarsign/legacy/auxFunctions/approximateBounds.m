function [rmin,rmax] = approximateBounds(restr,phi,q0,K,opt)
% Approximate bounds of the identified set via simulation using a Gibbs
% sampler to draw values of q satisfying the identifying retrictions.
% Inputs:
% - restr: structure containing information about restrictions
% - phi: structure containing reduced-form VAR parameters
% - q0: initial value of q satisfying sign restrictions in
% linear subspace.
% - K: change of basis matrix
% - opt: structure containing model information and options

H = opt.H;
ivar = opt.ivar;
L = opt.L;
vma = phi.vma;
cumIR = opt.cumIR;

% Draw q satisfying the restrictions L times.
q = drawqGibbs(restr.Sbar,size(restr.F,1),K,q0,L);

% Compute impulse responses for each draw of q.
etaDraw = zeros(H+1,length(ivar),L);

for hh = 1:H+1 % For each horizon

    % Multiply by each draw of q to obtain impulse response.
    etaDraw(hh,:,:) = vma(ivar,:,hh)*q;

end
 
if ~isempty(cumIR)
    
    % Compute cumulative impulse responses for relevant variables.
    etaDraw(:,cumIR,:,:) = cumsum(etaDraw(:,cumIR,:,:),1);

end

% Compute minimum and maximum impulse response over draws of Q.
rmin = min(etaDraw,[],3);
rmax = max(etaDraw,[],3);

end

