function [credlb,credub] = credibleRegion(rmin,rmax,opt)
% Compute robustified credible region with credibility level alpha.
% Inputs:
% - rmin: posterior draws of lower bound
% - rmax: posterior draws of upper bound
% - opt: structure containing options

aalpha = opt.aalpha; % Credibility level
gridLength = opt.gridLength; % Number of point on discrete grid

[K,H,ni] = size(rmin);

% Storage.
credlb = zeros(H,ni);
credub = credlb;

for ii = 1:ni % For each variable
    
    Cent = zeros(H,1);
    Rad = zeros(H,1);

    for hh = 1:H % For each horizon

        % Construct discrete grid.
        r = linspace(min(rmin(:,hh,ii)),max(rmax(:,hh,ii)),gridLength);
        gridr = kron(r,ones(K,1));

        % d(r,phi) = max{|r-l(phi)|,|r-u(phi)|}.
        d = max(abs(gridr-kron(rmin(:,hh,ii),ones(1,gridLength))),...
            abs(gridr-kron(rmax(:,hh,ii),ones(1,gridLength))));
        zhat = quantile(d,aalpha,1);
        [Rad(hh),ind] = min(zhat,[],2);
        Cent(hh) = gridr(1,ind);

    end

    % Approximate robustified credible region is an interval centered at 
    % Cent with radius Rad.
    credlb(:,ii) = Cent - Rad;
    credub(:,ii) = Cent + Rad;

end
    
end
