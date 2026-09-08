function [hpdlb,hpdub] = highestPosteriorDensity(rPost,opt)
% Compute highest posterior density interval with credibility level alpha.
% Inputs:
% - rPost: posterior draws
% - opt: structure containing options

aalpha = opt.aalpha; % Credibility level
gridLength = opt.gridLength; % Number of point on discrete grid

[K,H,ni] = size(rPost);

hpdlb = zeros(H,ni);
hpdub = hpdlb;

for ii = 1:ni % For each variable

    % Storage matrices
    Cent = zeros(H,1);
    Rad = zeros(H,1);

    for hh = 1:H % For each horizon

        % Construct grid over which to search for bounds.
        r = linspace(min(rPost(:,hh,ii)),max(rPost(:,hh,ii)),gridLength);
        gridr = kron(r,ones(K,1));

        % d(r,phi) = |r-r(phi)|
        d = abs(gridr-kron(rPost(:,hh,ii),ones(1,gridLength)));
        zhat = quantile(d,aalpha,1);
        [Rad(hh),ind] = min(zhat,[],2);
        Cent(hh) = gridr(1,ind);

    end

    % Approximate highest posterior density region is an interval centered 
    % at Cent with radius Rad.
    hpdlb(:,ii) = Cent - Rad;
    hpdub(:,ii) = Cent + Rad;
    
end
