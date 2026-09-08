% This script applies three different algorithms to check whether the
% identified set is nonempty at provided draws of reduced-form parameters.

%% Use Algorithm 4.1 to determine whether identified set is nonempty.

Qempty1 = zeros(opt.phiDraws,1);

tic
for kk = 1:opt.phiDraws
    
    % Put draw of reduced-form parameters into structure.
    phi.B = B(:,kk);
    phi.Sigmatr = Sigmatr(:,:,kk);
    phi.Sigmatrinv = Sigmatrinv(:,:,kk);
    phi.vma = vma(:,:,:,kk);
    
    [~,~,~,~,~,Qempty1(kk)] = chebyCheck(restr,phi,opt);
    
end
AlgoRunTime1 = toc;

%% Use rejection sampler to determine whether identified set is nonempty.

Qempty2 = zeros(opt.phiDraws,1);

tic
for kk = 1:opt.phiDraws
    
    % Put draw of reduced-form parameters into structure.
    phi.B = B(:,kk);
    phi.Sigmatr = Sigmatr(:,:,kk);
    phi.Sigmatrinv = Sigmatrinv(:,:,kk);
    phi.vma = vma(:,:,:,kk);
    
    [~,~,~,Qempty2(kk)] = drawq(restr,phi,opt);
    
end
AlgoRunTime2 = toc;

%% Use algorithm from GKV to determine whether identified set is nonempty.

Qempty3 = zeros(opt.phiDraws,1);

tic
for kk = 1:opt.phiDraws
    
    % Put draw of reduced-form parameters into structure.
    phi.B = B(:,kk);
    phi.Sigmatr = Sigmatr(:,:,kk);
    phi.Sigmatrinv = Sigmatrinv(:,:,kk);
    phi.vma = vma(:,:,:,kk);
    
    Qempty3(kk) = checkEmptyIS_GKV(restr,phi,opt);
    
end
AlgoRunTime3 = toc;