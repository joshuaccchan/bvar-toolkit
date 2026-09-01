% this function simulate a standard SV model with mean muh, AR(1)
% coefficient phih and varianc omegah2

function h = genSV(muh,phih,omegah2,T)
h = zeros(T,1); 
h(1) = muh + sqrt(omegah2/(1-phih^2))*rand; 
for t=2:T                
    h(t) = muh + phih*(h(t-1)-muh) + sqrt(omegah2)*randn;
end
end