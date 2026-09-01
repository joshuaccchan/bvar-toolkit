%% Implementation of the Devroye (2014) algorithm for sampling from 
% the generalized inverse Gaussian (GIG) distribution
%
% function X = gigrnd(p, a, b, sampleSize)
%
% The generalized inverse Gaussian (GIG) distribution is a continuous
% probability distribution with probability density function:
%
% p(x | p,a,b) = (a/b)^(p/2)/2/besselk(p,sqrt(a*b))*x^(p-1)*exp(-(a*x + b/x)/2)
%
% Parameters:
%   p \in Real, a > 0, b > 0
%
% See Wikipedia page for properties:
% https://en.wikipedia.org/wiki/Generalized_inverse_Gaussian_distribution 
%
% This is an implementation of the Devroye (2014) algorithm for GIG sampling.
% 
% Returns:
%   X      = random variates [sampleSize x 1] from the GIG(p, a, b)
%  
% References:
% L. Devroye
% Random variate generation for the generalized inverse Gaussian distribution 
% Statistics and Computing, Vol. 24, pp. 239-246, 2014.
%
% (c) Copyright Enes Makalic and Daniel F. Schmidt, 2015
function X = gigrnd(P, a, b, sampleSize)

%% Setup -- we sample from the two parameter version of the GIG(alpha,omega)
lambda = P;
omega = sqrt(a*b);

swap = false;
if(lambda < 0)
     lambda = lambda * -1;
     swap = true;
end

alpha = sqrt(omega^2 + lambda^2) - lambda;

%% Find t
x = -psi(1, alpha, lambda);
if((x >= 0.5) && (x <= 2))
    t = 1;
elseif(x > 2)
    t = sqrt(2 / (alpha + lambda));
elseif(x < 0.5)
    t = log(4/(alpha + 2*lambda));     
end

%% Find s
x = -psi(-1, alpha, lambda);
if((x >= 0.5) && (x <= 2))
    s = 1;
elseif(x > 2)
    s = sqrt(4/(alpha*cosh(1) + lambda));
elseif(x < 0.5)
    s = min(1/lambda, log(1 + 1/alpha + sqrt(1/alpha^2+2/alpha)));
end

%% Generation
eta = -psi(t, alpha, lambda);
zeta = -dpsi(t, alpha, lambda);
theta = -psi(-s, alpha, lambda);
xi = dpsi(-s, alpha, lambda);
p = 1/xi;
r = 1/zeta;
td = t - r*eta;
sd = s - p*theta;
q = td + sd;

X = zeros(sampleSize, 1);
for sample = 1:sampleSize
    done = false;
    while(~done)
        U = rand(1); 
        V = rand(1); 
        W = rand(1);
        if(U < (q / (p + q + r)))
            X(sample) = -sd + q*V;
        elseif(U < ((q + r) / (p + q + r)))
            X(sample) = td - r*log(V);
        else
            X(sample) = -sd + p*log(V);
        end

        %% Are we done?
        f1 = exp(-eta - zeta*(X(sample)-t));
        f2 = exp(-theta + xi*(X(sample)+s));
        if((W*g(X(sample), sd, td, f1, f2)) <= exp(psi(X(sample), alpha, lambda)))
            done = true;
        end
    end
end

%% Transform X back to the three parameter GIG(p,a,b)
X = exp(X) * (lambda / omega + sqrt(1 + (lambda/omega)^2));
if(swap)
    X = 1 ./ X;
end
X = X ./ sqrt(a/b);

end

function f = psi(x, alpha, lambda)
    f = -alpha*(cosh(x) - 1) - lambda*(exp(x) - x - 1);
end

function f = dpsi(x, alpha, lambda)
    f = -alpha*sinh(x) - lambda*(exp(x) - 1);
end

function f = g(x, sd, td, f1, f2)

a = 0;
b = 0;
c = 0;
if((x >= -sd) && (x <= td))
    a = 1;
elseif(x > td)
    b = f1;
elseif(x < -sd)
    c = f2;   
end

f = a + b + c;

end