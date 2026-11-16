% bvar.util.tnormrnd - N draws from the normal truncated to (a, b), by inverse
% transform.
%
%   t = bvar.util.tnormrnd(mu, sigma2, a, b)
%   t = bvar.util.tnormrnd(mu, sigma2, a, b, N)
%
%   mu, sigma2 : mean and VARIANCE of the untruncated normal, either scalars,
%                which are expanded to N identical values, or length-N vectors
%   a, b       : truncation bounds, scalar or length-N
%   N          : number of draws (default length(mu))
%   t          : N x 1 vector of draws
%
% Far in the upper tail, where both normal cdf values round to 1, a draw comes from the
% mirrored lower tail with the same uniform.
%
% See:
% Chan, J.C.C. (2023). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, 235(2): 1419-1446.

function t = tnormrnd(mu, sigma2, a, b, N)

if ( nargin < 4  )
    error( 'wrong # of arguments' );
end

K = length( mu );

if ( nargin < 5  )
     N = K;
end

if ( ( K ~= N ) | ( length( sigma2 ) ~= N ) ) & ( ( K ~= 1 ) )
    error( 'dimensions of mu and sigma must equal N')
end
    
if K == 1
    mu = ones( N, 1 ) * mu;
    sigma2 = ones( N, 1 ) * sigma2;
end

sigma = sqrt( sigma2 );
u = rand(N,1);
p1 = normcdf( ( a - mu ) ./ sigma );
p2 = normcdf( ( b - mu ) ./ sigma );
C = norminv( p1 + ( p2 - p1 ) .* u );
if any(~isfinite(C(:)))   % far in the upper tail p1 and p2 round to 1: -t is in (-b, -a)
    q1 = normcdf( ( mu - b ) ./ sigma );
    q2 = normcdf( ( mu - a ) ./ sigma );
    Cm = -norminv( q2 + ( q1 - q2 ) .* u );
    C(~isfinite(C)) = Cm(~isfinite(C));
end
t = mu + sigma .* C;
end