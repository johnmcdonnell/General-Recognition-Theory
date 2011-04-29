function Fxvec = normalcdf(xvec,mu,sigma)
%Fxvec = normalcdf(xvec,mu,sigma)
%  returns the normal cumulative distribution function with
%  mean mu and standard deviation sigma for each number in xvec.
%
%  Parameters:
%    xvec is a vector of real numbers
%    mu is a real number (optional:  default value is 0)
%    sigma is a non-negative real number (optional:
%      default value is 1)

% Created by Leola A. Alfonso-Reese / 8-August-95
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------


if nargin < 2;
	mu = 0;
elseif isempty(mu)
  mu = 0;
end

if nargin < 3, 
	sigma = 1;
elseif isempty(sigma)
  sigma = 1;
end

if sigma <= 0
	error('Sigma must be greater than zero.');
end

% Compute CDFs (use relation between erf and normal cdf)
Fxvec = 0.5*(1 + erf( (xvec-mu)./(sigma*sqrt(2)) ) );

% Check that values in Fxvec are not greater than 1 (due to round-off
% errors)
exceedone = find(Fxvec > 1);
if ~isempty(exceedone)
	Fxvec(exceedone) = ones(size(exceedone));
end

