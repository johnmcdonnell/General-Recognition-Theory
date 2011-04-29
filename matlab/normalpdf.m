function fxvec = normalpdf(xvec,mu,sigma)
%fxvec = normalpdf(xvec,mu,sigma)
% returns the normal probability distribution function value
% with mean mu and standard deviation sigma for each number
% in xvec
%
%  Parameters:
%    xvec is a vector of real numbers
%    mu is a real number (optional:  default value is 0)
%    sigma is a non-negative real number (optional:
%      default value is 1)

% Created by Leola A. Alfonso-Reese / 22-January-95
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

fxvec = 1/(sigma*sqrt(2*pi))*exp(-.5*((xvec-mu)./sigma).^2);

