function xvec = normalcdfinv(Fxvec,mu,sigma)
%xvec = normalcdfinv(Fxvec,mu,sigma)
%  returns the inverse of the normal cumulative distribution
%  function with mean mu and standard deviation sigma for
%  each probability in Fxvec.
%
%  Parameters:
%    Fxvec is a vector of probabilities
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

% Check that the Fxvec values are all probabilities
notprobs = find(Fxvec < 0 | Fxvec > 1);
if ~isempty(notprobs)
	error('Not all Fxvec values are probabilities:  between 0 and 1');
end

% Initialize xvec
xvec = zeros(size(Fxvec));

% Find infinity values
zero_indices = find(Fxvec == 0);
if ~isempty(zero_indices)
	xvec(zero_indices) = -inf*ones(size(zero_indices));
end

one_indices = find(Fxvec == 1);
if ~isempty(one_indices)
	xvec(one_indices) = inf*ones(size(one_indices));
end

% Compute remaining inv cdfs (use relation between erf and normal cdf)
% remember:  x = z*sigma + mu
%            erf(z/sqrt(2)) = 2*prob-1
%            erfinv((z/sqrt(2)) = z/sqrt(2)
%      So,   z/sqrt(2) * (sqrt(2)) * sigma + mu = x
others = find(Fxvec > 0  &  Fxvec < 1);
if ~isempty(others)
    xvec(others) = erfinv(2*Fxvec(others)-1).*sqrt(2)*sigma + mu;
end

