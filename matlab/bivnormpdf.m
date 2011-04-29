function f = bivnormpdf(XYs,mu,covmat)
% f = bivnormpdf(XYs,mu,covmat)
% returns the likelihood of the XY coordinates
% of a bivariate normal distribution specified
% with mean uvec and covariance matrix covmat
% If uvec and covmat are not specified, then the
% default values mu = [0 0]' and covmat = [1 0;0 1]
% are used.
%
% Parameters:
%  XYs is a matrix of (x,y) coordinates where each
%           row is one x-y pair
%

% by Leola
% 7/11/96

if nargin < 3
	covmat = [1 0;0 1];
end

if nargin < 2
	mu = [0 0]';
end

% Check the covariance matrix
if (covmat(1,1) == 0) | (covmat(2,2) == 0) | (covmat(1,2)~=covmat(2,1))
	error('The covariance matrix is not valid.');
else
  s1 = sqrt(covmat(1,1));
	s2 = sqrt(covmat(2,2));
	p = covmat(1,2)/(s1*s2);
end

X = XYs(:,1);
Y = XYs(:,2);
mu1 = mu(1)*ones(length(X),1);
mu2 = mu(2)*ones(length(Y),1);

% Initialize f
f = zeros(size(XYs,1),1);

c1 = 1/(2*pi*s1*s2*sqrt(1-p^2));
c2 = -1/(2*(1-p^2));
c3 = -2*p;

Z1 = (X-mu1)/s1;
Z2 = (Y-mu2)/s2;

f = c1 * exp( c2 * (Z1.^2  + c3*Z1.*Z2 + Z2.^2) );

