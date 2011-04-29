function Y = gensample(u,K,n,catnum)
%Y = gensample(u,K,n,catnum)
%  generates 'n' normally distributed data points.
%
%  Parameters:
%    u is a column vector containing the population mean
%    K is a column vector containing the population covariance
%      matrix
%    n is the number of data points
%    catnum is an integer category label
%
%  Row format of output:  [catnum x y ...]

% Created by Leola A. Alfonso-Reese & David H. Brainard/ 16-April-94
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------

% Step 1: Generate X, a matrix containing N normally
% distributed mean 0 covariance I random numbers.
dim = length(u);
X = randn(dim,n);

% Step 2: Transform X to have covariance K.  Recall that
% if we form A = A*X, then the covariance of A will be
% given by K = A*A'.  So we need to find A such that
% A*A' = K.  The cholesky decomposition does this.
A = chol(K)';
%K_check = A*A';

% Then form Y as simply Y = (A*X)'.
Y = (A*X)';

% Check the covariance of Y to see that it is close to
% what we expected.
%Sample_K = cov(Y);

% Update values in Y to reflect the desired mean of the distribution.
% Also, add a row to Y to indicate the category from which the points
% were generated.
o = ones(n,dim);
uvec = o*diag(u);
Y = [catnum*ones(n,1) Y+uvec];

