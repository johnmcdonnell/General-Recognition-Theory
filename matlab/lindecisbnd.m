function [a,b] = lindecisbnd(K,meanA,meanB)
%[a,b] = lindecisbnd(K,meanA,meanB)
%  finds the ideal linear decision boundary coefficients
%  such that  h(x,y,...) = a*(x,y,...)' + b.  The boundary,
%  [a b], is used to categorize points coming from one of
%  two distributions with the same covariance matrix and
%  different means.
%
%  Parameters:
%    K  is the common covariance matrix
%    meanA is a column vector containing the mean of category A
%    meanB is a column vector containing the mean of category B
%
%  Recommended use of this rule:
%    if h(x) < 0, then put x in category A;
%    otherwise, put x in category B.

% Created by Leola A. Alfonso-Reese / 24-February-94
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%   2-August-94     Normalize decision bound vector a
%                   by David H. Brainard
%   13-Sept-02      Add 'not full rank' error message and check
%                   by Alfonso-Reese & Spiering

% Check that covariance matrix is full rank.
if rank(K) ~= size(K,1)
	error('Covariance matrix not full rank.');
end

Kinv=inv(K);
udiff = meanB-meanA;
a = (udiff)'*Kinv;
normConstant = norm(a);
a = a/normConstant;
b = -1/2*udiff'*Kinv*(meanB+meanA)/normConstant;


