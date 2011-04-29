function [A,b,c] = sphdecisbnd(mean1,mean2)
%[A,b,c] = sphdecisbnd(mean1,mean2)
%  finds the spherical decision boundary coefficients
%  such that h(x,y,...) = (x,y,...)A(x,y,...) + b(x,y,...) + c.
%  The boundary, [A,b,c], is used to categorize points coming
%  from one of the two distributions.  The sphere surrounds mean1
%  and bisects the two means.

%  Parameters:
%    meanA is a column vector containing the mean of category A
%    meanB is a column vector containing the mean of category B
%
%  Recommended use of this rule:
%   if h(x) < 0, then put x in catagory A;
%   otherwise, put x in category B.

% Created by Leola A. Alfonso-Reese / 7-Sept-99
% Copyright (c) 1999
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------


r = sqrt(sum((mean2-mean1).^2))/2;
A = diag([1 1 1]);
b = -2*mean1;
c = sum(mean1.^2)-r^2;

