function [A,b,c] = quaddecisbnd(KA,KB,meanA,meanB)
%[A,b,c] = quaddecisbnd(KA,KB,meanA,meanB)
%  finds the ideal quadratic decision boundary coefficients
%  such that h(x,y,...) = (x,y,...)A(x,y,...) + b(x,y,...) + c.
%  The boundary, [A,b,c], is used to categorize points coming
%  from one of the two distributions having different covariance
%  matrices and different means.
%
%  Parameters:
%    KA  is the covariance matrix for category A
%    KB  is the covariance matrix for category B
%    meanA is a column vector containing the mean of category A
%    meanB is a column vector containing the mean of category B
%
%  Recommended use of this rule:
%   if h(x) < 0, then put x in catagory A;
%   otherwise, put x in category B.

% Created by Leola A. Alfonso-Reese / 25-March-96
% Copyright (c) 1996
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------


KAinv=inv(KA);
KBinv=inv(KB);
A = -.5*(KBinv-KAinv);
b = (meanB'*KBinv - meanA'*KAinv)';
c = .5*( -log(det(KB)/det(KA)) - meanB'*KBinv*meanB + meanA'*KAinv*meanA);

