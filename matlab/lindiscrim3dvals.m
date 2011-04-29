function discrimvals = lindiscrim3dvals(data,bnd)
%discrimvals = lindiscrim3dvals(data,bnd)
%  computes the linear discriminant values that a
%  categorizer uses to determine a response.  An example
%  of a response rule would be:
%
%    If a discriminant value is less than zero, then
%    respond category A; otherwise, respond category B.
%
%  Parameters:
%    data row format:  [x y z]
%    bnd format:  [a1 a2 a3 bias]
%          where a1*x + a2*y a3*z + bias = 0.

% Created by Leola A. Alfonso-Reese / 31-March-94
% Copyright (c) 1995
% $Revisions$
%   Date           Modification and Name
%   ----           ---------------------

discrimvals = data*bnd(1:3)' + bnd(4);

