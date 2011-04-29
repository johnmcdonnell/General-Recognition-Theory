function [discrimvals] = lindiscrim1dvals(data,bnd)
%[discrimvals] = lindiscrim1dvals(data,bnd)
%  computes the linear discriminant values that a
%  categorizer uses to determine a response.  An example
%  of a response rule would be:
%
%    If a discriminant value is less than zero, then
%    respond category A; otherwise, respond category B.
%
%  Parameters:
%    data row format:  [x]
%    bnd format:  [a bias]
%          where a*x + bias = 0.

% Created by Leola A. Alfonso-Reese & Brian Spiering/ 1-Nov-02
% Copyright (c) 2002
% $Revisions$
%   Date           Modification and Name
%   ----           ---------------------


discrimvals = data*bnd(1)' + bnd(2);

