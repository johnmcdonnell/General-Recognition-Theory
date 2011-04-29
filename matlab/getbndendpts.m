function [pt1,pt2] = getbndendpts(bnd,axisvals)
%[pt1,pt2] = getbndendpts(bnd,axisvals)
%  returns two points on the boundary that will extend
%  across a two dimensional graph given the axis values
%  for a graph.  If no axis values are specified, then
%  this routine gets the axis of the current figure.
%
%  Parameters:
%    bnd format: [a1 a2 b], where a1*x+ a2*y + b = 0
%    axisvals is [minX maxX minY maxY] (optional)

% Created by Leola A. Alfonso-Reese / 18-January-95
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------

if nargin < 2
	a = axis;
elseif isempty(axisvals)
	a = axis;
else
	a = axisvals;
end

if bnd(2)==0
	x = -bnd(3)/bnd(1);
	pt1 = [x a(3)];
	pt2 = [x a(4)];
else
  x = a(1);
	y = (-bnd(1)*x-bnd(3))/bnd(2);
	pt1 = [x y];
  x = a(2);
	y = (-bnd(1)*x-bnd(3))/bnd(2);
	pt2 = [x y];
end

