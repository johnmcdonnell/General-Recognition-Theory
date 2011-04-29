function [] = plot2dlinbnd(bnd,plotstr,xyaxes)
%plot2dlinbnd(bnd,plotstr,xyaxes)
%  plots a two-dimensional linear boundary on the
%  current figure.
%
%  Parameters:
%    bnd format:  [a1 a2 b]
%       where a1*x + a2*y + b = 0 is a boundary
%    plotstr is a character string indicating plot
%       symbols and colors (see MATLAB's built-in
%       plot function)
%    xyaxes format:  [minX maxX minY maxY] (optional)

% Created by Leola A. Alfonso-Reese / 7-August-95
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------


if nargin < 3
  xyaxes = [];
end

if ~isempty(xyaxes)
	axis(xyaxes);
	axis(axis);
	hold on;
end

% Display the boundary
[pt1,pt2] = getbndendpts(bnd);
plot([pt1(1) pt2(1)],[pt1(2) pt2(2)],plotstr);

