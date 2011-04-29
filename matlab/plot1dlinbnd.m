function [] = plot1dlinbnd(bnd,plotstr,xaxes)
%plot1dlinbnd(bnd,plotstr,xaxes)
%  plots a one-dimensional linear boundary on the
%  current figure.
%
%  Parameters:
%    bnd format:  [a b]
%       where a*x + b = 0 is a boundary
%    plotstr is a character string indicating plot
%       symbols and colors (see MATLAB's built-in
%       plot function)
%    xaxes format:  [minX maxX] (optional)

% Created by Leola A. Alfonso-Reese / 10-March-05
% Copyright (c) 2005
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------


if nargin < 3
  xaxes = [];
end

if ~isempty(xaxes)
	curaxes = axis;
	newaxes = [xaxes(1) xaxes(2) curaxes(3:4)];
	axis(newaxes);
	axis(axis);
	hold on;
end

% Display the boundary
plot([bnd(2) bnd(2)],[0 newaxes(4)],plotstr)

