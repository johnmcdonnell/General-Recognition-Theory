function [] = plot2dquadbnd(bnd,range,plotstr,xyaxes)
%plot2dquadbnd(bnd,range,plotstr,xyaxes)
%  plots a two-dimensional quadratic boundary on the
%  current figure.
%
%  Parameters:
%    bnd format:  [A b c]
%                 where (x,y)*A*(x,y)' + b'*(x,y) + c
%                 is a quadratic decision bound
%    range format:  [minX maxX
%                    minY maxY]
%                   where the last row is a clipping range for Y
%                   and it is optional
%    plotstr:  a character string indicating plot
%              symbols and colors (see MATLAB's built-in
%              plot function)
%    xyaxes format:  [minX maxX minY maxY] (optional)

% Created by Leola A. Alfonso-Reese/ 23-April-96
% (Thanks to Patricia Berretty for suggesting I use MATLAB's contour
% command.)
% Copyright (c) 1996
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
% 12/16/96         by Leola Alfonso-Reese
%                  fixed routine so that it behaves conventionally
%                  when "hold on" or "hold off" is set


if nargin < 4
  xyaxes = [];
end

handle = newplot;
next = lower(get(handle,'NextPlot'));
hold_state = ishold;

hold on;

if ~isempty(xyaxes)
	axis(xyaxes);
	axis(axis);
end

% Display the boundary
xinc = (range(1,2) - range(1,1))/99;
if size(range,1) < 2
	range = [range(1,:); range(1,:)];
end
yinc = (range(2,2) - range(2,1))/99;

[x,y] = meshgrid(range(1,1):xinc:range(1,2),range(2,1):yinc:range(2,2));
z = bnd(1)*x.^2 + bnd(2)*y.^2 + bnd(3)*x.*y + bnd(4)*x + bnd(5)*y + bnd(6);
contour(x,y,z,[0 0],plotstr);

if ~hold_state, set(handle,'NextPlot',next); end


