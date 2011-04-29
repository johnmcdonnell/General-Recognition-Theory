function [] = plot2dstim(data,xyaxes,newfig)
%plot2dstim(data,xyaxes,newfig)
%  plots a figure of two-dimensional stimulus data points.
%
%  Parameters:
%    data row format:  [cat x y resp]
%    xyaxes format:  [minX maxX minY maxY] (optional)
%    newfig indicates whether to plot the data onto a new
%       figure (newfig = 1) or not (newfig = 0). (optional)

% Created by Leola A. Alfonso-Reese / 7-August-95
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------


if nargin < 2
  xyaxes = [];
end

if nargin < 3
	newfig = 1;
end

if newfig == 1
	figure;
end

if ~isempty(xyaxes)
	axis(xyaxes);
	axis(axis);
	hold on;
end

% Find A and B response subsets
A_indices = find( data(:,1) == 1 );
B_indices = find( data(:,1) ~= 1 );

% Plot A and B responses
plot(data(A_indices,2),data(A_indices,3),'g+',...
	data(B_indices,2),data(B_indices,3),'bo');
hold off;

