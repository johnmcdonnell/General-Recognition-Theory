function [] = plot3dresp(data,xyaxes,newfig)
% plot3dresp(data,xyaxes,newfig)
%  plots a figure of three-dimensional response data points.
%
%  Parameters:
%    data row format:  [cat x y z resp]
%    xyaxes format:  [minX maxX minY maxY minZ maxZ] (optional)
%    newfig indicates whether to plot the data onto a new
%       figure (newfig = 1) or not (newfig = 0). (optional)

% Created by Leola A. Alfonso-Reese / 13-March-97
% Copyright (c) 1997
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
A_indices = find( data(:,5) == 1 );
B_indices = find( data(:,5) ~= 1 );

% Plot A and B responses
plot3(data(A_indices,2),data(A_indices,3),data(A_indices,4),'g+');
hold on;
plot3(data(B_indices,2),data(B_indices,3),data(B_indices,4),'bo');
hold off;

