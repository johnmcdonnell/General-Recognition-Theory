function [] = plot3dstim(data,xyzaxes,newfig)
%plot3dstim(data,xyzaxes,newfig)
%  plots a figure of three-dimensional stimulus data points.
%
%  Parameters:
%    data row format:  [cat x y z resp]
%    xyzaxes format:  [minX maxX minY maxY minZ maxZ] (optional)
%    newfig indicates whether to plot the data onto a new
%       figure (newfig = 1) or not (newfig = 0). (optional)

% Created by Leola A. Alfonso-Reese / 13-March-97
% Copyright (c) 1997
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
% 1/29/03          Fix bug when xyzaxes parameter is missing
% 2/21/03          Fix bug of generating two figures

if nargin < 3
	newfig = 1;
end

if newfig == 1
	figure;
end

if nargin < 2
  xyaxes = [];
else
  if ~isempty(xyzaxes)
    axis(xyzaxes);
    axis(axis);
    hold on;
  end
end

% Find A and B response subsets
A_indices = find( data(:,1) == 1 );
B_indices = find( data(:,1) ~= 1 );

% Plot A and B responses
plot3(data(A_indices,2),data(A_indices,3),data(A_indices,4),'g+');
hold on;
plot3(data(B_indices,2),data(B_indices,3),data(B_indices,4),'bo');
hold off;

