function [] = plot3drespcustom(data,s1,s2,xyaxes,newfig)
% plot3dresp(data,s1,s2,xyaxes,newfig)
%  plots a customized figure of three-dimensional response
%data points.
%
%  Parameters:
%    data row format:  [cat x y z resp]
%    str1:  a character string indicating plot symbol and color
%           for category 1
%    str2:  a character string indicating plot symbol and color
%           for category 2
%    xyaxes format:  [minX maxX minY maxY minZ maxZ] (optional)
%    newfig indicates whether to plot the data onto a new
%       figure (newfig = 1) or not (newfig = 0). (optional)

% Created by Leola A. Alfonso-Reese / 11-July-98
% Copyright (c) 1998
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------

if nargin < 4
  xyaxes = [];
end

if nargin < 5
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
plot3(data(A_indices,2),data(A_indices,3),data(A_indices,4),s1,...
      data(B_indices,2),data(B_indices,3),data(B_indices,4),s2);
hold off;

