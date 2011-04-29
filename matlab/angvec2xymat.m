function xymat = angvec2xymat(anglevec)
%xy = angvec2xymat(anglevec)
%  converts a vector of angles to a matrix of (x,y) coordinates
%  on a circle of radius one.

% Created by Leola A. Alfonso-Reese / 6-June-94
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%   28-July-95     Modified to work for vectors of angles
%                  by Leola A. Alfonso-Reese


xymat = cos(anglevec)';
xymat(:,2) = sin(anglevec)';

