function anglevec = xymat2angvec(xymat);
%anglevec = xymat2angvec(xymat);
%  converts a matrix of (x,y) coordinates on a
%  circle of radius one to a vector of angles.

% Created by Leola A. Alfonso-Reese / 6-June-94
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%   9/4/95         Fixed bug
%                  by Leola


% Treat xy pairs with y-coordinates greater than zero
% differently.
pos_ys = find(xymat(:,2) > 0);
neg_ys = find(xymat(:,2) <= 0);

anglevec = zeros(size(xymat,1),1);
anglevec(pos_ys) = acos(xymat(pos_ys,1));
anglevec(neg_ys) = 2*pi - acos(xymat(neg_ys,1));

