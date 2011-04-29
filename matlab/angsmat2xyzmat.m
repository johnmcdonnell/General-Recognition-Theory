function xyz = angsmat2xyzmat(angles)
%xyz = angsmat2xyzmat(angles)
%  converts a matrix of angle pairs to a matrix of
%  (x,y,z) coordinates on a circle of radious one.

% Created by Leola A. Alfonso-Reese / 1-September-92
% Copyright (c) 1992
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  3/12/97         by Leola A. Alfonso-Reese
%                  made minor changes so that it fits the mold
%                  of other routines in the GRT toolbox

xyz = zeros(size(angles,1),3);

xyz(:,3) = sin(angles(:,1));
tmp = cos(angles(:,1));
xyz(:,1) = tmp.*cos(angles(:,2));
xyz(:,2) = tmp.*sin(angles(:,2));

