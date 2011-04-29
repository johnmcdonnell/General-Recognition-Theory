function angles = xyzmat2angmat(xyz);
%angles = xyzmat2angmat(xymat);
%  converts a matrix of (x,y,z) coordinates on a
%  circle of radius one to a matrix of angle pairs.

% Created by Leola A. Alfonso-Reese / 1-September-92
% Copyright (c) 1992
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  3/12/97         by Leola A. Alfonso-Reese
%                  made minor changes so that it fits the mold
%                  of other routines in the GRT toolbox

angles = zeros(size(xyz,1),2);

angles(:,1) = asin(xyz(:,3));
if (xyz(:,2) > 0)
  angles(:,2) = real( acos( xyz(:,1)./cos(angles(:,1)) ) );
else
  angles(:,2) = 2*pi - real( acos( xyz(:,1)./cos(angles(:,1)) ) );
end

