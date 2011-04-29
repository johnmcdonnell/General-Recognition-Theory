function [slvec,intvec] = bnd2slint(bndmat,flag)
%[slvec intvec] = bnd2slint(bndmat,flag)
%  converts a matrix of two-dimensional linear boundaries
%  to slope and y_intercept vectors.
%
%  When the bound is a vertical line, bnd2slint returns
%  -infinity for the slope and y_intercept values.  To obtain
%  the x_intercept for these cases, set the flag to 1.
%
%  Parameters:
%    bnd row format: [a1 a2 b], where a1*x+ a2*y + b = 0
%    flag:  1 => return x-intercept when bound is vertical line
%           0 => (default) return -infinity when bound is vertical line

% Created by Leola A. Alfonso-Reese / 12-August-95
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  6/27/96         Added flag parameter

slvec = -bndmat(:,1)./bndmat(:,2);
intvec = -bndmat(:,3)./bndmat(:,2);

if nargin == 2
	if flag == 1
		x_intindices = find(bndmat(:,2) == 0);
		if ~isempty(x_intindices)
			intvec(x_intindices) = bndmat(x_intindices,3)./bndmat(x_intindices,1);
		end
	end
end

