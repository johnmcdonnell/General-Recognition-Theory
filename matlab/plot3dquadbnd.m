function [xx1,xx2,yy1,yy2,zz1,zz2] = plot3dquadbnd(bnd,range,linespec);
%[xx1,xx2,yy1,yy2,zz1,zz2] = plot3dquadbnd(bnd,range,linespec);
%  plots a three-dimensional quadratic boundary on the
%  current figure.
%
%  Parameters:
%    bnd format:  [a(1:6) b(1:3) c]
%        where h(x,y,z) = a1*xx + a2*yy + a3*zz
%                         + a4*xy + a5*xz + a6*yz
%                         + b1*x + b2*y + b3*z + c = 0
%        is a quadratic decision bound
%    range format:  [minX maxX
%                    minY maxY]
%                    minZ maxZ]
%       where the last row is a clipping range for Z
%       and it is optional
%    linespec is a line specification string (Default is '-r', red lines.)

% Created by Leola A. Alfonso-Reese / 20-October-94
% Copyright (c) 1994
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  3/13/97         by Leola A. Alfonso-Reese
%                  made changes so that it fits the mold
%                  of other routines in the GRT toolbox
%  11/24/04        by Leola A. Alfonso-Reese
%                  made more efficient and
%                  fixed problem of holes in the mesh

a = bnd(1:6);
b = bnd(7:9);
c = bnd(10);

xx1 = [];
xx2 = [];
yy1 = [];
yy2 = [];
zz1 = [];
zz2 = [];

%ncontourlines = 35;

if nargin < 3
	linespec = '-r';
elseif isempty(linespec)
	linespec = '-r';
end

if sum(abs(a)) == 0
	% The bound is linear
	colorval = colorstr2colorval(linespec);
	if nargout == 0
		plot3dlinbnd([b c],range,colorval);
	else
		[xx1,yy1,zz1] = plot3dlinbnd([b c],range,colorval);
	end
	return;
	
else
	% The bound is non-linear

	if nargout == 0
		handle = newplot;
		next = lower(get(handle,'NextPlot'));
		hold_state = ishold;
		hold on;
	end
	
	a = bnd(1:6);
	b = bnd(7:9);
	c = bnd(10);
	
	% Compute the quadratic mesh.
	Xvec = range(1,1):10:range(1,2);
	Yvec = range(2,1):10:range(2,2);
	if size(range,1) > 2
  		Zrange = range(3,:);
	else
		Zrange = range(1,:);
	end
	clevels = Zrange(1):(Zrange(2)-Zrange(1))/40:Zrange(2);
	
	Z = [];
	Z1 = [];
	Z2 = [];
	
	if a(3) == 0
		% There is only one solution to the quadratic equation.	
		[X,Y] = meshgrid(Xvec,Yvec);
		
		t1 = a(5)*X + a(6)*Y + b(3);
		t2 = a(1)*X.^2 + a(2)*Y.^2 + a(4)*X.*Y + b(1)*X + b(2)*Y + c;
		Z = -t2./t1;
%		t1_not_zero = find(t1 ~= 0);
%		Z(t1_not_zero) = -t2(t1_not_zero)./t1(t1_not_zero);
		contour3(X,Y,Z,clevels,linespec)		
	else
		% There are two solutions to the quadratic equation.
		% Generate two meshes.

		[X,Y] = meshgrid(Xvec,Yvec);
		t1 = a(5)*X + a(6)*Y + b(3);
		t2 = a(1)*X.^2 + a(2)*Y.^2 + a(4)*X.*Y + b(1)*X + b(2)*Y + c;
		Z1 = (-t1 - sqrt(t1.^2 - 4*a(3)*t2))/(2*a(3));  % Zminus
		Z2 = (-t1 + sqrt(t1.^2 - 4*a(3)*t2))/(2*a(3));  % Zplus
		contour3(X,Y,Z1,clevels,linespec)
		contour3(X,Y,Z2,clevels,linespec)
	end

% 	% If the Z range is specified, then clip values outside of
% 	% the Z range.
% 	if size(range,1) > 2
% 		if isempty(Z)
% 			clipindices = find(Z1 < Zrange(1) | Z1 > Zrange(2));
% 			Z1(clipindices) = NaN*ones(length(clipindices),1);
% 			clipindices = find(Z2 < Zrange(1) | Z2 > Zrange(2));
% 			Z2(clipindices) = NaN*ones(length(clipindices),1);
% 		else
% 			clipindices = find(Z < Zrange(1) | Z > Zrange(2));
% 			Z(clipindices) = NaN*ones(length(clipindices),1);
% 		end
% 	end

	xx1 = X;
	yy1 = Y;
	if isempty(Z)
		xx2 = xx1;
		yy2 = yy1;
		zz1 = Z1;
		zz2 = Z2;
	else
		zz1 = Z;
	end		

end

% Plot the mesh(es).
% if nargout == 0
% 	colormap(hsv);
% 	caxis([0,50]);
% 	mesh(xx1,yy1,zz1,color*ones(size(zz1)));
% 	if ~isempty(zz2)
% 		mesh(xx2,yy2,zz2,color*ones(size(zz2)));
% 	end
% end

if nargout == 0
	if ~hold_state, set(handle,'NextPlot',next); end
end

