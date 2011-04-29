function [xx,yy,zz] = plot3dlinbnd(bnd,range,linespec)
%[xx,yy,zz] = plot3dlinbnd(bnd,range,linespec)
%  plots a three-dimensional linear boundary on the
%  current figure.
%
%  Parameters:
%    bnd format:  [a1 a2 a3 b]
%       where a1*x + a2*y + a3*z + b = 0 is a plane
%    range format:  [minX maxX
%                    minY maxY]
%                    minZ maxZ]
%                   where the last row is a clipping range for Z
%                   and it is optional
%    linespec is a line specification string (Default is '-r', red lines.)

% Created by Leola A. Alfonso-Reese / 19-October-94
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  3/13/97         by Leola A. Alfonso-Reese
%                  made changes so that it fits the mold
%                  of other routines in the GRT toolbox
%  12/14/04        by Leola A. Alfonso-Reese
%                  changed color parameter to linespec

xx = [];
xx = [];
yy = [];

if nargout == 0
	handle = newplot;
	next = lower(get(handle,'NextPlot'));
	hold_state = ishold;
	hold on;
end

if nargin < 3
	linespec = '-r';
elseif isempty(linespec)
	linespec = '-r';
end

b = bnd(1:3);
c = bnd(4);

% Compute linear boundary
%Xvec = range(1,1):10:range(1,2);
%Yvec = range(2,1):10:range(2,2);
increment = (range(1,2) - range(1,1))/40;
Xvec = range(1,1):increment:range(1,2);
Yvec = range(2,1):increment:range(2,2);
if size(range,1) > 2
	Zrange = range(3,:);
else
	Zrange = range(1,:);
end

if b(3) == 0
	if b(2) ~= 0
	% Plot vertical lines for b1*x + b2*y + c = 0
	X = [];
	Y = [];
	Z = [];
	Yvec = (-b(1)*Xvec - c)./b(2);
	for i = 1:length(Xvec)
		X = [X [Xvec(i); Xvec(i)] ];
	  	Y = [Y [Yvec(i); Yvec(i)] ];
	  	Z = [Z Zrange'];
		end
		lastX = length(Xvec);
		lastY = length(Yvec);
		for i = 1:length(Xvec)
	  	X = [X [Xvec(1); Xvec(lastX)] ];
	  	Y = [Y [Yvec(1); Yvec(lastY)] ];
	  	Z = [Z [Xvec(i); Xvec(i)]];
		end
		
		if nargout == 0
			colormap(hsv);
			caxis([0,50]);  
			plot3(X,Y,Z,linespec);
		else
		   xx = X; yy = Y; zz = Z;
		end
		return;
		
	else
		% Make mesh for b1*x + c = 0
		Xvec = (-c/b(1))*ones(length(Xvec));
		Z = [];
		for i = 1:length(Yvec)
			Z = [Z; Zrange(1):10:Zrange(2)];
		end
		fprintf('\nWarning:\nThere are many possible boundaries where X = %6.2f.\n',b(1));
  end

else
	% Make mesh for b1*x + b2*y + b3*z + c = 0
	Z = [];
	for i = 1:length(Yvec)
		Y = Yvec(i)*ones(1,length(Xvec));
		Z = [Z; (-b(1)*Xvec - b(2)*Y - c)./b(3)];
	end
end

% Clip values outside of Z range if specified
if size(range,1) > 2
	clipindices = find(Z < Zrange(1) | Z > Zrange(2));
	Z(clipindices) = NaN*ones(length(clipindices),1);
end

% Graph linear bound mesh.
if nargout == 0
	colormap(hsv);
	caxis([0,50]); 
	color = colorstr2colorval(linespec);
	mesh(Xvec,Yvec,Z,color*ones(size(Z)));
	if size(range,1) == 3
		axis([range(1,1:2) range(2,1:2) range(3,1:2)]);
		axis(axis);
	end
else
	xx = Xvec; yy = Yvec; zz = Z;
end

if nargout == 0
	if ~hold_state, set(handle,'NextPlot',next); end
end

