function pc = quadbndpercorr(pts,bnd,refpnt)
%pc quadbndpercorr(pts,bnd,refpnt) 
%  computes the percent correct obtained by categorizing data points
%  using the quadratic decision boundary
%       h(x,y,...) = (x,y,...)A(x,y,...) + b(x,y,...) + c
%
%  This routine assumes that the points all belong to just one category.
%  It uses the reference point 'refpnt' to determine the correct side of
%  the bound for classifying 'pts'.  If no reference point is specified,
%  then the mean of the points becomes the reference point.
% 
%  Parameters:
%    pts row format:  [x y ...]
%    bnd is a row vector of the form [coeffx coeffy ... const]
%    refpnt is a column vector of the form [x y ...]

%	2/24/97	lar	Created
%       7/13/98 lar     Modified to work for 3 dimensions as well

% Default parameters
if nargin == 2
  refpnt = mean(pts);
end

[m,n] = size(pts);

if n == 2
	A = [   bnd(1) .5*bnd(3)
	     .5*bnd(3)    bnd(2)];
	b = bnd(4:5)';
	c = bnd(6);
else
	A = [   bnd(1) .5*bnd(4) .5*bnd(5)
	     .5*bnd(4)    bnd(2) .5*bnd(6)
	     .5*bnd(5) .5*bnd(6)    bnd(3)];
	b = bnd(7:9)';
	c = bnd(10);
end

% To determine whether a response is correct or not, I need to know on which
% side of the boundary the points lie.  This is determined by checking on
% which side of the boundary the reference point lies.
h = refpnt*A*refpnt' + b'*refpnt' + c;

% determine h([x y ...]) for each point
temp = pts*A;

% This is what used to work for 2 dimensions only
%temp(:,1).*pts(:,1) + temp(:,2).*pts(:,2) + pts*b + c*ones(m,1);

hvec = zeros(m,1);
for j = 1:n
	hvec = hvec + temp(:,j).*pts(:,j);
end
hvec = hvec + pts*b + c*ones(m,1);

if h < 0
  % any other h's less than 0 are categorized correctly
  numcorr = length(find(hvec < 0));
  else
  % any other h's greater than 0 are categorized correctly
  numcorr = length(find(hvec > 0));
end

pc = numcorr/m * 100;


