function [] = plot2bivnorms(muA,KA,muB,KB,style,smoothvals)
%plot2bivnorms(muA,KA,muB,KB,style,smoothvals)
%  generates a plot containing two bivariate normal
%  distributions
%
%  Parameters:
%    muA, muB are means of the two distributions
%    KA, KB are covariance matrices of the two distributions
%    style is type of plot (optional)
%               0 => mesh plot
%               1 => surface plot
%    smoothvals is a matrix of real numbers affecting number of grid points used
%              for generating the pdf plot (optional).
%              The greater the number, the smoother the plot.
%              However, larger numbers are computationally
%              more expensive. 
%
%              format for smoothvals is [CatA_npts_x CatA_npts_y
%                                        CatB_npts_x CatB_npts_y]

% Created by Leola A. Alfonso-Reese / 20-February-97
% Copyright (c) 1997
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------

handle = newplot;
next = lower(get(handle,'NextPlot'));
hold_state = ishold;
hold on;

sigmaAx = sqrt(KA(1,1));
sigmaAy = sqrt(KA(2,2));

sigmaBx = sqrt(KB(1,1));
sigmaBy = sqrt(KB(2,2));

midpt = mean([muA'; muB']);
% Determine range of X-Y space
range_x = [ min([muA(1)-5*sigmaAx muB(1)-5*sigmaBx]),...
            max([muA(1)+5*sigmaAx muB(1)+5*sigmaBx]) ];
range_y = [ min([muA(2)-5*sigmaAy muB(2)-5*sigmaBy]),...
            max([muA(2)+5*sigmaAy muB(2)+5*sigmaBy]) ];

lightsource = [0,20];

if nargin < 5
	style = 0;
elseif isempty(style)
	style = 0;
end

if nargin > 5
	if size(smoothvals,1) ~= 2 | size(smoothvals,2) ~= 2
		error('The parameter ''smoothvals'' must be a 2x2 matrix.');
	else
		smoothness = smoothvals;
	end
end

az = 35;
el = 50;

axisvals = [];
% Plot pdfs
for j = 1:2
	if j == 1
		mu = muA;
		K = KA;
	else
		mu = muB;
		K = KB;
	end
	sigmax = sqrt(K(1,1));
	sigmay = sqrt(K(2,2));

	if nargin < 6
		smoothness = [40 40];
	else
		smoothness = smoothvals(j,:);
	end
	fprintf('Smoothness values used for pdf plot #%i: %g, %g\n',j,smoothness);
	
	% Generate equally spaced grid of points
%	temp = midpt(1):(range_x(2)-midpt(1))/smoothness(1):range_x(2);
%	X = [range_x(1):(midpt(1)-range_x(1))/smoothness(1):midpt(1) temp(2:length(temp))]';
%	temp = midpt(2):(range_y(2)-midpt(2))/smoothness(2):range_y(2);
%	Y = [range_y(1):(midpt(2)-range_y(1))/smoothness(2):midpt(2) temp(2:length(temp))]';
	X = linspace(range_x(1), range_x(2), 2*smoothness(1))';
	Y = linspace(range_y(1), range_y(2), 2*smoothness(2))';
	
	lenX = length(X);
	
	Z = [];
	for k = 1:length(Y)
		Ys = Y(k)*ones(lenX,1);
		Z(:,k) = bivnormpdf([X Ys],mu,K);
	end

if 0 == 1
	X = linspace(range_x(1), range_x(2), 2*smoothness(1));
	Y = linspace(range_y(1), range_y(2), 2*smoothness(2));
	
	corrxy = K(1,2)/(sigmax*sigmay);
	[zx, zy] = meshgrid( linspace(-3.5, 3.5, 2*smoothness(1)), linspace(-3.5, 3.5, 2*smoothness(1)) );
	Q = ( zx.^2 - 2*corrxy*zx.*zy + zy.^2 ) ./ ( 1 - corrxy^2 );
    Z = exp( -0.5*Q ) ./ ( 2*pi*sigmax*sigmay*sqrt(1-corrxy^2) );
end

	if isempty(axisvals)
		axisvals = [range_x(1:2) range_y(1:2) 0 max(max(Z))];
	else
		axisvals = [range_x(1:2) range_y(1:2) 0 max([max(max(Z)) axisvals(6)]) ];
	end
	maxZ = max(max(Z));
	minZ = min(min(Z));
	
	% Dump nearly zero values
	smallones = find(Z < .0001*axisvals(6));
	if ~isempty(smallones)
		Z(smallones) = NaN*ones(length(smallones),1);
	end
	
	% Draw a pdf
	if style
		surfl(X,Y,Z',lightsource);
		axis(axisvals);
		view(az,el);
		shading interp
		colormap(gray);
	else
		axis(axisvals);
		view(az,el);
		colormap(gray);
		caxis([0 1]);
%		mesh(X,Y,Z',ones(size(Z')));
	mesh(X,Y,Z',zeros(size(Z')));
	end
	hold on;
	
end

% Label graph
xlabel('X');
ylabel('Y');
zlabel('f(X,Y)');

if ~hold_state, set(handle,'NextPlot',next); end

