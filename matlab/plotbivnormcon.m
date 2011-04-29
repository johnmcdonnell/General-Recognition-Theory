function plotbivnormcon(mu,K,likelihood,smoothvals)
%plotbivnormcon(mu,K,likelihood,smoothvals)
%  generates an equivocality contour plot of a bivariate normal
%  distribution
%
%  Parameters:
%    mu is a mean of the distribution
%    K is a covariance matrix of the distribution
%    likelihood is a real number specifying the likelihood of
%              the contour to be plotted
%    smoothvals is a vector of real numbers affecting number of grid 
%              points used for generating the pdf plot (optional).
%              The greater the number, the smoother the plot.
%              However, larger numbers are computationally more
%              expensive. 
%              (format for smoothval is [npts_x npts_y])

% Created by Leola A. Alfonso-Reese / 20-February-97
% Copyright (c) 1997
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
% 2/16/98          Added the likelihood parameter

handle = newplot;
next = lower(get(handle,'NextPlot'));
hold_state = ishold;
hold on;

sigmax = sqrt(K(1,1));
sigmay = sqrt(K(2,2));

% Determine range of X-Y space
range_x = [ mu(1)-5*sigmax mu(1)+5*sigmax ];
range_y = [ mu(2)-5*sigmay mu(2)+5*sigmay ];

lightsource = [0,20];

if nargin > 2
	level = [1 likelihood];
else
	level = 1;
end

if nargin > 3
	if ~(size(smoothvals,1) == 2 | size(smoothvals,2) == 2)
		error('The parameter ''smoothvals'' must be a 2x1 or 1x2 vector.');
	else
		smoothness = smoothvals;
	end
end

az = 35;
el = 50;

axisvals = [];

% Plot contours
if nargin < 5
	smoothness = 2*[5+30*(1-sigmax/(range_x(2)-range_x(1)))^2 5+30*(1-sigmay/(range_y(2)-range_y(1)))^2];
else
	smoothness = smoothvals;
end
fprintf('Smoothness values used for pdf plot #%i: %g, %g\n',j,smoothness);


% Generate equally spaced grid of points
%midpt = mu;
%temp = midpt(1):(range_x(2)-midpt(1))/smoothness(1):range_x(2);
%X = [range_x(1):(midpt(1)-range_x(1))/smoothness(1):midpt(1) temp(2:length(temp))]';
%temp = midpt(2):(range_y(2)-midpt(2))/smoothness(2):range_y(2);
%Y = [range_y(1):(midpt(2)-range_y(1))/smoothness(2):midpt(2) temp(2:length(temp))]';

X = linspace(range_x(1), range_x(2), 2*smoothness(1))';
Y = linspace(range_y(1), range_y(2), 2*smoothness(2))';

lenX = length(X);

Z = [];
for k = 1:length(Y)
	Ys = Y(k)*ones(lenX,1);
	Z(:,k) = bivnormpdf([X Ys],mu,K);
end
if isempty(axisvals)
	axisvals = [range_x(1:2) range_y(1:2) 0 max(max(Z))];
else
	axisvals = [range_x(1:2) range_y(1:2) 0 max([max(max(Z)) axisvals(6)]) ];
end

% Draw a contour plot
hold on;
contour(X,Y,Z,level,'k');

% Label graph
xlabel('X');
ylabel('Y');

if ~hold_state, set(handle,'NextPlot',next); end

