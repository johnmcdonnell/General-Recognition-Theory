%Plot2BVNpdfs
%
%  plots two bivariate normal probability density functions
%  and their equivocality contours

% Created by Leola A. Alfonso-Reese / 20-February-97
% Copyright (c) 1997
% $Revisions$
%   Date           Modification and Name
%   ----           ---------------------

format compact;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make one mesh plot of pdfs and its corresponding contour plots %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify means, variances, and correlations of the two distributions
muA = [0 0]';
sigmaAx = 1;
sigmaAy = 2;
corrA = .1;

muB = [5 5]';
sigmaBx = 2;
sigmaBy = 1;
corrB = .3;

% Covariance matrices
KA = [sigmaAx^2  corrA*sigmaAx*sigmaAy
      corrA*sigmaAx*sigmaAy  sigmaAy^2];
KB = [sigmaBx^2  corrB*sigmaBx*sigmaBy
      corrB*sigmaBx*sigmaBy  sigmaBy^2];

% Type of plot
%   0 => mesh plot
%   1 => surface plot
plotsurf = 0;

%plot2bivnorms(muA,KA,muB,KB);
plot2bivnorms(muA,KA,muB,KB);
grid on;
figure;
plotbivnormcon(muA,KA);
hold on;
plotbivnormcon(muB,KB);
axis([-3 8 -3 8]);
axis('square');
grid on;
hold off;

