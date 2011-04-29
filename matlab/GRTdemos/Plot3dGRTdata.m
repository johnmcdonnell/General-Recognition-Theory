%Plot3dGRTdata
%
%  plots a stimulus set and the optimal quadratic bound.
%
%  Then, this program plots a subject's response data from the
%  file 'subjdemo_3d.dat,' as well as the optimal boundary and
%  the best fitting linear and quadratic boundaries.


% Created by Leola A. Alfonso-Reese / 14-March-97
% Copyright (c) 1997
% $Revisions$
%   Date           Modification and Name
%   ----           ---------------------
%  12/21/04        by Leola A. Alfonso-Reese
%                  modified call to plot3d bound routines
%                  to match changes in the plot functions


format compact;
echo on;
clc;
% This demo plots a stimulus set and the optimal quadratic bound.
%
% Then, it plots a subject's response data from the file
% 'subjdemo_3d.dat,' as well as the optimal boundary and the best
% fitting linear and quadratic boundaries.

pause  % Press any key to continue.
echo off;

% Load subject's data file
disp('...Loading a data file');
load subjdemo_3d.dat
data = subjdemo_3d;

% Specify optimal boundary, subject's boundary
% any other boundaries, and xyz axis limits.
optquadbnd = [-0.0017   -0.0013   -0.0017    0  0.0042  0  0.5146,...
              0.7330   -0.9032 -109.8578];
subjlinbnd = [0.3414    0.6166   -0.7094 -102.7599];
subjquadbnd = [-0.0621 -0.0714 -0.0414  0.1178  0.1185 -0.0474  -9.8158,...
               13.4258   -9.1668  -93.3853];           
%optxbnd = [1.0000 0.0000 -200.0000];
%optybnd = [0.0000 1.0000 -120.0000];
xyzaxes = [0 400 0 400 0 400];
range = [xyzaxes(1) xyzaxes(2); xyzaxes(3) xyzaxes(4); xyzaxes(5) xyzaxes(6)];

% Plot stimuli, optimal bound and other interesting bounds
disp('...Plotting stimulus set and optimal bound');
plot3dstim(data,xyzaxes,1);
view(110,35);
hold on;
%plot3dquadbnd(optquadbnd,range);
plot3dquadbnd(optquadbnd,range,'r');
title('Stimuli');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
hold off;

% Plot responses, optimal bound and subject's bound
disp('...Plotting responses, optimal bound, best fitting GLC and best fitting GQC');
plot3dresp(data,xyzaxes,1);
view(110,35);
hold on;
plot3dquadbnd(subjquadbnd,range);
plot3dlinbnd(subjlinbnd,range(1:2,1:2),'m');
title('Subject Responses');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
hold off;

echo on;

% Created by Leola A. Alfonso-Reese / 14-March-97
% Copyright (c) 1997

echo off;


