%Plot2dGRTdata
%
%  plots a stimulus set, an optimal boundary and the best
%  dimensional boundaries.
%
%  Then, this program plots a subject's response data from the
%  file 'subjdemo_2d.dat,' as well as the optimal boundary and
%  the best fitting linear and quadratic boundaries.


% Created by Leola A. Alfonso-Reese / 16-May-95
% Copyright (c) 1995
% $Revisions$
%   Date           Modification and Name
%   ----           ---------------------
%  7/18/96         Added quadratic boundary
%                  by Leola


format compact;
echo on;
clc;
% This demo plots a stimulus set, an optimal boundary and the best
% dimensional boundaries.
%
% Then, it plots a subject's response data from the file
% 'subjdemo_2d.dat,' as well as the optimal boundary and the best
% fitting linear and quadratic boundaries.

pause  % Press any key to continue.
echo off;

% Load subject's data file
disp('...Loading a data file');
load subjdemo_2d.dat
data = subjdemo_2d;

% Specify optimal boundary, subject's boundary
% any other boundaries, and xy axis limits.
optbnd = [.6 -1 0];
subjlinbnd = [0.4441   -0.8960   11.5237];
subjquadbnd = [0.0003    0.0041   -0.0005    0.4099   -1.7901   65.4630];
optxbnd = [1.0000 0.0000 -200.0000];
optybnd = [0.0000 1.0000 -120.0000];
xyaxes = [0 350 0 350];

% Plot stimuli, optimal bound and other interesting bounds
disp('...Plotting stimulus set, optimal bound and best dimensional bounds');
plot2dstim(data,xyaxes,1);
hold on;
plot2dlinbnd(optbnd,'r-',xyaxes);
plot2dlinbnd(optxbnd,'r--',xyaxes);
plot2dlinbnd(optybnd,'r--',xyaxes);
title('Stimuli');
xlabel('Length');
ylabel('Orientation');
grid on;
hold off;

% Plot responses, optimal bound and subject's bound
disp('...Plotting responses, optimal bound, best fitting GLC and best fitting GQC');
plot2dresp(data,xyaxes,1);
hold on;
plot2dlinbnd(optbnd,'r-',xyaxes);
plot2dlinbnd(subjlinbnd,'c-',xyaxes);
plot2dquadbnd(subjquadbnd,xyaxes,'m-',xyaxes);
title('Subject Responses');
xlabel('Length');
ylabel('Orientation');
grid on;
hold off;


echo on;

% Created by Leola A. Alfonso-Reese / 18-April-95
% Copyright (c) 1995

echo off;


