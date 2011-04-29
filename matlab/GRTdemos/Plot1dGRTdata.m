%Plot1dGRTdata
%
%  plots a stimulus set, an optimal boundary and the best
%  fitting linear boundary.
%
%  Then, this program plots a subject's response data from the
%  file 'subjdemo_1d.dat,' as well as the optimal boundary and
%  the best fitting linear.


% Created by Leola A. Alfonso-Reese / 10-March-05
% Copyright (c) 2005
% $Revisions$
%   Date           Modification and Name
%   ----           ---------------------


format compact;
echo on;
clc;
% This demo plots a stimulus set, an optimal boundary and the best
% dimensional boundaries.
%
% Then, it plots a subject's response data from the file
% 'subjdemo_1d.dat,' as well as the optimal boundary and the best
% fitting linear boundary.

pause  % Press any key to continue.
echo off;

% Load subject's data file
disp('...Loading a data file');
load subjdemo_1d.dat
data = subjdemo_1d;

% Specify optimal boundary, subject's boundary
% any other boundaries, and xy axis limits.
optbnd = [-1 200];
subjlinbnd = [-1.0000  199.4435];
xaxes = [0 400];

% Plot stimuli, optimal bound and other interesting bounds
disp('...Plotting stimulus set, optimal bound and best fitting linear boundary');
plot1dstim(data,xaxes,1);
hold on;
plot1dlinbnd(optbnd,'r-',xaxes);
title('Stimuli');
xlabel('Length');
grid on;
hold off;

% Plot responses, optimal bound and subject's bound
disp('...Plotting responses, optimal bound, best fitting GLC and best fitting GQC');
plot1dresp(data,xaxes,1);
hold on;
plot1dlinbnd(optbnd,'r-',xaxes);
plot1dlinbnd(subjlinbnd,'c-',xaxes);
title('Subject Responses');
xlabel('Length');
grid on;
hold off;


echo on;

% Created by Leola A. Alfonso-Reese / 10-March-05
% Copyright (c) 2005

echo off;


