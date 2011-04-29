%Sim3dSubj
%
%  loads a 3D stimulus set and simulates response data of a subject
%  participating in a two-category classification experiment.

% Created by Leola A. Alfonso-Reese / 11-March-94
% Copyright (c) 1994
% $Revisions$
%   Date           Modification and Name
%   ----           ---------------------
%  3/13/97         by Leola A. Alfonso-Reese
%                  made changes so that it fits the mold
%                  of other routines in the GRT toolbox


format compact;
echo on;
clc;
% This demo loads a 3D stimulus set and simulates response data of a
% subject participating in a two-category classification experiment.
%
% The row format of the stimulus data file is:  [category x y z]
%
% After the responses are generated, they are saved into a data
% file called 'subjdemo_3d.dat' and finally plotted.

pause  % Press any key to continue.
echo off;

clc;
% Load stimulus data file
disp('...Loading a stimulus set');
load day1_3d.dat
stimuli = day1_3d;

% Specify noise for the hypothetical subject and the boundary that
% this subject uses to classify stimuli.
% [noise a1 a2 a3 b] where 0 = a1*x + a2*y + a3*z + b
subject_params = [10  0.5584    0.3657   -0.7446  -96.2391];


% Generate subject's responses
disp('...Generating hypothetical subject responses');
respdata = sim3dlin(stimuli,subject_params);

%fprintf('\n...Saving data to a file called day1.dat\n');
echo on;

%
% NOTE:  The following code that writes data to a file is
% commented out so that the data file called 'subjdemo_3d.dat'
% is not created each time someone runs this demo.  The data
% file that would be created if this code was run,
% 'subjdemo_3d.dat', is in the GRTdemos directory in the GRT
% Toolbox.
%
%fid = fopen('subjdemo_3d.dat','w');
%fprintf(fid,'%i %f %f %f %i\n',respdata');
%fclose(fid);
echo off;

fprintf('\nThe stimulus set is in the variable ''stimuli'' and\n');
fprintf('the response data set is in the varable ''respdata''.\n');

fprintf('\nFirst 10 category labels and stimulus points:\n');
stimuli(1:10,:)

fprintf('\nFirst 10 responses:\n');
respdata(1:10,5)

% Plot hypothetical subject responses
disp('...Plotting the hypothetical subject responses');
plot3dresp(respdata,[100 300 100 300 100 300],1);
title('Hypothetical Responses');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;

echo on;

% Created by Leola A. Alfonso-Reese / 29-October-96
% Copyright (c) 1996

echo off;

