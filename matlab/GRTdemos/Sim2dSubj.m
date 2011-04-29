%Sim2dSubj
%
%  loads a 2D stimulus set and simulates response data of a subject
%  participating in a two-category classification experiment.

% Created by Leola A. Alfonso-Reese / 11-March-94
% Copyright (c) 1995
% $Revisions$
%   Date           Modification and Name
%   ----           ---------------------


format compact;
echo on;
clc;
% This demo loads a 2D stimulus set and simulates response data of a
% subject participating in a two-category classification experiment.
%
% The row format of the stimulus data file is:  [category x y]
%
% After the responses are generated, they are saved into a data
% file called 'subjdemo_2d.dat' and finally plotted.

pause  % Press any key to continue.
echo off;

clc;
% Load stimulus data file
disp('...Loading a stimulus set');
load day1_2d.dat
stimuli = day1_2d;

% Specify noise for the hypothetical subject and the boundary that
% this subject uses to classify stimuli.
subject_params = [10 .50 -1 10];		% [noise a1 a2 b] where 0 = a1*x + a2*y +b

% Generate subject's responses
disp('...Generating hypothetical subject responses');
respdata = sim2dlin(stimuli,subject_params);

%fprintf('\n...Saving data to a file called day1.dat\n');
echo on;

%
% NOTE:  The following code that writes data to a file is
% commented out so that the data file called 'subjdemo_2d.dat'
% is not created each time someone runs this demo.  The data
% file that would be created if this code was run,
% 'subjdemo_2d.dat', is in the GRTdemos directory in the GRT
% Toolbox.
%
%fid = fopen('subjdemo_2d.dat','w');
%fprintf(fid,'%i %f %f %i\n',respdata');
%fclose(fid);
echo off;

fprintf('\nThe stimulus set is in the variable ''stimuli'' and\n');
fprintf('the response data set is in the varable ''respdata''.\n');

fprintf('\nFirst 10 category labels and stimulus points:\n');
stimuli(1:10,:)

fprintf('\nFirst 10 responses:\n');
respdata(1:10,4)

% Plot hypothetical subject responses
disp('...Plotting the hypothetical subject responses');
plot2dresp(respdata,[0 300 0 300],1);
title('Hypothetical Responses');
grid on;

echo on;

% Created by Leola A. Alfonso-Reese / 18-April-95
% Copyright (c) 1995

echo off;

