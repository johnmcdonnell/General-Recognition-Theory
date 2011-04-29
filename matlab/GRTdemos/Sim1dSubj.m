%Sim1dSubj
%
%  loads a 1D stimulus set and simulates response data of a subject
%  participating in a two-category classification experiment.

% Created by Leola A. Alfonso-Reese & Brian Spiering/ 1-Nov-02
% Copyright (c) 2002
% $Revisions$
%   Date           Modification and Name
%   ----           ---------------------
%  3/10/05         Fixed plot function of response data
%                  by LAR

format compact;
echo on;
clc;
% This demo loads a 1D stimulus set and simulates response data of a
% subject participating in a two-category classification experiment.
%
% The row format of the stimulus data file is:  [category x]
%
% After the responses are generated, they are saved into a data
% file called 'subjdemo_1d.dat' and finally plotted.

pause  % Press any key to continue.
echo off;

clc;
% Load stimulus data file
disp('...Loading a stimulus set');
load day1_1d.dat
stimuli = day1_1d;

% Specify noise for the hypothetical subject and the boundary that
% this subject uses to classify stimuli.
subject_params = [10 -1 200];		% [noise a b] where 0 = a*x + b

% Generate subject's responses
disp('...Generating hypothetical subject responses');
respdata = sim1dlin(stimuli,subject_params);

%fprintf('\n...Saving data to a file called day1.dat\n');
echo on;

%
% NOTE:  The following code that writes data to a file is
% commented out so that the data file called 'subjdemo_1d.dat'
% is not created each time someone runs this demo.  The data
% file that would be created if this code was run,
% 'subjdemo_1d.dat', is in the GRTdemos directory in the GRT
% Toolbox.
%
%fid = fopen('subjdemo_1d.dat','w');
%fprintf(fid,'%i %f %i\n',respdata');
%fclose(fid);
echo off;

fprintf('\nThe stimulus set is in the variable ''stimuli'' and\n');
fprintf('the response data set is in the varable ''respdata''.\n');

fprintf('\nFirst 10 category labels and stimulus points:\n');
stimuli(1:10,:)

fprintf('\nFirst 10 responses:\n');
respdata(1:10,3)

% Plot hypothetical subject responses
disp('...Plotting the hypothetical subject responses');
plot1dresp(respdata,[0 400],1);
title('Hypothetical Responses');

echo on;

% Created by Leola A. Alfonso-Reese / 1-Nov-02
% Copyright (c) 2002

echo off;

