%Find1dGLC
%
%  finds the best fitting General Linear Classifier for the subject's
%  responses in data file 'subjdemo_1d.dat'.  Specifications for the
%  General Linear Classifier may be found in Ashby, 1992.

% Created by Leola A. Alfonso-Reese / 10-March-05
% Copyright (c) 2005
% $Revisions$
%   Date           Modification and Name
%   ----           ---------------------


format compact;
echo on;
clc;
% This demo finds the best fitting General Linear Classifier for the
% subject's responses in data file 'subjdemo_1d.dat'.  Specifications
% for the General Linear Classifier may be found in Ashby, 1992.

pause  % Press any key to continue.
echo off;

clc;
%
disp('Fit the General Linear Classifier.');
disp('Use Fisher''s Discriminant for initial boundary values.');
fprintf('Parameter format: [noise xcoeff bias]\n\n');

sinit = 10;
fprintf('Initial sigma = %5.2f.\n',sinit);

% Load raw data
load subjdemo_1d.dat
data = subjdemo_1d;

% Massage raw data format
clear data1
data1(:,1) = data(:,3);
data1(:,2) = data(:,2);
data1(:,3) = ones(length(data),1);

% Use Fisher's linear discriminant for Initial search values.
% Parameter format: [noise a bias]
sinit = 10;
fishers_coeffs = fisherdiscrim1d(data,3);
raw_params = [sinit, fishers_coeffs]
start_params = norm_old_1dparams(raw_params)
fprintf('...Searching for best fit\n');
[final_params neglikelihood] = fit_1dGLC(start_params,data1,7);

GLC_results = [final_params neglikelihood];

% FInd the AIC score so that this GLC model may be compared
% with others.
% AIC = 2(-logL + r)
% where r = 2 (1 X coordinate + 1 bias
%              + 1 noise - 1 since the X & bias params
%              are normalized)
% Note:  neglikelihood from the Fit routine is -logL
negloglikeidx = size(GLC_results,2);
r = 2;
aicGLC = 2*(GLC_results(:,negloglikeidx)+r)';

fprintf('\n\nSEARCH RESULTS:\n');

fprintf('\nNoise\n');
disp(final_params(1));

fprintf('\nBoundary Parameters\n');
disp(final_params(2:3));

fprintf('\nNegative Loglikelihood\n');
disp(neglikelihood);

fprintf('\nAIC score\n');
disp(aicGLC)


echo on;

% Created by Leola A. Alfonso-Reese / 10-March-05
% Copyright (c) 2005

echo off;

