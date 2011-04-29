%Find2dGLC
%
%  finds the best fitting General Linear Classifier for the subject's
%  responses in data file 'subjdemo_2d.dat'.  Specifications for the
%  General Linear Classifier may be found in Ashby, 1992.

% Created by Leola A. Alfonso-Reese / 18-April-95
% Copyright (c) 1995
% $Revisions$
%   Date           Modification and Name
%   ----           ---------------------


format compact;
echo on;
clc;
% This demo finds the best fitting General Linear Classifier for the
% subject's responses in data file 'subjdemo_2d.dat'.  Specifications
% for the General Linear Classifier may be found in Ashby, 1992.

pause  % Press any key to continue.
echo off;

clc;
%
disp('Fit the General Linear Classifier.');
disp('Use Fisher''s Discriminant for initial boundary values.');
fprintf('Parameter format: [noise xycoeffs(1:2) bias]\n\n');

sinit = 10;
fprintf('Initial sigma = %5.2f.\n',sinit);

% Load raw data
load subjdemo_2d.dat
data = subjdemo_2d;

% Massage raw data format
data1(:,1) = data(:,4);
data1(:,2:3) = data(:,2:3);
data1(:,4) = ones(length(data),1);

% Use Fisher's linear discriminant for Initial search values.
% Parameter format: [noise unit_normal(1:2) bias]
fishers_coeffs = fisherdiscrim2d(data,3);
raw_params = [sinit, fishers_coeffs]
start_params = norm_old_2dparams(raw_params)
fprintf('...Searching for best fit\n');
[final_params neglikelihood] = fit_2dGLC(start_params,data1,7);
start_params
final_params
neglikelihood

GLC_results = [final_params neglikelihood];

% FInd the AIC score so that this GLC model may be compared
% with others.
% AIC = 2(-logL + r)
% where r = 3 (2 of the XY coordinates + 1 bias
%              + 1 noise - 1 since the XY & bias params
%              are normalized)
% Note:  neglikelihood from the Fit routine is -logL
negloglikeidx = size(GLC_results,2);
r = 3;
aicGLC = 2*(GLC_results(:,negloglikeidx)+r)';

% Find the slopes and y-intercepts for this bound
[slint(:,1) slint(:,2)] = bnd2slint(GLC_results(2:4));

fprintf('\n\nSEARCH RESULTS:\n');

fprintf('\nNoise\n');
disp(final_params(1));

fprintf('\nBoundary Parameters\n');
disp(final_params(2:4));

fprintf('\nNegative Loglikelihood\n');
disp(neglikelihood);

fprintf('\nAIC score\n');
disp(aicGLC)

fprintf('\nSlope and y-intercept\n');
disp(slint)


echo on;

% Created by Leola A. Alfonso-Reese / 18-April-95
% Copyright (c) 1995

echo off;

