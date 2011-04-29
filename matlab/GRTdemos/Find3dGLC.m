%Find3dGLC
%
%  finds the best fitting General Linear Classifier for the subject's
%  responses in data file 'subjdemo_3d.dat'.  Specifications for the
%  General Linear Classifier may be found in Ashby, 1992.

% Created by Leola A. Alfonso-Reese / 13-March-97
% Copyright (c) 1997
% $Revisions$
%   Date           Modification and Name
%   ----           ---------------------


format compact;
echo on;
clc;
% This demo finds the best fitting General Linear Classifier for the
% subject's responses in data file 'subjdemo_3d.dat'.  Specifications
% for the General Linear Classifier may be found in Ashby, 1992.

pause  % Press any key to continue.
echo off;

clc;
%
disp('Fit the General Linear Classifier.');
disp('Use Fisher''s Discriminant for initial boundary values.');
fprintf('Parameter format: [noise xyzcoeffs(1:3) bias]\n\n');

sinit = 10;
fprintf('Initial sigma = %5.2f.\n',sinit);

% Load raw data
load subjdemo_3d.dat
data = subjdemo_3d;

% Massage raw data format
data1(:,1) = data(:,5);
data1(:,2:4) = data(:,2:4);
data1(:,5) = ones(length(data),1);

% Use Fisher's linear discriminant for Initial search values.
% Parameter format: [noise unit_normal(1:3) bias]
fishers_coeffs = fisherdiscrim3d(data1,9);
raw_params = [sinit, fishers_coeffs]
start_params = norm_old_3dparams(raw_params)
fprintf('...Searching for best fit\n');
[final_params neglikelihood] = fit_3dGLC(start_params,data1,7);
start_params
final_params
neglikelihood

GLC_results = [final_params neglikelihood];

% FInd the AIC score so that this GLC model may be compared
% with others.
% AIC = 2(-logL + r)
% where r = 4 (3 of the XYZ coordinates + 1 bias
%              + 1 noise - 1 since the XYZ & bias params
%              are normalized)
% Note:  neglikelihood from the Fit routine is -logL
negloglikeidx = size(GLC_results,2);
r = 4;
aicGLC = 2*(GLC_results(:,negloglikeidx)+r)';

fprintf('\n\nSEARCH RESULTS:\n');

fprintf('\nNoise\n');
disp(final_params(1));

fprintf('\nBoundary Parameters\n');
disp(final_params(2:5));

fprintf('\nNegative Loglikelihood\n');
disp(neglikelihood);

fprintf('\nAIC score\n');
disp(aicGLC)


echo on;

% Created by Leola A. Alfonso-Reese / 13-March-97
% Copyright (c) 1997

echo off;

