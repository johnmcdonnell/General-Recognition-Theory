%Find3dGQC
%
%  finds the best fitting General Quadratic Classifier for the subject's
%  responses in data file 'subjdemo_3d.dat'.  Specifications for the
%  General Quadratic Classifier may be found in Ashby, 1992.

% Created by Leola A. Alfonso-Reese / 13-March-97
% Copyright (c) 1997
% $Revisions$
%   Date           Modification and Name
%   ----           ---------------------


format compact;
echo on;
clc;
% This demo finds the best fitting General Quadratic Classifier for the
% subject's responses in data file 'subjdemo_3d.dat'.  Specifications
% for the General Quadratic Classifier may be found in Ashby, 1992.

pause  % Press any key to continue.
echo off;

clc;
%
disp('Fit the General Quadratic Classifier.');
disp('Use Fisher''s Discriminant for initial boundary values.');
fprintf('Parameter format: [pnoise cnoise Amat_entries(1:6) b_entries(1:3) c_bias]\n\n');

pnoiseinit = 10;     % perceptual noise variance
cnoiseinit = 100;    % criterial noise variance
fprintf('Initial pnoisenit = %5.2f.\n',pnoiseinit);
fprintf('Initial cnoiseinit = %5.2f.\n',cnoiseinit);

% Load raw data
load subjdemo_3d.dat
data = subjdemo_3d;

% Massage raw data format
data1(:,1) = data(:,5);
data1(:,2:4) = data(:,2:4);
data1(:,5) = ones(length(data),1);

% Use Fisher's quadratic discriminant for Initial search values.
% Parameter format: [noise unit_normal(1:2) bias]
fishers_coeffs = fisherdiscrim3d(data,22);
raw_params = [pnoiseinit, cnoiseinit, fishers_coeffs]
start_params = raw_params
fprintf('...Searching for best fit\n');
[final_params neglikelihood] = fit_3dGQC(start_params,data1,7);
start_params
final_params
neglikelihood

GQC_results = [final_params neglikelihood];

% FInd the AIC score so that this GQC model may be compared
% with others.
% AIC = 2(-logL + r)
% where r = 11 (6 matrix entries + 3 XYZ coordinates + 1 bias
%              + 2 noise - 1 since we could set one parameter to equal 1))
negloglikeidx = size(GQC_results,2);
r = 11;
aicGQC = 2*(GQC_results(:,negloglikeidx)+r)';

fprintf('\n\nSEARCH RESULTS:\n');

fprintf('\nPerceptual Noise\n');
disp(final_params(1));

fprintf('\nCriterial Noise\n');
disp(final_params(2));

fprintf('\nBoundary Parameters\n');
disp(final_params(3:12));

fprintf('\nNegative Loglikelihood\n');
disp(neglikelihood);

fprintf('\nAIC score\n');
disp(aicGQC)


echo on;

% Created by Leola A. Alfonso-Reese / 13-March-97
% Copyright (c) 1997

echo off;

