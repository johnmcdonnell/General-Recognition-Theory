%Comp2dStats
%
%  finds statistics of interest for the GRT data set, 'subjdemo_2d.dat'.
%  Specifically, this program computes percent correct, dprime, and
%  efficiency.

% Created by Leola A. Alfonso-Reese / 18-April-95
% Copyright (c) 1995
% $Revisions$
%   Date           Modification and Name
%   ----           ---------------------

format compact;
echo on;
clc;
% This demo finds statistics of interest for the GRT data set,
% 'subjdemo_2d.dat'.  Specifically, it computes percent correct,
% dprime, and efficiency.

pause  % Press any key to continue.
echo off;

% These are the population parameters (means, and common
% covariance matrix) for the experiment.
muA = [186.6231
       142.2948];
muB = [213.3769
       97.7052];
covmat = [25^2  0
          0    25^2];

dprimeidealpop = dprimef(muA,muB,covmat);

load subjdemo_2d.dat
data = subjdemo_2d;

% Find percent correct
pc = percorr(data,2);

% Find d-primes and efficiency
dprimeidealsample = dprime(data(:,1:4),'SampleIdeal');
dprimeobs = dprime(data(:,1:4),'Observer');
efficiencypop = (dprimeobs/dprimeidealpop)^2;
efficiencysample = (dprimeobs/dprimeidealsample)^2;

fprintf('\npercent correct = %f\n',pc);
fprintf('d prime ideal (pop.) = %f\n',dprimeidealpop);
fprintf('d prime ideal (sample) = %f\n',dprimeidealsample);
fprintf('d prime observer = %f\n',dprimeobs);
fprintf('efficiency (using pop. d'') = %f\n',efficiencypop);
fprintf('efficiency (using sample d'') = %f\n',efficiencysample);


echo on;

% Created by Leola A. Alfonso-Reese / 18-April-95
% Copyright (c) 1995

echo off;

