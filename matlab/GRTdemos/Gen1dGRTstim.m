%Gen1dGRTstim
%
%  generates 250 points from a Univariate Normal Distribution with mean
%  uA and variance KA.  Then 250 more points are generated with mean
%  uB and variance KB = KA.

% Created by Leola A. Alfonso-Reese / 10-March-05
% Copyright (c) 2005
% $Revisions$
%   Date           Modification and Name
%   ----           ---------------------

format compact;
echo on;

clc;
% This demo generates 250 points from a Univariate Normal Distribution with
% mean uA and variance KA.  Then 250 more points are generated with mean
% uB and covariance KB = KA.  This program makes sure that the points all
% fall within a specified clipping range.
%
% The row format of the resulting data files is:  [category length]
%
% After one stimulus set is generated, it is randomized twice and saved
% each time into a data file so that the experiment may be run for two
% days.  The resulting data files are 'day1_1d.dat' and 'day2_1d.dat'.
% Finally, the set of stimulus points and decision boundaries are displayed.
%
% Note:
% The generated stimuli are meant to be presented as lines varying in
% length.  

pause  % Press any key to continue after pauses.

echo off;

clc;
% Generate 'n' points from a Univariate Normal Distribution
n = 250;  % number of points to be generated for each category

% Category parameter values and boundary specs
uA = [247];
uB = [153];
KA = [40^2];
KB = KA;
optbnd = [-1 200];

% Pick number of standard deviations for clipping range
numstd = 4;

% GENERATE CATEGORY A AND B SAMPLE
% Keep generating stimuli until the range of the angle dimension
% is within 180 points.
ok = 0;
while ok == 0
	fprintf('...Generating data for Category A\n');
	Xa = gensample(uA,KA,n,1);							% row format:  [1 xval]
	fprintf('...Generating data for Category B\n');
	Xb = gensample(uB,KB,n,2);							% row format:  [2 xval]

	% Transform the sample so that the sample means and covariance
	% matrices are equal to the population means and covariance matrices.
	fprintf('\n...Transforming sample\n');
	Xat = transample(Xa(:,2),uA,KA);
	Xbt = transample(Xb(:,2),uB,KB);

	fprintf('\n...Checking Dim1 range\n');
	minyA = min(Xat(:,1));
	maxyA = max(Xat(:,1));
	minyB = min(Xbt(:,1));
	maxyB = max(Xbt(:,1));
	miny = min([minyA,minyB]);
	maxy = max([maxyA,maxyB]);
	if and(miny>=0,maxy<401)
		ok = 1;
	end
	fprintf('Dim1 range = %6.4f\n', maxy-miny);
end
fprintf('Final range of dimension 1 = %6.4f\n', maxy-miny);

% Clip values outside of range in both dimensions for category A.
fprintf('\n...Clipping values outside of range for category A\n');
rangeA = [uA-numstd*sqrt(KA) uA+numstd*sqrt(KA)]
Xat = clip_vals(Xat,[max(rangeA(1),0) min(rangeA(2),400)]);

% Final category A data points
Ya = [Xa(:,1) Xat(:,1)];

% Clip values outside of range in both dimensions for category B.
fprintf('...Clipping values outside of range for category B\n');
rangeB = [uB-numstd*sqrt(KB) uB+numstd*sqrt(KB)]
Xbt = clip_vals(Xbt,[max(rangeB(1),0) min(rangeB(2),400)]);

% Final category B data points
Yb = [Xb(:,1) Xbt(:,1)];

% Report the following:
%	1) minimum and maximum values of X, Y and Z for both categories
%	2) standard deviations of X, Y and Z for Category A and Category B

fprintf('\n...Reporting statistics\n');

minA = min(Ya(:,2));
disp('Min x from category A');
disp(minA);

maxA = max(Ya(:,2));
disp('Max x category A');
disp(maxA);

disp('Sample mean of category A:');
disp(mean(Ya(:,2)));
disp('Sample standard dev of category A:');
disp(std(Ya(:,2)));

minB = min(Ya(:,2));
disp('Min x from category B');
disp(minB);

maxB = max(Yb(:,2));
disp('Max x from category B');
disp(maxB);

disp('Sample mean of category B:');
disp(mean(Yb(:,2)));
disp('Sample standard dev of category B:');
disp(std(Yb(:,2)));

% Randomize the order of the stimuli and
% save the generated data points into a file.
fprintf('\n...Randomizing order of stimuli\n');
rand('seed',sum(100*clock));
sample1 = randrows([Ya;Yb]);

fprintf('\n...Saving data to a file called day1_1d.dat\n');
echo on;
%
% NOTE:  The following code that writes data to a file is
% commented out so that data files called 'day1_1d.dat' and
% 'day2_1d.dat' are not created each time someone runs this
% demo.  Sample files that would be created if this code was
% run, 'day1_1d.dat' and 'day2_1d.dat', are in the GRTdemos
% directory in the GRT Toolbox.
%
%fid = fopen('day1_1d.dat','w');
%fprintf(fid,'%i %f \n',sample1');
%fclose(fid);
echo off;

fprintf('\n...Randomizing order of stimuli for day 2\n');
sample2 = randrows([Ya;Yb]);
fprintf('\n...Saving data to a file called day2_1d.dat\n');
echo on;
%
% See NOTE above.
%
%fid = fopen('day2_1d.dat','w');
%fprintf(fid,'%i %f \n',sample2');
%fclose(fid);
echo off;

% Display the stimuli
fprintf('\n...Displaying the stimuli, optimal bound,');
fprintf('\nand best dimensional bounds.\n');
minX = min([minA,minB]);
maxX = max([maxA,maxB]);
minY = minX;
maxY = maxX;
xaxes = [0  400];

plot1dstim(sample1,xaxes,1);
hold on;
plot1dlinbnd(optbnd,'r-',xaxes);
title('Stimuli');
hold off;


echo on;

% Created by Leola A. Alfonso-Reese / 10-March-05
% Copyright (c) 2005

echo off;

