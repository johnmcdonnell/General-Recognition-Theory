%Gen2dGRTstim
%
%  generates 250 points from a Bivariate Normal Distribution with mean
%  uA and covariance KA.  Then 250 more points are generated with mean
%  uB and covariance KB = KA.

% Created by Leola A. Alfonso-Reese / 4-December-94
% Copyright (c) 1995
% $Revisions$
%   Date           Modification and Name
%   ----           ---------------------
%  8/21/95         Fixed comment to refer to file day1_2d.dat,day2_2d.dat
%                  by Leola
%  8/31/95         Transform the stimuli before checking the 180 degree
%                  range restriction on angle
%                  by Leola 

format compact;
echo on;

clc;
% This demo generates 250 points from a Bivariate Normal Distribution with
% mean uA and covariance KA.  Then 250 more points are generated with mean
% uB and covariance KB = KA.  This program makes sure that the points all
% fall within a specified clipping range.
%
% The row format of the resulting data files is:  [category length angle]
%
% After one stimulus set is generated, it is randomized twice and saved
% each time into a data file so that the experiment may be run for two
% days.  The resulting data files are 'day1_2d.dat' and 'day2_2d.dat'.
% Finally, the set of stimulus points and decision boundaries are displayed.
%
% Note:
% The generated stimuli are meant to be presented as lines varying in
% length and orientation.  We want orientation to vary within a 65 degree
% range when presented on the screen; so, this program repeatedly generates
% samples until the range of the orientation dimension is less than
% 65 degrees.  Each point value corresponds to pi/500 radians according to
% the way we plan to present the stimuli.  Thus, a range of 65 degrees
% corresponds to a range of approximately 180 points.


pause  % Press any key to continue after pauses.

echo off;

clc;
% Generate 'n' points from a Bivariate Normal Distribution
n = 250;  % number of points to be generated for each category

% Category parameter values and boundary specs
% Note:  These parameters are set when you run DesignExp.m.
% If you have just run DesignExp, and you would rather use those
% specifications, you may comment these out.
uA = [186.6231
      142.2948];
uB = [213.3769
      97.7052];
KA = [25^2  0
      0     25^2];
KB = KA;
optbnd = [.6 -1 0];
optxbnd = [1.0000 0.0000 -200.0000];
optybnd = [0.0000 1.0000 -120.0000];

% Pick number of standard deviations for clipping range
numstd = 4;

% GENERATE CATEGORY A AND B SAMPLE
% Keep generating stimuli until the range of the angle dimension
% is within 180 points.
ok = 0;
while ok == 0
	fprintf('...Generating data for Category A\n');
	Xa = gensample(uA,KA,n,1);							% row format:  [1 xval yval]
	fprintf('...Generating data for Category B\n');
	Xb = gensample(uB,KB,n,2);							% row format:  [2 xval yval]

	% Transform the sample so that the sample means and covariance
	% matrices are equal to the population means and covariance matrices.
	fprintf('\n...Transforming sample\n');
	Xat = transample(Xa(:,2:3),uA,KA);
	Xbt = transample(Xb(:,2:3),uB,KB);

	fprintf('\n...Checking angle range\n');
	minangA = min(Xat(:,2));
	maxangA = max(Xat(:,2));
	minangB = min(Xbt(:,2));
	maxangB = max(Xbt(:,2));
	minang = min([minangA,minangB]);
	maxang = max([maxangA,maxangB]);
	if maxang-minang < 180
		ok = 1;
	end
	fprintf('angle range = %6.4f\n', maxang-minang);
end
fprintf('Final range of angle dimension = %6.4f\n', maxang-minang);

% Clip values outside of range in both dimensions for category A.
fprintf('\n...Clipping values outside of range for category A\n');
rangeAx = [uA(1)-numstd*sqrt(KA(1,1)) uA(1)+numstd*sqrt(KA(1,1))]
rangeAy = [uA(2)-numstd*sqrt(KA(2,2)) uA(2)+numstd*sqrt(KA(2,2))]
Xat = clip_vals(Xat,[rangeAx;rangeAy]);

% Final category A data points
Ya = [Xa(:,1) Xat(:,1) Xat(:,2)];

% Clip values outside of range in both dimensions for category B.
fprintf('...Clipping values outside of range for category B\n');
rangeBx = [uB(1)-numstd*sqrt(KB(1,1)) uB(1)+numstd*sqrt(KB(1,1))]
rangeBy = [uB(2)-numstd*sqrt(KB(2,2)) uB(2)+numstd*sqrt(KB(2,2))]
Xbt = clip_vals(Xbt,[rangeBx;rangeBy]);

% Final category B data points
Yb = [Xb(:,1) Xbt(:,1) Xbt(:,2)];

% Report the following:
%	1) minimum and maximum values of X, Y and Z for both categories
%	2) standard deviations of X, Y and Z for Category A and Category B

fprintf('\n...Reporting statistics\n');
minA = min(Ya(:,2:3));
disp('Min x, and y from category A');
disp(minA);

maxA = max(Ya(:,2:3));
disp('Max x, and y from category A');
disp(maxA);

disp('Sample mean of category A:');
disp(mean(Ya(:,2:3)));

disp('Sample covariance of category A:');
disp(cov(Ya(:,2:3)));

minB = min(Yb(:,2:3));
disp('Min x, and y from category B');
disp(minB);

maxB = max(Yb(:,2:3));
disp('Max x, and y from category B');
disp(maxB);

disp('Sample mean of category B:');
disp(mean(Yb(:,2:3)));

disp('Sample covariance of category B:');
disp(cov(Yb(:,2:3)));

% Randomize the order of the stimuli and
% save the generated data points into a file.
fprintf('\n...Randomizing order of stimuli\n');
rand('seed',sum(100*clock));
sample1 = randrows([Ya;Yb]);

fprintf('\n...Saving data to a file called day1_2d.dat\n');
echo on;
%
% NOTE:  The following code that writes data to a file is
% commented out so that data files called 'day1_2d.dat' and
% 'day2_2d.dat' are not created each time someone runs this
% demo.  Sample files that would be created if this code was
% run, 'day1_2d.dat' and 'day2_2d.dat', are in the GRTdemos
% directory in the GRT Toolbox.
%
%fid = fopen('day1_2d.dat','w');
%fprintf(fid,'%i %f %f \n',sample1');
%fclose(fid);
echo off;

fprintf('\n...Randomizing order of stimuli for day 2\n');
sample2 = randrows([Ya;Yb]);
fprintf('\n...Saving data to a file called day2_2d.dat\n');
echo on;
%
% See NOTE above.
%
%fid = fopen('day2_2d.dat','w');
%fprintf(fid,'%i %f %f \n',sample2');
%fclose(fid);
echo off;

% Display the stimuli
fprintf('\n...Displaying the stimuli, optimal bound,');
fprintf('\nand best dimensional bounds.\n');
minX = min([minA,minB]);
maxX = max([maxA,maxB]);
minY = minX;
maxY = maxX;
xyaxes = [minX  maxX minY maxY];

plot2dstim(sample1,xyaxes,1);
hold on;
plot2dlinbnd(optbnd,'r-',xyaxes);
plot2dlinbnd(optxbnd,'r--',xyaxes);
plot2dlinbnd(optybnd,'r--',xyaxes);
title('Stimuli');
hold off;


echo on;

% Created by Leola A. Alfonso-Reese / 4-December-94
% Copyright (c) 1995

echo off;

