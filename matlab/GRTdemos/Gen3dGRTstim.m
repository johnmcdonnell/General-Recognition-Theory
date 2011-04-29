%Gen3dGRTstim
%
%  generates 250 points from a Trivariate Normal Distribution with mean
%  uA and covariance KA.  Then 250 more points are generated with mean
%  uB and covariance KB = KA.

% Created by Leola A. Alfonso-Reese / 13-March-97
% Copyright (c) 1997
% $Revisions$
%   Date           Modification and Name
%   ----           ---------------------
%  3/13/97         by Leola A. Alfonso-Reese
%                  made changes so that it fits the mold
%                  of other routines in the GRT toolbox
%  12/14/04        by Leola A. Alfonso-Reese
%                  modified call to plot3d bound routines
%                  to match changes in the plot functions

format compact;
echo on;

clc;
% This demo generates 250 points from a Bivariate Normal Distribution with
% mean uA and covariance KA.  Then 250 more points are generated with mean
% uB and covariance KB = KA.  This program makes sure that the points all
% fall within a specified clipping range.
%
% The row format of the resulting data files is:  [category bar1 bar2 bar3]
%
% After one stimulus set is generated, it is randomized twice and saved
% each time into a data file so that the experiment may be run for two
% days.  The resulting data files are 'day1_2d.dat' and 'day2_2d.dat'.
% Finally, the set of stimulus points and decision boundaries are displayed.

pause  % Press any key to continue after pauses.

echo off;

clc;
% Generate 'n' points from a Trivariate Normal Distribution
n = 250;  % number of points to be generated for each category

% Category parameter values and boundary specs
% Note:  These parameters are set when you run DesignExp.m.
% If you have just run DesignExp, and you would rather use those
% specifications, you may comment these out.
uA = [200
      200
			200];
uB = [250
      250
			100];
KA = [23^2   0       0
      0     23^2    0
			0      0     23^2];
KB = [ 23^2      0     23*23*.8
         0      15^2      0
			23*23*.8   0     23^2  ];

% Get optimal quadratic bound, a reasonable linear bound,
% and the best rule that ignores dimension X
[A b c] = quaddecisbnd(KA,KB,uA,uB);
optquadbnd = [diag(A)' 2*A(1,2:3) 2*A(2,3) b' c];

Kpooled = ((n-1)*KA+(n-1)*KB)/(n+n-2);
[a b] = lindecisbnd(Kpooled,uA,uB);
genlinbnd = [a b];

newKA = getcovs(3,KA);
newKB = getcovs(3,KB);
[newuA,newuB] = getmeans(3,uA,uB);
newKpooled = ((n-1)*newKA+(n-1)*newKB)/(n+n-2);
[a b] = lindecisbnd(newKpooled,newuA,newuB);
ignlinbnd = [a(1) 0 a(2) b];

% Pick number of standard deviations for clipping range
numstd = 4;

% GENERATE CATEGORY A AND B SAMPLE
fprintf('...Generating data for Category A\n');
Xa = gensample(uA,KA,n,1);							% row format:  [1 xval yval zval]
fprintf('...Generating data for Category B\n');
Xb = gensample(uB,KB,n,2);							% row format:  [2 xval yval zval]]

% Transform the sample so that the sample means and covariance
% matrices are equal to the population means and covariance matrices.
fprintf('\n...Transforming sample\n');
Xat = transample(Xa(:,2:4),uA,KA);
Xbt = transample(Xb(:,2:4),uB,KB);

% Clip values outside of range in both dimensions for category A.
fprintf('\n...Clipping values outside of range for category A\n');
rangeAx = [uA(1)-numstd*sqrt(KA(1,1)) uA(1)+numstd*sqrt(KA(1,1))]
rangeAy = [uA(2)-numstd*sqrt(KA(2,2)) uA(2)+numstd*sqrt(KA(2,2))]
rangeAz = [uA(3)-numstd*sqrt(KA(3,3)) uA(3)+numstd*sqrt(KA(3,3))]
Xat = clip_vals(Xat,[rangeAx;rangeAy;rangeAz]);

% Final category A data points
Ya = [Xa(:,1) Xat(:,1) Xat(:,2) Xat(:,3)];

% Clip values outside of range in both dimensions for category B.
fprintf('...Clipping values outside of range for category B\n');
rangeBx = [uB(1)-numstd*sqrt(KB(1,1)) uB(1)+numstd*sqrt(KB(1,1))]
rangeBy = [uB(2)-numstd*sqrt(KB(2,2)) uB(2)+numstd*sqrt(KB(2,2))]
rangeBz = [uB(3)-numstd*sqrt(KB(3,3)) uB(3)+numstd*sqrt(KB(3,3))]
Xbt = clip_vals(Xbt,[rangeBx;rangeBy;rangeBz]);

% Final category B data points
Yb = [Xb(:,1) Xbt(:,1) Xbt(:,2) Xbt(:,3)];

% Report the following:
%	1) minimum and maximum values of X, Y and Z for both categories
%	2) standard deviations of X, Y and Z for Category A and Category B

fprintf('\n...Reporting statistics\n');
minA = min(Ya(:,2:4));
disp('Min x, y and z from category A');
disp(minA);

maxA = max(Ya(:,2:3));
disp('Max x, y and z from category A');
disp(maxA);

disp('Sample mean of category A:');
disp(mean(Ya(:,2:4)));

disp('Sample covariance of category A:');
disp(cov(Ya(:,2:4)));

minB = min(Yb(:,2:4));
disp('Min x, y and z from category B');
disp(minB);

maxB = max(Yb(:,2:4));
disp('Max x, y and z from category B');
disp(maxB);

disp('Sample mean of category B:');
disp(mean(Yb(:,2:4)));

disp('Sample covariance of category B:');
disp(cov(Yb(:,2:4)));

% Randomize the order of the stimuli and
% save the generated data points into a file.
fprintf('\n...Randomizing order of stimuli\n');
rand('seed',sum(100*clock));
sample1 = randrows([Ya;Yb]);

fprintf('\n...Saving data to a file called day1_3d.dat\n');
echo on;
%
% NOTE:  The following code that writes data to a file is
% commented out so that data files called 'day1_2d.dat' and
% 'day2_3d.dat' are not created each time someone runs this
% demo.  Sample files that would be created if this code was
% run, 'day1_3d.dat' and 'day2_3d.dat', are in the GRTdemos
% directory in the GRT Toolbox.
%
%fid = fopen('day1_3d.dat','w');
%fprintf(fid,'%i %f %f %f\n',sample1');
%fclose(fid);
echo off;

fprintf('\n...Randomizing order of stimuli for day 2\n');
sample2 = randrows([Ya;Yb]);
fprintf('\n...Saving data to a file called day2_3d.dat\n');
echo on;
%
% See NOTE above.
%
%fid = fopen('day2_3d.dat','w');
%fprintf(fid,'%i %f %f %f\n',sample2');
%fclose(fid);
echo off;

% Display the stimuli
fprintf('\n...Displaying the stimuli');
minX = min([minA,minB]);
maxX = max([maxA,maxB]);
minY = minX;
maxY = maxX;
minZ = minX;
maxZ = maxX;
xyzaxes = [minX-50 maxX+50 minY-50 maxY+50 minZ-50 maxZ+50];
range = [xyzaxes(1)  xyzaxes(2); xyzaxes(3) xyzaxes(4); xyzaxes(5) xyzaxes(6)];

plot3dstim(sample1,xyzaxes,1);
view(110,35);
grid;
echo on;

pause  % Press any key to continue after pauses.

echo off;
%fprintf('\n...Displaying the optimal bound, a good linear');
%fprintf('\nbound, and a bound that ignores the Y dimension.\n');
fprintf('\n...Displaying the optimal bound and a good linear bound.');
hold on;
plot3dquadbnd(optquadbnd,range,'r-');
plot3dlinbnd(genlinbnd,range,'m-');
%plot3dlinbnd(ignlinbnd,range,'m-');
title('Stimuli');
xlabel('X');
ylabel('Y');
zlabel('Z');

hold off;

echo on;

% Created by Leola A. Alfonso-Reese / 13-March-97
% Copyright (c) 1997

echo off;

