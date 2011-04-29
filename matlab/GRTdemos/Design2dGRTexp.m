%Design2dGRTexp
%
%  finds the category population means for an experiment that requires ...
%
%    1. two normal category distributions with equal covariance matrices,
%    2. an optimal linear bound [a1 a2 b] where a1*x + a2*y + b = 0, and
%    3. optimal performance in terms of probability correct specified by pcorr.

% Created by Leola A. Alfonso-Reese / 17-January-95
% Copyright (c) 1995
% $Revisions$
%   Date           Modification and Name
%   ----           ---------------------


format compact;
echo on;

clc;
% This demo finds the category population means for an experiment that
% requires ...
%
%   1. two normal category distributions with equal covariance matrices,
%   2. an optimal linear bound [a1 a2 b] where a1*x + a2*y + b = 0, and
%   3. a specified optimal performance in terms of probability correct.
%
% To find the category means, this program takes the midpoint of the means
% (i.e., some point on the bound) and increments the distance of the
% category means to the midpoint until the desired probability correct is
% obtained.

pause  % Press any key to continue after pauses.


clc;
% Experiment design goals:

K = [25^2 0              % covariance matrix for both categories
     0    25^2];
		 
optbnd = [.6 -1 0];      % optimal boundary [a1 a2 b] such that
                         %   ax + by + c = 0
													
midpt = [200 200*.6]';   % point on bound between A and B prototype

pcorrgoal = .85;         % desired probability correct
echo off;

% Initialize variables
dist_to_u = 5;           % initial distance from the boundary to the category A mean
stepsize = 1;            % size of step to take while searching for
                         %   the best distance from the bound

% Find means for categories A and B that yield 90% correct
%   1. Start with a point dist_to_u values away from the bound and
%      towards the category A mean.
%   2. Increase the distance of the point from the bound (by the amount stepsize)
%      until the maximum probability correct using the specified means and covariance
%      matrix is 90% correct.
echo off;
fprintf('\nSearching for category means...\n');
normal_to_bnd = optbnd(1:2)'/norm(optbnd(1:2));	% normal vector (of length 1) to the bound
pcorr = 0;
dist_to_u = dist_to_u - stepsize;
while pcorr < pcorrgoal
  dist_to_u = dist_to_u + stepsize;
  uA = midpt-dist_to_u*normal_to_bnd;
  uB = midpt+dist_to_u*normal_to_bnd;
  pcorr = linprobcorr(K,uA,uB);
  fprintf('Probability correct = %6.4f\n', pcorr);	
end


% Final values for this experimental design
fprintf('\nFinal values\n');
disp('Category A mean:');
disp(uA);
disp('Category B mean:');
disp(uB);
disp('Category A covariance matrix:');
KA = K;
disp(KA);
disp('Category B covariance matrix:');
KB = K;
disp(KB);

% Find probability correct of the optimal dimensional rules
pcorrx = linprobcorr(K(1,1),uA(1),uB(1));
fprintf('\nProbability correct using dimensional rule in x = %6.4f\n',pcorrx);
[a,b] = lindecisbnd(K(1,1),uA(1),uB(1));
optxbnd = [a,b]/a;
fprintf('Dimensional bound = [%6.4f %6.4f %6.4f]\n',[optxbnd(1) 0 optxbnd(2)]);

pcorry = linprobcorr(K(2,2),uA(2),uB(2));
fprintf('\nProbability correct using dimensional rule in y = %6.4f\n',pcorry);
[a,b] = lindecisbnd(K(2,2),uA(2),uB(2));
optybnd = [a,b]/a;
fprintf('Dimensional bound = [%6.4f %6.4f %6.4f]\n',[0 optybnd]);


echo on;

% Created by Leola A. Alfonso-Reese / 17-January-95
% Copyright (c) 1995

echo off;

