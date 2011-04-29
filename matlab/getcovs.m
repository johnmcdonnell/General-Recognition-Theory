function [Kest,Kideal] = getcovs(modelnum,K)
%[Kest,Kideal] = getcovs(modelnum,K)
%  returns the estimated and ideal covariance matrices
%  'Kest' and 'Kideal' that correspond with the model
%  indicated by 'modelnum'.
%
%  Parameters:
%    modelnum may be an integer from 1 to 10 (See fisherdiscrim3d)
%    K is a covariance matrix, i.e., the presumed structure is
%       K = [sdx^2   covxy   covxz
%            covxy   sdy^2   covyz
%            covxz  covyz    sdz^2]

% Created by Leola A. Alfonso-Reese / 10-March-94
% Copyright (c) 1994
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  3/12/97         by Leola A. Alfonso-Reese
%                  made minor changes so that it fits the mold
%                  of other routines in the GRT toolbox

if modelnum == 1
   % ignore dimension X
   Kest = reducemat(K,2,3);
   Kideal = Kest;
elseif modelnum == 2
   % ignore dimension Y
   Kest = reducemat(K,1,3);
   Kideal = Kest;
elseif modelnum == 3
   % ignore dimension Z
   Kest = reducemat(K,1,2);
   Kideal = Kest;
elseif modelnum == 4
   % use one common variance
   comvar = (K(1,1)+K(2,2)+K(3,3))/3;
   Kest = [comvar 0      0
           0      comvar 0
           0      0      comvar];
   Kideal = K;
elseif modelnum == 5
   % use all three variances
   Kest = [K(1,1) 0      0
           0      K(2,2) 0
           0      0      K(3,3)];
   Kideal = K;
elseif modelnum == 6
   % all variances and covxy
   Kest = [K(1,1) K(1,2) 0
           K(2,1) K(2,2) 0
           0      0      K(3,3)];
   Kideal = K;
elseif modelnum == 7
   % all variances and covyz
   Kest = [K(1,1) 0      0
           0      K(2,2) K(2,3)
           0      K(3,2) K(3,3)];
   Kideal = K;
elseif modelnum == 8
   % all variances and covzx
   Kest = [K(1,1) 0      K(1,3)
           0      K(2,2) 0
           K(3,1) 0      K(3,3)];
   Kideal = K;
elseif modelnum == 10
   % K is the identity matrix
   Kest = [1 0 0
           0 1 0
           0 0 1];
   Kideal = K;
else
   % all variances and covariances
   Kest = K;
   Kideal = K;
end

