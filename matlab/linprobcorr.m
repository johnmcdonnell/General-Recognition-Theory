function prc = linprobcorr(K,uA,uB,noise)
%prc = linprobcorr(K,uA,uB)
%  returns the probability correct using the optimal boundary to
%  categorize the data from categories having means uA, uB and
%  a common covariance matrix, K, and noise.
%
%  Parameters:
%    K  is the common covariance matrix
%    uA is a column vector containing the mean of category A
%    uB is a column vector containing the mean of category B
%    noise is a real number (optional)
%
%  If noise is not specified, it is assumed to be zero.

% Created by Leola A. Alfonso-Reese / 12-April-94
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  20-Nov-99       added noise parameter
%                  by Leola Alfonso-Reese 
%  02-Feb-02       added check for equal means
%                  by Leola Alfonso-Reese 

if nargin < 4
	noise = 0;
end

if uA == uB
	prc = .5;
else
	% Set up for computing probability correct
	invK = inv(K);
	uhx = (uB-uA)'*invK*uA + .5*(uA'*invK*uA-uB'*invK*uB);
	varhx = (uB-uA)'*invK*(uB-uA);
	
	% Correct for unit length discriminant function
	normConstant = norm((uB-uA)'*invK);
	uhx = uhx/normConstant;
	varhx = varhx/(normConstant^2);
	z = (0-uhx)/sqrt(varhx + noise*noise);
	prc = normalcdf(z,0,1);

end

