function d = dprimef(uA,uB,covmat,noise)
%d = dprimef(uA,uB,covmat)
%  computes d-prime from the population parameters.
%
%  D-prime is the standardized distance between means
%  of two distributions which are assumed to be normal.
%
%  Parameters:
%    uA is a column vector containing the population mean
%       of category A
%    uB is a column vector containing the population mean
%       of category B
%    covmat is the covariance matrix for each category
%    noise represents perceptual and criterial noise (expressed
%          as stdev).  This is an optional parameter.

% Created by Leola A. Alfonso-Reese / 3-March-94
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%   6-Sept-99      by Leola Alfonso-Reese
%                  Added noise parameter

if nargin > 3
	if isempty(noise)		
		 error('Input argument noise must be specified.');
	else
		if noise >= 0
			covmat = covmat +  diag(noise^2*(ones(length(uA),1)));
		else
			error('Noise cannot be negative.');
		end
	end
end
d = sqrt((uA - uB)'*inv(covmat)*(uA-uB));
	

