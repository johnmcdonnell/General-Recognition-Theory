function negloglike = negloglike_2dGLC(params,data,z_limit)
%negloglike = negloglike_2dGLC(params,data,z_limit)
%  returns the negative loglikelihood of the 2d data for the
%  General Linear Classifier.
%
%  Parameters:
%    params format:  [noise angle bias]
%    data row format:  [response x y 1]
%    z_limit is the z-score value beyond which one should truncate

% Created by Leola A. Alfonso-Reese / 3-March-94
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  4/29/96         Checked for cases of no A or no B responses
%                  by LAR
%  5/4/96          Modified help text
%                  by LAR
%  3/10/05         Modified to work with Optimization Toolbox v3.0
%                  by LAR

% Normalize plane coefficients so that orthogonal
% is expressed in units of standard deviation
% Note:  the distance and likelihood models yield the same
%        sigma parameter estimate because the norm of the xyz
%        coordinates to the bound is 1
xy = angvec2xymat(params(2));
z_coefs = [xy, params(3)]' ./ params(1);

% Compute z-scores for each data point
zscores = data(:,2:4) * z_coefs;
   
% Truncate the large z-scores
low_indices = find(zscores < -z_limit);
zscores(low_indices) = -z_limit*ones(length(low_indices),1);
high_indices = find(zscores > z_limit);
zscores(high_indices) = z_limit*ones(length(high_indices),1);
   
% Find A and B response subsets
A_indices = find( data(:,1) == 1 );
B_indices = find( data(:,1) ~= 1 );

% Find log of cumulative probability
% Note:  Since "h(x) < 0 => cat A," I use 1-CDF for A points
% and CDF for B points.
if ~isempty(A_indices)
	log_A_probs = log( 1-fastnormalcdf(zscores(A_indices)) );
else
	log_A_probs = zeros(length(A_indices),1);
end
if ~isempty(B_indices)
	log_B_probs = log( fastnormalcdf(zscores(B_indices)) );
else
	log_B_probs = zeros(length(B_indices),1);
end

% Sum them up and return the negative
negloglike = -(sum(log_A_probs)+sum(log_B_probs));

end
