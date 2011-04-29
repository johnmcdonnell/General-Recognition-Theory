function negloglike = negloglike_1dGLC(params,data,z_limit)
%negloglike = negloglike_1dGLC(params,data,z_limit)
%  returns the negative loglikelihood of the 1d data for the
%  General Linear Classifier.
%
%  Parameters:
%    params format:  [noise a bias]
%    data row format:  [response x 1]
%    z_limit is the z-score value beyond which one should truncate

% Created by Leola A. Alfonso-Reese / 23-August-04
% Copyright (c) 2004
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  3/10/05         Modified to work with Optimization Toolbox v3.0
%                  by LAR

z_coefs = [params(2) params(3)]' ./ params(1);
%params

% Compute z-scores for each data point
zscores = data(:,2:3) * z_coefs;
   
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
