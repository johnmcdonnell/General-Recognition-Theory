function negloglike = negloglike_2dGQC(params,data,z_limit)
%negloglike = negloglike_2dGQC(params,data,z_limit)
%  returns the negative loglikelihood of the 2d data for the
%  General Quadratic Classifier.
%  Note:  These computations are based on the likelihood model
%
%  Parameters:
%    params format:  [pnoise cnoise Amat_entries(1:3) b_entries(1:2) c_bias]
%    data row format:  [correct_response x y 1]
%    z_limit is the z-score value beyond which one should truncate

% Created by Leola A. Alfonso-Reese / 25-March-96
% Copyright (c) 1996
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  3/10/05         Modified to work with Optimization Toolbox v3.0
%                  by LAR

pnoise = params(1);
cnoise = params(2);
A = [   params(3) .5*params(5)
     .5*params(5)    params(4)];
b = params(6:7)';
c = params(8);

% Compute mean of h(x) for each x (as shown in Eqn. 11, Ashby p.462).
temp1 = trace(A*pnoise);
xAxs = params(3)*data(:,2).^2 + params(4)*data(:,3).^2 + ( params(5)*prod(data(:,2:3)') )';
meanhxs = temp1 + xAxs + data(:,2:3)*b + c;

% Compute variance of h(x) for each x (as shown in Eqn. 12, Ashby p.462).
% Take advantage of fact that perceptual noise matrix is diagonal with
% pnoise the same in all three dimensions
ndatapts = size(data,1);
bvec = [b(1)*ones(1,ndatapts);b(2)*ones(1,ndatapts)];
temp2 = bvec+2*A*data(:,2:3)';
varhxs = 2*temp1*temp1 + pnoise*(sum(temp2.^2))';

% Compute z-scores for each data point
zscores = meanhxs./(sqrt(varhxs+cnoise));
   
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
end
if ~isempty(B_indices)
	log_B_probs = log( fastnormalcdf(zscores(B_indices)) );
end

% Sum them up and return the negative
negloglike = -(sum(log_A_probs)+sum(log_B_probs));

   
   
   


