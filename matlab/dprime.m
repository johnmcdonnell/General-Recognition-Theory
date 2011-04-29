function d = dprime(data,kind,noise,bnd)
%d = dprime(data,kind,noise,bnd)
%  computes d-prime from the data.
%
%  D-prime is the standardized distance between means
%  of two distributions which are assumed to be normal,
%  unless a bound is specified.
%
%  Parameters:
%    data row format:  [category x y ... response]
%    kind  is a string indicating the kind of dprime desired
%          Valid strings are 'SampleIdeal' or 'Observer'.
%    noise represents perceptual and criterial noise (expressed
%          as stdev).  This is an optional parameter.
%    bnd   is the category boundary.  Format is [a1 ... an bias]
%          where a1*x + a2*y ... + bias = 0.  This is an optional
%          parameter.

% Created by Leola A. Alfonso-Reese / 7-March-94
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
% 14-Sept-99       by Leola Alfonso-Reese
%                  put noise in stdev units rather than var

 
% Get number of dimensions
[m,n] = size(data);
dim = n - 2;

if strcmp(kind,'SampleIdeal')

  % Find the A and B stimulus sets
  Aindices = find(data(:,1) == 1);
  Alength = length(Aindices);
  Bindices = find(data(:,1) ~= 1);
  Blength = length(Bindices);

  if nargin > 3
	
	  % Use method based on a specified noise and boundary

		if isempty(noise)		
		  error('Input argument noise must be specified.');
			
		elseif noise > 0
	    z_limit = 7;
	
	    % Compute h(x) for linear boundary
	    bnd_norm = norm(bnd(1:dim));
	    h_coefs = [bnd(1:dim)';bnd(dim+1)]/bnd_norm;
	
	    % Compute h-scores for each data point
	    % and divide each h-score by noise
	    h_vals = -[data(:,2:dim+1) ones(size(data,1),1)] * (h_coefs/noise);
	
	    % Truncate the large z-scores
	    low_indices = find(h_vals < -z_limit);
	    h_vals(low_indices) = -z_limit*ones(length(low_indices),1);  
	    high_indices = find(h_vals > z_limit);
	    h_vals(high_indices) = z_limit*ones(length(high_indices),1);
	
	    % Find hit rate
			hit_freq = fastnormalcdf(h_vals(Aindices));
			hitrate = sum(hit_freq)/Alength;
		
			% Find false alarm rate
			falsealarm_freq = fastnormalcdf(h_vals(Bindices));
			falsealarmrate = sum(falsealarm_freq)/Blength;
	
			% Compute d'.
			d = normalcdfinv(1-falsealarmrate)-normalcdfinv(1-hitrate);
			
		else
		    error('Input argument noise must be greater than zero.');
		end
		
  else
	  % Use method based on sample statistics
		
    % Find the means and variances of A and B stimulus sets
    meanA = mean(data(Aindices,2:dim+1))';
    meanB = mean(data(Bindices,2:dim+1))';
    KA = cov(data(Aindices,2:dim+1));
    KB = cov(data(Bindices,2:dim+1));
    Kpooled = ((Alength-1)*KA+(Blength-1)*KB)/(Alength+Blength-2);

    if (nargin > 2)
		  if noise >= 0
				Kpooled = Kpooled + diag(noise^2*(ones(length(meanA),1)));
			else
		    error('Input argument noise must be greater than zero.');
			end
		end
	
		% Compute the ideal d' based on the sample stimulus set
    d = sqrt((meanA - meanB)'*inv(Kpooled)*(meanA-meanB));
  end	

elseif strcmp(kind,'Observer')
  
  % Find the hit and falsealarm rates
  hits = length( find((data(:,1) == 1)&(data(:,dim+2) == 1)) );
  falsealarms = length( find((data(:,1) == 2)&(data(:,dim+2) == 1)) );
  misses = length( find((data(:,1) == 1)&(data(:,dim+2) == 2)) );
  correjects = length( find((data(:,1) == 2)&(data(:,dim+2) == 2)) );
	
  falsealarmrate = falsealarms/(falsealarms+correjects);
  hitrate = hits/(hits+misses);
	
  % Compute the observer's d'.
  d = normalcdfinv(1-falsealarmrate)-normalcdfinv(1-hitrate);
	
else
  
  error('Invalid string:  kind');
	
end

