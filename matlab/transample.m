function newsample = transample(sample,newmean,newcov)
%newsample = transample(sample,newmean,newcov)
%  transforms the data in sample so that it has a mean,
%  newmean, and a covariance matrix, newcov.
%
%  Parameters:
%    row format of sample: [var1 var2 ... varn]
%    newmean is a column vector containing the desired
%       category mean
%    newcov is the desired covariance matrix
%
%  Note: If the sample is one-dimensional, then newmean
%  and newcov are simply scalars, where newcov is the
%  new variance.

% Created by Leola A. Alfonso-Reese / 15-April-94
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%   28-July-95     Fixed bug
%                  by Leola A. Alfonso-Reese
%   4-Sept-99      Made it work for 1-dimensional sample
%                  by Leola A. Alfonso-Reese


% Get sample statistics
oldmean = mean(sample);
oldcov = cov(sample);

if length(newmean) == 1
	% Sample is one-dimensional
	if oldcov ~= 0
		newsample = sqrt(newcov/oldcov)*(sample-oldmean)+newmean;
	else
		newsample = sample-oldmean+newmean;
	end
else
	% Sample is multi-dimensional
	
	% Get Cholesky factorization of the sample and covariance matrices
	% where Aold'*Aold = oldcov and Anew'*Anew = newcov.
	Aold = chol(oldcov)';
	Anew = chol(newcov)';
	
	% Transform Sample
	o = ones(size(sample,1),1);
	oldmeanvec = o*oldmean;
	newmeanvec = o*newmean';
	newsample = (Anew*inv(Aold)*(sample-oldmeanvec)' + newmeanvec')';
end

