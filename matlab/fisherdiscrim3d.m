function fishers_coeffs = fisherdiscrim3d(data,modelnum)
%fishers_coeffs = fisherdiscrim3d(data,modelnum)
%  computes coefficients of Fisher's Discriminant function
%  for the model indicated.
%
%  Warnings:
%    a) The data must contain at least 4 data points.
%    b) If all the data points are assigned a single
%       category label, then the returned coefficients
%       are arbitrary.
%
%  Parameters:
%    row format of data:  [category x y z response]
%    modelnum may be an integer from 1 to 10 indicating...
%     01. estimate best 2-D linear classifier
%         based on response means and covariance matrix
%         ignoring var and mean of dimension X
%         (i.e., detect yz covariances)
%     02. estimate best 2-D linear classifier
%         based on response means and covariance matrix
%         ignoring var and mean of dimension Y
%         (i.e., detect xz covariances)
%     03. estimate best 2-D linear classifier
%         based on response means and covariance matrix
%         ignoring var and mean of dimension Z
%         (i.e., detect xy covariances)
%     04. estimate best 3-D linear classifier
%         based on response means and covariance matrix
%         using a diagonal covariance matrix with one common
%         variance (i.e., no correlations)
%     05. estimate best 3-D linear classifier
%         based on response means and covariance matrix
%         using a diagonal covariance matrix with three different
%         variances (i.e., no correlations)
%     06. estimate best 3-D linear classifier
%         based on response means and covariance matrix
%         using a diagonal covariance matrix with three different
%         variances and a nonzero covxy (i.e., other covariances
%         are zero)
%     07. estimate best 3-D linear classifier
%         based on response means and covariance matrix
%         using a diagonal covariance matrix with three different
%         variances and a nonzero covyz (i.e., other covariances
%         are zero)
%     08. estimate best 3-D linear classifier
%         based on response means and covariance matrix
%         using a diagonal covariance matrix with three different
%         variances and a nonzero covzx (i.e., other covariances
%         are zero)
%     09. estimate the optimal linear classifier
%         based on sample means and covariance matrix
%     10. estimate minimum distance classifier
%         based on sample means and covariance matrix
%     11. estimate the optimal linear classifier
%         based on response means and response covariance matrix
%
%     21. estimate optimal quadratic classifier
%         based on sample means and sample covariance matrices
%     22. estimate optimal quadratic classifier
%         based on response means and response covariance matrices
%     23. estimate sphere around category A that bisects midpoints
%         based on sample means
%     24. estimate sphere around category B that bisects midpoints
%         based on sample means
%     25. estimate sphere around category A that bisects midpoints
%         based on response means
%     26. estimate sphere around category B that bisects midpoints
%         based on response means
%
%  Note that models 1-8, 11, 22 and 23 are based on the response categories
%  while models 9, 10, 21 and 24 are based on the sample categories

% Created by Leola A. Alfonso-Reese / 24-February-94
% Copyright (c) 1994
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  12-Mar-97       by Leola A. Alfonso-Reese
%                  made minor changes so that it fits the mold
%                  of other routines in the GRT toolbox
%  6-Sept-99       by Leola Alfonso-Reese
%                  added models 23-26
%  13_Sept-02      by Alfonso-Reese & Spiering
%                  Added check for full rank cov matrix
%                  If not full, add noise to data.

if (modelnum == 9 | modelnum == 10 | modelnum == 21 | modelnum == 23 | modelnum == 24) 
	Aindices = find(data(:,1) == 1);
	Alength = length(Aindices);
	Bindices = find(data(:,1) ~= 1);
	Blength = length(Bindices);
elseif ( (modelnum >= 1 & modelnum <= 8) | modelnum == 11 | modelnum == 22 | modelnum == 25 | modelnum == 26) 
	Aindices = find(data(:,5) == 1);
	Alength = length(Aindices);
	Bindices = find(data(:,5) ~= 1);
	Blength = length(Bindices);
else
	error('Model number not available.');
end

% Pick arbitary bound if all responses are assigned one category label
if (Alength == 0 | Blength == 0)
	% take half of the points to be A's, the other half B's
	mididx = floor(size(data,1)/2);
	Aindices = 1:mididx;
  Alength = length(Aindices);
	Bindices = mididx+1:size(data,1);
  Blength = length(Bindices);
end

meanA = mean(data(Aindices,2:4))';
meanB = mean(data(Bindices,2:4))';
KA = cov(data(Aindices,2:4));
KB = cov(data(Bindices,2:4));
Kpooled = ((Alength-1)*KA+(Blength-1)*KB)/(Alength+Blength-2);

[Kest,Kact] = getcovs(modelnum,Kpooled);
[meanAest,meanBest] = getmeans(modelnum,meanA,meanB);

if modelnum < 21
	% Decision boundary for the linear model is
	% h(x,y,z) = a1*x + a2*y + a3*z + b
	if rank(Kest) == size(Kest,1)
		[a,b] = lindecisbnd(Kest,meanAest,meanBest);
	else
		% Add noise to the data and calculate new estimates
		noise = randn(size(data,1),3).*10;
		noisy_data = data(:,2:4)+noise;
		meanA = mean(noisy_data(Aindices,:))';
		meanB = mean(noisy_data(Bindices,:))';
		KA = cov(noisy_data(Aindices,:));
		KB = cov(noisy_data(Bindices,:));
		Kpooled = ((Alength-1)*KA+(Blength-1)*KB)/(Alength+Blength-2);
		
		[Kest,Kact] = getcovs(modelnum,Kpooled);
		[meanAest,meanBest] = getmeans(modelnum,meanA,meanB);

		[a,b] = lindecisbnd(Kest,meanAest,meanBest);

	end
elseif (modelnum == 23 | modelnum == 25)
	% Decision boundary for spherical model around A
	[A,b,c] = sphdecisbnd(meanA,meanB);
elseif (modelnum == 24 | modelnum == 26)
	% Decision boundary for spherical model around B
	[A,b,c] = sphdecisbnd(meanB,meanA);
else
	% Decision boundary for the quadratic model is
	% h(x,y,z) = x'Ax + b'x + c
	% where A is a 3x3 matrix, b is a 3x1 vector and c is a constant.
	[A,b,c] = quaddecisbnd(KA,KB,meanA,meanB);
end

if modelnum == 1
   % ignore dimension X
   fishers_coeffs = [0, a(1), a(2), b];
elseif modelnum == 2
   % ignore dimension Y
   fishers_coeffs = [a(1), 0, a(2), b];
elseif modelnum == 3
   % ignore dimension Z
   fishers_coeffs = [a(1), a(2), 0, b];
elseif modelnum == 10
   % the covariance matrix is the identity
   fishers_coeffs = [a, b];
elseif modelnum < 21
   % use all three dimensions
   fishers_coeffs = [a, b];
else % modelnum == 21..26
	fishers_coeffs = [diag(A)' 2*A(1,2:3) 2*A(2,3) b' c];
end

