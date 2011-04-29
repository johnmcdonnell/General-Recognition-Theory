function fishers_coeffs = fisherdiscrim1d(data,modelnum)
%fishers_coeffs = fisherdiscrim2d(data,modelnum)
%  computes coefficients of Fisher's Discriminant function
%  for the model indicated.
%
%  Warnings:
%    a) The data must contain at least 4 data points.
%    b) If all the data points are assigned a single
%       category label, then the returned coefficients
%       are arbitrary.
%
%
%  Parameters:
%    row format of data:  [category x response]
%    modelnum may be an integer from 1 to 5, 13 or 14 indicating...
%     01. estimate best 1-D linear classifier
%         based on sample means and sample covariance matrix
%     03. estimate optimal linear classifier
%         based on sample means and sample covariance matrix
%     05. estimate optimal linear classifier
%         based on response means and response covariance matrix
%
%  Note that model 3 is based on the sample categories
%  while models 1 and 5 are based on the subject's responses.

% Created by Leola A. Alfonso-Reese / 23-August-2004
% Copyright (c) 2004
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------


if (modelnum == 3 )
	% Get data using sample category labels
	Aindices = find(data(:,1) == 1);
	Alength = length(Aindices);
	Bindices = find(data(:,1) ~= 1);
	Blength = length(Bindices);
elseif (modelnum == 1 | modelnum == 5 )
	% Get data using response category labels
	Aindices = find(data(:,3) == 1);
	Alength = length(Aindices);
	Bindices = find(data(:,3) ~= 1);
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

meanA = mean(data(Aindices,2))';
meanB = mean(data(Bindices,2))';
KA = var(data(Aindices,2));
KB = var(data(Bindices,2));
Kpooled = ((Alength-1)*KA+(Blength-1)*KB)/(Alength+Blength-2);

if modelnum == 1
	meanAest = meanA;
	meanBest = meanB;
	Kest = Kpooled;
elseif (modelnum == 3 | modelnum == 5)
	meanAest = meanA;
	meanBest = meanB;
	Kest = Kpooled;
elseif modelnum == 4
	meanAest = meanA;
	meanBest = meanB;
	% K is the identity matrix
	Kest = eye( 2 )
end

if modelnum < 10
	% Decision boundary for this model is h(x,y)=a(x,y)+b where a and b are
	[a,b] = lindecisbnd(Kest,meanAest,meanBest);
end

if modelnum == 1
	fishers_coeffs = [a, b];
elseif (modelnum == 3 | modelnum == 5)
   % use both dimensions
   fishers_coeffs = [a, b];
elseif modelnum == 4
	% the covariance matrix is the identity
	fishers_coeffs = [a, b];
end

