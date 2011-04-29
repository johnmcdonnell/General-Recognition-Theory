function fishers_coeffs = fisherdiscrim2d(data,modelnum)
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
%    row format of data:  [category x y response]
%    modelnum may be an integer from 1 to 5, 13 or 14 indicating...
%     01. estimate best 1-D linear classifier
%         based on sample means and sample covariance matrix
%         ignoring var and mean of dimension X
%     02. estimate best 1-D linear classifier
%         based on sample means and sample covariance matrix
%         ignoring var and mean of dimension Y
%     03. estimate optimal linear classifier
%         based on sample means and sample covariance matrix
%     04. estimate minimum distance classifier
%         based on sample means
%     05. estimate optimal linear classifier
%         based on response means and response covariance matrix
%
%     13. estimate optimal quadratic classifier
%         based on sample means and sample covariance matrix
%     15. estimate optimal quadratic classifier
%         based on response means and response covariance matrix
%
%  Note that models 3, 4 and 13 are based on the sample categories
%  while models 1, 2, 5 and 15 are based on the subject's responses.

% Created by Leola A. Alfonso-Reese / 18-January-95
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  9/4/95          Access correct column when model number is 1 or 2
%                  by Leola
%  1/18/96         Added model 5
%                  by Leola and Elliot
%  3/25/96         Added models 13 and 15
%                  by Leola
%  4/29/96         Checked for case where all stimuli are assigned one
%                  category label
%                  by Leola


if (modelnum == 3 | modelnum == 4 | modelnum == 13)
	% Get data using sample category labels
  Aindices = find(data(:,1) == 1);
  Alength = length(Aindices);
  Bindices = find(data(:,1) ~= 1);
  Blength = length(Bindices);
elseif (modelnum == 1 | modelnum == 2 | modelnum == 5 | modelnum == 15)
	% Get data using response category labels
  Aindices = find(data(:,4) == 1);
  Alength = length(Aindices);
  Bindices = find(data(:,4) ~= 1);
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

meanA = mean(data(Aindices,2:3))';
meanB = mean(data(Bindices,2:3))';
KA = cov(data(Aindices,2:3));
KB = cov(data(Bindices,2:3));
Kpooled = ((Alength-1)*KA+(Blength-1)*KB)/(Alength+Blength-2);

if modelnum == 1
	% ignore dimension X
	meanAest = meanA(2);
	meanBest = meanB(2);
	Kest = Kpooled(2,2);
elseif modelnum == 2
	% ignore dimension Y
	meanAest = meanA(1);
	meanBest = meanB(1);
	Kest = Kpooled(1,1);
elseif (modelnum == 3 | modelnum == 5)
	meanAest = meanA;
	meanBest = meanB;
	Kest = Kpooled;
elseif modelnum == 4
	meanAest = meanA;
	meanBest = meanB;
	% K is the identity matrix
	Kest = eye( 2 );
end

if modelnum < 10
	% Decision boundary for this model is h(x,y)=a(x,y)+b where a and b are
	[a,b]   = lindecisbnd(Kest,meanAest,meanBest);
else % modelnum is in the 20's
	% Decision boundary for the quadratic model is h(x,y) = x'Ax + b'x + c
	% where A is a 2x2 matrix, b is a 2x1 vector and c is a constant.
	[A,b,c] = quaddecisbnd(KA,KB,meanA,meanB);
end

if modelnum == 1
	% ignore dimension X
	fishers_coeffs = [0, a, b];
elseif modelnum == 2
	% ignore dimension Y
	fishers_coeffs = [a, 0, b];
elseif (modelnum == 3 | modelnum == 5)
   % use both dimensions
   fishers_coeffs = [a, b];
elseif modelnum == 4
	% the covariance matrix is the identity
	fishers_coeffs = [a, b];
else % modelnum == 13 or 15
	fishers_coeffs = [diag(A)' 2*A(1,2) b' c];
end

