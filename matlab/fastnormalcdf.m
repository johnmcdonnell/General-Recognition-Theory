function Fx = fastnormalcdf(xvec)
%Fx = fastnormalcdf(xvec)
%  returns the normal cumulative distribution function with
%  mean 0 and standard deviation 1 for each number in xvec.
%
%  This function is optimized for speed when processing vectors.
%
%  Parameters:
%    xvec is a vector of real numbers
%
%  To speed things up, this routine uses a rational approximation
%  for the normal cdf taken from Milton Abramowitz and Irene A.
%  Stegun's "Handbook of mathematical functions"

% Created by Leola A. Alfonso-Reese / 10-March-94
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------


% Allocate space
Fx = ones(length(xvec),1);

% Handle positive and negative cases separately
neg_indices = find(xvec < 0);
tval = 1 ./ (1+ 0.33267.*(-xvec(neg_indices)) );
tval2 = tval.*tval;
tval3 = tval2.*tval;
Fx(neg_indices) = (exp(-(xvec(neg_indices).^2)/2)/2.50662827463100).*(0.4361836.*tval-0.1201676.*tval2+0.9372980.*tval3);

pos_indices = find(xvec>= 0);
tval = 1 ./ (1 + 0.33267.*xvec(pos_indices));
tval2 = tval.*tval;
tval3 = tval2.*tval;
Fx(pos_indices) = 1-(exp(-(xvec(pos_indices).^2)/2)/2.50662827463100).*(0.4361836.*tval-0.1201676.*tval2+0.9372980.*tval3);

