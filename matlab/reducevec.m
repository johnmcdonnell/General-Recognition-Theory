function outvec = reducevec(invec,dim1,dim2)
% outvec reducevec(invec,dim1,dim2)
%   reduces an n-dimensional vector 'invec' to a 2-dimensional
%   vector.  The parameters 'dim1' and 'dim2' are the indices
%   of elements that are to be extracted from the original vector.
%
%  Parameters:
%    invec is a vector of length n
%		 dimA and dimB are indices of invec

% Created by Leola A. Alfonso-Reese / 10-March-94
% Copyright (c) 1994
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  3/12/97         by Leola A. Alfonso-Reese
%                  made minor changes so that it fits the mold
%                  of other routines in the GRT toolbox

temp = invec;
temp(1) = invec(dim1);
temp(2) = invec(dim2);
outvec = temp(1:2);

