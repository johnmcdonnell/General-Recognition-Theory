function outmat = reducemat(inmat,dim1,dim2)
% outmat reducemat(inmat,dim1,dim2)
%   reduces an n-dimensional matrix 'inmat' to a 2-dimensional 
%   matrix.  The parameters 'dim1' and 'dim2' are the indices
%   of elements that are to be extracted from the original matrix.
%
% The parameters:
% 	inmat is an nxn matrix
%		dimA and dimB are indices of 'inmat'

% Created by Leola A. Alfonso-Reese / 10-March-94
% Copyright (c) 1994
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  3/12/97         by Leola A. Alfonso-Reese
%                  made minor changes so that it fits the mold
%                  of other routines in the GRT toolbox

outmat(1,1) = inmat(dim1,dim1);
outmat(1,2) = inmat(dim1,dim2);
outmat(2,1) = inmat(dim2,dim1);
outmat(2,2) = inmat(dim2,dim2);

