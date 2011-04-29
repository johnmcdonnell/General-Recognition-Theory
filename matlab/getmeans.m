function [meanAest,meanBest] = getmeans(modelnum,meanA,meanB)
%[meanAest,meanBest] = getmeans(modelnum,meanA,meanB)
%  returns the mean vectors that correspond with the
%  model indicated by 'modelnum'.
%
%  Parameters:
%    modelnum may be an integer from 1 to 10 (See fisherdiscrim3d)
%    meanA and meanB are 1x3 vectors, i.e., the presumed structure is
%       mean = [meanX
%               meanY
%               meanZ]

% Created by Leola A. Alfonso-Reese / 10-March-94
% Copyright (c) 1994
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  3/12/97         by Leola A. Alfonso-Reese
%                  made minor changes so that it fits the mold
%                  of other routines in the GRT toolbox

if modelnum == 1
   % ignore dimension X
   meanAest = reducevec(meanA,2,3);
   meanBest = reducevec(meanB,2,3);
elseif modelnum == 2
   % ignore dimension Y
   meanAest = reducevec(meanA,1,3);
   meanBest = reducevec(meanB,1,3);
elseif modelnum == 3
   % ignore dimension Z
   meanAest = reducevec(meanA,1,2);
   meanBest = reducevec(meanB,1,2);
else
   % use original three-dimensional vectors
   meanAest = meanA;
   meanBest = meanB;
end

