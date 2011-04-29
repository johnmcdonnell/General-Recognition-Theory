function pc = percorr(data,d)
%pc = percorr(data,d)
%  determines the percent of correct responses
%  in the categorization data.
%
%  Parameters:
%    data row format:  [correct_response x y ... subject_response]
%    d is number of dimensions
%
% Note: If data row format is [correct_resp subject_response]
%       then, enter zero for d

% Created by Leola A. Alfonso-Reese / 3-September-92
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%   5-February-95  Modified for any dimensional data
%                  by Leola A. Alfonso-Reese
%   27-May-96      Added parameter to size function call


if nargin > 1
	indices = find(data(:,1) == data(:,d+2));
	pc = length(indices)/size(data,1)*100;
else
  error('Required Parameter Missing:  d');
end


