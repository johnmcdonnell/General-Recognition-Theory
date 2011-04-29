function out_params = norm_old_1dparams(old_params)
%out_params = norm_old_1dparams(old_params)
%  normalizes the old parameters [noise a b],
%  where a*x + b = 0.

% Created by Leola A. Alfonso-Reese/ 23-August-04
% Copyright (c) 2004
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------

out_params = old_params / abs(old_params(2));


