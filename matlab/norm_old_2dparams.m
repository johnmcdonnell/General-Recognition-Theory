function out_params = norm_old_2dparams(old_params)
%out_params = norm_old_2dparams(old_params)
%  normalizes the old parameters [noise a1 a2 b],
%  where a1*x + a2*y + b = 0.

% Created by Leola A. Alfonso-Reese & David H. Brainard/ 16-May-94
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------


out_params(1) = old_params(1) / norm(old_params(2:3));
out_params(2:3) = old_params(2:3) ./ norm(old_params(2:3));
out_params(4) = old_params(4) ./ norm(old_params(2:3));


