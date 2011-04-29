function new_params = old2new_2dparams(old_params)
%new_params = old2new_2dparams(old_params)
%  converts the format of the search parameters from
%  [noise a1 a2 bias] to [noise angle bias].

% Created by Leola A. Alfonso-Reese & David H. Brainard/ 14-March-94
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------


norm_params = norm_old_2dparams(old_params);
noise = norm_params(1);
xy = norm_params(2:3);
bias = norm_params(4);
angle = xymat2angvec(xy);
new_params = [noise angle bias];

