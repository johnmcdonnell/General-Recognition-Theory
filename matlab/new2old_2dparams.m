function old_params = new2old_2dparams(new_params)
%old_params = new2old_2dparams(new_params)
%  converts the format of the search parameters from
%  [noise angle bias] to [noise a1 a2 bias].

% Created by Leola A. Alfonso-Reese & David H. Brainard/ 14-March-94
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------


noise = new_params(1);
xy = angvec2xymat( new_params(2) );
bias = new_params(3);
old_params = [noise xy bias];


