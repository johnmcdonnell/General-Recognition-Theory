function [out_params, neglikelihood] = fit_3dGLC(in_params, data, z_limit)
%[out_params, neglikelihood] = fit_3dGLC(in_params, data, z_limit)
%  sets up the parameter search range, sets up the search options, and
%  calls the search function, "constr", found in the optimization toolbox
%  for MATLAB.  The function being minimized is the negative log
%  likelihood of the 3-dimensional data, negloglike_3dGLC.
%
%  Parameters:
%    in_params format:  [noise a1 a2 a3 bias]
%    data row format:  [response x y z 1]
%    z_limit is the z-score value beyond which one should truncate

% Created by Leola A. Alfonso-Reese / 10-March-94
% Copyright (c) 1994
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  3/12/97         by Leola A. Alfonso-Reese
%                  made minor changes so that it fits the mold
%                  of other routines in the GRT toolbox
%  3/10/05         Modified options and search function to work 
%                  with Optimization Toolbox v 3.0
%                  by LAR

% Set contraints on parameters
%vlb = [ 0    -pi/2    0   ];
%vub = [ 500   pi/2    2*pi];

vlb = [ .001 ];
vub = [ 500 ];

% Convert params to search format.  
% Search format: [noise angles(1:2) bias]
start_params = old2new_3dparams(in_params);

% Set search options
options = optimset(...
    'Display', 'iter',...  % Spit out progress info
    'TolFun', .001',...    % Termination tolerance on f
    'MaxIter', 2000);      % Max iterations

% Search!
final_params = fmincon('negloglike_3dGLC',start_params,[],[],[],[],vlb,vub,[],options,data,z_limit);

% Report results
out_params = new2old_3dparams(final_params);
neglikelihood = negloglike_3dGLC(final_params,data,z_limit);

