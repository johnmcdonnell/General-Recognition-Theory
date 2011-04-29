function [out_params, neglikelihood] = fit_2dGLC(in_params, data, z_limit)
%[out_params, neglikelihood] = fit_2dGLC(in_params, data, z_limit)
%  sets up the parameter search range, sets up the search options, and
%  calls the search function, "constr", found in the optimization toolbox
%  for MATLAB.  The function being minimized is the negative log
%  likelihood of the 2 dimensional data, negloglike_2dGLC.
%
%  Parameters:
%    in_params format:  [noise a1 a2 bias]
%    data row format:  [response x y ... 1]
%    z_limit is the z-score value beyond which one should truncate

% Created by Leola A. Alfonso-Reese / 18-April-95
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  5/4/96          Modified help text);
%                  by LAR
%  3/10/05         Modified options and search function to work 
%                  with Optimization Toolbox v 3.0
%                  by LAR

% Set contraints on parameters
vlb = [ .001 ];
vub = [ 500 ];

% Convert params to search format.  
% Search format: [noise angle bias]
start_params = old2new_2dparams(in_params);

% Set search options
options = optimset(...
    'Display', 'iter',...  % Spit out progress info
    'TolFun', .001');      % Termination tolerance on f    

% Search!
final_params = fmincon('negloglike_2dGLC',start_params,[],[],[],[],vlb,vub,[],options,data,z_limit);

% Report results
out_params = new2old_2dparams(final_params);
neglikelihood = negloglike_2dGLC(final_params,data,z_limit);

