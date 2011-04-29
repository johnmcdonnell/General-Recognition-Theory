function [out_params, neglikelihood] = fit_1dGLC(in_params, data, z_limit)
%[out_params, neglikelihood] = fit_1dGLC(in_params, data, z_limit)
%  sets up the parameter search range, sets up the search options, and
%  calls the search function, "constr", found in the optimization toolbox
%  for MATLAB.  The function being minimized is the negative log
%  likelihood of the 1 dimensional data, negloglike_1dGLC.
%
%  Parameters:
%    in_params format:  [noise a bias]
%    Assuming these parameters were normalized before calling this subroutine,
%    then the 'a' parameter should simply be one.
%    data row format:  [response x 1]
%    z_limit is the z-score value beyond which one should truncate

% Created by Leola A. Alfonso-Reese / 23-August-04
% Copyright (c) 2004
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  3/10/05         Modified options and search function to work 
%                  with Optimization Toolbox v 3.0
%                  by LAR


% Set contraints on parameters
vlb = [ .001 in_params(2)];
vub = [ 500 in_params(2)];

% Convert params to search format.  
% Search format: [noise bias]
start_params = in_params;

% Set search options
options = optimset(...
    'Display', 'iter',...  % Spit out progress info
    'TolFun', .001');      % Termination tolerance on f    

% Search!
[out_params, neglikelihood] = fmincon('negloglike_1dGLC',start_params,[],[],[],[],vlb,vub,[],options,data,z_limit);


