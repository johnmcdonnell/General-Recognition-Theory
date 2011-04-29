function [out_params, neglikelihood] = fit_3dGQC(in_params,data,z_limit)
%[out_params, neglikelihood] = fit_3dGQC(in_params,data,z_limit)
%  sets up the parameter search range, sets up the search options, and
%  calls the search function, "constr", found in the optimization toolbox
%  for MATLAB.  The function being minimized is the negative log
%  likelihood of the 3-dimensional data, negloglike_3dGQC.
%
%  Parameters:
%    in_params format: [pnoise cnoise Amat_entries(1:6) b_entries(1:3) c_bias]
%    data row format:  [correct_response x y z 1]
%    z_limit is the z-score value beyond which one should truncate

% Created by Leola A. Alfonso-Reese / 10-July-94
% Copyright (c) 1994
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  3/10/05         Modified to work with Optimization Toolbox v3.0
%                  by LAR

% Set contraints on parameters
%vlb = [ .01   .01];
vlb = [ .1   .1];
vub = [ 5000 5000];

% Convert params to search format.  
% Search format: [pnoise cnoise Amat_entries(1:6) b_entries(1:3) c_bias]
start_params = in_params;

% Set search options
options = optimset(...
    'Display', 'iter',...  % Spit out progress info
    'TolFun', .001',...    % Termination tolerance on f
    'MaxIter', 4000);      % Max iterations

% Search!
final_params = fmincon('negloglike_3dGQC',start_params,[],[],[],[],vlb,vub,[],options,data,z_limit);


% Report results
out_params = final_params;
neglikelihood = negloglike_3dGQC(final_params,data,z_limit);

