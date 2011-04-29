function [out_params, neglikelihood] = fit_2dGQC(in_params, data, z_limit)
%[out_params, neglikelihood] = fit_2dGQC(in_params, data, z_limit)
%  sets up the parameter search range, sets up the search options, and
%  calls the search function, "constr", found in the optimization toolbox
%  for MATLAB.  The function being minimized is the negative log
%  likelihood of the 2 dimensional data, negloglike_2dGQC.
%
%  Parameters:
%    in_params format: [pnoise cnoise Amat_entries(1:3) b_entries(1:2) c_bias]
%    data row format:  [correct_response x y ... 1]
%    z_limit is the z-score value beyond which one should truncate

% Created by Leola A. Alfonso-Reese / 25-March-96
% Copyright (c) 1996
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  3/10/05         Modified options and search function to work 
%                  with Optimization Toolbox v 3.0
%                  by LAR

% Set contraints on parameters
vlb = [ .001   .001];
vub = [ 5000  5000];

% No conversion for search format.  
% Search format: [pnoise cnoise Amat_entries(1:3) b_entries(1:2) c_bias]
start_params = in_params;

% Set search options
options = optimset(...
    'Display', 'iter',...  % Spit out progress info
    'TolFun', .001');      % Termination tolerance on f    

% Search!
final_params = fmincon('negloglike_2dGQC',start_params,[],[],[],[],vlb,vub,[],options,data,z_limit);

% Report results
out_params = final_params;
neglikelihood = negloglike_2dGQC(final_params,data,z_limit);

