function out_params = norm_old_3dparams(old_params,model)
%out_params = norm_old_3dparams(old_params,model)
%  normalizes the old parameters [noise a1 a2 a3 b],
%  where a1*x + a2*y + a3*z + b = 0.
%
%  Parameters:
%    old_params is a vector, [noise a1 a2 a3 b]
%    model is a string (optional)
%      'Bound' (default) =>     
%                'Trans' =>
%            'TransLong' =>

% Created by Leola A. Alfonso-Reese & David H. Brainard/ 16-May-94
% Copyright (c) 1994
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------

% Default parameters
if nargin == 1
  model = 'Bound';
end

if strcmp(model,'Bound')
  out_params(1) = old_params(1) / norm(old_params(2:4));
  out_params(2:4) = old_params(2:4) ./ norm(old_params(2:4));
  out_params(5) = old_params(5) ./ norm(old_params(2:4));
elseif strcmp(model,'Trans')
  out_params(1) = old_params(1) / norm(old_params(2:4));
  out_params(2:4) = old_params(2:4) ./ norm(old_params(2:4));
  out_params(5:7) = old_params(5:7) ./ norm(old_params(5:7));
elseif strcmp(model,'TransLong')
  out_params(1) = old_params(1) / norm(old_params(2:4));
	% two rows in the transformation matrix A
  out_params(2:4) = old_params(2:4) ./ norm(old_params(2:4));
  out_params(5:7) = old_params(5:7) ./ norm(old_params(5:7));
	% the bound coefficients
	out_params(8:9) = old_params(8:9) ./ norm(old_params(8:9));
	out_params(10) = old_params(10) ./ norm(old_params(8:9));
else  
  error('Invalid string');
end  

