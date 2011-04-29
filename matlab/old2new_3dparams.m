function new_params = old2new_3dparams(old_params,model)
%new_params = old2new_3dparams(old_params,model)
%  converts the format of the search parameters from
%  [noise a1 a2 a3 bias] to [noise angle1 angle2 bias].
%
%  Parameters:
%    old_params is a vector, [noise a1 a2 a3 b]
%    model is a string (optional)
%      'Bound' (default) =>     
%                'Trans' =>
%            'TransLong' =>

% Created by Leola A. Alfonso-Reese & David H. Brainard/ 14-March-94
% Copyright (c) 1994
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  3/12/97         by Leola A. Alfonso-Reese
%                  made minor changes so that it fits the mold
%                  of other routines in the GRT toolbox

% Default parameters
if nargin == 1
  model = 'Bound';
end

if strcmp(model,'Bound')
    norm_params = norm_old_3dparams(old_params,'Bound');
    noise = norm_params(1);
    xyz = norm_params(2:4);
    bias = norm_params(5);
    angles = xyzmat2angmat(xyz);
    new_params = [noise angles bias];
  elseif strcmp(model,'Trans')
    norm_params = norm_old_3dparams(old_params,'Trans');
    noise = norm_params(1);
    row1 = norm_params(2:4);
    row2 = norm_params(5:7);
    anglesRow1 = xyzmat2angmat(row1);
		anglesRow2 = xyzmat2angmat(row2);
    new_params = [noise anglesRow1 anglesRow2];
  elseif strcmp(model,'TransLong')
    norm_params = norm_old_3dparams(old_params,'TransLong');
    noise = norm_params(1);
    row1 = norm_params(2:4);
    row2 = norm_params(5:7);
		xy = norm_params(8:9);
		bias = norm_params(10);
    anglesRow1 = xyzmat2angmat(row1);
		anglesRow2 = xyzmat2angmat(row2);
		angle = xymat2angvec(xy);
    new_params = [noise anglesRow1 anglesRow2 angle bias];
  else  
    error('Invalid string');
end  

