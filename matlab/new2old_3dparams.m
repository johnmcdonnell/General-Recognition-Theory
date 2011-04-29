function [old_params] = new2old_3dparams(new_params,model)
%old_params = new2old_3dparams(new_params,model)
%  converts the format of the search parameters from
%  [noise angle1 angle2 bias] to [noise a1 a2 a3 bias].
%
%  Parameters:
%    new_params is a vector, [noise angle1 angle2 b]
%    model is a string (optional)
%      'Bound' (default) =>     
%                'Trans' =>
%            'TransLong' =>

% Created by Leola A. Alfonso-Reese & David H. Brainard/ 16-May-94
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
    noise = new_params(1);
    xyz = angsmat2xyzmat( new_params(2:3) );
    bias = new_params(4);
    old_params = [noise xyz bias];
  elseif strcmp(model,'Trans')
    noise = new_params(1);
    row1 = angsmat2xyzmat( new_params(2:3) );
    row2 = angsmat2xyzmat( new_params(4:5) );
    old_params = [noise row1 row2];
  elseif strcmp(model,'TransLong')
    noise = new_params(1);
    row1 = angsmat2xyzmat( new_params(2:3) );
    row2 = angsmat2xyzmat( new_params(4:5) );
	xy = angvec2xymat( new_params(6) );
	bias = new_params(7);
    old_params = [noise row1 row2 xy bias];	
  else  
    error('Invalid string');
end  

