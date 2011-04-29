function s = colorval2colorstr(c)
%color = colorval2colorstr(c)
%  returns a valid string from the plot function that is closest to
%  the specified color c

% Created by Leola A. Alfonso-Reese / 7-October-94
% Copyright (c) 1994
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------

if c <= 5
  color = 'r';
elseif c <= 10
  color = 'y';
elseif c <= 21
  color = 'g';
elseif c <= 28
  color = 'c';
elseif c <= 37
  color = 'b';
else
  color = 'm';
end

