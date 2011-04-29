function c = colorstr2colorval(s)
%colorval = colorstr2colorval(c)
%  returns a valid color number that resembles the color
%  specified in the string s
%
%  Parameters:
%    s is a 1, 2, or 3 character string made from the characters
%      listed in the plot command

% Created by Leola A. Alfonso-Reese / 24-November-04
% Copyright (c) 2004
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------

strarray = cellstr(s');
colorlist = 'ymcrgbwk';
keeplooking = 1;
colorindex = 1;
while keeplooking
	indices = strmatch(colorlist(colorindex),strarray,'exact');
	if isempty(indices)
		if colorindex == length(colorlist)
			keeplooking = 0;
		else			
			colorindex = colorindex+1;
		end		
	else
		keeplooking = 0;
	end
end

if isempty(indices)
	error('Error in color/linetype argument');
else
	switch colorindex
		case 1 % yellow
			c = 9;
		case 2 % magenta
			c = 40;
		case 3 % cyan
			c = 25;
		case 4 % red
			c = 1;
		case 5 % green
			c = 20;
		case 6 % blue
			c = 35;
		case 7 % white => yellow;
			c =  9
		otherwise % black => blue
			c = 35;
	end

end

