function [newdata,nclip] = clip_vals(data,range,msg)
%[newdata,nclip] = clip_vals(data,rangemsg,msg)
%  clips values outside of the range for each dimension
%  that has a range specified.
%
%  Parameters:
%    data row format: [var1 var2 ... varn]
%    range = [minvar1 maxvar1
%             minvar2 maxvar2
%             ...
%             minvarn maxvarn]
%    To avoid clipping along any of the dimensions, specify
%    NaN entries for the minimum or maximum range values.
%    msg:  flag indicating whether to display clipping message
%          0 => do not display message
%          1 => display message

% Created by Leola A. Alfonso-Reese / 7-August-95
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------
%  20-October-97   Added check for msg argument
%                  by Leola A. Alfonso-Reese
%  6-April-01      Added nclip output parameter
%  5-Dec-02        Fix bug where only large values were being reported

if nargin < 3
	msg = 1;
elseif isempty(msg)
	msg = 1;
end

[dims,m] = size(range);
newdata = data;

nclip = 0;

for i = 1:dims
	minval = range(i,1);
	maxval = range(i,2);
	if ~isnan(minval)
		lessthanmin = find(data(:,i) < minval);
		if ~isempty(lessthanmin)
			newdata(lessthanmin,i) = minval*ones(length(lessthanmin),1);
			printmsg = 'fprintf(''Clipped %i values less than %6.4f.\n'',length(lessthanmin),minval);';
			nclip = nclip + length(lessthanmin);
		else
			printmsg = 'fprintf(''No values found less than %6.4f.\n'',minval);';  
		end
	end
	if msg == 1
		eval(printmsg);
	end
	if ~isnan(maxval)
		exceedmax = find(data(:,i) > maxval);
		if ~isempty(exceedmax)
			newdata(exceedmax,i) = maxval*ones(length(exceedmax),1);
			printmsg = 'fprintf(''Clipped %i values greater than %6.4f.\n'',length(exceedmax),maxval);';
			nclip = nclip + length(exceedmax);
		else
			printmsg = 'fprintf(''No values found greater than %6.4f.\n'',maxval);';	  
		end
	end
	if msg == 1
		eval(printmsg);
	end
end

