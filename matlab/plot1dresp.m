function [freqs,xbinctrs] = plot1dresp(data,xaxis,newfig)
%plot1dstim(data,xaxis,newfig)
%  plots a figure of one-dimensional stimulus data points.
%
%  Parameters:
%    data row format:  [cat x resp]
%    xaxis format:  [minX maxX] (optional)
%    newfig indicates whether to plot the data onto a new
%       figure (newfig = 1) or not (newfig = 0). (optional)
%
%[freqs,xbinctrs] = plot1dresp(...) returns the frequencies for
%  each bin and the position of the bin centers in X.

% Created by Leola Alfonso-Reese / 27-August-04
% Copyright (c) 2004
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------

if nargin < 2
  xaxis = [];
end

if nargin < 3
	newfig = 1;
end

if newfig == 1
	figure;
end

% Find A and B response subsets
aind= find(data(:,3)==1);
bind= find(data(:,3)==2);
dataa=data(aind,:);
datab=data(bind,:);

lena = size(dataa,1);
lenb = size(datab,1);

% Equalize length of A and B 
if lena > lenb
	len = lena;
	xtra = len-lenb;
	datab = [datab; 2*ones(xtra,1) ones(xtra,1)*nan nan*ones(xtra,1)];	
elseif lena < lenb
	len = lenb;
	xtra = len-lena;
	dataa = [dataa;  ones(xtra,1) ones(xtra,1)*nan nan*ones(xtra,1)]; 	
else len = lena;
end

% Plot A and B responses
[freqs,xbinctrs]=hist([dataa(1:len,2),datab(1:len,2)]);
hist([dataa(1:len,2),datab(1:len,2)]);
xlabel('X');
ylabel('Frequency');

curaxes = axis;

if ~isempty(xaxis)
	axis([xaxis curaxes(3:4)]);
%	axis(axis);
end

hold off;

