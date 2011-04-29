function [] = plot1dstim(data,xaxis,newfig)
%plot1dstim(data,xaxis,newfig)
%  plots a figure of one-dimensional stimulus data points.
%
%  Parameters:
%    data row format:  [cat x y resp]
%    xaxis format:  [minX maxX] (optional)
%    newfig indicates whether to plot the data onto a new
%       figure (newfig = 1) or not (newfig = 0). (optional)

% Created by Brian Spiering / 12/14/03
% from plot2dstim (Leola A. Alfonso-Reese / 7-August-95)
% Copyright (c) 2003
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
aind= find(data(:,1)==1);
bind= find(data(:,1)==2);
dataa=data(aind,:);
datab=data(bind,:);

lena = size(dataa,1);
lenb = size(datab,1);

% Equalize length of A and B 
if lena > lenb
	len = lena;
	xtra = len-lenb;
	datab = [datab; 2*ones(xtra,1) ones(xtra,1)*nan ones(xtra,1)*nan];	
elseif lena < lenb
	len = lenb;
	xtra = len-lena;
	dataa = [dataa; ones(xtra,1) ones(xtra,1)*nan ones(xtra,1)*nan]; 	
else len = lena;
end

% Plot A and B responses
hist([dataa(1:len,2),datab(1:len,2)]);
xlabel('X');
ylabel('Frequency');

curaxes = axis;

if ~isempty(xaxis)
	axis([xaxis curaxes(3:4)]);
%	axis(axis);
end

hold off;

