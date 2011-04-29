function [responsedata] = sim3dlin(stimuli,decision_params)
%[responsedata] = sim3dlin(stimuli,decision_params)
%  simulates an observer's response data according to the
%  noise and linear decision bound specified in
%  decision_params.
%
%  Parameters:
%    stimuli row format:  [cat x y]
%    decision_params format:  [noise a1 a2 a3 b]

% Created by Leola A. Alfonso-Reese / 29-October-96
% Copyright (c) 1996
% $Revisions$
%   Date           Modification and Name
%   ----           ---------------------


[m n] = size(stimuli);

noise = decision_params(1);

responsedata(:,1:4) = stimuli;

% Generate random noise
if noise == 0
  noise_sample = zeros(m,1);
else
	noise_sample = randn(m,1).*noise;
end

% Compute discriminant values and generate responses
decision_var = lindiscrim3dvals(stimuli(:,2:4),decision_params(2:5)) + noise_sample;
Aindices = find(decision_var <= 0);
Bindices = find(decision_var > 0);
if (~isempty(Aindices))
	responsedata(Aindices,5) = ones(length(Aindices),1);
end
if (~isempty(Bindices))
	responsedata(Bindices,5) = 2*ones(length(Bindices),1);
end

