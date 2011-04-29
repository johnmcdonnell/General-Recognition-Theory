function [responsedata] = sim1dlin(stimuli,decision_params)
%[responsedata] = sim1dlin(stimuli,decision_params)
%  simulates an observer's response data according to the
%  noise and linear decision bound specified in
%  decision_params.
%
%  Parameters:
%    stimuli row format:  [cat x]
%    decision_params format:  [noise a b]

% Created by Leola A. Alfonso-Reese & Brian Spiering/ 1-Nov-02
% Copyright (c) 2002
% $Revisions$
%   Date           Modification and Name
%   ----           ---------------------


[m n] = size(stimuli);

noise = decision_params(1);

responsedata(:,1:2) = stimuli;

% Generate random noise
if noise == 0
  noise_sample = zeros(m,1);
else
	noise_sample = randn(m,1).*noise;
end

% Compute discriminant values and generate responses
decision_var = lindiscrim1dvals(stimuli(:,2),decision_params(2:3)) + noise_sample;
Aindices = find(decision_var <= 0);
Bindices = find(decision_var > 0);
if (~isempty(Aindices))
	responsedata(Aindices,3) = ones(length(Aindices),1);
end
if (~isempty(Bindices))
	responsedata(Bindices,3) = 2*ones(length(Bindices),1);
end

