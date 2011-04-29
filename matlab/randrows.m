function Y = randrows(X)
%Y = randrows(X)
%  randomizes the order of the rows in the matrix X.
%
%  Warning:  It is recommended that you set the 'seed'
%  of the uniform generator before calling this function
%  within a MATLAB session.  To do this, type
%        rand('seed',sum(100*clock))

% Created by Leola A. Alfonso-Reese / 17-April-94
% Copyright (c) 1995
% $Revisions:$
%   Date           Modification and Name
%   ----           ---------------------

num_rows = size(X,1);
upd_num_rows = num_rows;
Y = [];
for i = 1:num_rows
  % Pick a row from available rows in X and append it to the new matrix Y.
	rand_idx = ceil(rand*upd_num_rows);

  Y = [Y; X(rand_idx,:)];
   
  % If the row that was copied from X was not the last row in X, then
  % move the last row in X into the place of the one that was copied.
  if rand_idx ~= upd_num_rows
    X(rand_idx,:) = X(upd_num_rows,:);
  end
  upd_num_rows = upd_num_rows-1;
end

