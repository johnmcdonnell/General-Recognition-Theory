function [ neglikelihood , bic, params ] = fitsubj( fittype, this ) 
% fittype is as follows:
%  0 : 0d with fixed p=0.5
%  1 : 0d with floating p
%  2 : 1d: x axis
%  3 : 1d: y axis
%  4 : 2d
% datafile should be of form data/41.dat' (unix readable)

sinit = 10;

if fittype == 2
    data = [this(:,1:2), this(:,4)];
elseif fittype == 3
    data = [this(:,1), this(:,3:4)];
else
    data = this
end

n = size( data, 1 )

switch fittype
	case 0
		neglikelihood = - log( .5 ) * n
		params = [.5]
	case 1
		numa = sum( data(:, 1) )
		numb = n - numa
		proba = mean( data(:,1) )
		neglikelihood = - ( log( proba ) * numa + log( 1-proba ) * numb )
		params = [proba]
	case 4
		% Use Fisher's linear discriminant for Initial search values.
		% Parameter format: [noise unit_normal(1:2) bias]
		fishers_coeffs = fisherdiscrim2d(data, 3);
		raw_params = [sinit, fishers_coeffs];
		start_params = norm_old_2dparams(raw_params);
		fprintf('...Searching for best fit\n'); 
		[final_params neglikelihood] = fit_2dGLC(start_params,data, 7);
		params = final_params;
	otherwise % 1d case.
		% Use Fisher's linear discriminant for Initial search values.
		% Parameter format: [noise unit_normal(1:2) bias]
		fishers_coeffs = fisherdiscrim1d(data, 3);
		raw_params = [sinit, fishers_coeffs];
		start_params = norm_old_1dparams(raw_params);
		fprintf('...Searching for best fit\n'); 
		[final_params neglikelihood] = fit_1dGLC(start_params,data, 7);
		params = final_params
end


% FInd the AIC score so that this GLC model may be compared
% with others.
% AIC = 2(-logL + r)
% where r = 3 (2 of the XY coordinates + 1 bias
%              + 1 noise - 1 since the XY & bias params
%              are normalized)
% Note:  neglikelihood from the Fit routine is -logL
if fittype == 0
	r = 0;
elseif fittype == 1
	r = 1;
elseif fittype == 4
	[slint(:,1) slint(:,2)] = bnd2slint(params);
	r = 3;
else  % case of 1d fits.
	r = 2;
end

aic = 2 * neglikelihood + r * 2;
bic = 2 * neglikelihood + r * log(n);

% Find the slopes and y-intercepts for this bound

% plot_fit( data, slint )

end
