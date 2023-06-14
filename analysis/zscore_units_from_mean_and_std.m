function [z_units] = zscore_units_from_mean_and_std(units, u_means, u_stds)
%--------------------------------------
% z scores firing rates based on activity during a different period
% (whatever went into u_means and u_stds
%--------------------------------------
% INPUTS
% units - n_trialls x n_times x n_units array that contains raw firing rates
% u_means - an n_units x 1 array with the means of each unit from a some time point
% u_stds - an n_units x 1 array with the stds of each unit from some time point

% OUTPUTS
% z_units - n_trials x n_times x n_units array (same as the raw firing rates) but z-scored to the input means and stds
%--------------------------------------

% pre-allocate output
z_units = NaN(size(units));

[n_trials, n_times, n_units] = size(units);


for u = 1:n_units
    
    z_units(:,:,u) = (units(:,:,u) - u_means(u)) / u_stds(u);
        
end % of looping over units


end % of function