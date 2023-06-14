function [z_out] = zscore_tensor_v01(indata)


[n_trials, n_times, n_units] = size(indata);

z_out = NaN(size(indata));

for u = 1:n_units
    
    u_mean = nanmean(indata(:,:,u),'all');
    u_std = nanstd(indata(:,:,u),[],'all');
    
    z_out(:,:,u) = (indata(:,:,u) - u_mean) / u_std;
    
    
    
end % of looping over units



end % of function 