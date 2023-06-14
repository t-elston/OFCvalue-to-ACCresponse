function [ch_states] = get_ch_val_states(mean_ch_pps, magnitude_thresh, time_thresh)

% initialize output
ch_states = zeros(size(mean_ch_pps));

[n_trials, n_times] = size(mean_ch_pps);

% loop over individual trials
for t = 1:n_trials
    
    t_data_over_mag_thresh = mean_ch_pps(t, :) >= magnitude_thresh;
    
    % now see how long this run was
    [run_starts, run_lens, ~] = ZeroOnesCount(t_data_over_mag_thresh);
    
    valid_runs = find(run_lens > time_thresh);
    
    
    for r = 1:numel(valid_runs)
        
        ch_states(t, run_starts(r) : run_starts(r) + run_lens(r)) = 1;
        
    end % of looping over valid runs
    
end % of looping over trials

end % of function