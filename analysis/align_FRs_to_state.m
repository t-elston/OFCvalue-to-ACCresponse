function [state_aligned_FRs, state_sides] = align_FRs_to_state(firing_rates, state_times, ts, offsets, trial_side)


state_aligned_FRs=[];
state_sides=[];
[n_trials, n_times, n_units] = size(firing_rates);


for t = 1:n_trials
    
    t_state_FRs=[];
    t_states = state_times(t,:);
    t_states(isnan(t_states))=[];

    % now loop over each state and find window aligned on state onset
    for ch_s = 1:numel(t_states)
        [~,state_start] = min(abs(ts - (t_states(ch_s)-offsets(1))));
        [~,state_end] = min(abs(ts - (t_states(ch_s)+offsets(2))));
        
        t_state_FRs(ch_s,:,:) = firing_rates(t,state_start:state_end,:);
    end
    
    state_aligned_FRs = [state_aligned_FRs; t_state_FRs];
    state_sides = [state_sides; ones(ch_s,1)*trial_side(t)];

end % of looping over trials



end % of function 