function [rand_state_ix] = get_random_states_from_trial(n_states_per_trial)


n_trials = numel(n_states_per_trial);

rand_state_ix=NaN(size(n_states_per_trial));



for t = 1:n_trials
    
    rand_state_ix(t,1) = randi(n_states_per_trial(t));

end % of looping over trials


end % of function