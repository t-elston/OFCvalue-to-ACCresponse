function [PP_baselines] = get_value_PPs_when_not_available(PPs, ch_val, unch_val)

values = unique(ch_val);

for v = 1:numel(values)
    
    % find trials where this value was not present
    this_val_na_trials = ch_val~=v & unch_val~=v;
    
    % now get the mean of these trials for this state/value
    PP_baselines(v,:) = nanmean(PPs(this_val_na_trials,:,v));
    
    
end % of looping over values




end % of function