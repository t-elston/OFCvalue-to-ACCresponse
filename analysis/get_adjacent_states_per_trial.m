function [states2use] = get_adjacent_states_per_trial(ch_state_times, unch_state_times, trials2use, trial_ch_details, trial_unch_details)

state_times = [];
state_in_trial = [];
states2use=[];

for t = 1:numel(trials2use)
    
    % get the indices of the states associated with this trial
    ch_ix = ismember(trial_ch_details(:,1), trials2use(t));
    unch_ix = ismember(trial_unch_details(:,1), trials2use(t));
    
    % grab those trials 
    ch_times = ch_state_times(ch_ix);
    unch_times = unch_state_times(unch_ix);
    
    % accumulate the times of the chosen and unchosen states for this trial
    t_state_times = NaN(numel(unch_times),3);
    t_state_n =  NaN(numel(unch_times),4);
    
    for u = 1:numel(unch_times)
        
        this_unch_state_time = unch_times(u);        
        
        % get the times of the chosen states before/after the unchosen one
        ch_after_unch = ch_times(ch_times > this_unch_state_time);
        ch_before_unch = ch_times(ch_times < this_unch_state_time);
        
        % create vectors of state_in_trial
        
        ch_before_state_n = 1:numel(ch_before_unch);
        ch_after_state_n = numel(ch_before_unch)+1 : numel(ch_before_unch) + numel(ch_after_unch);
        
        
        % are there any chosen states before this unchosen one?
        if numel(ch_before_unch) > 0   
            
            % find the nearest one
            [~, nearest_before_state] = min(abs(this_unch_state_time - ch_before_unch)); 

            t_state_times(u,1) = ch_before_unch(nearest_before_state);   
            t_state_n(u,1) = ch_before_state_n(nearest_before_state);
            
        end
        
        % are there any chosen states after this unchosen one?
        if numel(ch_after_unch) > 0   
            
           [~, nearest_after_state] = min(abs(this_unch_state_time - ch_after_unch)); 
            t_state_times(u,2) = ch_after_unch(nearest_after_state);   
            t_state_n(u,2) = ch_after_state_n(nearest_after_state);
            
        end
        
        t_state_times(u,3) = this_unch_state_time;
        t_state_n(u,3) = u;
        t_state_n(u,4) = t;
        
    end % of looping over unchosen states
    
    state_times = [state_times ; t_state_times];
    state_in_trial = [state_in_trial; t_state_n];
    
    % now randomly select a chosen and unchosen state for use in this bootstrap
    
    has_adjacent_ch_state = false; 
    while ~has_adjacent_ch_state
        unch_state2_use = randi(numel(unch_times));
        has_adjacent_ch_state = nansum(t_state_n(unch_state2_use,1:2)) > 0;
    end
    
    
    ch_state2_use = NaN; 
    
    while isnan(ch_state2_use)
        ch_state2_use = t_state_n(unch_state2_use,randi(2));
    end
    
    states2use(t,:) = [unch_state2_use , ch_state2_use];
    
    used_state_times(t,:) = [unch_times(unch_state2_use) , ch_times(ch_state2_use)];
    
    
end % of looping over trials


mean_ch = nanmean(state_times(:,1:2),2);

% [h, p] = ttest2(used_state_times(:,1), used_state_times(:,2));

% figure; 
% hold on
% histogram(mean_ch);
% histogram(state_times(:,3));


end % of function