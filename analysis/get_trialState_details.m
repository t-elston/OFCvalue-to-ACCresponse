function [ch_val, unch_val, ch_nStates, unch_nStates, ch_times, unch_times] = get_trialState_details(t_values,choice_dirs, animal_id, ch_states, unch_states, t_mids)

[n_trials, n_times] = size(ch_states);

chosen_side = zeros(size(choice_dirs));
unchosen_side = zeros(size(choice_dirs));
chosen_side(choice_dirs== 1) = 2; 
chosen_side(choice_dirs== -1)= 1; 
chosen_side(chosen_side==0) = 1;

unchosen_side(choice_dirs== -1) = 2; 
unchosen_side(choice_dirs== 1)= 1; 
unchosen_side(unchosen_side==0) = 1;



if animal_id == 0
    [~,state_win_start_time_ix] = min(abs(t_mids - -440));
    [~,max_time_ix] = min(abs(t_mids - -50));
else
    [~,state_win_start_time_ix] = min(abs(t_mids - -420));
    [~,max_time_ix] = min(abs(t_mids - -110));
end



%     [~,state_win_start_time_ix] = min(abs(t_mids - -500));
%     [~,max_time_ix] = min(abs(t_mids - -50));



for t = 1:n_trials
        
    % what were the values for this trial?
    ch_val(t,1) = t_values(t,chosen_side(t));
    unch_val(t,1) = t_values(t,unchosen_side(t));
    
    t_ch_state = ch_states(t,:);
    
    % get this trial's unchosen state
    t_unch_state = unch_states(t,:);

    % find when the chosen state was 'up'
    [~, ch_state_starts, ch_state_widths] = findpeaks(t_ch_state);

    
    % find when the unchosen state was 'up'
    [~, unch_state_starts, unch_state_widths] = findpeaks(t_unch_state);
    
    ch_durations = mean(diff(t_mids))*ch_state_widths;
    unch_durations = mean(diff(t_mids))*unch_state_widths;

    % get the indices of the state starts that happened after the pics turned on
    ch_state_ix = ch_state_starts((ch_state_starts>state_win_start_time_ix) & (ch_state_starts<max_time_ix));
    unch_state_ix = unch_state_starts((unch_state_starts>state_win_start_time_ix) & (unch_state_starts<max_time_ix));
       
    ch_durations = ch_durations((ch_state_starts>state_win_start_time_ix) & (ch_state_starts<max_time_ix));
    unch_durations = unch_durations((unch_state_starts>state_win_start_time_ix) & (unch_state_starts<max_time_ix));
    
    % now get the times in ms of when the states started (will use these to align the units later)
    ch_state_start_times = t_mids(ch_state_ix);
    unch_state_start_times = t_mids(unch_state_ix);
    
    % now save the chosen and unchosen state times
    ch_times(t,1:numel(ch_state_start_times)) = ch_state_start_times; 
    unch_times(t,1:numel(unch_state_start_times)) = unch_state_start_times; 
    
    ch_nStates(t,1) = numel(ch_state_ix);
    unch_nStates(t,1) = numel(unch_state_ix);

end % of cycling over trials


ch_times(ch_times==0) = NaN;
unch_times(unch_times==0) = NaN;







end % of function