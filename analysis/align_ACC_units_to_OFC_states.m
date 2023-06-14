function [win_ch, win_unch, ch_FRs, unch_FRs, xt] = align_ACC_units_to_OFC_states(units, unit_ts, ch_details, unch_details)

% format of ch_details and unch_details is:
% column 1 = trial state occured in
% column 2 = the state within the trial (first, second, third, etc)
% column 3 = the value of the state
% column 4 = the screen side of the state (-1 = left, 1 = right)
% column 5 = the time of state onset                                                                                             


win_ch=[];
win_unch=[];
ch_FRs=[];
unch_FRs=[];
xt=[];


%-------------
% cycle through the chosen states
[n_ch_states,~] = size(ch_details);

for i = 1:n_ch_states
    
    trial_num = ch_details(i,1); 
    onset_time = ch_details(i,5);
    
    % now find the nearest timestep to this state's onset and get the indices for the windows
    [~,onset_ix] = min(abs(unit_ts - onset_time));
    [~,small_after_end] = min(abs(unit_ts - (onset_time+100)));
    [~,small_before_start] = min(abs(unit_ts - (onset_time -200)));
    [~,small_before_end] = min(abs(unit_ts - (onset_time -50)));
    [~,big_window_start] = min(abs(unit_ts - (onset_time-700)));
    [~,big_window_end] = min(abs(unit_ts - (onset_time +700)));
    
    % get the mean firing rates in the windows just before and after this state onset for all units
    win_ch(i,1,:) = squeeze(nanmean(units(trial_num,small_before_start:small_before_end,:),2));
    win_ch(i,2,:) = squeeze(nanmean(units(trial_num,onset_ix:small_after_end,:),2));
    
    ch_FRs(i,:,:) = units(trial_num,big_window_start:big_window_end,:);
    


end % of cycling over chosen states

%-------------
% cycle through the unchosen states
[n_unch_states,~] = size(unch_details);

for j = 1:n_unch_states
    
    trial_num = unch_details(j,1); 
    onset_time = unch_details(j,5);
    
    % now find the nearest timestep to this state's onset and get the indices for the windows
    [~,onset_ix] = min(abs(unit_ts - onset_time));
    [~,small_after_end] = min(abs(unit_ts - (onset_time+100)));
    [~,small_before_start] = min(abs(unit_ts - (onset_time -200)));
    [~,small_before_end] = min(abs(unit_ts - (onset_time -50)));
    [~,big_window_start] = min(abs(unit_ts - (onset_time-700)));
    [~,big_window_end] = min(abs(unit_ts - (onset_time +700)));
    
    % get the mean firing rates in the windows just before and after this state onset for all units
    win_unch(j,1,:) = squeeze(nanmean(units(trial_num,small_before_start:small_before_end,:),2));
    win_unch(j,2,:) = squeeze(nanmean(units(trial_num,onset_ix:small_after_end,:),2));

    unch_FRs(j,:,:) = units(trial_num,big_window_start:big_window_end,:);


end % of cycling over chosen states

[~, n_times, n_units] = size(ch_FRs);


xt = [1:n_times]*mean(diff(unit_ts)) - 700;




end % of function