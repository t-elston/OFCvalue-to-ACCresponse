function [ch_FRs, unch_FRs, xt] = align_unit_to_states(u_data, unit_ts, ch_details, unch_details)

% format of ch_details and unch_details is:
% column 1 = trial state occured in
% column 2 = the state within the trial (first, second, third, etc)
% column 3 = the value of the state
% column 4 = the screen side of the state (-1 = left, 1 = right)
% column 5 = the time of state onset                                                                                             

ch_FRs=[];
unch_FRs=[];
xt=[];

offset = 500;

%-------------
% cycle through the chosen states
[n_ch_states,~] = size(ch_details);

for i = 1:n_ch_states
    
    trial_num = ch_details(i,1); 
    onset_time = ch_details(i,5);
    
    % now find the nearest timestep to this state's onset and get the indices for the windows
    [~,onset_ix] = min(abs(unit_ts - onset_time));
    [~,big_window_start] = min(abs(unit_ts - (onset_time-offset)));
    [~,big_window_end] = min(abs(unit_ts - (onset_time +offset)));
      
    ch_FRs(i,:) = u_data(trial_num,big_window_start:big_window_end);
   
end % of cycling over chosen states

%-------------
% cycle through the unchosen states
[n_unch_states,~] = size(unch_details);

for j = 1:n_unch_states
    
    trial_num = unch_details(j,1); 
    onset_time = unch_details(j,5);
    
    % now find the nearest timestep to this state's onset and get the indices for the windows
    [~,onset_ix] = min(abs(unit_ts - onset_time));
    [~,big_window_start] = min(abs(unit_ts - (onset_time-offset)));
    [~,big_window_end] = min(abs(unit_ts - (onset_time +offset)));

    unch_FRs(j,:) = u_data(trial_num,big_window_start:big_window_end);


end % of cycling over chosen states

[~, n_times, n_units] = size(ch_FRs);


xt = [1:n_times]*mean(diff(unit_ts)) - offset;




end % of function