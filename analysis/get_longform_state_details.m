function [ch_details, unch_details] = get_longform_state_details(ch_nStates, unch_nStates, ch_times,...
                                                                 unch_times,ch_val, unch_val, choice_dirs, trialtype)

% format of ch_details and unch_details is:
% column 1 = trial state occured in
% column 2 = the state within the trial (first, second, third, etc)
% column 3 = the value of the state
% column 4 = the screen side of the state (-1 = left, 1 = right)                                                                                            
% column 5 = the time of state onset                                                                                             
                                                                                             
ch_details=[];
unch_details=[];


chosen_dir = choice_dirs; 
unchosen_dir = choice_dirs*-1;

[n_trials, ~] = size(ch_nStates);

for t = 1:n_trials
    
    if  ch_nStates(t) > 0  
        ch_trial_num = ones(ch_nStates(t),1)*t;
        ch_state_in_trial = (1:ch_nStates(t))';   
        ch_state_value = ones(ch_nStates(t),1)*ch_val(t);
        ch_state_alt_value = ones(ch_nStates(t),1)*unch_val(t);
        ch_state_dir = ones(ch_nStates(t),1)*chosen_dir(t);
        ch_state_trialtype = ones(ch_nStates(t),1)*trialtype(t);
        
        t_ch_times = ch_times(t, 1:ch_nStates(t))';
        
        ch_details = [ch_details ; [ch_trial_num , ch_state_in_trial, ch_state_value, ch_state_dir,...
                                    t_ch_times, ch_state_alt_value, ch_state_trialtype] ];
        
    end
    
    if  unch_nStates(t) > 0
        unch_trial_num = ones(sum(unch_nStates(t)),1)*t;
        unch_state_in_trial = (1:unch_nStates(t))'; 
        unch_state_value = ones(unch_nStates(t),1)*unch_val(t);
        unch_state_alt_value = ones(unch_nStates(t),1)*ch_val(t);
        unch_state_dir = ones(unch_nStates(t),1)*unchosen_dir(t);
        unch_state_trialtype = ones(unch_nStates(t),1)*trialtype(t);

        t_unch_times = unch_times(t, 1:unch_nStates(t))';


        unch_details = [unch_details ; [unch_trial_num , unch_state_in_trial, unch_state_value, unch_state_dir,...
                                        t_unch_times, unch_state_alt_value, unch_state_trialtype] ];
        
    end
    
    
    

end % of looping over trials


end % of function