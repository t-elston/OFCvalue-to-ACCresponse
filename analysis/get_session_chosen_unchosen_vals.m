function [t_ch_val, t_unch_val] = get_session_chosen_unchosen_vals(trialinfo)

n_trials = numel(trialinfo.TrialNumber);

t_ch_val=[];
t_unch_val=[];

for t = 1:n_trials
    
    if trialinfo.lever(t) == -1
        
        t_ch_val(t,1) = trialinfo.valbin_expval(t,1);
        t_unch_val(t,1) = trialinfo.valbin_expval(t,2);
    end
    
    if trialinfo.lever(t) == 1
        
        t_ch_val(t,1) = trialinfo.valbin_expval(t,2);
        t_unch_val(t,1) = trialinfo.valbin_expval(t,1);

    end
    
    if trialinfo.lever(t) == 0
        
        t_ch_val(t,1) = NaN;
        t_unch_val(t,1) = NaN;
    end
    
    
end % of looping over trials

end % of function