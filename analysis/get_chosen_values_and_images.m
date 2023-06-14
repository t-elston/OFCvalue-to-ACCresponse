function [ch_val, unch_val, ch_img, unch_img] = get_chosen_values_and_images(trialinfo)

n_trials = numel(trialinfo.TrialNumber);

% initialize output
ch_val = NaN(n_trials,1);
unch_val = NaN(n_trials,1);
ch_img = NaN(n_trials,1);
unch_img = NaN(n_trials,1);

for t = 1:n_trials
    
    if trialinfo.lever(t) == -1
        ch_val(t,1) = trialinfo.valbin_expval(t,1);
        ch_img(t,1) = trialinfo.jpg(t,1);
        unch_val(t,1) = trialinfo.valbin_expval(t,2);
        unch_img(t,1) = trialinfo.jpg(t,2);
    end
    
    if trialinfo.lever(t) == 1
        ch_val(t,1) = trialinfo.valbin_expval(t,2);
        ch_img(t,1) = trialinfo.jpg(t,2);
        unch_val(t,1) = trialinfo.valbin_expval(t,1);
        unch_img(t,1) = trialinfo.jpg(t,1);
    end
    
end % of looping over trials


end % of function