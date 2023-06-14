function [ch_states, unch_states] = run_value_LDA_and_get_states(unit_data, trialinfo, ts)

[n_trials, n_times, n_units] = size(unit_data);

% initialize output
ch_states = NaN(n_trials, n_times);

% get indices of window for firing rates used to train decoder
[~,win_start] = min(abs(ts - 100));
[~,win_end] = min(abs(ts - 400));

% get the mean firing rates for the window used for training the decoder
win_FRs = squeeze(nanmean(unit_data(:,win_start:win_end,:),2));

% now find the chosen/unchosen value and images for each trial
[ch_val, unch_val, ch_img, unch_img] = get_chosen_values_and_images(trialinfo);

% find the trials which could be used for training the decoder: forced choice with a response
candidate_train_ix = trialinfo.trialtype==1 & trialinfo.lever~= 0 & ~isnan(ch_val);
free_trials = trialinfo.trialtype==2 & trialinfo.lever~= 0 & ~isnan(ch_val);

% now let's do some bootstrapping and get mean posteriors

n_boots = 50;

for b = 1:50
    
    % get this bootstrap's training trials 
    [train_trials, ~]  = BalanceTrainSet_v01(candidate_train_ix,{ch_img, trialinfo.lever});
    
    train_labels = ch_val(train_trials);
    
    % define the PC space for the training trials
    [train_pcs,coeff,PC95,mu] = getBestPCs_v01(win_FRs(train_trials,:));
    
    % train the decoder
    LDA_mdl = fitcdiscr(train_pcs,train_labels,'FillCoeffs','off');
    
    % now loop over each timepoint and test the decoder
    for t = 1:n_times
        
        % project the test data into that PC space
        [test_pcs] = ProjectIntoPCspace_v02(squeeze(unit_data(:,t,:)),mu,coeff,PC95);
        
        [~, PPs(:,t,:), ~] = predict(LDA_mdl,test_pcs);
        
    end % of looping over times
      
    % set the training trials to NaNs for this bootstrap
    PPs(train_trials,:,:) = NaN;
    
    % calculate baselines for the value states when they were unavailble on a trial
    [PP_baselines] = get_value_PPs_when_not_available(PPs(free_trials,:,:), ch_val(free_trials), unch_val(free_trials));
    
    % rescale the data to the baselines
    [rescaled_PPs] = rescale_PPs_to_baseline(PPs, PP_baselines); 
        
    [ch_pps(:,:,b), unch_pps(:,:,b)] = GetChosenUnChosenPPs_v01(rescaled_PPs,ch_val,unch_val);

end % of looping over bootstraps

% now take the mean over the bootstraps
mean_ch_pps = nanmean(ch_pps,3);
mean_unch_pps = nanmean(unch_pps,3);

% now find the states on each trial based on the chosen posteriors
[ch_states] = get_ch_val_states(mean_ch_pps, 200, 3);
[unch_states] = get_ch_val_states(mean_unch_pps, 200, 3);

ch_states = ch_states(:,1:n_times);
unch_states = unch_states(:,1:n_times);


end % of function