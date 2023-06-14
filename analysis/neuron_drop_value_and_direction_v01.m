function [val_acc, dir_acc, n_neurons] = neuron_drop_value_and_direction_v01(val_win_FRs, dir_win_FRs, trialinfo)

[n_trials, n_units] = size(val_win_FRs);
unit_shuffle_ix = shuffle([1:n_units]);

[t_ch_val, t_unch_val] = get_session_chosen_unchosen_vals(trialinfo);

free_ix = trialinfo.trialtype==2 & trialinfo.lever ~=0 & t_ch_val ~= t_unch_val & ~isnan(t_ch_val);
forced_ix = trialinfo.trialtype==1 & trialinfo.lever ~=0 & ~isnan(t_ch_val);

%----
% partition the data into training and testing sets for the direction decoder
%----
% grab 80% of trials for the test-set
[candidate_train_ix, dir_test_trials]  = GrabRandomProportionOfArray(find(free_ix),.8);

candidate_train_trials = zeros(numel(t_ch_val),1);
candidate_train_trials(candidate_train_ix) = 1; candidate_train_trials = logical(candidate_train_trials);

% now pull a balanced set of left/right trials from these candidate trials
[dir_train_trials, ~]  = BalanceTrainSet_v01(candidate_train_trials,{t_ch_val,t_unch_val, trialinfo.lever});

dir_train_labels = trialinfo.lever(dir_train_trials);
dir_test_labels = trialinfo.lever(dir_test_trials);

%----
% partition the data into training and testing sets for the value decoder
%----
[val_train_trials, ~]  = BalanceTrainSet_v01(forced_ix,{t_ch_val, trialinfo.lever});
val_test_trials = free_ix;

val_train_labels = t_ch_val(val_train_trials);
val_test_labels = t_ch_val(val_test_trials);

%---------------------------
% *****start neuron dropping
ix =0;
for i = 5:5:130
    
    ix = ix+1;
    n_neurons(1,ix) = i;
    
    % be sure there are enough neurons
    if i <= n_units
        
        % grab the random neurons for this iteration
        val_i_units = val_win_FRs(:,unit_shuffle_ix(1:i));
        dir_i_units = dir_win_FRs(:,unit_shuffle_ix(1:i));
        
        % define the PC space on the training data
        [val_train_pcs,val_coeff,val_PC95,val_mu] = getBestPCs_v01(val_i_units(val_train_trials,:));
        [dir_train_pcs,dir_coeff,dir_PC95,dir_mu] = getBestPCs_v01(dir_i_units(dir_train_trials,:));
        
        % project the test data into that PC space
        [val_test_pcs] = ProjectIntoPCspace_v02(val_i_units(val_test_trials,:),val_mu,val_coeff,val_PC95);
        [dir_test_pcs] = ProjectIntoPCspace_v02(dir_i_units(dir_test_trials,:),dir_mu,dir_coeff,dir_PC95);
        
        % train the decoders
        val_LDA = fitcdiscr(val_train_pcs,val_train_labels);
        dir_LDA = fitcdiscr(dir_train_pcs,dir_train_labels);
        
        % test the decoders
        val_predicted_labels = predict(val_LDA,val_test_pcs);
        dir_predicted_labels = predict(dir_LDA,dir_test_pcs);
        
        % get the decoding accuracies
        val_acc(1,ix) = mean(val_predicted_labels == val_test_labels);
        dir_acc(1,ix) = mean(dir_predicted_labels == dir_test_labels);
          
    else
        
        val_acc(1,ix) = NaN;
        dir_acc(1,ix) = NaN;
        
    end % of ensuring there are enough neurons
    
end % of looping over numbers of neurons


end % of function