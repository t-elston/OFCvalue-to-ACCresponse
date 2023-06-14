function [no_val_acc, acc, n_neurons] = compare_val_neuron_dropping(win_FRs, pos_val_units, train_trials, test_trials,...
                                                                                                train_labels, test_labels)

all_units = win_FRs;
no_val_units = win_FRs(:,~pos_val_units);

% now figure out how many neurons there are
[n_trials, n_all_units] = size(all_units);
[~, n_no_val_units] = size(no_val_units);

% now generate some random indices for grabbing the neurons
all_unit_shuffle_ix = shuffle([1:n_all_units]);
no_val_unit_shuffle_ix = shuffle([1:n_no_val_units]);

ix = 0;
for i = 5:5:200
    ix = ix+1;
    n_neurons(1,ix) = i;
    
    % get the neurons for this iteration
    if i <= n_all_units
        i_all_units = all_units(:,all_unit_shuffle_ix(1:i));
        
        % define the PC space on the training data
        [train_pcs,coeff,PC95,mu] = getBestPCs_v01(i_all_units(train_trials,:));
        
        % project the test data into that PC space
        [test_pcs] = ProjectIntoPCspace_v02(i_all_units(test_trials,:),mu,coeff,PC95);
        
        % fit the decoder
        LDA_mdl = fitcdiscr(train_pcs,train_labels);
        
        % test the decoder
        predicted_labels = predict(LDA_mdl,test_pcs);
        
        acc(1,ix) = mean(predicted_labels == test_labels);
        
    else
        acc(1,ix) = NaN;
        
        
    end
    

   %---------------------------------
   % look at neurons without positive value
    if i <= n_no_val_units
        i_no_val_units = no_val_units(:,no_val_unit_shuffle_ix(1:i));
        
        % define the PC space on the training data
        [no_val_train_pcs,no_val_coeff,no_val_PC95,no_val_mu] = getBestPCs_v01(i_no_val_units(train_trials,:));
        
        % project the test data into that PC space
        [no_val_test_pcs] = ProjectIntoPCspace_v02(i_no_val_units(test_trials,:),no_val_mu,no_val_coeff,no_val_PC95);
        
        % fit the decoder
        no_val_LDA_mdl = fitcdiscr(no_val_train_pcs,train_labels);
        
        % test the decoder
        no_val_predicted_labels = predict(no_val_LDA_mdl,no_val_test_pcs);
        
        no_val_acc(1,ix) = mean(no_val_predicted_labels == test_labels);
    else
        no_val_acc(1,ix) = NaN;
          
    end
   
end % of looping over including different numbers of neurons

% figure; 
% hold on
% plot(acc); 
% plot(no_val_acc);


            
end % of function