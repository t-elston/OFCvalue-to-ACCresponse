% aaa_DirDecoding_withoutValue_units_v01
%-------------------------------
bhv_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\bhv\';
rec_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\recdata\';
alignment = 'choice';
brain_area = 'ACC';
offsets=[-1500, 1000];

load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\ACC_dirdec_output.mat');
files_for_decoding = unique(dirdata_choice.session);
%-------------------------------

% find names of bhv files and rec files
bhv_folder_names = get_subfolder_names(bhv_dir);

n_files = numel(files_for_decoding);

% now loop over each file
for f_ix = 1:n_files
    
    current_file = files_for_decoding(f_ix);
    monkey_id = contains(current_file,'George'); % 0 = Cap; 1 = George
    
    fprintf('\n') 
    disp(['***' current_file{1}]);
    
    % pull out this file's behavior
    [trialinfo] = load_behavior(bhv_dir, current_file);
    
    % pull out this file's neurons
    [units, u_hemi_ix, choice_ts] = load_neurons(rec_dir, current_file, alignment, brain_area, offsets);
    
    if ~isempty(units)
    
    [t_ch_val, t_unch_val] = get_session_chosen_unchosen_vals(trialinfo);
    
    % now loop over each unit and check out direction encoding
    reg_tbl = table;
    reg_tbl.trial_type = trialinfo.trialtype;
    reg_tbl.trial_type(reg_tbl.trial_type==2) = -1;
    reg_tbl.max_val = nanmax(trialinfo.subjval_expval ,[],2);
    reg_tbl.t_num = trialinfo.TrialNumber;
    
    % now loop over times and units
    [n_trials, n_times, n_units] = size(units);
    
    % initialize output of GLMs
    sig_times = NaN(n_units, n_times,2);
    t_betas = NaN(n_units, n_times,2);
    t_CPD = NaN(n_units, n_times,2);
    
    u_sig=NaN(n_units,2);
    u_beta=NaN(n_units,3);
    
    % now loop over the units
    pw = PoolWaitbar(n_units, 'Assessing file units...');
    parfor u = 1:n_units
        increment(pw);
        
        [sig_times(u,:,:), t_betas(u,:,:), ~, ~, u_sig(u,:), u_beta(u,:), t_CPD(u,:,:)] =...
            parallel_sliding_GLM_v01(reg_tbl, trialinfo.lever*u_hemi_ix(u), units(:,:,u), choice_ts);
        
    end % of looping over units
    delete(pw);
    
    % now find the units with a positive relationship to value
    pos_val_units = u_beta(:,3) > 0  & u_sig(:,2)==1;
    
    % OK, now we've found the units which significantly and positively code for value in this session
    % Now, let's train a decoder and specifically hold out those units, according to the paper, I'll do PCA and then 25 bootstraps
    % on just the FREE choice trials
   
    
    n_boots = 1000;

    free_ix = trialinfo.trialtype==2 & trialinfo.lever ~=0 & t_ch_val > t_unch_val;
    no_val_acc = NaN(n_boots, n_times);
    acc = NaN(n_boots, n_times);  
    no_val_drop_acc=[];
    all_drop_acc=[];
    
    parfor b = 1:n_boots
        
        % print a little message saying which bootstrap we're on
        disp(['b = ' num2str(b) ' / ' num2str(n_boots)]);
        
        
        % now grab 80% of trials for the test-set
        [candidate_train_ix, test_trials]  = GrabRandomProportionOfArray(find(free_ix),.8);
        
        candidate_train_trials = zeros(n_trials,1);
        candidate_train_trials(candidate_train_ix) = 1; candidate_train_trials = logical(candidate_train_trials);
        
        % now pull a balanced set of left/right trials from these candidate trials
        [train_trials, ~]  = BalanceTrainSet_v01(candidate_train_trials,{t_ch_val,t_unch_val, trialinfo.lever});
        
        train_labels = trialinfo.lever(train_trials);
        test_labels = trialinfo.lever(test_trials);
        
        % now, determine which monkey this is and pull a fixed 200ms window around the decoder peak for the neuron-dropping
        % analysis
        if monkey_id ==0
            win_FRs = squeeze(mean(units(:,50:58,:),2));
        else
            win_FRs = squeeze(mean(units(:,54:62,:),2));
        end
        
        % do the neuron-dropping analysis
        [no_val_drop_acc(b,:), all_drop_acc(b,:), n_neurons] = compare_val_neuron_dropping(win_FRs, pos_val_units, train_trials, test_trials,...
            train_labels, test_labels);
        
        % now loop over time steps
        for t = 1:n_times
            
            % first do the decoding with positive-value coding cells held out
            % define the PC space on the training data
            [no_val_train_pcs,no_val_coeff,no_val_PC95,no_val_mu] = getBestPCs_v01(squeeze(units(train_trials,t,~pos_val_units)));
            
            % project the test data into that PC space
            [no_val_test_pcs] = ProjectIntoPCspace_v02(squeeze(units(test_trials,t,~pos_val_units)),no_val_mu,no_val_coeff,no_val_PC95);
            
            % fit the decoder
            no_val_LDA_mdl = fitcdiscr(no_val_train_pcs,train_labels);
            
            % test the decoder
            no_val_predicted_labels = predict(no_val_LDA_mdl,no_val_test_pcs);
            
            no_val_acc(b,t) = mean(no_val_predicted_labels == test_labels);
            
            %----------------
            % now do the decoding with all cells
            % define the PC space on the training data
            [train_pcs,coeff,PC95,mu] = getBestPCs_v01(squeeze(units(train_trials,t,:)));
            
            % project the test data into that PC space
            [test_pcs] = ProjectIntoPCspace_v02(squeeze(units(test_trials,t,:)),mu,coeff,PC95);
            
            % fit the decoder
            LDA_mdl = fitcdiscr(train_pcs,train_labels);
            
            % test the decoder
            predicted_labels = predict(LDA_mdl,test_pcs);
            
            acc(b,t) = mean(predicted_labels == test_labels);
            
        end % of looping over times
        
    end % of looping over bootstraps
    
            % do the neuron-dropping analysis
    n_neurons = 5:5:200;
    
    f_no_val_acc(f_ix,:) = mean(no_val_acc);
    f_acc(f_ix,:) = mean(acc);
    f_animal(f_ix,1) = monkey_id;
    
    f_no_val_drop_acc(f_ix,:) = nanmean(no_val_drop_acc);
    f_drop_acc(f_ix,:) = nanmean(all_drop_acc);
    
    f_no_val_reg = regress(nanmean(no_val_drop_acc)', [ones(size(n_neurons')), n_neurons']);
    f_reg = regress(nanmean(all_drop_acc)', [ones(size(n_neurons')), n_neurons']);
    
    f_no_val_beta(f_ix) = f_no_val_reg(2);
    f_all_beta(f_ix) = f_reg(2);
    
    else
        
        f_no_val_acc(f_ix,:) = NaN;
        f_acc(f_ix,:) = NaN;
        f_animal(f_ix,1) = NaN;
        
        f_no_val_drop_acc(f_ix,:) = NaN;
        f_drop_acc(f_ix,:) = NaN;
       
        f_no_val_beta(f_ix) = NaN;
        f_all_beta(f_ix) = NaN;
        
    end % of making sure we only using things with units
    
    
end % of looping over files


[chap_no_val_mean, chap_no_val_sem] = GetMeanCI(f_no_val_acc(f_animal==0,:),'sem');
[chap_all_mean, chap_all_sem] = GetMeanCI(f_acc(f_animal==0,:),'sem');

[chap_drop_acc_mean, chap_drop_acc_sem] = GetMeanCI(f_drop_acc(f_animal==0,:),'sem');
[chap_no_val_drop_acc_mean, chap_no_val_drop_acc_sem] = GetMeanCI(f_no_val_drop_acc(f_animal==0,:),'sem');


[george_no_val_mean, george_no_val_sem] = GetMeanCI(f_no_val_acc(f_animal==1,:),'sem');
[george_all_mean, george_all_sem] = GetMeanCI(f_acc(f_animal==1,:),'sem');

[george_drop_acc_mean, george_drop_acc_sem] = GetMeanCI(f_drop_acc(f_animal==1,:),'sem');
[george_no_val_drop_acc_mean, george_no_val_drop_acc_sem] = GetMeanCI(f_no_val_drop_acc(f_animal==1,:),'sem');




CT = cbrewer('qual','Dark2',8);


figure; 
subplot(3,2,1); 
hold on
% plot(choice_ts, chap_all_mean,'LineWidth',2,'color',CT(1,:));
% plot(choice_ts, chap_no_val_mean,'LineWidth',2,'color',CT(3,:));
shadedErrorBar(choice_ts, chap_all_mean, chap_all_sem,'LineProps',{'LineWidth',2,'color',CT(1,:)});
shadedErrorBar(choice_ts, chap_no_val_mean, chap_no_val_sem,'LineProps',{'LineWidth',2,'color',CT(3,:)});
xlabel('Time from Choice (ms)');
ylabel('Decoder Accuracy');
title('Animal C');
xlim([-1000, 500]);
% legend('All Units', 'No +Val Units');
axis square


subplot(3,2,2); 
hold on
shadedErrorBar(choice_ts, george_all_mean, george_all_sem,'LineProps',{'LineWidth',2,'color',CT(1,:)});
shadedErrorBar(choice_ts, george_no_val_mean, george_no_val_sem,'LineProps',{'LineWidth',2,'color',CT(3,:)});
title('Animal G');
xlim([-1000, 500]);
axis square


subplot(3,2,3); 
hold on
plot(n_neurons, f_drop_acc(f_animal==0,:)','color',CT(1,:));
plot(n_neurons, f_no_val_drop_acc(f_animal==0,:)','color',CT(3,:));

shadedErrorBar(n_neurons, chap_drop_acc_mean, chap_drop_acc_sem,'LineProps',{'LineWidth',2,'color',CT(1,:)});
shadedErrorBar(n_neurons, chap_no_val_drop_acc_mean, chap_no_val_drop_acc_sem,'LineProps',{'LineWidth',2,'color',CT(3,:)});
ylabel('Decoder Accuracy');
xlabel('Number of Neurons');
axis square


subplot(3,2,4); 
hold on
plot(n_neurons, f_drop_acc(f_animal==1,:)','color',CT(1,:));
plot(n_neurons, f_no_val_drop_acc(f_animal==1,:)','color',CT(3,:));

shadedErrorBar(n_neurons, george_drop_acc_mean, george_drop_acc_sem,'LineProps',{'LineWidth',2,'color',CT(1,:)});
shadedErrorBar(n_neurons, george_no_val_drop_acc_mean, george_no_val_drop_acc_sem,'LineProps',{'LineWidth',2,'color',CT(3,:)});
ylabel('Decoder Accuracy');
xlabel('Number of Neurons');
axis square


subplot(3,2,5)
hold on
plot(ones(size(f_all_beta(f_animal==0))),f_all_beta(f_animal==0),'.','color',CT(1,:),'MarkerSize',20);
plot(ones(size(f_no_val_beta(f_animal==0)))*1.5,f_no_val_beta(f_animal==0),'.','color',CT(3,:),'MarkerSize',20);

chap_beta_pairs  = [ f_all_beta(f_animal==0)', f_no_val_beta(f_animal==0)'];
for p = 1:numel(chap_beta_pairs(:,1))
    
plot([1,1.5], chap_beta_pairs(p,:),'color',CT(8,:),'LineWidth',1);

end
xlim([.8,1.7]);
xticks([1, 1.5]);
xticklabels({'All Units', 'No +Val Units'});
ylabel('Info. Gain per Unit (Accuracy / Neuron)');
axis square
[~, c_pval, ~, c_stats] = ttest(chap_beta_pairs(:,1), chap_beta_pairs(:,2));



subplot(3,2,6)
hold on
plot(ones(size(f_all_beta(f_animal==1))),f_all_beta(f_animal==1),'.','color',CT(1,:),'MarkerSize',20);
plot(ones(size(f_no_val_beta(f_animal==1)))*1.5,f_no_val_beta(f_animal==1),'.','color',CT(3,:),'MarkerSize',20);

george_beta_pairs  = [ f_all_beta(f_animal==1)', f_no_val_beta(f_animal==1)'];
for p = 1:numel(george_beta_pairs(:,1))
    
plot([1,1.5], george_beta_pairs(p,:),'color',CT(8,:),'LineWidth',1);

end
xlim([.8,1.7]);
xticks([1, 1.5]);
xticklabels({'All Units', 'No +Val Units'});
axis square

[~, g_pval, ~, g_stats] = ttest(george_beta_pairs(:,1), george_beta_pairs(:,2));






xx=[];