%%% aaa_Assess_OFC_DIRdecoding_UnitEncoding

%-------------------------------
bhv_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\bhv\';
rec_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\recdata\';
alignment = 'choice';
brain_area = 'OFC';
offsets=[-1000,600];
%-------------------------------

% find names of bhv files and rec files
bhv_folder_names = get_subfolder_names(bhv_dir);

% load the OFC value decoder output and see which files were used. We'll use those same files for direction-decoding.
load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\OFC_valdec_output.mat');

files_for_decoding = unique(bhvdata.session);

n_files = numel(bhv_folder_names);

all_u_sig=[];
all_u_beta=[];
all_u_animal=[];
all_u_hemi=[];
%-------
% initialize decoder output
all_ch_dir = [];
all_dir_animal = [];
ctr=0;

% now loop over each file
for f_ix = 1:n_files
    
    current_file = bhv_folder_names(f_ix);
    monkey_id = contains(current_file,'George'); % 0 = Chap; 1 = George
    
    fprintf('\n')
    disp(['***' current_file{1}]);
    
    % pull out this file's behavior
    [trialinfo] = load_behavior(bhv_dir, current_file);
    
    % pull out this file's neurons
    [units, u_hemi_ix, choice_ts] = load_neurons(rec_dir, current_file, alignment, brain_area, offsets);
    
    % make sure there are in fact neurons for this file
    if ~isempty(units)
        
        % now loop over times and units
        [n_trials, n_times, n_units] = size(units);
                
        [t_ch_val, t_unch_val] = get_session_chosen_unchosen_vals(trialinfo);
        
        % now loop over each unit and check out direction encoding
        reg_tbl = table;
        reg_tbl.trial_type = trialinfo.trialtype;
        reg_tbl.trial_type(reg_tbl.trial_type==2) = -1;
        reg_tbl.max_val = nanmax(trialinfo.subjval_expval ,[],2);
        reg_tbl.t_num = trialinfo.TrialNumber;
        

        
%         % initialize output of GLMs
%         u_sig=NaN(n_units,3);
%         u_beta=NaN(n_units,3);
%         
%         % now loop over the units
%         pw = PoolWaitbar(n_units, 'Assessing file units...');
%         parfor u = 1:n_units
%             increment(pw);
%             
%             [u_sig(u,:), u_beta(u,:)] =parallel_sliding_GLM_v02(reg_tbl, trialinfo.lever*u_hemi_ix(u), units(:,:,u), choice_ts);
%             
%         end % of looping over units
%         delete(pw);
%         
%         % accumulate results
%         all_u_sig = [  all_u_sig ; u_sig];
%         all_u_beta = [  all_u_beta ; u_beta];
%         all_u_animal = [all_u_animal ; ones(n_units, 1)*monkey_id];
%         all_u_hemi = [all_u_hemi ; u_hemi_ix];
% %         
        ch_dir = zeros(size(trialinfo.lever));
        ch_dir(trialinfo.lever == 1) = 2;
        ch_dir(trialinfo.lever == -1) = 1;
        
        unch_dir = zeros(size(trialinfo.lever));
        unch_dir(trialinfo.lever == 1) = 1;
        unch_dir(trialinfo.lever == -1) = 2;

        %--------------------------------
        % now check whether we should do direction decoding for this file
        if any(contains(files_for_decoding, current_file))
            ctr = ctr+1;
            
            disp('...running decoder');
            
            n_boots = 100;
            
            s_ch_dir_pps = NaN(n_trials, n_times, n_boots);
            
            % only do direction decoding on free-choice trials
            free_ix = trialinfo.trialtype==2 & trialinfo.lever ~=0;
            
            acc=[];
            % do bootstraps
            for b = 1:n_boots
                
                % now grab 80% of trials for the test-set
                [candidate_train_ix, candidate_test_ix]  = GrabRandomProportionOfArray(find(free_ix),.8);
                
                candidate_test_ix(ch_dir(candidate_test_ix)==0) = [];
                
                candidate_train_trials = zeros(n_trials,1);
                candidate_train_trials(candidate_train_ix) = 1; 
                candidate_train_trials = logical(candidate_train_trials);
                
                candidate_test_trials = zeros(n_trials,1);
                candidate_test_trials(candidate_test_ix) = 1; 
                candidate_test_trials = logical(candidate_test_trials);
                
                % now pull a balanced set of left/right trials from these candidate trials
                [train_trials, ~]  = BalanceTrainSet_v01(candidate_train_trials,{t_ch_val, t_unch_val, trialinfo.lever});
                [test_trials, ~]  = BalanceTrainSet_v01(candidate_test_trials,{trialinfo.lever});

                
                train_labels = ch_dir(train_trials);
                test_labels = ch_dir(test_trials);
                
                % now loop over timesteps for training/testing the decoder
                predicted_labels=[];
                PPs=[];
                for t = 1:n_times
                    
                    % define the PC space on the training data
                    [train_pcs,coeff,PC95,mu] = getBestPCs_v01(squeeze(units(train_trials,t,:)));
                    
                    % project the test data into that PC space
                    [test_pcs] = ProjectIntoPCspace_v02(squeeze(units(test_trials,t,:)),mu,coeff,PC95);
                    
                    % fit the decoder
                    dir_LDA_mdl = fitcdiscr(train_pcs,train_labels,'FillCoeffs','off');
                    
                    % test the decoder
                    [predicted_labels(:,t), PPs(:,t,:), ~] = predict(dir_LDA_mdl,test_pcs);
                    
                end % of looping over times
                
                % get the posteriors associated with chosen and unchosen directions
                b_ch_dir_pps=[];
                [b_ch_dir_pps, ~] = GetChosenUnChosenPPs_v01(PPs,ch_dir(test_trials),unch_dir(test_trials));
                
                s_ch_dir_pps(test_trials, :,b) = b_ch_dir_pps;
                
                % get accuracy
                acc(b,:) = nanmean(predicted_labels == test_labels);
                
            end % of looping over bootstraps
            
            % take the mean over the bootstraps
            s_ch_dir_pps = nanmean(s_ch_dir_pps,3);
            
            all_ch_dir = [all_ch_dir; s_ch_dir_pps];
            all_dir_animal = [all_dir_animal; ones(n_trials,1)*monkey_id];
            
            session_acc(ctr,:) = nanmean(acc);
            session_animal(ctr, 1) = monkey_id; 
            
        end % of doing direction decoding for this file

    end % of ensuring there are units for this file
    
end % of looping over files

% c_ix = all_u_animal ==0;
% g_ix = all_u_animal ==1; 
% 
% dir_sig = all_u_sig(:,1)==1;
% val_sig = all_u_sig(:,2)==1;
% type_sig = all_u_sig(:,3)==1;
% 
% g_any_sig = sum(dir_sig(g_ix) | val_sig(g_ix) | type_sig(g_ix));
% c_any_sig = sum(dir_sig(c_ix) | val_sig(c_ix) | type_sig(c_ix));
% 
% 
% sum(~val_sig(g_ix) & dir_sig(g_ix) & type_sig(g_ix))
% 
% 
% g_n_value = sum(all_u_sig(g_ix, 2));
% g_n_dir = sum(all_u_sig(g_ix, 1));
% 
% g_pos_value = sum(all_u_sig(g_ix, 2) & all_u_beta(g_ix,2) > 0);
% 
% g_n_dir = sum(all_u_sig(g_ix, 1));
% g_contra_dir = sum(all_u_sig(g_ix, 1) & all_u_beta(g_ix,1) < 0);
% g_left_dir = sum(all_u_sig(g_ix, 1) & (all_u_beta(g_ix,1).*all_u_hemi(g_ix)) < 0);
% 
% g_n_trialtype = sum(all_u_sig(g_ix, 3));
% 
% 
% c_n_value = sum(all_u_sig(c_ix, 2));
% c_pos_value = sum(all_u_sig(c_ix, 2) & all_u_beta(c_ix,2) > 0);
% 
% c_n_dir = sum(all_u_sig(c_ix, 1));
% c_contra_dir = sum(all_u_sig(c_ix, 1) & all_u_beta(c_ix,1) < 0);
% c_left_dir = sum(all_u_sig(c_ix, 1) & (all_u_beta(c_ix,1).*all_u_hemi(c_ix)) < 0);
% 
% c_n_trialtype = sum(all_u_sig(c_ix, 3));

%------------------


OFC_val_PPs = valdata_choice.ch_ppd;
OFC_val_ts = valdata_choice.t_mids; 
OFC_val_animal = double(contains(valdata_choice.subject,'George'));

load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\ACC_valdec_output.mat');
ACC_val_PPs = valdata_choice.ch_ppd;
ACC_val_ts = valdata_choice.t_mids; 
ACC_val_animal = double(contains(valdata_choice.subject,'George'));

load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\ACC_dirdec_output.mat');
ACC_dir_PPs = dirdata_choice.postprob_ch;
ACC_dir_ts = dirdata_choice.t_mids; 
ACC_dir_free = bhvdata.trialtype == 2; 
ACC_dir_animal = double(contains(dirdata_choice.subject,'George'));


% plot direction decoder posteriors
CT = cbrewer('qual','Set1',9);
[c_OFC_dir_mean, c_OFC_dir_sem] = GetMeanCI(all_ch_dir(all_dir_animal==0,:),'sem');
[g_OFC_dir_mean, g_OFC_dir_sem] = GetMeanCI(all_ch_dir(all_dir_animal==1,:),'sem');

[c_ACC_dir_mean, c_ACC_dir_sem] = GetMeanCI(ACC_dir_PPs(ACC_dir_animal==0 & ACC_dir_free,:),'sem');
[g_ACC_dir_mean, g_ACC_dir_sem] = GetMeanCI(ACC_dir_PPs(ACC_dir_animal==1 & ACC_dir_free,:),'sem');


[c_OFC_val_mean, c_OFC_val_sem] = GetMeanCI(OFC_val_PPs(OFC_val_animal==0,:),'sem');
[g_OFC_val_mean, g_OFC_val_sem] = GetMeanCI(OFC_val_PPs(OFC_val_animal==1,:),'sem');

[c_ACC_val_mean, c_ACC_val_sem] = GetMeanCI(ACC_val_PPs(ACC_val_animal==0,:),'sem');
[g_ACC_val_mean, g_ACC_val_sem] = GetMeanCI(ACC_val_PPs(ACC_val_animal==1,:),'sem');







CT = cbrewer('qual','Dark2',8);
figure; 
subplot(1,2,1); 
hold on
xlim([-1000, 500]);
ylim([.49, .6]);
yticks([.5, .55, .6])
shadedErrorBar(choice_ts,c_OFC_dir_mean-.01, c_OFC_dir_sem,'lineprops',{'color',CT(1,:),'LineWidth',2});     
shadedErrorBar(ACC_dir_ts,c_ACC_dir_mean, c_ACC_dir_sem,'lineprops',{'color',CT(3,:),'LineWidth',2});     
plot(xlim,[.5,.5],'k','LineStyle',':');
plot([0, 0], ylim,'k','LineStyle',':');
xlabel('Time from Choice (ms)');
ylabel('Direction PPs');
title('Animal C');
axis square


subplot(1,2,2); 
hold on
xlim([-1000, 500]);
ylim([.49, .6]);
yticks([.5, .55, .6])
shadedErrorBar(choice_ts,g_OFC_dir_mean-.01, g_OFC_dir_sem,'lineprops',{'color',CT(1,:),'LineWidth',2}); 
shadedErrorBar(ACC_dir_ts,g_ACC_dir_mean, g_ACC_dir_sem,'lineprops',{'color',CT(3,:),'LineWidth',2});     
plot(xlim,[.5,.5],'k','LineStyle',':');
plot([0, 0], ylim,'k','LineStyle',':');
title('Animal G');
axis square




