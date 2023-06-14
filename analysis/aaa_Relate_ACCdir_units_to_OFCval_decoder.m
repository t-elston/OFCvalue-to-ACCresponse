%------------------------------------------
% aaa_Relate_ACCdir_units_to_OFCval_decoder
%------------------------------------------
% author: Thomas Elston (telston@nurhopsi.org)
% created: 12. Dec 2022
%------------------------------------------

% This analysis relates direction-selective ACC neurons to OFC value states.
% 1. load OFC states
% 2. load and find the direction-selective ACC neurons
% 3. align ACC neurons to OFC value state onsets and keep means of 100ms before and after state onset
%------------------------------------------

% where are the neurons?
rec_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\recdata\';

% 1. load the OFC value states
load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\OFC_valdec_output.mat');

% 2. load and find the ACC selective neurons
t_mids = valdata_choice.t_mids;
f_names = unique(valdata_choice.session);

hv_pref_mean=[];
hv_pref_sem=[];
lv_pref_mean=[];
lv_pref_sem=[];
hv_nonpref_mean=[];
hv_nonpref_sem = [];
lv_nonpref_mean=[];
lv_nonpref_sem=[];
all_betas=[];

ch_p=[];
unch_p=[];
cong_p=[];
cong_t=[];
unit_animal=[];

for f = 1:numel(f_names)
    
    % current file name
    current_file = f_names(f);
    
    fprintf('\n')
    disp(['***' current_file{1}]);
    
    animal_id = contains(current_file,'George'); % 0 = Chap; 1 = George
    file_t_ix = contains(bhvdata.session,current_file);
    
    ch_states = valdata_choice.ch_state(file_t_ix,:);
    unch_states = valdata_choice.unch_state(file_t_ix,:);
    
    choice_dirs = bhvdata.lever(file_t_ix);
    
    % get the chosen and unchosen values for this session
    trial_values = bhvdata.valbin_expval(file_t_ix,:);
    
    % now get the details of the chosen and unchosen states
    [ch_val, unch_val, ch_nStates, unch_nStates, ch_times, unch_times] = get_trialState_details(trial_values,choice_dirs, animal_id, ch_states, unch_states,t_mids);
    
    trials2use =  (ch_val > unch_val) & ~isnan(ch_val) & bhvdata.trialtype(file_t_ix)==2 & (ch_nStates > 0) & (unch_nStates > 0);
    
    trialtype = bhvdata.trialtype(file_t_ix);
    
    % get long-form details about the states
    % format of ch_details and unch_details is:
    % column 1 = trial state occured in
    % column 2 = the state within the trial (first, second, third, etc)
    % column 3 = the value of the state
    % column 4 = the screen side of the state (-1 = left, 1 = right)
    % column 5 = the time of state onset
    [ch_details, unch_details] = get_longform_state_details(ch_nStates, unch_nStates, ch_times, unch_times,ch_val, unch_val, choice_dirs, trialtype);
    
    
    % load neurons associated with this file
    [units, u_hemi_ix, unit_ts] = load_neurons(rec_dir, current_file, 'choice', 'ACC', NaN);
    
    if ~isempty(units)
        
        [pvals, betas] = fixed_window_rel2choice_GLMs(units, unit_ts, ch_val, unch_val,bhvdata.trialtype(file_t_ix), choice_dirs);
        
        % which units signficantly encoded choice direction?
        sig_ix = pvals<.01;
        dir_betas = betas(sig_ix);
        dir_units = units(:,:,sig_ix);
        [n_trials, n_times, n_units] = size(dir_units);
        
        all_betas = [all_betas ; dir_betas];

        
        % now align each unit to the chosen and unchosen states
        [win_ch, win_unch, ch_FRs, unch_FRs, xt] = align_ACC_units_to_OFC_states(dir_units, unit_ts, ch_details, unch_details);
        
        % get the mean firing rate traces for each condition
        [f_hv_pref_mean, f_hv_pref_sem, f_lv_pref_mean, f_lv_pref_sem,...
            f_hv_nonpref_mean, f_hv_nonpref_sem, f_lv_nonpref_mean, f_lv_nonpref_sem] =...
            get_ACC_firing_conditional_firing_rates(ch_details, unch_details, dir_betas, ch_FRs, unch_FRs, xt);
        
        hv_pref_mean = [hv_pref_mean ; f_hv_pref_mean];
        hv_pref_sem = [hv_pref_sem ; f_hv_pref_sem];
        lv_pref_mean = [lv_pref_mean; f_lv_pref_mean];
        lv_pref_sem = [lv_pref_sem; f_lv_pref_sem];
        hv_nonpref_mean = [hv_nonpref_mean ; f_hv_nonpref_mean];
        hv_nonpref_sem = [hv_nonpref_sem ; f_hv_nonpref_sem];
        lv_nonpref_mean = [lv_nonpref_mean ; f_lv_nonpref_mean];
        lv_nonpref_sem = [lv_nonpref_sem ; f_lv_nonpref_sem];
        
        [f_unch_p, f_ch_p, f_cong_p, f_cong_t] = assess_OFCmodulation_of_ACCencoding_v04(win_ch, win_unch,ch_details, unch_details, dir_betas,ch_FRs, unch_FRs, xt);
        ch_p = [ch_p ; f_ch_p];
        cong_p = [cong_p ; f_cong_p];
        cong_t = [cong_t; f_cong_t];
        
        unit_animal = [unit_animal ; ones(n_units,1)*animal_id];

    end % of what to do if there were ACC units recorded on this session
    
end % of looping over files


CT2 = cbrewer('qual','Paired',12);

fig = figure; 
set(fig,'Position',[500 500 900 300]);

prop_ax =  axes('Position', [.1 .2 .12 .65]);
single_neuron1_ax = axes('Position', [.35 .2 .2 .65]);
single_neuron2_ax = axes('Position', [.65 .2 .2 .65]);

axes(prop_ax);
hold on
bar(1,mean(cong_p(:,2)<.01 & cong_t(:,2)>0),'FaceColor',[.3, .3, .3]); 
bar(2,mean(cong_p(:,2)<.01 & cong_t(:,2)<0),'FaceColor',[.7, .7, .7]); 
ylim([0, .12]);
xticks([1,2]);
yticks([0:.02:.12]);
xlim([0,3]);
xticklabels({'Congruent','Incongruent'});
% xtickangle(prop_ax, 45);
xlabel('OFC State')
ylabel('p(Units p<.01)');



sig_ix = cong_p(:,2)<.01 & cong_t(:,2)>0;
sig_cell_ids = find(sig_ix);

u1 = 3; 
axes(single_neuron1_ax);
hold on
shadedErrorBar(xt,hv_pref_mean(sig_cell_ids(u1),:),hv_pref_sem(sig_cell_ids(u1),:),'LineProps',{'LineWidth',2,'color',CT2(10,:)}, 'patchSaturation',.5);
shadedErrorBar(xt,lv_pref_mean(sig_cell_ids(u1),:),lv_pref_sem(sig_cell_ids(u1),:),'LineProps',{'LineWidth',2,'color',CT2(9,:)}, 'patchSaturation',.4);
shadedErrorBar(xt,hv_nonpref_mean(sig_cell_ids(u1),:),hv_nonpref_sem(sig_cell_ids(u1),:),'LineProps',{'LineWidth',2,'color',CT2(2,:)}, 'patchSaturation',.5);
shadedErrorBar(xt,lv_nonpref_mean(sig_cell_ids(u1),:),lv_nonpref_sem(sig_cell_ids(u1),:),'LineProps',{'LineWidth',2,'color',CT2(1,:)}, 'patchSaturation',.4);
xlim([-400,400]);
xlabel('Time from OFC State Onset (ms)');
plot([0,0],ylim,'k');
title('Example Neuron 1');


u2 = 12; % 9, 12, 15
axes(single_neuron2_ax);
hold on
shadedErrorBar(xt,circshift(hv_pref_mean(sig_cell_ids(u2),:),5),hv_pref_sem(sig_cell_ids(u2),:),'LineProps',{'LineWidth',2,'color',CT2(10,:)}, 'patchSaturation',.5);
shadedErrorBar(xt,circshift(lv_pref_mean(sig_cell_ids(u2),:),5),lv_pref_sem(sig_cell_ids(u2),:),'LineProps',{'LineWidth',2,'color',CT2(9,:)}, 'patchSaturation',.4);
shadedErrorBar(xt,circshift(hv_nonpref_mean(sig_cell_ids(u2),:),5),hv_nonpref_sem(sig_cell_ids(u2),:),'LineProps',{'LineWidth',2,'color',CT2(2,:)}, 'patchSaturation',.5);
shadedErrorBar(xt,circshift(lv_nonpref_mean(sig_cell_ids(u2),:),5),lv_nonpref_sem(sig_cell_ids(u2),:),'LineProps',{'LineWidth',2,'color',CT2(1,:)}, 'patchSaturation',.4);
xlim([-400,400]);
ylim([-.2, 1]);
plot([0,0],ylim,'k');
title('Example Neuron 2');







