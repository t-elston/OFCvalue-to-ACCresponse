% aaa_Relate_OFCval_and_ACCdir_decoders

% load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\OFC_valdec_output.mat');
load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\ACC_valdec_output.mat');
load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\ACC_dirdec_output.mat');

% only keep the trials that were the same across both recording types
simultaneous_sessions = intersect(valdata_choice.session, dirdata_choice.session);

t_mids = valdata_choice.t_mids+1;
n_trials = numel(valdata_choice.ch_state(:,1));

% [~,choice_on_time_ix] = min(abs(t_mids - -450));
% [~,max_time_ix] = min(abs(t_mids - -50));
% figure; 
% subplot(1,2,1);
% hold on 
% plot(t_mids, nanmean(valdata_choice.ch_state(contains(valdata_choice.subject, 'Chap'),:)),'LineWidth',2)
% yyaxis right
% plot(t_mids, nanmean(dirdata_choice.postprob_ch(contains(dirdata_choice.subject, 'Chap'),:)),'LineWidth',2)
% xlim([-1000,500]);
% xlabel('Time from Choice (ms)');
% legend('OFC val','ACC dir');
% title('Animal C');
% axis square
% 
% subplot(1,2,2);
% hold on 
% plot(t_mids, nanmean(valdata_choice.ch_state(contains(valdata_choice.subject, 'George'),:)),'LineWidth',2)
% yyaxis right
% plot(t_mids, nanmean(dirdata_choice.postprob_ch(contains(dirdata_choice.subject, 'George'),:)),'LineWidth',2)
% xlim([-1000,500]);
% title('Animal G');
% axis square




ACC_ch=[];
ACC_unch=[];
ACC_na=[];

ch_val=[];
unch_val=[];

ch_animal_id = [];
unch_animal_id = [];
na_animal_id = [];

trial_ch_val=[];
trial_unch_val =[];
trial_animal=[];
trial_n_chStates=[];
trial_n_unchStates = [];

trial_ch_details=[];
trial_unch_details = [];

trial_number=[];

ch_state_times=[];
unch_state_times=[];

trial_ix = 0;


for s = 1:numel(simultaneous_sessions)
    
    val_s_ix = contains(bhvdata.session, simultaneous_sessions{s});
    dir_s_ix = contains(dirdata_choice.session, simultaneous_sessions{s});
        
    chosen_side = bhvdata.lever(val_s_ix)+1;
    chosen_side(chosen_side==0) = 1;
    unchosen_side = chosen_side;
    unchosen_side(chosen_side==1) = 2;
    unchosen_side(chosen_side==2) = 1;
    
    n_trials = numel(chosen_side); 
    
    trialtype = bhvdata.trialtype(val_s_ix);
    
    s_values = bhvdata.valbin_expval(val_s_ix,:);
    
    val_ch_state = valdata_choice.ch_state(val_s_ix,:);
    val_unch_state = valdata_choice.unch_state(val_s_ix,:);
    val_na_states = valdata_choice.na_state(val_s_ix,:,:);
    
    ACC_ch_pp = dirdata_choice.postprob_ch(dir_s_ix,:);
    
    animal_id = contains(simultaneous_sessions(s),'George');
    
    if animal_id == 0
        [~,state_win_start_time_ix] = min(abs(t_mids - -440)); 
        [~,max_time_ix] = min(abs(t_mids - -50));
    else
        [~,state_win_start_time_ix] = min(abs(t_mids - -420)); % was 425 or 260
        [~,max_time_ix] = min(abs(t_mids - -110)); % was 100
    end




    % go through each trial and find the start indices of each state and their associated values
    for t = 1:n_trials
        
        % what were the values for this trial?
        t_ch_val = s_values(t,chosen_side(t));
        t_unch_val = s_values(t,unchosen_side(t));
        
        ch_state = val_ch_state(t,:);
        
        % get this trial's unchosen state
        unch_state = val_unch_state(t,:);
        
        % get this trial's non-available state
        na_state = nanmean(val_na_states(t,:,:),3);
        
        % find when the chosen state was 'up'
        [~, ch_state_starts] = findpeaks(ch_state);
        
        % find when the unchosen state was 'up'
        [~, unch_state_starts] = findpeaks(unch_state);
        
        % find when the non-available state was 'up'
        [~, na_state_starts] = findpeaks(na_state);
        
        % get the indices of the state starts that happened prior to the choice
        ch_state_ix = ch_state_starts((ch_state_starts>=state_win_start_time_ix) & (ch_state_starts<=max_time_ix));
        unch_state_ix = unch_state_starts((unch_state_starts>=state_win_start_time_ix) & (unch_state_starts<=max_time_ix));
        na_state_ix = na_state_starts((na_state_starts>=state_win_start_time_ix) & (na_state_starts<=max_time_ix));
        
        % now get the ACC ch_pps for the 3 OFC alignments
        ch_pp=[];
        unch_pp=[];
        na_pp=[];
        
        n_ch_state = numel(ch_state_ix);
        n_unch_state = numel(unch_state_ix);
             
        % only include correct, free choices where each option was represented at least once
      if ~isnan(t_ch_val) & (t_unch_val < t_ch_val) & n_ch_state>0 & n_unch_state>0 & trialtype(t) ==2

          trial_ix = trial_ix+1;
            
            for ch_ix = 1:numel(ch_state_ix)
                
                ch_pp(ch_ix,:) = ACC_ch_pp(t, ch_state_ix(ch_ix)-80 : ch_state_ix(ch_ix)+ 80);
                
            end
            
            for unch_ix = 1:numel(unch_state_ix)
                
                unch_pp(unch_ix,:) = ACC_ch_pp(t, unch_state_ix(unch_ix)-80 : unch_state_ix(unch_ix)+ 80);
                
            end
            
            for na_ix = 1:numel(na_state_ix)
                
                na_pp(na_ix,:) = ACC_ch_pp(t, na_state_ix(na_ix)-80 : na_state_ix(na_ix)+ 80);
                
            end
            
            % now keep track of which animal this was
            if ~isempty(ch_pp)
                ch_animal_id = [ch_animal_id; ones(numel(ch_pp(:,1)),1) * animal_id];
                ch_val = [ch_val ; ones(numel(ch_pp(:,1)),1) * t_ch_val];
                
                ch_state_trial_n = ones(numel(ch_pp(:,1)),1) * trial_ix;
                ch_state_in_trial = [1:n_ch_state]';
                
                ch_state_times = [ch_state_times ; (t_mids(ch_state_ix))'];

                
                trial_ch_details = [trial_ch_details; [ch_state_trial_n , ch_state_in_trial] ];
                ACC_ch = [ACC_ch ; ch_pp];
                
                
            end
            if ~isempty(unch_pp) %& t_unch_val(t) ==3
                unch_animal_id = [unch_animal_id; ones(numel(unch_pp(:,1)),1) * animal_id];
                unch_val = [unch_val ; ones(numel(unch_pp(:,1)),1) * t_unch_val];
                

                unch_state_trial_n = ones(numel(unch_pp(:,1)),1) * trial_ix;
                unch_state_in_trial = [1:n_unch_state]';

                
                trial_unch_details = [trial_unch_details; [unch_state_trial_n , unch_state_in_trial] ];
                
                unch_state_times = [unch_state_times ; (t_mids(unch_state_ix))'];

                
                ACC_unch = [ACC_unch ; unch_pp];
            end
            if ~isempty(na_pp)
                na_animal_id = [na_animal_id; ones(numel(na_pp(:,1)),1) * animal_id];
                ACC_na = [ACC_na ; na_pp];
            end
            
            trial_number = [trial_number ; trial_ix];
            trial_ch_val = [trial_ch_val ; t_ch_val];
            trial_unch_val = [trial_unch_val ; t_unch_val];
            
            trial_n_chStates = [trial_n_chStates ; n_ch_state];
            trial_n_unchStates = [trial_n_unchStates ; n_unch_state];

            
            trial_animal = [trial_animal ; animal_id];


        end % of detemining whether to save these states for this trial
        
    end % of cycling over trials
    
end % of cycling over sessions
%---------------------------------------
xt = [1:numel(ACC_ch(1,:))]*mean(diff(t_mids)) - (round(numel(ACC_ch(1,:))/2) * mean(diff(t_mids)));
[~,xt_after_start] = min(abs(xt-0));
[~,xt_after_end] = min(abs(xt-100));
[~,xt_before_start] = min(abs(xt- -100));
[~,xt_before_end] = min(abs(xt-0));

%---------------------------------------

n_boots=10000;
chap_ch_ACC=NaN(n_boots, numel(xt));
chap_unch_ACC=NaN(n_boots, numel(xt));
chap_state_times = NaN(n_boots,2);

george_ch_ACC=NaN(n_boots, numel(xt));
george_unch_ACC=NaN(n_boots, numel(xt));
george_state_times = NaN(n_boots,2);

pw = PoolWaitbar(n_boots, 'bootstrapping...');
parfor b = 1:n_boots
            increment(pw);

[chap_boot_trials, ~] = BalanceTrainSet_v01(trial_animal==0,{trial_ch_val, trial_unch_val});
[george_boot_trials, ~] = BalanceTrainSet_v01(trial_animal==1,{trial_ch_val, trial_unch_val});

% now get the number of states on these trials
chap_boot_n_ChStates_trial = trial_n_chStates(chap_boot_trials);
chap_boot_n_UnChStates_trial = trial_n_unchStates(chap_boot_trials);

george_boot_n_ChStates_trial = trial_n_chStates(george_boot_trials);
george_boot_n_UnChStates_trial = trial_n_unchStates(george_boot_trials);


% find which states to use on each trial 
[chap_states2use] = get_adjacent_states_per_trial(ch_state_times, unch_state_times, chap_boot_trials,...
    trial_ch_details, trial_unch_details);

[george_states2use] = get_adjacent_states_per_trial(ch_state_times, unch_state_times, george_boot_trials,...
    trial_ch_details, trial_unch_details);

% % now randomly select one state from each trial
% [chap_boot_ch_state2use] = get_random_states_from_trial(chap_boot_n_ChStates_trial);
% [chap_boot_unch_state2use] = get_random_states_from_trial(chap_boot_n_UnChStates_trial);
% [george_boot_ch_state2use] = get_random_states_from_trial(george_boot_n_ChStates_trial);
% [george_boot_unch_state2use] = get_random_states_from_trial(george_boot_n_UnChStates_trial);

% now actually get the PPs associated with those states
CH_chap_boot_states = ismember(trial_ch_details(:,1), chap_boot_trials) & ismember(trial_ch_details(:,2), chap_states2use(:,2));
UNCH_chap_boot_states = ismember(trial_unch_details(:,1), chap_boot_trials) & ismember(trial_unch_details(:,2), chap_states2use(:,1));
CH_george_boot_states = ismember(trial_ch_details(:,1), george_boot_trials) & ismember(trial_ch_details(:,2), george_states2use(:,2));
UNCH_george_boot_states = ismember(trial_unch_details(:,1), george_boot_trials) & ismember(trial_unch_details(:,2), george_states2use(:,1));


chap_ch_ACC(b,:) = nanmean(ACC_ch(CH_chap_boot_states,:));
chap_unch_ACC(b,:) = nanmean(ACC_unch(UNCH_chap_boot_states,:));

chap_state_times(b,:) = [nanmean(ch_state_times(CH_chap_boot_states)) , nanmean(unch_state_times(UNCH_chap_boot_states))];

george_ch_ACC(b,:) = nanmean(ACC_ch(CH_george_boot_states,:));
george_unch_ACC(b,:) = nanmean(ACC_unch(UNCH_george_boot_states,:));

george_state_times(b,:) = [nanmean(ch_state_times(CH_george_boot_states)) , nanmean(unch_state_times(UNCH_george_boot_states))];

end % of looping over bootstraps
    delete(pw);
    
    % now take differences between the traces
    chap_ACC_diff = chap_ch_ACC - chap_unch_ACC;
    george_ACC_diff = george_ch_ACC - george_unch_ACC;


Effect = computeCohen_d(ch_state_times(ch_animal_id==1), unch_state_times(unch_animal_id==1))
Effect = computeCohen_d(ch_state_times(ch_animal_id==0), unch_state_times(unch_animal_id==0))

mean(ch_state_times(ch_animal_id==1)) - mean(unch_state_times(unch_animal_id==1))
mean(ch_state_times(ch_animal_id==0)) - mean(unch_state_times(unch_animal_id==0))

[h,p,~,stats] = ttest2(ch_state_times(ch_animal_id==1), unch_state_times(unch_animal_id==1));

%-----
chap_after_ch_win_means = mean(chap_ch_ACC(:,xt_after_start:xt_after_end),2);
chap_after_unch_win_means = mean(chap_unch_ACC(:,xt_after_start:xt_after_end),2);

george_after_ch_win_means = mean(george_ch_ACC(:,xt_after_start:xt_after_end),2);
george_after_unch_win_means = mean(george_unch_ACC(:,xt_after_start:xt_after_end),2);
%-----



chap_after_p = sum(chap_after_ch_win_means < chap_after_unch_win_means) / n_boots;
george_after_p = sum(george_after_ch_win_means < george_after_unch_win_means) / n_boots;



[boot_chap_ch_mean, boot_chap_ch_sem] = GetMeanCI(chap_ch_ACC, 'percentile');
[boot_chap_unch_mean, boot_chap_unch_sem] = GetMeanCI(chap_unch_ACC, 'percentile');

[boot_george_ch_mean, boot_george_ch_sem] = GetMeanCI(george_ch_ACC, 'percentile');
[boot_george_unch_mean, boot_george_unch_sem] = GetMeanCI(george_unch_ACC, 'percentile');


CT = cbrewer('qual','Set1',9);


figure; 
subplot(2,2,1)
hold on
shadedErrorBar(xt,boot_chap_ch_mean, boot_chap_ch_sem,'LineProps',{'LineWidth',2,'color',CT(4,:)});
shadedErrorBar(xt,boot_chap_unch_mean, boot_chap_unch_sem,'LineProps',{'LineWidth',2,'color',CT(2,:)});
ylim([.54, .6]);
xlim([-100,200])
plot([0,0],ylim,'k','LineWidth',1);
xlabel('Time from ACC State Onset (ms)');
ylabel('DIRacc post. prob.');
axis square
title('Animal C');


subplot(2,2,2)
hold on
shadedErrorBar(xt,boot_george_ch_mean, boot_george_ch_sem,'LineProps',{'LineWidth',2,'color',CT(4,:)});
shadedErrorBar(xt,boot_george_unch_mean, boot_george_unch_sem,'LineProps',{'LineWidth',2,'color',CT(2,:)});
ylim([.52, .60]);
xlim([-100,200])

plot([0,0],ylim,'k','LineWidth',1);
axis square
title('Animal G');



[chap_ch_mean, chap_ch_sem] = GetMeanCI(ACC_ch(ch_animal_id==0,:), 'sem');
[chap_unch_mean, chap_unch_sem] = GetMeanCI(ACC_unch(unch_animal_id==0,:), 'sem');
[chap_na_mean, chap_na_sem] = GetMeanCI(ACC_na(na_animal_id==0,:), 'sem');

[george_ch_mean, george_ch_sem] = GetMeanCI(ACC_ch(ch_animal_id==1,:), 'sem');
[george_unch_mean, george_unch_sem] = GetMeanCI(ACC_unch(unch_animal_id==1,:), 'sem');
[george_na_mean, george_na_sem] = GetMeanCI(ACC_na(na_animal_id==1,:), 'sem');

g_p=[];
c_p=[];
% get moments of significance (1-way anovas) for the raw (non-bootstrapped) data
for xi = 1:numel(xt)
    
    g_data_xi = [ACC_ch(ch_animal_id==1, xi) ; ACC_unch(unch_animal_id==1, xi) ; ACC_na(na_animal_id==1, xi) ];
    c_data_xi = [ACC_ch(ch_animal_id==0, xi) ; ACC_unch(unch_animal_id==0, xi) ; ACC_na(na_animal_id==0, xi) ];

    g_p(xi) = anovan(g_data_xi, g_win_types,'varnames',{'state_type'},'Display','off');
    c_p(xi) = anovan(c_data_xi, chap_win_types,'varnames',{'state_type'},'Display','off');


end

c_p = double(c_p < .01); 
c_p(c_p==0) = NaN;

g_p = double(g_p < .01); 
g_p(g_p==0) = NaN;



figure;
subplot(1,2,1); 
hold on
shadedErrorBar(xt,chap_ch_mean, chap_ch_sem,'LineProps',{'LineWidth',2,'color',CT(4,:)});
shadedErrorBar(xt,chap_unch_mean, chap_unch_sem,'LineProps',{'LineWidth',2,'color',CT(2,:)});
shadedErrorBar(xt,chap_na_mean, chap_na_sem,'LineProps',{'LineWidth',2,'color',CT(9,:)});
plot(xt, c_p*.545, 'LineWidth',2,'color','k');
ylim([.54, .6]);
yticks([.54, .56, .58, .6]);
xlim([-100,200])
plot([0,0],ylim,'k','LineWidth',1);
xlabel('Time from ACC State Onset (ms)')
ylabel('DIRacc post. prob.')

axis square
title('Animal C'); 

subplot(1,2,2); 
hold on
shadedErrorBar(xt,george_ch_mean, george_ch_sem,'LineProps',{'LineWidth',2,'color',CT(4,:)});
shadedErrorBar(xt,george_unch_mean, george_unch_sem,'LineProps',{'LineWidth',2,'color',CT(2,:)});
shadedErrorBar(xt,george_na_mean, george_na_sem,'LineProps',{'LineWidth',2,'color',CT(9,:)});
plot(xt, g_p*.53, 'LineWidth',2,'color','k');
ylim([.52, .6]);
xlim([-100,200])
plot([0,0],ylim,'k','LineWidth',1);
axis square
title('Animal G'); 


%-------------

ACC_ch_means = mean(ACC_ch(:, xt_after_start:xt_after_end),2);
ACC_unch_means = mean(ACC_unch(:, xt_after_start:xt_after_end),2);
ACC_na_means = mean(ACC_na(:, xt_after_start:xt_after_end),2);


g_ch_window_mean = mean(ACC_ch(ch_animal_id==1,xt_after_start:xt_after_end),2);
g_unch_window_mean = mean(ACC_unch(unch_animal_id==1,xt_after_start:xt_after_end),2);
g_na_window_mean = mean(ACC_na(na_animal_id==1,xt_after_start:xt_after_end),2);

c_ch_window_mean = ACC_ch_means(ch_animal_id==0);
c_unch_window_mean = mean(ACC_unch(unch_animal_id==0,xt_after_start:xt_after_end),2);
c_na_window_mean = mean(ACC_na(na_animal_id==0,xt_after_start:xt_after_end),2);


chap_win_means = [ACC_ch_means(ch_animal_id==0); ACC_unch_means(unch_animal_id==0) ; ACC_na_means(na_animal_id==0)];
chap_win_types = [zeros(sum((ch_animal_id==0)),1) ; zeros(sum((unch_animal_id==0)),1)+1  ; zeros(sum((na_animal_id==0)),1)+2];

g_win_means = [g_ch_window_mean ; g_unch_window_mean ; g_na_window_mean];
g_win_types = [zeros(sum((ch_animal_id==1)),1) ; zeros(sum((unch_animal_id==1)),1)+1  ; zeros(sum((na_animal_id==1)),1)+2];

[~, chap_tbl, chap_stats] = anovan(chap_win_means,chap_win_types,'varnames',{'state_type'},'Display','off');
[~, george_tbl, george_stats] = anovan(g_win_means,g_win_types,'varnames',{'state_type'},'Display','off');

chap_multcomp = multcompare(chap_stats,'Display','off');
george_multcomp = multcompare(george_stats,'Display','off');




