function [pref_pp_pval, nonpref_pp_pval, hv_pref_FRs, lv_pref_FRs, hv_nonpref_FRs, lv_nonpref_FRs, pref_beta, nonpref_beta, hv_pref_sem, lv_pref_sem, hv_nonpref_sem, lv_nonpref_sem] =...
    regress_Dir_selective_ACCunits_against_OFCstates_v03(units, betas, ts, t_ch_state_times, t_ch_val, t_unch_state_times, t_unch_val, trial_types, chosen_side,...
    ch_pp, unch_pp, full_ch_pp, full_unch_pp, decoder_ts)


hv_pref_FRs=[];
lv_pref_FRs=[]; 
hv_nonpref_FRs=[]; 
lv_nonpref_FRs=[];
hv_pref_sem=[];
lv_pref_sem=[];
hv_nonpref_sem=[];
lv_nonpref_sem=[];

trials_with_ch_state = sum(~isnan(t_ch_state_times),2)>0;
trials_with_unch_state = sum(~isnan(t_unch_state_times),2)>0;

xx=[];
%---------------------------------------------
% identify trials to use
trials2keep = ~isnan(t_ch_val) & trial_types ==2 & (t_ch_val > t_unch_val) & t_unch_val > 1;


units = units(trials2keep,:,:);
t_ch_state_times = t_ch_state_times(trials2keep,:);
t_ch_val = t_ch_val(trials2keep);
t_unch_state_times = t_unch_state_times(trials2keep,:);
t_unch_val = t_unch_val(trials2keep);
chosen_side = chosen_side(trials2keep);
ch_pp = ch_pp(trials2keep,:);
unch_pp = unch_pp(trials2keep,:);
%---------------------------------------------

n_trials = numel(t_ch_val);
all_chState_FRs=[];
all_ch_pp=[];
all_full_ch_pp=[];
all_chState_chVal=[];
all_chState_unchVal=[];
all_chState_chSide=[];
all_chState_unchSide=[];
all_unchState_chSide=[];
all_unchState_unchSide=[];

all_unchState_FRs=[];
all_unch_pp=[];
all_unchState_unchVal=[];
all_unchState_chVal=[];
% loop through each trial and accumulate the firing rates


for t = 1:n_trials
    
    t_ch_FRs=[];
    t_ch_pp=[];
    t_unch_FRs=[];
    t_unch_pp=[];
    t_full_ch_pps=[];
    % get the times of the chosen states on this trial 
    t_ch_side = chosen_side(t);
    t_unch_side = t_ch_side*-1;

    t_ch_states = t_ch_state_times(t,:);
    t_ch_states(isnan(t_ch_states))=[];
    
    t_unch_states = t_unch_state_times(t,:);
    t_unch_states(isnan(t_unch_states))=[];
    
    % now loop over each state and find window aligned on state onset
    for ch_s = 1:numel(t_ch_states)
        [~,ch_state_start] = min(abs(ts - (t_ch_states(ch_s)-700)));
        [~,ch_state_end] = min(abs(ts - (t_ch_states(ch_s)+700)));
        
        t_ch_FRs(ch_s,:,:) = units(t,ch_state_start:ch_state_end,:);
        t_ch_pp(ch_s,:) = ch_pp(t, ch_state_start:ch_state_end);
        
        % for full frequency decoder data
        [~,full_ch_state_start] = min(abs(decoder_ts - (t_ch_states(ch_s)-700)));
        [~,full_ch_state_end] = min(abs(decoder_ts - (t_ch_states(ch_s)+700)));
        t_full_ch_pps(ch_s,:) = full_ch_pp(t,full_ch_state_start:full_ch_state_end);
    end
    
    % now loop over each state and find window aligned on state onset
    for unch_s = 1:numel(t_unch_states)
        [~,unch_state_start] = min(abs(ts - (t_unch_states(unch_s)-700)));
        [~,unch_state_end] = min(abs(ts - (t_unch_states(unch_s)+700)));
        
        t_unch_FRs(unch_s,:,:) = units(t,unch_state_start:unch_state_end,:);
        t_unch_pp(unch_s,:) = unch_pp(t, unch_state_start:unch_state_end);

    end
    
    all_chState_FRs = [all_chState_FRs; t_ch_FRs];
    all_ch_pp = [all_ch_pp ; t_ch_pp];
    all_full_ch_pp = [all_full_ch_pp ; t_full_ch_pps];
    all_chState_chVal = [all_chState_chVal; ones(ch_s,1)*t_ch_val(t)];
    all_chState_unchVal = [all_chState_unchVal; ones(ch_s,1)*t_unch_val(t)];
    
    all_chState_chSide = [all_chState_chSide; ones(ch_s,1)*t_ch_side];
    all_chState_unchSide = [all_chState_unchSide; ones(ch_s,1)*t_unch_side];
    
    
    all_unchState_FRs = [all_unchState_FRs; t_unch_FRs];
    all_unch_pp = [all_unch_pp ; t_unch_pp];

    all_unchState_unchVal = [all_unchState_unchVal; ones(unch_s,1)*t_unch_val(t)];
    all_unchState_chVal = [all_unchState_chVal; ones(unch_s,1)*t_ch_val(t)];
    all_unchState_chSide = [all_unchState_chSide; ones(unch_s,1)*t_ch_side];
    all_unchState_unchSide = [all_unchState_unchSide; ones(unch_s,1)*t_unch_side];

end % of looping over trials

xt = [1:numel(all_ch_pp(1,:))]*mean(diff(ts)) - 700;
xt_full = [1:numel(all_full_ch_pp(1,:))]*mean(diff(decoder_ts)) - 700;
% 
% figure; 
% hold on
% plot(xt_full,nanmean(all_full_ch_pp));
% plot(xt,nanmean(all_ch_pp));



% now loop over the individual cells and check whether units are differentiating states
        left_Chstate_ix = all_chState_chSide ==-1; 
        right_Chstate_ix = all_chState_chSide ==1; 
        left_unChstate_ix = all_unchState_unchSide ==-1; 
        right_unChstate_ix = all_unchState_unchSide ==1; 
        
%         dir_factor =[ones(sum(left_Chstate_ix),1) ; ones(sum(right_Chstate_ix),1)*-1 ; ones(sum(left_unChstate_ix),1) ; ones(sum(right_unChstate_ix),1)*-1 ] ;
%         Ch_factor =[ones(sum(left_Chstate_ix),1) ; ones(sum(right_Chstate_ix),1) ; ones(sum(left_unChstate_ix),1)*-1 ; ones(sum(right_unChstate_ix),1)*-1 ] ;
        
        dir_factor =[ones(sum(left_Chstate_ix),1)*-1 ; ones(sum(right_Chstate_ix),1)];

for u = 1:numel(all_chState_FRs(1,1,:))
    
    if betas(u) < 0
        
        % pref side is ultimately chosen
        % OFC represents that side
        hv_pref_fr = all_chState_FRs(left_Chstate_ix,:,u);
        hv_pref_pp = all_ch_pp(left_Chstate_ix,:);
        
        % OFC represents the other side
        lv_pref_fr = all_unchState_FRs(right_unChstate_ix,:,u);
        lv_pref_pp = all_unch_pp(right_unChstate_ix,:);

        
        %nonpref side is ultimately chosen
        % OFC represents the nonpref side
        hv_nonpref_fr = all_chState_FRs(right_Chstate_ix,:,u);
        hv_nonpref_pp = all_ch_pp(right_Chstate_ix,:);

        % OFC represents the preferred side
        lv_nonpref_fr = all_unchState_FRs(left_unChstate_ix,:,u);
        lv_nonpref_pp = all_unch_pp(left_unChstate_ix,:);

        
    else
        
        % pref side is ultimately chosen
        % OFC represents that side
        hv_pref_fr = all_chState_FRs(right_Chstate_ix,:,u);
        hv_pref_pp = all_ch_pp(right_Chstate_ix,:);

        
        % OFC represents the other side
        lv_pref_fr = all_unchState_FRs(left_unChstate_ix,:,u);
        lv_pref_pp = all_unch_pp(left_unChstate_ix,:);

        
        %nonpref side is ultimately chosen
        % OFC represents the nonpref side
        hv_nonpref_fr = all_chState_FRs(left_Chstate_ix,:,u);
        hv_nonpref_pp = all_ch_pp(left_Chstate_ix,:);

        % OFC represents the preferred side
        lv_nonpref_fr = all_unchState_FRs(right_unChstate_ix,:,u);
        lv_nonpref_pp = all_unch_pp(right_unChstate_ix,:);

    end
    
    
    for time = 1: numel(hv_pref_fr(1,:,1))
        
       
%         pref_pp_factor = [hv_pref_pp(:,time) ; lv_pref_pp(:,time)];
        pref_pp_factor = [ones(size(hv_pref_pp(:,time))) ; ones(size(lv_pref_pp(:,time)))*-1];
        pref_uFR = [hv_pref_fr(:,time) ; lv_pref_fr(:,time)];
        
%         nonpref_pp_factor = [hv_nonpref_pp(:,time) ; lv_nonpref_pp(:,time)];
        nonpref_pp_factor = [ones(size(hv_nonpref_pp(:,time))) ; ones(size(lv_nonpref_pp(:,time)))*-1];
        nonpref_uFR = [hv_nonpref_fr(:,time) ; lv_nonpref_fr(:,time)];
        
        
        
        [pref_b,~,~,~,pref_stats] = regress(pref_uFR,[ones(size(pref_pp_factor)), pref_pp_factor]); 
        
        pref_pp_pval(u,time) = pref_stats(3);
        pref_beta(u,time) = pref_b(2);
        
        [nonpref_b,~,~,~,nonpref_stats] = regress(nonpref_uFR,[ones(size(nonpref_pp_factor)), nonpref_pp_factor]); 
        
        nonpref_pp_pval(u,time) = nonpref_stats(3);
        nonpref_beta(u,time) = nonpref_b(2);




    end % of looping over times
    
    [hv_pref_FRs(u,:), hv_pref_sem(u,:)] = GetMeanCI(hv_pref_fr,'sem');
    [lv_pref_FRs(u,:), lv_pref_sem(u,:)] = GetMeanCI(lv_pref_fr,'sem');
    [hv_nonpref_FRs(u,:), hv_nonpref_sem(u,:)] = GetMeanCI(hv_nonpref_fr,'sem');
    [lv_nonpref_FRs(u,:), lv_nonpref_sem(u,:)] = GetMeanCI(lv_nonpref_fr,'sem');


    
    
end % of looping over units

% figure; 
% hold on
% plot(xt,nanmean(pref_pp_pval<.01))
% plot(xt,nanmean(nonpref_pp_pval<.01))


end % of function