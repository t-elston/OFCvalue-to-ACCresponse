function [ChDir_pvals, Dir_pvals  Ch_pvals, hv_Left_FRs, hv_Right_FRs, lv_Left_FRs, lv_Right_FRs] =...
    regress_Dir_selective_ACCunits_against_OFCstates_v01(units, betas, ts, t_ch_state_times, t_ch_val, t_unch_state_times, t_unch_val, trial_types, chosen_side)


hv_Left_FRs=[];
hv_Right_FRs=[]; 
lv_Left_FRs=[]; 
lv_Right_FRs=[];

xx=[];
%---------------------------------------------
% only include free choice trials and those in which he made a choice
% trials2keep = trial_types==2  & ((t_ch_val ==4 & t_unch_val ==3) | (t_ch_val ==3 & t_unch_val ==2));
% trials2keep = trial_types==2  & ((t_ch_val ==2 & t_unch_val ==3) | (t_ch_val ==3 & t_unch_val ==2));
% trials2keep = trial_types==2  & ((t_ch_val ==3 & t_unch_val ==2) | (t_ch_val ==4 & t_unch_val ==3));
trials2keep = trial_types==2  & (t_ch_val ==3  |  t_unch_val ==3);



units = units(trials2keep,:,:);
t_ch_state_times = t_ch_state_times(trials2keep,:);
t_ch_val = t_ch_val(trials2keep);
t_unch_state_times = t_unch_state_times(trials2keep,:);
t_unch_val = t_unch_val(trials2keep);
chosen_side = chosen_side(trials2keep);
%---------------------------------------------

n_trials = numel(t_ch_val);
all_chState_FRs=[];
all_chState_chVal=[];
all_chState_unchVal=[];
all_chState_chSide=[];
all_chState_unchSide=[];
all_unchState_chSide=[];
all_unchState_unchSide=[];

all_unchState_FRs=[];
all_unchState_unchVal=[];
all_unchState_chVal=[];
% loop through each trial and accumulate the firing rates


for t = 1:n_trials
    
    t_ch_FRs=[];
    t_unch_FRs=[];
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
    end
    
    % now loop over each state and find window aligned on state onset
    for unch_s = 1:numel(t_unch_states)
        [~,unch_state_start] = min(abs(ts - (t_unch_states(unch_s)-700)));
        [~,unch_state_end] = min(abs(ts - (t_unch_states(unch_s)+700)));
        
        t_unch_FRs(unch_s,:,:) = units(t,unch_state_start:unch_state_end,:);
    end
    
    all_chState_FRs = [all_chState_FRs; t_ch_FRs];
    all_chState_chVal = [all_chState_chVal; ones(ch_s,1)*t_ch_val(t)];
    all_chState_unchVal = [all_chState_unchVal; ones(ch_s,1)*t_unch_val(t)];
    
    all_chState_chSide = [all_chState_chSide; ones(ch_s,1)*t_ch_side];
    all_chState_unchSide = [all_chState_unchSide; ones(ch_s,1)*t_unch_side];
    
    
    all_unchState_FRs = [all_unchState_FRs; t_unch_FRs];
    all_unchState_unchVal = [all_unchState_unchVal; ones(unch_s,1)*t_unch_val(t)];
    all_unchState_chVal = [all_unchState_chVal; ones(unch_s,1)*t_ch_val(t)];
    all_unchState_chSide = [all_unchState_chSide; ones(unch_s,1)*t_ch_side];
    all_unchState_unchSide = [all_unchState_unchSide; ones(unch_s,1)*t_unch_side];

end % of looping over trials


% now loop over the individual cells and check whether units are differentiating states
        left_Chstate_ix = all_chState_chSide ==-1; 
        right_Chstate_ix = all_chState_chSide ==1; 
        left_unChstate_ix = all_unchState_unchSide ==-1; 
        right_unChstate_ix = all_unchState_unchSide ==1; 
        
        dir_factor =[ones(sum(left_Chstate_ix),1) ; zeros(sum(right_Chstate_ix),1) ; ones(sum(left_unChstate_ix),1) ; zeros(sum(right_unChstate_ix),1) ] ;
        Ch_factor =[ones(sum(left_Chstate_ix),1) ; ones(sum(right_Chstate_ix),1) ; zeros(sum(left_unChstate_ix),1) ; zeros(sum(right_unChstate_ix),1) ] ;

for u = 1:numel(all_chState_FRs(1,1,:))
    
    for time = 1: numel(all_chState_FRs(1,:,1))
        
        
        % make factors for an anova model
        uFR = [all_chState_FRs(left_Chstate_ix,time,u) ; all_chState_FRs(right_Chstate_ix,time,u) ; all_unchState_FRs(left_unChstate_ix,time,u) ; all_unchState_FRs(right_unChstate_ix,time,u)];
        

        
        [p,tbl] = anovan(uFR,{Ch_factor dir_factor},'model','interaction','varnames',{'chosen','dir'},'Display','off');

        ChDir_pvals(u, time) = p(3);
        Dir_pvals(u, time) = p(2);
        Ch_pvals(u, time) = p(1);

            
    end % of looping over times
    
    hv_Left_FRs(u,:) = mean(all_chState_FRs(left_Chstate_ix,:,u),1);
    hv_Right_FRs(u,:) = mean(all_chState_FRs(right_Chstate_ix,:,u),1);
    
    lv_Left_FRs(u,:) = mean(all_unchState_FRs(left_unChstate_ix,:,u),1);
    lv_Right_FRs(u,:) = mean(all_unchState_FRs(right_unChstate_ix,:,u),1);

    
    
end % of looping over units

% adjust the timecourses to account for smoothing window bias
hv_Left_FRs = circshift(hv_Left_FRs,[0,5]);
hv_Right_FRs = circshift(hv_Right_FRs,[0,5]);
lv_Left_FRs = circshift(lv_Left_FRs,[0,5]);
lv_Right_FRs = circshift(lv_Right_FRs,[0,5]);



% 
% figure; 
% subplot(1,2,1);
% hold on
% plot(mean(hv_Right_FRs(betas>0,:)), 'LineWidth',2)
% plot(mean(lv_Right_FRs(betas>0,:)), 'LineWidth',2)
% plot(mean(hv_Left_FRs(betas>0,:)), 'LineWidth',2)
% plot(mean(lv_Left_FRs(betas>0,:)), 'LineWidth',2)
% title('Right Coding cells');
% legend('HV Pref','LV Pref','HV Non-Pref','LV Non-Pref')
% 
% 
% subplot(1,2,2);
% hold on
% plot(mean(hv_Left_FRs(betas<0,:)), 'LineWidth',2)
% plot(mean(lv_Left_FRs(betas<0,:)), 'LineWidth',2)
% plot(mean(hv_Right_FRs(betas<0,:)), 'LineWidth',2)
% plot(mean(lv_Right_FRs(betas<0,:)), 'LineWidth',2)
% title('Left Coding cells');
% 
% legend('HV Pref','LV Pref','HV Non-Pref','LV Non-Pref')





end % of function