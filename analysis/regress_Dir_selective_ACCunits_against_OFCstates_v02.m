function [dir_pp_pvals, dir_pvals  pp_pvals, hv_Left_FRs, hv_Right_FRs, lv_Left_FRs, lv_Right_FRs] =...
    regress_Dir_selective_ACCunits_against_OFCstates_v02(units, betas, ts, t_ch_state_times, t_ch_val, t_unch_state_times, t_unch_val, trial_types, chosen_side,...
    ch_pp, unch_pp, full_ch_pp, full_unch_pp, decoder_ts)


pref_side(betas<0)=-1;
pref_side(betas>0)=1;


hv_Left_FRs=[];
hv_Right_FRs=[]; 
lv_Left_FRs=[]; 
lv_Right_FRs=[];

xx=[];
%---------------------------------------------
% identify trials to use and pull those out
trials2keep = ((t_ch_val ==4 & t_unch_val ==3) | (t_ch_val ==3 & t_unch_val ==2));
trials2keep = ~isnan(t_ch_val) & trial_types == 2; 


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

% xt = [1:numel(all_ch_pp(1,:))]*mean(diff(ts)) - 700;
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
all_ch_pp = circshift(all_ch_pp,[0,7]);

for u = 1:numel(all_chState_FRs(1,1,:))
    
    for time = 1: numel(all_chState_FRs(1,:,1))
        
       
%        pp_factor = [all_ch_pp(left_Chstate_ix,time) ; all_ch_pp(right_Chstate_ix,time) ; all_unch_pp(left_unChstate_ix,time) ; all_unch_pp(right_unChstate_ix,time)];
% 
%         % make factors for an anova model
%         uFR = [all_chState_FRs(left_Chstate_ix,time,u) ; all_chState_FRs(right_Chstate_ix,time,u) ; all_unchState_FRs(left_unChstate_ix,time,u) ; all_unchState_FRs(right_unChstate_ix,time,u)];
        
        pp_factor = [all_ch_pp(left_Chstate_ix,time) ; all_ch_pp(right_Chstate_ix,time)];

        % make factors for an anova model
        uFR = [all_chState_FRs(left_Chstate_ix,time,u) ; all_chState_FRs(right_Chstate_ix,time,u)];
        
        reg_tbl =table;
        reg_tbl.fr = uFR; 
        reg_tbl.dir = dir_factor*pref_side(u);
        reg_tbl.pp = pp_factor; 
%         reg_tbl.ch = Ch_factor;
%         
%         mdl = fitglm(reg_tbl,'fr ~ dir*pp*ch');
% 
%         dir_pvals(u, time) = mdl.Coefficients.pValue(2);
%         pp_pvals(u, time) = mdl.Coefficients.pValue(3);
%         ch_pvals(u, time) = mdl.Coefficients.pValue(4);
%         dir_pp_pvals(u, time) = mdl.Coefficients.pValue(5);
%         dir_ch_pvals(u, time) = mdl.Coefficients.pValue(6);
%         pp_ch_pvals(u, time) = mdl.Coefficients.pValue(7);
%         dir_pp_ch_pvals(u, time) = mdl.Coefficients.pValue(8);

        mdl = fitglm(reg_tbl,'fr ~ pp*dir');

        dir_pvals(u, time) = mdl.Coefficients.pValue(2);
        pp_pvals(u, time) = mdl.Coefficients.pValue(3);
        dir_pp_pvals(u, time) = mdl.Coefficients.pValue(4);

    end % of looping over times
    
    hv_Left_FRs(u,:) = mean(all_chState_FRs(left_Chstate_ix,:,u),1);
    hv_Right_FRs(u,:) = mean(all_chState_FRs(right_Chstate_ix,:,u),1);
    
    lv_Left_FRs(u,:) = mean(all_unchState_FRs(left_unChstate_ix,:,u),1);
    lv_Right_FRs(u,:) = mean(all_unchState_FRs(right_unChstate_ix,:,u),1);

    
    
end % of looping over units


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

xt = [1:numel(dir_pvals(1,:))]*mean(diff(ts)) - 700;

% 
% figure; 
% hold on
% plot(xt, mean(dir_pvals<.01), 'LineWidth',2);
% plot(xt, mean(pp_pvals<.01), 'LineWidth',2);
% plot(xt, mean(ch_pvals<.01), 'LineWidth',2);
% plot(xt, mean(dir_pp_pvals<.01), 'LineWidth',2);
% plot(xt, mean(dir_ch_pvals<.01), 'LineWidth',2);
% plot(xt, mean(pp_ch_pvals<.01), 'LineWidth',2);
% plot(xt, mean(dir_pp_ch_pvals<.01), 'LineWidth',2);
% 
% legend({'dir','pp','ch','dir x pp','dir x ch','pp x ch','dir x pp x ch'});

% figure; 
% hold on
% plot(xt, mean(dir_pvals<.01), 'LineWidth',2);
% plot(xt, mean(pp_pvals<.01), 'LineWidth',2);
% plot(xt, mean(dir_pp_pvals<.01), 'LineWidth',2);




end % of function