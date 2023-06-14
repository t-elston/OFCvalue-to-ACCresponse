% aaa_Relate_valencodingACCunits_to_OFC_val_decoder_v01


CT = cbrewer('qual','Set1',9);
CT2 = cbrewer('qual','Dark2',8);
% how do OFC single-unit value representations map onto population value representations?
decoder_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\';
rec_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\recdata\';
alignment = 'pics';

brain_area = 'OFC';

% load the decoder output
% load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\OFC_valdec_output.mat');
load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\ACC_valdec_output.mat');

t_mids = valdata_pics.t_mids;

f_names = unique(valdata_pics.session);

n_trials = numel(valdata_pics.ch_state(:,1));

chosen_side = bhvdata.lever+1;
chosen_side(chosen_side==0) = 1;
unchosen_side = chosen_side;
unchosen_side(chosen_side==1) = 2;
unchosen_side(chosen_side==2) = 1;

% initialize and array to save the indices of the state times
t_ch_state_times = NaN(n_trials, 15);
t_unch_state_times = NaN(n_trials, 15);

% set a threshold of when to look at the states (i.e states that occur after the pics come on)
[~,choice_on_time_ix] = min(abs(t_mids - 0));
[~,max_time_ix] = min(abs(t_mids - 500));

% go through each trial and find the start indices of each state and their associated values
for t = 1:n_trials
        
    % what were the values for this trial?
    t_ch_val(t,1) = bhvdata.valbin_expval(t,chosen_side(t));
    t_unch_val(t,1) = bhvdata.valbin_expval(t,unchosen_side(t));
    
    % get this trial's chosen state
    ch_state = valdata_pics.ch_state(t,:);
    
    % get this trial's unchosen state
    unch_state = valdata_pics.unch_state(t,:);
    
    % find when the chosen state was 'up'
    [~, ch_state_starts] = findpeaks(ch_state);
    
    % find when the unchosen state was 'up'
    [~, unch_state_starts] = findpeaks(unch_state);
    
    % get the indices of the state starts that happened after the pics turned on
    ch_state_ix = ch_state_starts((ch_state_starts>choice_on_time_ix) & (ch_state_starts<max_time_ix));
    unch_state_ix = unch_state_starts((unch_state_starts>choice_on_time_ix) & (unch_state_starts<max_time_ix));
    
    % now get the times in ms of when the states started (will use these to align the units later)
    ch_state_start_times = t_mids(ch_state_ix);
    unch_state_start_times = t_mids(unch_state_ix);
    
    % now save the chosen and unchosen state times
    t_ch_state_times(t,1:numel(ch_state_start_times)) = ch_state_start_times; 
    t_unch_state_times(t,1:numel(unch_state_start_times)) = unch_state_start_times; 

end % of cycling over trials
%---------------------------------------

% now go file-by-file and look at its cells
% initialize some variables 
match_betas=[];
match_pvals=[];
non_match_betas=[];
non_match_pvals=[];
LR_match_betas=[];
LR_match_pvals=[];
chState_FRs=[];
chState_chVal=[];
chState_unchVal=[];
monkey_id=[];
for f = 1:numel(f_names)
    
    % current file name
    current_file = f_names(f);
    
    fprintf('\n')
    disp(['***' current_file{1}]);
    animal_id =contains(current_file,'George');
    
    % load neurons associated with this file
    [units, u_hemi_ix, ts] = load_neurons(rec_dir, current_file, alignment, brain_area, NaN);
    
    
    if ~isempty(units)
        % isolate the parts of the decoder output that corresponds to this file
        file_t_ix = contains(bhvdata.session,current_file);
        
        % now loop over cells and trials for this file
        [f_match_betas, f_match_pvals, f_non_match_betas, f_non_match_pvals, f_LR_match_betas, f_LR_match_pvals, f_all_chState_FRs, f_all_chState_chVal, f_all_chState_unchVal] =...
            regress_corresponding_and_noncorresponding_PPs_with_FRs_v01(units, ts, t_ch_state_times(file_t_ix,:),t_ch_val(file_t_ix), t_unch_state_times(file_t_ix,:), t_unch_val(file_t_ix), bhvdata.trialtype(file_t_ix), bhvdata.lever(file_t_ix));
        
        % now save these results into some bigger arrays
        match_betas = [match_betas; f_match_betas];
        match_pvals = [match_pvals; f_match_pvals];
        non_match_betas = [non_match_betas; f_non_match_betas];
        non_match_pvals = [non_match_pvals; f_non_match_pvals];
        LR_match_betas = [LR_match_betas; f_LR_match_betas];
        LR_match_pvals = [LR_match_pvals; f_LR_match_pvals];
        chState_FRs = cat(3, f_all_chState_FRs);
        chState_chVal = [chState_chVal; f_all_chState_chVal];
        chState_unchVal = [chState_unchVal; f_all_chState_unchVal];
        
        monkey_id = [monkey_id ; ones(size(u_hemi_ix))*animal_id];
        
    end
    
    
end % of cycling over files

[ch_rvals, unch_rvals, corr_animal_id] = get_sessionwiseACCOFC_valstate_correlations;

[c_ch_corr_mean, c_ch_corr_sem] = GetMeanCI(ch_rvals(corr_animal_id==0,:),'sem');
[c_unch_corr_mean, c_unch_corr_sem] = GetMeanCI(unch_rvals(corr_animal_id==0,:),'sem');
[g_ch_corr_mean, g_ch_corr_sem] = GetMeanCI(ch_rvals(corr_animal_id==1,:),'sem');
[g_unch_corr_mean, g_unch_corr_sem] = GetMeanCI(unch_rvals(corr_animal_id==1,:),'sem');

[chap_sig_ch] = get_significant_correlation_times(ch_rvals(corr_animal_id==0,:), 24, t_mids);
[chap_sig_unch] = get_significant_correlation_times(unch_rvals(corr_animal_id==0,:), 24, t_mids);
[george_sig_ch] = get_significant_correlation_times(ch_rvals(corr_animal_id==1,:), 24, t_mids);
[george_sig_unch] = get_significant_correlation_times(unch_rvals(corr_animal_id==1,:), 24, t_mids);






g_ix = monkey_id==1;
c_ix = monkey_id==0;

xt = [1:numel(match_pvals(1,:))]*mean(diff(ts)) - 500;
[~,pre_win_start] = min(abs(xt - -125));
[~,pre_win_end] = min(abs(xt - -25));
[~,post_win_start] = min(abs(xt - 25));
[~,post_win_end] = min(abs(xt - 125));


n_g_units = numel(match_pvals(g_ix,1));
n_c_units = numel(match_pvals(c_ix,1));


Fig = figure;
set(Fig,'renderer','Painters');
subplot(3,2,1);
hold on
shadedErrorBar(t_mids,c_unch_corr_mean, c_unch_corr_sem,'LineProps',{'LineWidth',1,'color',CT(2,:)});
shadedErrorBar(t_mids,c_ch_corr_mean, c_ch_corr_sem,'LineProps',{'LineWidth',1,'color',CT(4,:)});
plot(t_mids(chap_sig_ch), ones(size(chap_sig_ch))*-.02,'s', 'color',CT(4,:),'MarkerFaceColor',CT(4,:));
plot(t_mids(chap_sig_unch), ones(size(chap_sig_unch))*-.03,'s', 'color',CT(2,:),'MarkerFaceColor',CT(2,:));
ylim([-.05,.2]);

plot(xlim,[0 0],'k','LineWidth',1);
xlim([-400,600]);
ylabel('ACC-OFC PP Correlation');
axis square
title('Animal C');
xlabel('Time from Pictures On');
set(gca, 'FontSize',10, 'LineWidth',1);



subplot(3,2,2);
hold on
shadedErrorBar(t_mids,g_unch_corr_mean, g_unch_corr_sem,'LineProps',{'LineWidth',1,'color',CT(2,:)});
shadedErrorBar(t_mids,g_ch_corr_mean, g_ch_corr_sem,'LineProps',{'LineWidth',1,'color',CT(4,:)});
plot(t_mids(george_sig_ch), ones(size(george_sig_ch))*-.02,'s', 'color',CT(4,:),'MarkerFaceColor',CT(4,:));
plot(t_mids(george_sig_unch), ones(size(george_sig_unch))*-.03,'s', 'color',CT(2,:),'MarkerFaceColor',CT(2,:));
plot(xlim,[0 0],'k','LineWidth',1);
xlim([-400,600]);
ylim([-.05,.2]);

axis square
title('Animal G');
set(gca, 'FontSize',10, 'LineWidth',1);




subplot(3,2,3);

for t = 1:numel(xt)
    
    contab(1,1) = sum(match_pvals(c_ix,t)<.01);
    contab(1,2) = n_c_units;
    contab(2,1) = sum(non_match_pvals(c_ix,t)<.01);
    contab(2,2) = n_c_units;
    [pval,x2] = chisquarecont(contab);
    
    if pval<.01
        chap_sigX2(t) = 1;
    else
        chap_sigX2(t) = NaN;
    end
    
end % of cycling over times


hold on
bar(xt,mean(match_pvals(c_ix,:)<.01), 'EdgeColor',CT2(1,:), 'FaceColor',CT2(1,:),'BarWidth', 1);
bar(xt,mean(non_match_pvals(c_ix,:)<.01), 'EdgeColor',CT(9,:), 'FaceColor',CT(9,:),'BarWidth', 1);
% plot(xt,chap_sigX2*.45,'*','color','k')
ylim([0, .5]);
plot([0,0],ylim,'k','LineStyle',':');
xlim([-400, 400]);
xlabel('Time from State Onset (ms)');
ylabel('Prop. Units p<.01');
set(gca, 'FontSize',10, 'LineWidth',1);
axis square

n_c_total = sum(c_ix);
n_c_match_prewin = median(sum((match_pvals(c_ix,pre_win_start:pre_win_end)<.01)));
n_c_match_postwin = median(sum((match_pvals(c_ix,post_win_start:post_win_end)<.01)));
n_c_nonmatch_prewin = median(sum((non_match_pvals(c_ix,pre_win_start:pre_win_end)<.01)));
n_c_nonmatch_postwin = median(sum((non_match_pvals(c_ix,post_win_start:post_win_end)<.01)));




subplot(3,2,5);
hold on
plot(LR_match_betas(c_ix,1),LR_match_betas(c_ix,2),'.','color','k');
xlim([-.5,.5]);
ylim([-.5,.5]);
plot(xlim,[0,0],'k','LineWidth',1);
plot([0,0],ylim,'k','LineWidth',1);
l = lsline;
l.Color = [.5,.5,.5];
l.LineWidth=1;
axis square


% xlabel('\beta left | left state');
% ylabel('\beta right | right state');
set(gca, 'FontSize',10, 'LineWidth',1);
[chap_r, chap_p] = corr(LR_match_betas(c_ix,1),LR_match_betas(c_ix,2));

subplot(3,2,4);

for t = 1:numel(xt)
    
    contab(1,1) = sum(match_pvals(g_ix,t)<.01);
    contab(1,2) = n_g_units;
    contab(2,1) = sum(non_match_pvals(g_ix,t)<.01);
    contab(2,2) = n_g_units;
    [pval,x2] = chisquarecont(contab);
    
    
    if pval<.01
        george_sigX2(t) = 1;
    else
        george_sigX2(t) = NaN;
    end
    
end % of cycling over times


hold on
bar(xt,mean(match_pvals(g_ix,:)<.01), 'EdgeColor',CT2(1,:), 'FaceColor',CT2(1,:),'BarWidth', 1);
bar(xt,mean(non_match_pvals(g_ix,:)<.01), 'EdgeColor',CT(9,:), 'FaceColor',CT(9,:),'BarWidth', 1);
% plot(xt,george_sigX2*140,'*','color','k')
ylim([0, .3]);
plot([0,0],ylim,'k','LineStyle',':');
xlim([-400, 400]);
% xlabel('Time from State Onset (ms)');
% ylabel('Prop. Units p<.01');
set(gca, 'FontSize',10, 'LineWidth',1);
axis square


n_g_total = sum(g_ix);
n_g_match_prewin = median(sum((match_pvals(g_ix,pre_win_start:pre_win_end)<.01)));
n_g_match_postwin = median(sum((match_pvals(g_ix,post_win_start:post_win_end)<.01)));
n_g_nonmatch_prewin = median(sum((non_match_pvals(g_ix,pre_win_start:pre_win_end)<.01)));
n_g_nonmatch_postwin = median(sum((non_match_pvals(g_ix,post_win_start:post_win_end)<.01)));


% legend({'State = Value', 'State ~= Value', 'X^2 p<.01'},'FontSize',10);

subplot(3,2,6);
hold on
plot(LR_match_betas(g_ix,1),LR_match_betas(g_ix,2),'.','color','k');
xlim([-.5,.5]);
ylim([-.5,.5]);
plot(xlim,[0,0],'k','LineWidth',1);
plot([0,0],ylim,'k','LineWidth',1);
l = lsline;
l.Color = [.5,.5,.5];
l.LineWidth=1;
axis square

% xlabel('\beta left | left state');
% ylabel('\beta right | right state');
set(gca, 'FontSize',10, 'LineWidth',1);
[g_r, g_p] = corr(LR_match_betas(g_ix,1),LR_match_betas(g_ix,2));













