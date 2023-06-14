% aaa_Relate_OFCunits_to_OFC_val_decoder_v01


CT = cbrewer('qual','Set1',9);
CT2 = cbrewer('qual','Dark2',8);
% how do OFC single-unit value representations map onto population value representations?
decoder_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\';
rec_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\recdata\';
alignment = 'pics';

brain_area = 'ACC';

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
u_animal=[];
for f = 1:numel(f_names)
    
    % current file name
    current_file = f_names(f);
    
    fprintf('\n')
    disp(['***' current_file{1}]);
    monkey_id = contains(current_file,'George');
    
    % load neurons associated with this file
    [units, u_hemi_ix, ts] = load_neurons(rec_dir, current_file, alignment, brain_area, NaN);
    [n_trials, n_times, n_units] = size(units);
    
    
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
        
        u_animal = [u_animal ; ones(n_units,1)*monkey_id];
        
        
    end
    
    
end % of cycling over files

cmap_c = 1;

trials2keep = bhvdata.trialtype==2 & ~isnan(t_ch_val) & (t_ch_val ~= t_unch_val);

chap_ch_state = valdata_pics.ch_state(trials2keep & contains(bhvdata.subject,'Chap'),:);
george_ch_state = valdata_pics.ch_state(trials2keep & contains(bhvdata.subject,'George'),:);
[chap_ch_mean, chap_ch_sem] = GetMeanCI(chap_ch_state,'sem');
[george_ch_mean, george_ch_sem] = GetMeanCI(george_ch_state,'sem');

Fig = figure;
set(Fig,'renderer','Painters','Position',[200 200 600 800]);

xt = [1:numel(match_pvals(1,:))]*mean(diff(ts)) - 500;

subplot(3,2,1); 
hold on
shadedErrorBar(t_mids,chap_ch_mean,chap_ch_sem,'LineProps',{'LineWidth',2,'color',CT2(cmap_c,:)});
ylim([0, .45]);
xlim([-800,1500]);
plot([0,0],ylim,'k','LineStyle',':');
xlabel('Time from Pics On (ms)');
ylabel('p(State) Occurs');
set(gca, 'FontSize',10, 'LineWidth',1);
axis square

subplot(3,2,2); 
hold on
shadedErrorBar(t_mids,george_ch_mean, george_ch_sem,'LineProps',{'LineWidth',2,'color',CT2(cmap_c,:)});
xlim([-800,1500]);
ylim([0, .45]);
plot([0,0],ylim,'k','LineStyle',':');
% xlabel('Time from Pics On (ms)');
% ylabel('p(State) Occurs');
set(gca, 'FontSize',10, 'LineWidth',1);
axis square


subplot(3,2,3);
n_units = numel(match_pvals(:,1));

hold on
bar(xt,mean(match_pvals(u_animal==0,:)<.01), 'EdgeColor',CT2(cmap_c,:), 'FaceColor',CT2(cmap_c,:),'BarWidth', 1);
bar(xt,mean(non_match_pvals(u_animal==0,:)<.01), 'EdgeColor',CT(9,:), 'FaceColor',CT(9,:),'BarWidth', 1);
ylim([.05, .6]);
plot([0,0],ylim,'k','LineStyle',':');
xlim([-400, 400]);
xlabel('Time from State Onset (ms)');
ylabel('Prop. Units p<.01');
set(gca, 'FontSize',10, 'LineWidth',1);
axis square

subplot(3,2,4);
hold on
bar(xt,mean(match_pvals(u_animal==1,:)<.01), 'EdgeColor',CT2(cmap_c,:), 'FaceColor',CT2(cmap_c,:),'BarWidth', 1);
bar(xt,mean(non_match_pvals(u_animal==1,:)<.01), 'EdgeColor',CT(9,:), 'FaceColor',CT(9,:),'BarWidth', 1);
ylim([.05, .5]);
plot([0,0],ylim,'k','LineStyle',':');
xlim([-400, 400]);
xlabel('Time from State Onset (ms)');
ylabel('Prop. Units p<.01');
set(gca, 'FontSize',10, 'LineWidth',1);
axis square


subplot(3,2,5); 
hold on
plot(LR_match_betas(u_animal==0,1),LR_match_betas(u_animal==0,2),'.','color',CT2(cmap_c,:));
xlim([-.5,.5]);
ylim([-.5,.5]);
plot(xlim,[0,0],'k','LineWidth',1);
plot([0,0],ylim,'k','LineWidth',1);


l = lsline; 
l.Color = [.5,.5,.5];
l.LineWidth=1;

% xlabel('\beta left | left state');
% ylabel('\beta right | right state');
set(gca, 'FontSize',10, 'LineWidth',1);
[chap_r, chap_p] = corr(LR_match_betas(u_animal==0,1),LR_match_betas(u_animal==0,2));
axis square


subplot(3,2,6); 
hold on
plot(LR_match_betas(u_animal==1,1),LR_match_betas(u_animal==1,2),'.','color',CT2(cmap_c,:));
xlim([-.5,.5]);
ylim([-.5,.5]);
plot(xlim,[0,0],'k','LineWidth',1);
plot([0,0],ylim,'k','LineWidth',1);


l = lsline; 
l.Color = [.5,.5,.5];
l.LineWidth=1;

% xlabel('\beta left | left state');
% ylabel('\beta right | right state');
set(gca, 'FontSize',10, 'LineWidth',1);
[george_r, george_p] = corr(LR_match_betas(u_animal==1,1),LR_match_betas(u_animal==1,2));
axis square








c_before_median = median(mean((match_pvals(u_animal==0,14:18)<.01)));
g_before_median = median(mean((match_pvals(u_animal==1,14:18)<.01)));

c_after_median = median(mean((match_pvals(u_animal==0,21:24)<.01)));
g_after_median = median(mean((match_pvals(u_animal==1,21:24)<.01)));

c_after_alt_median = median(mean((non_match_pvals(u_animal==0,21:24)<.01)));
g_after_alt_median = median(mean((non_match_pvals(u_animal==1,21:24)<.01)));




