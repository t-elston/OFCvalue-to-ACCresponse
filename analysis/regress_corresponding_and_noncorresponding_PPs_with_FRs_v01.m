function [match_betas, match_pvals, non_match_betas, non_match_pvals, LR_match_betas, LR_match_pvals, all_chState_FRs, all_chState_chVal, all_chState_unchVal] =...
    regress_corresponding_and_noncorresponding_PPs_with_FRs_v01(units, ts, t_ch_state_times, t_ch_val, t_unch_state_times, t_unch_val, trial_types, chosen_side)

xx=[];
%---------------------------------------------
% only include free choice trials and those in which he made a choice
trials2keep = trial_types==2 & ~isnan(t_ch_val) & (t_ch_val > t_unch_val);
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
all_ch_side=[];

all_unchState_FRs=[];
all_unchState_unchVal=[];
all_unchState_chVal=[];
all_unch_side=[];
% loop through each trial and accumulate the firing rates

for t = 1:n_trials
    
    t_ch_FRs=[];
    t_unch_FRs=[];
    % get the times of the chosen states on this trial 
    t_ch_side = chosen_side(t);
    t_ch_states = t_ch_state_times(t,:);
    t_ch_states(isnan(t_ch_states))=[];
    
    t_unch_side =  chosen_side(t)*-1;
    t_unch_states = t_unch_state_times(t,:);
    t_unch_states(isnan(t_unch_states))=[];
    
    % now loop over each state and find window aligned on state onset
    for ch_s = 1:numel(t_ch_states)
        [~,ch_state_start] = min(abs(ts - (t_ch_states(ch_s)-500)));
        [~,ch_state_end] = min(abs(ts - (t_ch_states(ch_s)+500)));
        
        t_ch_FRs(ch_s,:,:) = units(t,ch_state_start:ch_state_end,:);
    end
    
    % now loop over each state and find window aligned on state onset
    for unch_s = 1:numel(t_unch_states)
        [~,unch_state_start] = min(abs(ts - (t_unch_states(unch_s)-500)));
        [~,unch_state_end] = min(abs(ts - (t_unch_states(unch_s)+500)));
        
        t_unch_FRs(unch_s,:,:) = units(t,unch_state_start:unch_state_end,:);
    end
    
    all_chState_FRs = [all_chState_FRs; t_ch_FRs];
    all_chState_chVal = [all_chState_chVal; ones(ch_s,1)*t_ch_val(t)];
    all_chState_unchVal = [all_chState_unchVal; ones(ch_s,1)*t_unch_val(t)];
    all_ch_side = [all_ch_side; ones(ch_s,1)*t_ch_side];
    
    
    all_unchState_FRs = [all_unchState_FRs; t_unch_FRs];
    all_unchState_unchVal = [all_unchState_unchVal; ones(unch_s,1)*t_unch_val(t)];
    all_unchState_chVal = [all_unchState_chVal; ones(unch_s,1)*t_ch_val(t)];
    all_unch_side = [all_unch_side ; ones(unch_s,1)*t_unch_side];

end % of looping over trials


% now loop over the individual cells and fit regressions

% let's first regress over the chosen vals
% go unit by unit, time by time

% create the match regressor
chosen_match_regressor = [ones(size(all_chState_chVal)), all_chState_chVal];
chosen_dir_regressor = [ones(size(all_ch_side)), all_ch_side];

unchosen_match_regressor = [ones(size(all_unchState_unchVal)), all_unchState_unchVal];
unchosen_dir_regressor = [ones(size(all_unch_side)), all_unch_side];


all_match_regressors = [chosen_match_regressor ; unchosen_match_regressor];

all_match_dir_regressors = [chosen_dir_regressor ; unchosen_dir_regressor];

all_nonmatch_dir_regressors = [unchosen_dir_regressor ; chosen_dir_regressor];


% create the non-match regressor
chosen_nonmatch_regressor = [ones(size(all_chState_unchVal)), all_chState_unchVal];
unchosen_nonmatch_regressor = [ones(size(all_unchState_chVal)), all_unchState_chVal];

all_nonmatch_regressors = [chosen_nonmatch_regressor; unchosen_nonmatch_regressor];


for u = 1:numel(all_chState_FRs(1,1,:))
    
    for time = 1: numel(all_chState_FRs(1,:,1))
        
      

        % create the match predictors for value
        reg_y = [all_chState_FRs(:,time,u) ; all_unchState_FRs(:,time,u)];
        
        [match_b,~,~,~,match_stats] = regress(reg_y, all_match_regressors);
        [nonmatch_b,~,~,~,nonmatch_stats] = regress(reg_y, all_nonmatch_regressors);
        
        % extract the beta weights
        match_betas(u, time) =  match_b(2);
        non_match_betas(u, time) =  nonmatch_b(2);
        
        % extract the pvalues
        match_pvals(u, time) =  match_stats(3);
        non_match_pvals(u, time) =  nonmatch_stats(3);

        
        % now look at direction
        [dir_match_b,~,~,~,dir_match_stats] = regress(reg_y, all_match_dir_regressors);
        [dir_nonmatch_b,~,~,~,dir_nonmatch_stats] = regress(reg_y, all_nonmatch_dir_regressors);
        
        % extract the beta weights
        dir_match_betas(u, time) =  dir_match_b(2);
        dir_non_match_betas(u, time) =  dir_nonmatch_b(2);
        
        % extract the pvalues
        dir_match_pvals(u, time) =  dir_match_stats(3);
        dir_non_match_pvals(u, time) =  dir_nonmatch_stats(3);

        
    end % of looping over times
    
    left_match_FRs = mean(all_chState_FRs(all_ch_side==-1,16:17,u),2);
    right_match_FRs = mean(all_chState_FRs(all_ch_side==1,16:17,u),2);
    
    % VALUE
    % now get the matching stats for left and right
    [LEFT_match_b,~,~,~,LEFT_match_stats] = regress(left_match_FRs,...
        [ones(size(left_match_FRs)), all_chState_chVal(all_ch_side==-1)]);
    
    [RIGHT_match_b,~,~,~,RIGHT_match_stats] = regress(right_match_FRs,...
        [ones(size(right_match_FRs)), all_chState_chVal(all_ch_side==1)]);
    
    
    LR_match_betas(u, 1) =  LEFT_match_b(2);
    LR_match_betas(u, 2) =  RIGHT_match_b(2);
    
    LR_match_pvals(u, 1) =  LEFT_match_stats(3);
    LR_match_pvals(u, 2) =  RIGHT_match_stats(3);
    
%     % DIRECTION
%     % now get the matching stats for left and right
%     [LEFT_dir_match_b,~,~,~,LEFT_dir_match_stats] = regress(left_match_FRs,...
%         [ones(size(left_match_FRs)), all_match_dir_regressors(all_ch_side==-1)]);
%     
%     [RIGHT_dir_match_b,~,~,~,RIGHT_dir_match_stats] = regress(right_match_FRs,...
%         [ones(size(right_match_FRs)), all_match_dir_regressors(all_ch_side==1)]);
%     
%     
%     LR_match_betas(u, 3) =  LEFT_dir_match_b(2);
%     LR_match_betas(u, 4) =  RIGHT_dir_match_b(2);
%     
%     LR_match_pvals(u, 3) =  LEFT_dir_match_stats(3);
%     LR_match_pvals(u, 4) =  RIGHT_dir_match_stats(3);
    

    
end % of looping over units







end % of function