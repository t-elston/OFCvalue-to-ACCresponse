function [xx] = make_state_likelihoods_over_time(t_mids, valdata_choice, bhvdata, t_ch_val, t_unch_val,...
    t_ch_state_times, t_unch_state_times, t_ch_state_dur, t_unch_state_dur, animal)
xx=[];

CT = cbrewer('qual','Set1',9);

n_ch_states = nansum(~isnan(t_ch_state_times),2);
n_unch_states = nansum(~isnan(t_unch_state_times),2);

animal_ix = contains(valdata_choice.subject,animal);

trials_with_ch_state = n_ch_states>1;
trials_with_unch_state = n_unch_states>1;
free_choice_ix =~isnan(t_ch_val);

candidate_trials = trials_with_ch_state & trials_with_unch_state & free_choice_ix & t_ch_val>t_unch_val;


% now get value combinations
[unchosen_vals, chosen_vals] = ndgrid([1:3],[2:4]);
combinations = [unchosen_vals(:), chosen_vals(:)];
combinations(diff(combinations,1,2)<=0,:)=[];

% now figure out the panel numbers for a 3x3 grid
panel_nums = [1, 4,5, 7,8,9];
figure; 
pvals=[];
sgtitle(animal)
for c = 1:numel(combinations(:,1))
    
    ch_ix = t_ch_val == combinations(c,2);
    unch_ix = t_unch_val == combinations(c,1);

      % NUM STATES
%     mean_ch(c) = mean(sum(n_ch_states(candidate_trials & ch_ix & animal_ix),2));
%     mean_unch(c) = mean(sum(n_unch_states(candidate_trials & unch_ix & animal_ix),2));
%     
%     ch_trial_sums = sum(n_ch_states(candidate_trials & ch_ix & animal_ix),2);
%     unch_trial_sums = sum(n_unch_states(candidate_trials & unch_ix & animal_ix),2);
%     
%     
%     [~,pvals(c)] = ttest2(ch_trial_sums, unch_trial_sums);

% DURATION OF STATES
    mean_ch(c) = mean(nanmean(t_ch_state_dur(candidate_trials & ch_ix & animal_ix,:),2));
    mean_unch(c) = mean(nanmean(t_unch_state_dur(candidate_trials & unch_ix & animal_ix,:),2));
    
    ch_trial_sums = nanmean(t_ch_state_dur(candidate_trials & ch_ix & animal_ix,:),2);
    unch_trial_sums = nanmean(t_unch_state_dur(candidate_trials & unch_ix & animal_ix,:),2);
    
    
    [~,pvals(c)] = ttest2(ch_trial_sums, unch_trial_sums);


%----------------------
    
    sig_level = 0;
    if pvals(c) < .01
        sig_level = 1;
    end
    if pvals(c) < .001
        sig_level = 2;
    end
    if pvals(c) < .0001
        sig_level = 3;
    end
    
    
    subplot(3,3,panel_nums(c))
    hold on
    bar(1,mean_ch(c),'FaceColor',CT(4,:))
    bar(2,mean_unch(c),'FaceColor',CT(2,:))
    xticks([1,2]);
    xticklabels({combinations(c,2), combinations(c,1)})
%     ylim([1, 1.7]);
    ylabel('mean state duration (ms)');
    xlabel('values');
    
    switch sig_level
        case 1
        plot(1.5, max(ylim),'*','MarkerSize',5,'color','k');
        case 2
        plot([1.3,1.6], [max(ylim), max(ylim)],'*','MarkerSize',5,'color','k');
        case 3
        plot([1.3,1.5,1.7], [max(ylim), max(ylim), max(ylim)],'*','MarkerSize',5,'color','k');
        end
    
    if c ==1
        legend('chosen','unchosen');
    end
    

  
end % of looping over combinations


end % of function