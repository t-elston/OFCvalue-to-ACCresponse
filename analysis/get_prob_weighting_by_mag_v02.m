function [distortion_tbl] = get_prob_weighting_by_mag_v02(bhv)


% loop over individual animals
monkey_ids = unique(bhv.monkey);

% initialize a counter and a table to save the data in
distortion_tbl =table;
ctr=0;
for a = 1:numel(monkey_ids)
    
    % get data for this animal
    a_bhv = bhv(bhv.monkey==monkey_ids(a),:);
    
    mag_ids = unique(a_bhv.left_amnt);
    
    mag_ids(mag_ids==0) = [];
    
    % now find all trials where the magnitudes were the same
    same_mag_trial_ix = a_bhv.left_amnt == a_bhv.right_amnt;
    same_mag_trials = a_bhv(same_mag_trial_ix,:);
        
        prob_ids = unique(same_mag_trials.right_prob);
        
        pChooseP_mean=[];
        pChooseP_isbest_mean=[];
        pChooseP_sem=[];
        pBias_mean=[];
        pChooseP_pval=[];
        
        % now loop over the probabilities
        for p = 1:numel(prob_ids)
            
            % find the trials where this probability was present
            prob_ix = (same_mag_trials.left_prob == prob_ids(p)) | (same_mag_trials.right_prob == prob_ids(p));
            
            prob_trials = same_mag_trials(prob_ix,:);
            
            % find out which side of the screen this option was on
            prob_side = NaN(size(prob_trials.monkey));
            prob_side(prob_trials.left_prob == prob_ids(p)) = -1;
            prob_side(prob_trials.right_prob == prob_ids(p)) = 1;

            % what's the likelihood that he picked this probability?
            pChooseP_mean(p) = nanmean(prob_trials.lever == prob_side);
            pChooseP_sem(p) = nanstd(prob_trials.lever == prob_side) / sqrt(numel(prob_side)); 
            pChooseP_isbest_mean(p) = nanmean(prob_trials.best_side == prob_side);
            pChooseP_isbest_sem(p) = nanstd(prob_trials.best_side == prob_side) / sqrt(numel(prob_side)); 
            pBias_mean(p) = pChooseP_mean(p) - pChooseP_isbest_mean(p);
            pBias_sem(p) = nanstd(prob_trials.lever == prob_side - prob_trials.best_side == prob_side) /sqrt(numel(prob_side));
            
            [~,pChooseP_pval(p)] = ttest(prob_trials.lever == prob_side, prob_trials.best_side == prob_side);

        end % of looping over probabilities
        
            %*** SAVE THE RESULTS
            distortion_tbl.monkey(a) = monkey_ids(a);
            distortion_tbl.pChooseP_mean(a,:) = pChooseP_mean;
            distortion_tbl.pChooseP_sem(a,:) = pChooseP_sem;
            distortion_tbl.pChooseP_isbest(a,:) = pChooseP_isbest_mean;
            distortion_tbl.pChooseP_isbest_sem(a,:) = pChooseP_isbest_sem;
            distortion_tbl.pBias_mean(a,:) =  pBias_mean;
            distortion_tbl.pBias_sem(a,:) =  pBias_sem;
            distortion_tbl.pval(a,:) =  pChooseP_pval;
            distortion_tbl.prob_ids(a,:) = prob_ids';
            

end % of looping over files








end % of function