function [distortion_tbl] = get_prob_weighting_by_mag_v01(bhv)


% loop over individual files
file_ids = unique(bhv.file_num);

% initialize a counter and a table to save the data in
ctr = 0;
distortion_tbl =table;

for f = 1:numel(file_ids)
    
    f_bhv = bhv(bhv.file_num==f,:);
    
    mag_ids = unique(f_bhv.left_amnt);
    
    mag_ids(mag_ids==0) = [];
    
    % now loop over possible reward amounts
    for m = 1:numel(mag_ids)
        
        % find trials where the magnitudes were the same between the left and right options
        same_mag_trial_ix = (f_bhv.left_amnt == mag_ids(m)) & (f_bhv.right_amnt == mag_ids(m));
        
        same_mag_trials = f_bhv(same_mag_trial_ix,:);
        
        prob_ids = unique(same_mag_trials.right_prob);
        
        pChooseP=[];
        pChooseP_isbest=[];
        
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
            pChooseP(p) = nanmean(prob_trials.lever == prob_side);
            pChooseP_isbest(p) = nanmean(prob_trials.best_side == prob_side)

        end % of looping over probabilities
        
            %*** SAVE THE RESULTS
            ctr = ctr+1; % increment the counter for saving
            distortion_tbl.monkey(ctr) = prob_trials.monkey(1);
            distortion_tbl.session(ctr) = f;
            distortion_tbl.amnt(ctr)   = mag_ids(m);
            distortion_tbl.pChooseP(ctr,:) = pChooseP;
            distortion_tbl.pChooseP_isbest(ctr,:) = pChooseP_isbest;
            distortion_tbl.pEfficiency(ctr,:) =  pChooseP- pChooseP_isbest;
            distortion_tbl.prob_ids(ctr,:) = prob_ids';
            
        
    end % of looping over magnitudes

end % of looping over files








end % of function