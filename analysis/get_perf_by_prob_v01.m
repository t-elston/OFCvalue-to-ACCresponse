function [prob_tbl] = get_perf_by_prob_v01(bhv)


% loop over individual files
file_ids = unique(bhv.file_num);

% initialize a counter and a table to save the data in
ctr = 0;
prob_tbl =table;


for f = 1:numel(file_ids)
    
    f_bhv = bhv(bhv.file_num==f,:);
    monkey_id = f_bhv.monkey(1);
    
    f_bhv = f_bhv(f_bhv.trialtype==2,:);
    % now find how well the animals picked when each probability was present
    prob_ids = unique(f_bhv.left_prob);
    prob_ids(prob_ids==0)=[];
    
    prob_acc=[];
    
    for p = 1:numel(prob_ids)
        
        prob_ix = (f_bhv.left_prob == prob_ids(p)) | (f_bhv.right_prob == prob_ids(p));
        
        prob_acc(p) = nanmean(f_bhv.picked_best(prob_ix));        
        
    end % of looping over probabilities
    
        ctr = ctr+1;
        prob_tbl.monkey(ctr) = monkey_id;
        prob_tbl.file_num(ctr) = f;
        prob_tbl.prob_acc(ctr,:) = prob_acc;
        prob_tbl.prob_ids(ctr,:) = prob_ids';

    end % of looping over files

end % of function

