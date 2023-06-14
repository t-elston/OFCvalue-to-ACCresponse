function [pvals, betas] = fixed_window_rel2choice_GLMs_v02(units, ts, trialinfo)

chosen_side = trialinfo.lever; 

chosen_val = NaN(size(trialinfo.lever)); 

chosen_val(chosen_side==-1) = trialinfo.subjval_expval(chosen_side==-1,1);
chosen_val(chosen_side==1) = trialinfo.subjval_expval(chosen_side==1,2);

% get indices of the fixed window
[~,w_start] = min(abs(ts- -400));
[~,w_end] = min(abs(ts- -0));

win_frs = squeeze(mean(units(:,w_start:w_end,:),2));

%---------------------------------------------
% only include free choice trials and those in which he made a choice
trials2keep = trialinfo.trialtype==2 & ~isnan(chosen_val); 
units = units(trials2keep,:,:);
chosen_val = chosen_val(trials2keep);
chosen_side = chosen_side(trials2keep);
win_frs = win_frs(trials2keep,:);
%---------------------------------------------


% make the table for the GLM
reg_tbl = table;
reg_tbl.ch_val = chosen_val; 
reg_tbl.side = chosen_side; 



% now do a GLM on each unit

for u = 1:numel(win_frs(1,:))
    
    reg_tbl.FRs = win_frs(:,u);
    
    u_mdl = fitglm(reg_tbl, 'FRs~ ch_val + side');
    
    % value terms
    pvals(u,1) = u_mdl.Coefficients.pValue(2);
    betas(u,1) = u_mdl.Coefficients.Estimate(2);

    % direction terms    
    pvals(u,2) = u_mdl.Coefficients.pValue(3);
    betas(u,2) = u_mdl.Coefficients.Estimate(3);
   
end % of looping over units


end % of function