function [pvals, betas] = fixed_window_rel2choice_GLMs(units, ts, t_ch_val, t_unch_val, trial_types, chosen_side)

% get indices of the fixed window
[~,w_start] = min(abs(ts- -450));
[~,w_end] = min(abs(ts- -0));

win_frs = squeeze(mean(units(:,w_start:w_end,:),2));

%---------------------------------------------
% only include free choice trials and those in which he made a choice
trials2keep = trial_types==2 & ~isnan(t_ch_val);
units = units(trials2keep,:,:);
t_ch_val = t_ch_val(trials2keep);
chosen_side = chosen_side(trials2keep);
win_frs = win_frs(trials2keep,:);
%---------------------------------------------


% make the table for the GLM
reg_tbl = table;
reg_tbl.ch_val = t_ch_val; 
reg_tbl.side = chosen_side; 



% now do a GLM on each unit

for u = 1:numel(win_frs(1,:))
    
    reg_tbl.FRs = win_frs(:,u);
    
    u_mdl = fitglm(reg_tbl, 'FRs~ ch_val + side');
    
    pvals(u,1) = u_mdl.Coefficients.pValue(3);
    betas(u,1) = u_mdl.Coefficients.Estimate(3);

    
end % of looping over units


end % of function