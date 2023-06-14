% aaa_ACC_subjective_value_analysis_v01

% top of stack for assessing choice behavior during value-based choice
% author: Thomas Elston 

bhv_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\bhv\';

% load all of the behavior
bhv = load_all_bhv_v01(bhv_dir);

make_bhv_response_fig_v01(bhv);


%-------------
% make RT plots
bhv2 = bhv(bhv.trialtype==2,:);
c_ix = bhv2.monkey==0;
g_ix = bhv2.monkey==1;

reg_table = table;
reg_table.max_val =  max([bhv2.left_expval, bhv2.right_expval],[],2);
reg_table.val_diff = abs(round(bhv2.left_expval - bhv2.right_expval,1));
reg_table.rt = bhv2.rt;

% fit two models per animal - one for max value, another for delta-value
c_max_val_mdl = fitglm(reg_table(c_ix,:), 'rt ~ max_val');
c_val_diff_mdl = fitglm(reg_table(c_ix,:), 'rt ~ val_diff');
g_max_val_mdl = fitglm(reg_table(g_ix,:), 'rt ~ max_val');
g_val_diff_mdl = fitglm(reg_table(g_ix,:), 'rt ~ val_diff');

% full models
chap_rt_mdl = fitglm(reg_table(c_ix,:), 'rt ~ max_val + val_diff');
george_rt_mdl = fitglm(reg_table(g_ix,:), 'rt ~ max_val + val_diff');


g_max_vals = unique(reg_table.max_val(g_ix));
c_max_vals = unique(reg_table.max_val(c_ix));
g_val_diffs = unique(reg_table.val_diff(g_ix));
c_val_diffs = unique(reg_table.val_diff(c_ix));

% now get the means and SEMs for the alternative models based on residuals
% max_val residuals grouped by value difference
for i = 1:numel(g_val_diffs)
    
    [g_max_val_mean(i,1), g_max_val_sem(i,1)] = GetMeanCI(g_max_val_mdl.Residuals.Raw(reg_table.val_diff(g_ix) == g_val_diffs(i)),'sem'); 
    
end % of looping over animal g

for i = 1:numel(c_val_diffs)
    
    [c_max_val_mean(i,1), c_max_val_sem(i,1)] = GetMeanCI(c_max_val_mdl.Residuals.Raw(reg_table.val_diff(c_ix) == c_val_diffs(i)),'sem'); 
    
end % of looping over animal c

for i = 1:numel(g_max_vals)
    
    [g_valdiff_mean(i,1), g_valdiff_sem(i,1)] = GetMeanCI(g_val_diff_mdl.Residuals.Raw(reg_table.max_val(g_ix) == g_max_vals(i)),'sem'); 
    
end % of looping over animal g

for i = 1:numel(c_max_vals)
    
    [c_valdiff_mean(i,1), c_valdiff_sem(i,1)] = GetMeanCI(c_val_diff_mdl.Residuals.Raw(reg_table.max_val(c_ix) == c_max_vals(i)),'sem'); 
    
end % of looping over animal c


CT = cbrewer('qual','Set1',9);

figure; 
subplot(2,2,1);
hold on
errorbar(c_val_diffs, c_max_val_mean, c_max_val_sem,'.','CapSize',0,'color',CT(3,:),'LineWidth',1,'MarkerSize',10)
title('Animal C');
xlabel('abs(delta value)');
ylabel({'Residual RT (ms)', 'after max value'});
axis square

subplot(2,2,2);
hold on
errorbar(g_val_diffs, g_max_val_mean, g_max_val_sem,'.','CapSize',0,'color',CT(3,:),'LineWidth',1,'MarkerSize',10)
title('Animal G');
xlabel('abs(delta value)');
axis square


subplot(2,2,3);
hold on
errorbar(c_max_vals, c_valdiff_mean, c_valdiff_sem,'.','CapSize',0,'color',CT(3,:),'LineWidth',1,'MarkerSize',10)
xlabel('Max Value');
ylabel({'Residual RT (ms)', 'after abs(value delta)'});
axis square

subplot(2,2,4);
hold on
errorbar(g_max_vals, g_valdiff_mean, g_valdiff_sem,'.','CapSize',0,'color',CT(3,:),'LineWidth',1,'MarkerSize',10)
xlabel('Max Value');
axis square





g_max_val_residuals = g_max_val_mdl.Residuals;
