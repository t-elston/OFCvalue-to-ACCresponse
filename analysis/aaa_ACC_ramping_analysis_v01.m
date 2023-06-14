% aaa_ACC_ramping_analysis_v01
% top of stack for assessing ACC ramping
% author: Thomas Elston 

bhv_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\bhv\';
rec_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\recdata\';
alignment = 'choice';
brain_area = 'ACC';

[all_times, all_betas, all_hemi_ix, all_monkey_ix, all_ipsi, all_contra, choice_ts, all_pval, all_beta, all_meanFRs, all_CPD] =...
                                                                                    assessing_ACC_ramping_v01(bhv_dir, rec_dir, alignment, brain_area);

plot_aligned_ramps_v03(all_pval, all_ipsi,all_contra, choice_ts, all_monkey_ix);


%----------------------
% is there a positive coding bias in ACC?

val_sig = logical(all_pval(:,2));

n_val_sig = sum(val_sig);

n_positive_val_sig = sum(all_beta(val_sig, 3)> 0);
n_negative_val_sig = sum(all_beta(val_sig, 3)< 0);

[tbl,chi2stat,pval] = chiSquareWithFrequencies_v02(n_positive_val_sig,n_val_sig,n_negative_val_sig,n_val_sig);
%----------------------
% is there a direction-encoding bias in ACC?
c_ix = all_monkey_ix ==0;
g_ix = all_monkey_ix ==1;
dir_sig = logical(all_pval(:,1));
side_beta = all_beta(:,1) .* all_hemi_ix;
CI_beta = all_beta(:,1);

c_n_dir_sig = sum(dir_sig & c_ix);
c_n_right = sum(side_beta(dir_sig & c_ix)> 0);
c_n_left = sum(side_beta(dir_sig & c_ix)< 0);
c_n_ipsi = sum(CI_beta(dir_sig & c_ix)> 0);
c_n_contra = sum(CI_beta(dir_sig & c_ix)< 0);


g_n_dir_sig = sum(dir_sig & g_ix);
g_n_right = sum(side_beta(dir_sig & g_ix)> 0);
g_n_left = sum(side_beta(dir_sig & g_ix)< 0);
g_n_ipsi = sum(CI_beta(dir_sig & g_ix)> 0);
g_n_contra = sum(CI_beta(dir_sig & g_ix)< 0);


figure; 
%---------
subplot(2,2,1);
plot(all_meanFRs(val_sig & c_ix,1), all_beta(val_sig & c_ix,2),'.','color',[.5,.5,.5],'MarkerSize',15)
ylim([0,.5]);
xlim([0,60]);
l1 = lsline; l1.Color='k'; l1.LineWidth=2;

[c_r_fix, c_p_fix] = corr(all_meanFRs(val_sig & c_ix,1), all_beta(val_sig & c_ix,2));

xlabel('Firing Rate at Fix (Hz)'); 
ylabel('CPD');
title('Animal C');
set(gca,'FontSize',12,'LineWidth',1,'Box','off');
axis square

%---------
subplot(2,2,2);
plot(all_meanFRs(val_sig & g_ix,1), all_beta(val_sig & g_ix,2),'.','color',[.5,.5,.5],'MarkerSize',15)
ylim([0,.5]);
xlim([0,60]);
l2 = lsline; l2.Color='k'; l2.LineWidth=2;

[g_r_fix, g_p_fix] = corr(all_meanFRs(val_sig & g_ix,1), all_beta(val_sig & g_ix,2));

xlabel('Firing Rate at Fix (Hz)'); 
title('Animal G');
set(gca,'FontSize',12,'LineWidth',1,'Box','off');
axis square


%---------
subplot(2,2,3);
plot(all_meanFRs(val_sig & c_ix,2), all_beta(val_sig & c_ix,2),'.','color',[.5,.5,.5],'MarkerSize',15)
ylim([0,.5]);
xlim([0,60]);
l3 = lsline; l3.Color='k'; l3.LineWidth=2;

[c_r_choice, c_p_choice] = corr(all_meanFRs(val_sig & c_ix,2), all_beta(val_sig & c_ix,2));

xlabel('Firing Rate at Choice (Hz)'); 
ylabel('CPD');
set(gca,'FontSize',12,'LineWidth',1,'Box','off');
axis square

%---------
subplot(2,2,4);
plot(all_meanFRs(val_sig & g_ix,2), all_beta(val_sig & g_ix,2),'.','color',[.5,.5,.5],'MarkerSize',15)
ylim([0,.5]);
xlim([0,60]);
l4 = lsline; l4.Color='k'; l4.LineWidth=2;

[g_r_choice, g_p_choice] = corr(all_meanFRs(val_sig & g_ix,2), all_beta(val_sig & g_ix,2));

xlabel('Firing Rate at Choice (Hz)'); 
set(gca,'FontSize',12,'LineWidth',1,'Box','off');
axis square

