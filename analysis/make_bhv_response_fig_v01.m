function make_bhv_response_fig_v01(bhv)
%-------------------------------------------------------------------------------
% 1. How many saccades before lever?
% 2. How often was the first saccade to the ultimately-chosen option?
% 3. What was the distribution of sacc. response times (first and second)?
% 4. How risk seeking or risk averse were the animals?
% 5. How large was the lever side bias?
% 6. How large was the first-saccade bias?
%-------------------------------------------------------------------------------

free_ix = bhv.trialtype==2;

% 1. how many saccades before the lever?
[sacc_means, sacc_sems] = grpstats(bhv.saccN(free_ix), bhv.monkey(free_ix),{'mean','sem'});
[sacc_modes] = grpstats(bhv.saccN(free_ix), bhv.monkey(free_ix),{'mode'});

% 2. How often was the first sacc. to the ultimately-chosen option?
first_sacc_to_lever = bhv.first_sacc_loc == bhv.lever;
[first_sacc_to_lever_means, first_sacc_to_lever_sems] = grpstats(first_sacc_to_lever(free_ix), bhv.monkey(free_ix),{'mean','sem'});

% 3. What was the distribtuion of sacc. response times?
bin_width = 20;
figure; 
subplot(1,2,1)
hold on
histogram(bhv.sacc1_rt(bhv.monkey==0 & free_ix),'BinWidth',bin_width,'Normalization','probability','EdgeColor','none');
histogram(bhv.sacc2_rt(bhv.monkey==0 & free_ix),'BinWidth',bin_width,'Normalization','probability','EdgeColor','none');
title('Animal C');
xlabel('Time from Pics On (ms)');
legend('1st saccade','2nd saccade');
xlim([1,1000]);
set(gca, 'FontSize',12,'LineWidth',1);
axis square


subplot(1,2,2)
hold on
histogram(bhv.sacc1_rt(bhv.monkey==1 & free_ix),'BinWidth',bin_width,'Normalization','probability','EdgeColor','none');
histogram(bhv.sacc2_rt(bhv.monkey==1 & free_ix),'BinWidth',bin_width,'Normalization','probability','EdgeColor','none');
title('Animal G');
xlabel('Time from Pics On (ms)');
xlim([1,1000]);
set(gca, 'FontSize',12,'LineWidth',1);
axis square


sacc1_medians = grpstats(bhv.sacc1_rt(free_ix), bhv.monkey(free_ix),{'median'});

% get details of side bias
trials2use = free_ix & ~isnan(bhv.first_sacc_loc);

n_chap_trials = sum(trials2use & bhv.monkey==0);
n_george_trials = sum(trials2use & bhv.monkey==1);

n_george_right_lever = sum(bhv.lever(trials2use &  bhv.monkey==1)==1);
n_chap_right_lever = sum(bhv.lever(trials2use &  bhv.monkey==0)==1);

n_george_right_sacc = sum(bhv.first_sacc_loc(trials2use &  bhv.monkey==1)==1);
n_chap_right_sacc = sum(bhv.first_sacc_loc(trials2use &  bhv.monkey==0)==1);





% 4. How risk seeking/averse were the animals?
[distortion_tbl] = get_prob_weighting_by_mag_v02(bhv);
m1_probs = distortion_tbl.prob_ids(1,:);
m2_probs = distortion_tbl.prob_ids(2,:);
p_means = distortion_tbl.pChooseP_mean;
p_sems = distortion_tbl.pChooseP_sem; 
p_best_means = distortion_tbl.pChooseP_isbest;
p_best_sems = distortion_tbl.pChooseP_isbest_sem;
p_distortion_mean = distortion_tbl.pBias_mean;
p_distortion_sem = distortion_tbl.pBias_sem;


CT = cbrewer('qual','Set1',9);


figure; 
subplot(1,3,1);
hold on
plot(m1_probs, p_means(1,:),'LineWidth',3,'color',CT(1,:));
plot(m1_probs, p_best_means(1,:),'LineWidth',3,'color',CT(9,:));
errorbar(m1_probs, p_means(1,:), p_sems(1,:),'LineWidth',3,'color',CT(1,:));
xticks([0,.5,1])
yticks([0,.5,1]);
xlabel('Probability in Trial');
ylabel('p(Choose Prob.)');
ylim([0,1]);
xlim([0,1]);
set(gca,'FontSize',12,'LineWidth',1);
legend({'Monkey C','Empirical Best'},'Location','northwest');
axis square


subplot(1,3,2);
hold on
plot(m2_probs, p_means(2,:),'LineWidth',3,'color',CT(2,:));
plot(m2_probs, p_best_means(2,:),'LineWidth',3,'color',CT(9,:));
errorbar(m2_probs, p_means(2,:), p_sems(2,:),'LineWidth',3,'color',CT(2,:));
xticks([0,.5,1])
yticks([0,.5,1]);
xlabel('Probability in Trial');
ylabel('p(Choose Prob.)');
ylim([0,1]);
xlim([0,1]);
xticks([0,.5,1])
yticks([0,.5,1]);
legend({'Monkey G','Empirical Best'},'Location','northwest');
set(gca,'FontSize',12,'LineWidth',1);
axis square


subplot(1,3,3);
hold on
plot(m1_probs, p_distortion_mean(1,:),'color',CT(1,:),'LineWidth',3);
plot(m2_probs, p_distortion_mean(2,:),'color',CT(2,:),'LineWidth',3);
errorbar(m1_probs, p_distortion_mean(1,:), p_distortion_sem(1,:), 'color',CT(1,:),'LineWidth',3);  
errorbar(m2_probs, p_distortion_mean(2,:), p_distortion_sem(2,:), 'color',CT(2,:),'LineWidth',3);  
plot([0,1],[0,0],'color',CT(9,:),'LineWidth',2);
xlim([0,1]);
ylim([-.2,.2]);
xticks([0,.5,1])
xlabel('Probability in Trial');
ylabel('Choice Bias Relative to Optimal');
set(gca,'FontSize',12,'LineWidth',1);
legend({'Monkey C','Monkey G'},'Location','northwest');
axis square


valid_eyes = ~isnan(bhv.first_sacc_loc);

% 5. How large was the lever side bias?
% express as mean direction during error trials (since that's most likely to reflect a bias).
total_lever_side_means =grpstats(bhv.lever(free_ix & valid_eyes), bhv.monkey(free_ix & valid_eyes));

% 6. How large was the first-saccade bias?
total_sacc_side_means =grpstats(bhv.first_sacc_loc(free_ix & valid_eyes), bhv.monkey(free_ix & valid_eyes));

mean(bhv.first_sacc_loc(free_ix & valid_eyes & bhv.monkey==0)==1) - .5



end % of function