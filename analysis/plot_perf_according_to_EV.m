function plot_perf_according_to_EV(prob_tbl, bhv)

[prob_means, prob_sems] = grpstats(prob_tbl.prob_acc,prob_tbl.monkey,{'mean','sem'});
m1_probs = prob_tbl.prob_ids(1,:);
m2_probs = prob_tbl.prob_ids(9,:);

figure; 
hold on
errorbar(m1_probs, prob_means(1,:), prob_sems(1,:),'LineWidth',2);
errorbar(m2_probs, prob_means(2,:), prob_sems(2,:),'LineWidth',2);
xlabel('Probability in Trial');
ylabel('p(Choose Best EV)');
ylim([0,1]);
xlim([0,1]);
set(gca,'FontSize',12,'LineWidth',1);
legend({'Monkey C','Monkey G'});




end % of function