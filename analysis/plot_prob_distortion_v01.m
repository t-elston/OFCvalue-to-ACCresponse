function plot_prob_distortion_v01(distortion_tbl)

monkey_ids = unique(distortion_tbl.monkey);

ctr = 0;

for m = 1:numel(monkey_ids)
    
    this_monkey_data = distortion_tbl(distortion_tbl.monkey == monkey_ids(m),:);
    probs = this_monkey_data.prob_ids(1,:);
    
    x_probs = [probs;probs;probs;probs];
    
    % comput summary stats
    [mag_means, mag_sems, mag_ids] = grpstats(this_monkey_data.pChooseP,this_monkey_data.amnt,{'mean','sem','gname'});
    [d_means, d_sems, mag_ids] = grpstats(this_monkey_data.pDistortion,this_monkey_data.amnt,{'mean','sem','gname'});
    
    ctr = ctr+1;
    subplot(2, 2, ctr);
    hold on
    errorbar(x_probs',mag_means', mag_sems','LineWidth',3);
    plot([0,1],[0,1],'k','LineWidth',1);
    legend(mag_ids,'Location','northwest')
    xlim([0,1])
    axis tight
    ylabel('p(Choose P)');
    xlabel('Probability');
    set(gca,'FontSize',12,'LineWidth',1);
    
    if m ==1
        title('p(Choose P) by Amount');   
    end
    
    ctr=ctr+1;
    subplot(2, 2, ctr);
    hold on
    errorbar(x_probs',d_means', d_sems','LineWidth',3);
    plot([0,1],[0,0],'k','LineWidth',1);
    axis tight
    ylabel('\Delta Unity');
    xlabel('Probability');
    set(gca,'FontSize',12,'LineWidth',1);
    if m ==1
        title('Prob. Distortion by Amount');
    end
    
    
end % of looping over animals





end % of function