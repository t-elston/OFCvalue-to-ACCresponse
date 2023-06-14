
[ch_rvals, unch_rvals, corr_animal_id, t_mids] = get_sessionwiseACCOFC_valstate_correlations;

[c_ch_corr_mean, c_ch_corr_sem] = GetMeanCI(ch_rvals(corr_animal_id==0,:),'sem');
[c_unch_corr_mean, c_unch_corr_sem] = GetMeanCI(unch_rvals(corr_animal_id==0,:),'sem');
[g_ch_corr_mean, g_ch_corr_sem] = GetMeanCI(ch_rvals(corr_animal_id==1,:),'sem');
[g_unch_corr_mean, g_unch_corr_sem] = GetMeanCI(unch_rvals(corr_animal_id==1,:),'sem');

[chap_sig_ch] = get_significant_correlation_times(ch_rvals(corr_animal_id==0,:), 24, t_mids);
[chap_sig_unch] = get_significant_correlation_times(unch_rvals(corr_animal_id==0,:), 24, t_mids);
[george_sig_ch] = get_significant_correlation_times(ch_rvals(corr_animal_id==1,:), 24, t_mids);
[george_sig_unch] = get_significant_correlation_times(unch_rvals(corr_animal_id==1,:), 24, t_mids);


CT = cbrewer('qual','Set1',9);

Fig = figure;
set(Fig,'renderer','Painters');
subplot(1,2,1);
hold on
shadedErrorBar(t_mids,c_unch_corr_mean, c_unch_corr_sem,'LineProps',{'LineWidth',1,'color',CT(2,:)});
shadedErrorBar(t_mids,c_ch_corr_mean, c_ch_corr_sem,'LineProps',{'LineWidth',1,'color',CT(1,:)});
plot(t_mids(chap_sig_ch), ones(size(chap_sig_ch))*-.02,'s', 'color',CT(1,:),'MarkerFaceColor',CT(1,:));
plot(t_mids(chap_sig_unch), ones(size(chap_sig_unch))*-.03,'s', 'color',CT(2,:),'MarkerFaceColor',CT(2,:));
ylim([-.05,.2]);

plot(xlim,[0 0],'k','LineWidth',1);
xlim([-400,600]);
ylabel('ACC-OFC PP Correlation');
axis square
title('Animal C');
xlabel('Time from Pictures On');
set(gca, 'FontSize',10, 'LineWidth',1);


subplot(1,2,2);
hold on
shadedErrorBar(t_mids,g_unch_corr_mean, g_unch_corr_sem,'LineProps',{'LineWidth',1,'color',CT(2,:)});
shadedErrorBar(t_mids,g_ch_corr_mean, g_ch_corr_sem,'LineProps',{'LineWidth',1,'color',CT(1,:)});
plot(t_mids(george_sig_ch), ones(size(george_sig_ch))*-.02,'s', 'color',CT(1,:),'MarkerFaceColor',CT(1,:));
plot(t_mids(george_sig_unch), ones(size(george_sig_unch))*-.03,'s', 'color',CT(2,:),'MarkerFaceColor',CT(2,:));
plot(xlim,[0 0],'k','LineWidth',1);
xlim([-400,600]);
ylim([-.05,.2]);

axis square
title('Animal G');
set(gca, 'FontSize',10, 'LineWidth',1);