
load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\OFC alignment differences.mat')
load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\ACC alignment differences.mat')
load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\ACC_valdec_output.mat');

t_mids = valdata_choice.t_mids+1;

xt = [1:numel(chap_OFC_diff(1,:))]*mean(diff(t_mids)) - (round(numel(chap_OFC_diff(1,:))/2) * mean(diff(t_mids)));
[~,xt_after_start] = min(abs(xt-0));
[~,xt_after_end] = min(abs(xt-100));
[~,xt_before_start] = min(abs(xt- -70));
[~,xt_before_end] = min(abs(xt-0));

chap_diffs = chap_ACC_diff - chap_OFC_diff;
george_diffs = george_ACC_diff - george_OFC_diff;


chap_before_diffs = mean(chap_diffs(:, xt_before_start: xt_before_end),2);
george_before_diffs = mean(george_diffs(:, xt_before_start: xt_before_end),2);

CT = cbrewer('qual','Set1',9);

chap_first = prctile(chap_before_diffs,1);
george_first = prctile(george_before_diffs,1);



figure; 
subplot(1,2,1); 
hold on
histogram(chap_before_diffs,'EdgeColor','None','FaceColor',CT(9,:), 'Normalization','probability');
plot([chap_first,chap_first],ylim,'k');
xlabel('ACC - OFC Delta DIRacc')
ylabel('Probability')
title('Animal C');   
axis square

subplot(1,2,2);
hold on
histogram(george_before_diffs,'EdgeColor','None','FaceColor',CT(9,:),'Normalization','probability');
plot([george_first,george_first],ylim,'k');

title('Animal G');   
axis square
