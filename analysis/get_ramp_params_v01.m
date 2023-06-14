function plot_aligned_ramps_v02(all_pval, all_ipsi,all_contra, choice_ts, all_monkey_ix)

% find the selective units and put them at the top
sel_ix = all_pval(:,1)==1;
sel_monkey_ix = all_monkey_ix(sel_ix);

[~, w_start] = min(abs(choice_ts - -1100));
[~, w_end] = min(abs(choice_ts - 1100));

[ipsi_norm, contra_norm] = normalize_by_row_v02(all_ipsi,all_contra, w_start, w_end);

[n_units, n_times] = size(ipsi_norm);

ipsi_norm = ipsi_norm(:,1:n_times-2);
contra_norm = contra_norm(:,1:n_times-2);

s_ipsi_norm = ipsi_norm(sel_ix,:);
s_contra_norm = contra_norm(sel_ix,:);

s_contra_norm_full = contra_norm(sel_ix,:);
s_ipsi_norm_full = ipsi_norm(sel_ix,:);

s_contra_z = all_contra(sel_ix,w_start:w_end);
s_ipsi_z = all_ipsi(sel_ix,w_start:w_end);


contra_thresh_time = NaN(numel(s_contra_norm(:,1)),1);
contra_peak = NaN(numel(s_contra_norm(:,1)),1);
contra_peak_time = NaN(numel(s_contra_norm(:,1)),1);
contra_ramp_z = NaN(numel(s_contra_norm(:,1)),1);

ipsi_thresh_time = NaN(numel(s_ipsi_norm(:,1)),1);
ipsi_peak = NaN(numel(s_ipsi_norm(:,1)),1);
ipsi_peak_time = NaN(numel(s_ipsi_norm(:,1)),1);

win_times = choice_ts(w_start+2:w_end);



% find the first time the firing rates cross some threshold
thresh = 1;

% find the first time the firing rates cross some threshold, the peak activity, and the time of the peak
for u = 1:numel(s_contra_norm(:,1))
    
    try
        
    contra_thresh_time(u) = min(find(s_contra_norm(u,:)>thresh));
    [contra_peak(u), contra_peak_time(u)] = max(s_contra_norm(u,:));
     
    end
    
end % of looping over units

% find the first time the firing rates cross some threshold, the peak activity, and the time of the peak
for u = 1:numel(s_ipsi_norm(:,1))
    
    try
        
    ipsi_thresh_time(u) = min(find(s_ipsi_norm(u,:)>thresh));
    [ipsi_peak(u), ipsi_peak_time(u)] = max(s_ipsi_norm(u,:));
    
    end
    
end % of looping over units

% find out which cells to keep 
units2_plot = (~isnan(contra_thresh_time) & ~isnan(ipsi_thresh_time));
monkey2keep = sel_monkey_ix(units2_plot);
contra_times = win_times(contra_thresh_time(units2_plot));
ipsi_times = win_times(ipsi_thresh_time(units2_plot));

contra_peaks = contra_peak(units2_plot);
ipsi_peaks = ipsi_peak(units2_plot);

contra_peak_times = win_times(contra_peak_time(units2_plot));
ipsi_peak_times = win_times(ipsi_peak_time(units2_plot));


% sort when the threshold was reached
[~, contra_ix] = sort(contra_thresh_time);
[~, ipsi_ix] = sort(ipsi_thresh_time);

C_sorted_contra_frs = s_contra_norm_full(contra_ix,:);
C_sorted_ipsi_frs = s_ipsi_norm_full(contra_ix,:);

I_sorted_contra_frs = s_contra_norm_full(ipsi_ix,:);
I_sorted_ipsi_frs = s_ipsi_norm_full(ipsi_ix,:);

monkey_ix = sel_monkey_ix == 1;

n_units = numel(ipsi_thresh_time(monkey_ix));

clim = [1,3];
fig = figure('units','inch','position',[2,2,5,7]);
subplot(2,2,1);
imagesc(choice_ts, 1:n_units, C_sorted_contra_frs(monkey_ix & units2_plot,:));
xlim([-1000,1000]);
caxis(clim);
xlabel('Time from Choice (ms)');
ylabel('unit #');
title('Contra Move, Contra Sort');
colormap('hot');
cb = colorbar;
cb.FontSize = 12;
cb.Position = [.19 .62 .02 .05];
cb.Color = 'w';
set(gca,'FontSize',12,'LineWidth',1,'Box','off');

subplot(2,2,2);
imagesc(choice_ts, 1:n_units, C_sorted_ipsi_frs(monkey_ix & units2_plot,:));
xlim([-1000,1000]);
caxis(clim);
colormap('hot');
set(gca,'FontSize',12,'LineWidth',1,'Box','off');
title('Ipsi Move, Contra Sort');

subplot(2,2,3);
imagesc(choice_ts, 1:n_units, I_sorted_contra_frs(monkey_ix & units2_plot,:));
xlim([-1000,1000]);
caxis(clim);
colormap('hot');
set(gca,'FontSize',12,'LineWidth',1,'Box','off');
title('Contra Move, Ipsi Sort');

subplot(2,2,4);
imagesc(choice_ts, 1:n_units, I_sorted_ipsi_frs(monkey_ix & units2_plot,:));
xlim([-1000,1000]);
caxis(clim);
colormap('hot');
set(gca,'FontSize',12,'LineWidth',1,'Box','off');
title('Ipsi Move, Ipsi Sort');

CT = cbrewer('qual','Set1',9);


% now plot the summary figure
fig2 = figure;
hold on
subplot(3,2,1);
plot(contra_times(monkey2keep == 0), ipsi_times(monkey2keep == 0),'.','color',CT(1,:),'MarkerSize',15);
% xlim([-600,0])
% ylim([-600,0])

t_line = lsline;
ylim([-600,0])

set(t_line(1),'color','k','LineWidth',1)
xlabel('Contra ramp onset (ms)');
ylabel('Ipsi ramp onset (ms)');
set(gca, 'FontSize',10, 'LineWidth',1,'Box','off');
title('Animal C');
axis square


[chap_ramp_start_r,chap_ramp_start_p] = corr(contra_times(monkey2keep == 0)', ipsi_times(monkey2keep == 0)');


subplot(3,2,2);
plot(contra_times(monkey2keep == 1), ipsi_times(monkey2keep == 1),'.','color',CT(2,:),'MarkerSize',15);
% xlim([-600,0])
% ylim([-600,0])
t_line = lsline;
set(t_line(1),'color','k','LineWidth',1)
xlabel('Contra ramp onset (ms)');
ylabel('Ipsi ramp onset (ms)');
set(gca, 'FontSize',10, 'LineWidth',1,'Box','off');
title('Animal G');
axis square

[george_ramp_start_r,george_ramp_start_p] = corr(contra_times(monkey2keep == 1)', ipsi_times(monkey2keep == 1)');


subplot(3,2,3);
hold on
plot(contra_peak_times(monkey2keep == 0), ipsi_peak_times(monkey2keep == 0),'.','color',CT(1,:),'MarkerSize',15);
t_line = lsline;
set(t_line(1),'color','k','LineWidth',1)
xlabel('Contra ramp peak time (ms)');
ylabel('Ipsi ramp peak time (ms)');
set(gca, 'FontSize',10, 'LineWidth',1,'Box','off');
axis square


[chap_peak_time_r,chap_peak_time_p] = corr(contra_peak_times(monkey2keep == 0)', ipsi_peak_times(monkey2keep == 0)');


subplot(3,2,4);
hold on
plot(contra_peak_times(monkey2keep == 1), ipsi_peak_times(monkey2keep == 1),'.','color',CT(2,:),'MarkerSize',15);
t_line = lsline;
set(t_line(1),'color','k','LineWidth',1)
xlabel('Contra ramp peak time (ms)');
ylabel('Ipsi ramp peak time (ms)');
set(gca, 'FontSize',10, 'LineWidth',1,'Box','off');
axis square


[george_peak_time_r,george_peak_time_p] = corr(contra_peak_times(monkey2keep == 1)', ipsi_peak_times(monkey2keep == 1)');


subplot(3,2,5);
hold on
plot(contra_peaks(monkey2keep == 0), ipsi_peaks(monkey2keep == 0),'.','color',CT(1,:),'MarkerSize',15);
xlim([1,5]);
ylim([1,5]);
t_line = lsline;
set(t_line(1),'color','k','LineWidth',1)
xlabel('Contra ramp peak (z)');
ylabel('Ipsi ramp peak (z)');
set(gca, 'FontSize',10, 'LineWidth',1,'Box','off');
axis square


[chap_peak_r,chap_peak_p] = corr(contra_peaks(monkey2keep == 0), ipsi_peaks(monkey2keep == 0));

subplot(3,2,6);
hold on
plot(contra_peaks(monkey2keep == 1), ipsi_peaks(monkey2keep == 1),'.','color',CT(2,:),'MarkerSize',15);
xlim([1,5]);
ylim([1,5]);
t_line = lsline;
set(t_line(1),'color','k','LineWidth',1)
xlabel('Contra ramp peak (z)');
ylabel('Ipsi ramp peak (z)');
set(gca, 'FontSize',10, 'LineWidth',1,'Box','off');
axis square


[george_peak_r, george_peak_p] = corr(contra_peaks(monkey2keep == 1), ipsi_peaks(monkey2keep == 1));



end % of function