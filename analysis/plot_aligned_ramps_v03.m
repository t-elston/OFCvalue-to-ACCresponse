function plot_aligned_ramps_v03(all_pval, all_ipsi,all_contra, choice_ts, all_monkey_ix)

% find the selective units and put them at the top
sel_ix = all_pval(:,1)==1;
sel_monkey_ix = all_monkey_ix(sel_ix);

[~, w_start] = min(abs(choice_ts - -1000));
[~, w_end] = min(abs(choice_ts - 1000));

[ipsi_norm, contra_norm] = normalize_by_row(all_ipsi,all_contra);


s_ipsi_norm = ipsi_norm(sel_ix,w_start:w_end);
s_contra_norm = contra_norm(sel_ix,w_start:w_end);

s_contra_norm_full = contra_norm(sel_ix,:);
s_ipsi_norm_full = ipsi_norm(sel_ix,:);




contra_thresh_ix = NaN(numel(s_contra_norm(:,1)),1);
contra_peak = NaN(numel(s_contra_norm(:,1)),1);
contra_peak_time = NaN(numel(s_contra_norm(:,1)),1);
contra_times = NaN(numel(s_contra_norm(:,1)),1);
contra_peak_times = NaN(numel(s_contra_norm(:,1)),1);

ipsi_thresh_ix = NaN(numel(s_ipsi_norm(:,1)),1);
ipsi_peak = NaN(numel(s_ipsi_norm(:,1)),1);
ipsi_peak_time = NaN(numel(s_ipsi_norm(:,1)),1);

ipsi_times = NaN(numel(s_contra_norm(:,1)),1);
ipsi_peak_times = NaN(numel(s_contra_norm(:,1)),1);

win_times = choice_ts(w_start:w_end);



% find the first time the firing rates cross some threshold
thresh = 1;

% find the first time the firing rates cross some threshold, the peak activity, and the time of the peak
for u = 1:numel(s_contra_norm(:,1))
    
    try
        
    contra_thresh_ix(u) = min(find(s_contra_norm(u,:)>thresh));
    [contra_peak(u), contra_peak_time(u)] = max(s_contra_norm(u,:));
    
    contra_times(u) = win_times(contra_thresh_ix(u));
    contra_peak_times(u) = win_times(contra_peak_time(u));
     
    end
    
end % of looping over units

% find the first time the firing rates cross some threshold, the peak activity, and the time of the peak
for u = 1:numel(s_ipsi_norm(:,1))
    
    try
        
    ipsi_thresh_ix(u) = min(find(s_ipsi_norm(u,:)>thresh));
    [ipsi_peak(u), ipsi_peak_time(u)] = max(s_ipsi_norm(u,:));
    
    ipsi_times(u) = win_times(ipsi_thresh_ix(u));
    ipsi_peak_times(u) = win_times(ipsi_peak_time(u));
    
    end
    
end % of looping over units

% sort when the threshold was reached
[~, contra_ix] = sort(contra_thresh_ix);
[~, ipsi_ix] = sort(ipsi_thresh_ix);

C_sorted_contra_frs = s_contra_norm_full(contra_ix,:);
C_sorted_ipsi_frs = s_ipsi_norm_full(contra_ix,:);

I_sorted_contra_frs = s_contra_norm_full(ipsi_ix,:);
I_sorted_ipsi_frs = s_ipsi_norm_full(ipsi_ix,:);


% find out which cells to keep 
good_units = (~isnan(contra_thresh_ix) | ~isnan(ipsi_thresh_ix)); 
bad_units = ~good_units;

monkey2keep = sel_monkey_ix(good_units);

contra_times = contra_times(good_units);
ipsi_times = ipsi_times(good_units);

contra_peaks = contra_peak(good_units);
ipsi_peaks = ipsi_peak(good_units);

contra_peak_times = contra_peak_times(good_units);
ipsi_peak_times = ipsi_peak_times(good_units);


contra_peak_times = contra_peak_times-50;
ipsi_peak_times = ipsi_peak_times-50;

contra_times = contra_times-50;
ipsi_times = ipsi_times-50;

sum(ipsi_times(monkey2keep==1) < 50 | ipsi_times(monkey2keep==1) < 50)
sum(ipsi_times(monkey2keep==0) < 50 | ipsi_times(monkey2keep==0) < 50)



% now windsorize the data 
contra_peaks(contra_peak_times > 50 | contra_times > 50) = NaN;
contra_peaks(ipsi_peak_times > 50 | ipsi_times > 50) = NaN;

contra_times(contra_times>50) = 0;
ipsi_times(ipsi_times>50) = 0;

contra_peak_times(contra_peak_times>50) = 0;
ipsi_peak_times(ipsi_peak_times>50) = 0;

n_c_units = sum(sel_monkey_ix == 0);
n_g_units = sum(sel_monkey_ix == 1);




clim = [1,3];
fig = figure('units','inch','position',[2,2,10,7]);
subplot(2,4,1);
imagesc(choice_ts, 1:n_c_units, C_sorted_contra_frs(sel_monkey_ix(contra_ix) == 0,:));
xlim([-1000,1000]);
caxis(clim);
xlabel('Time from Choice (ms)');
ylabel('unit #');
title('Contra Move, Contra Sort');
colormap('hot');
cb = colorbar;
cb.FontSize = 12;
cb.Position = [.17 .62 .01 .05];
cb.Color = 'w';
set(gca,'FontSize',12,'LineWidth',1,'Box','off');

subplot(2,4,2);
imagesc(choice_ts, 1:n_c_units, C_sorted_ipsi_frs(sel_monkey_ix(contra_ix) == 0,:));
xlim([-1000,1000]);
caxis(clim);
colormap('hot');
set(gca,'FontSize',12,'LineWidth',1,'Box','off');
title('Ipsi Move, Contra Sort');

subplot(2,4,5);
imagesc(choice_ts, 1:n_c_units, I_sorted_contra_frs(sel_monkey_ix(ipsi_ix) == 0,:));
xlim([-1000,1000]);
caxis(clim);
colormap('hot');
set(gca,'FontSize',12,'LineWidth',1,'Box','off');
title('Contra Move, Ipsi Sort');

subplot(2,4,6);
imagesc(choice_ts, 1:n_c_units, I_sorted_ipsi_frs(sel_monkey_ix(ipsi_ix) == 0,:));
xlim([-1000,1000]);
caxis(clim);
colormap('hot');
set(gca,'FontSize',12,'LineWidth',1,'Box','off');
title('Ipsi Move, Ipsi Sort');

%----
subplot(2,4,3);
imagesc(choice_ts, 1:n_g_units, C_sorted_contra_frs(sel_monkey_ix(contra_ix) == 1,:));
xlim([-1000,1000]);
caxis(clim);
xlabel('Time from Choice (ms)');
ylabel('unit #');
title('Contra Move, Contra Sort');
colormap('hot');
set(gca,'FontSize',12,'LineWidth',1,'Box','off');

subplot(2,4,4);
imagesc(choice_ts, 1:n_g_units, C_sorted_ipsi_frs(sel_monkey_ix(contra_ix) == 1,:));
xlim([-1000,1000]);
caxis(clim);
colormap('hot');
set(gca,'FontSize',12,'LineWidth',1,'Box','off');
title('Ipsi Move, Contra Sort');

subplot(2,4,7);
imagesc(choice_ts, 1:n_g_units, I_sorted_contra_frs(sel_monkey_ix(ipsi_ix) == 1,:));
xlim([-1000,1000]);
caxis(clim);
colormap('hot');
set(gca,'FontSize',12,'LineWidth',1,'Box','off');
title('Contra Move, Ipsi Sort');

subplot(2,4,8);
imagesc(choice_ts, 1:n_g_units, I_sorted_ipsi_frs(sel_monkey_ix(ipsi_ix) == 1,:));
xlim([-1000,1000]);
caxis(clim);
colormap('hot');
set(gca,'FontSize',12,'LineWidth',1,'Box','off');
title('Ipsi Move, Ipsi Sort');

%-----------


CT = cbrewer('qual','Set1',9);



% now plot the summary figure
fig2 = figure;
hold on
subplot(3,2,1);
plot(contra_times(monkey2keep == 0), ipsi_times(monkey2keep == 0),'.','color',CT(1,:),'MarkerSize',15);
xlim([-1000,0])
ylim([-1000,0])
xticks([-1000,-500,0]);
yticks([-1000,-500,0]);
t_line = lsline;
xlim([-1000,0])
ylim([-1000,0])
set(t_line(1),'color','k','LineWidth',1)
xlabel('Contra ramp onset (ms)');
ylabel('Ipsi ramp onset (ms)');
title('Animal C');
axis square


[chap_ramp_start_r, chap_ramp_start_p] = corr(contra_times(monkey2keep == 0), ipsi_times(monkey2keep == 0),'Rows','complete');


subplot(3,2,2);
plot(contra_times(monkey2keep == 1), ipsi_times(monkey2keep == 1),'.','color',CT(2,:),'MarkerSize',15);
xlim([-1000,0])
ylim([-1000,0])
xticks([-1000,-500,0]);
yticks([-1000,-500,0]);
t_line = lsline;
set(t_line(1),'color','k','LineWidth',1)
xlabel('Contra ramp onset (ms)');
ylabel('Ipsi ramp onset (ms)');
title('Animal G');
axis square

[george_ramp_start_r,george_ramp_start_p] = corr(contra_times(monkey2keep == 1), ipsi_times(monkey2keep == 1),'Rows','complete');


subplot(3,2,3);
hold on
plot(contra_peak_times(monkey2keep == 0), ipsi_peak_times(monkey2keep == 0),'.','color',CT(1,:),'MarkerSize',15);
xlim([-1000,0])
ylim([-1000,0])
xticks([-1000,-500,0]);
yticks([-1000,-500,0]);
t_line = lsline;
set(t_line(1),'color','k','LineWidth',1)
xlabel('Contra ramp peak time (ms)');
ylabel('Ipsi ramp peak time (ms)');
axis square


[chap_peak_time_r,chap_peak_time_p] = corr(contra_peak_times(monkey2keep == 0), ipsi_peak_times(monkey2keep == 0),'Rows','complete');


subplot(3,2,4);
hold on
plot(contra_peak_times(monkey2keep == 1), ipsi_peak_times(monkey2keep == 1),'.','color',CT(2,:),'MarkerSize',15);
xlim([-1000,0])
ylim([-1000,0])
xticks([-1000,-500,0]);
yticks([-1000,-500,0]);
t_line = lsline;
xlim([-1000,0])
ylim([-1000,0])
set(t_line(1),'color','k','LineWidth',1)
xlabel('Contra ramp peak time (ms)');
ylabel('Ipsi ramp peak time (ms)');
axis square


[george_peak_time_r,george_peak_time_p] = corr(contra_peak_times(monkey2keep == 1), ipsi_peak_times(monkey2keep == 1),'Rows','complete');


subplot(3,2,5);
hold on
plot(contra_peaks(monkey2keep == 0), ipsi_peaks(monkey2keep == 0),'.','color',CT(1,:),'MarkerSize',15);
xlim([1,5]);
ylim([1,5]);
xticks([1,3,5]);
yticks([1,3,5]);
t_line = lsline;
set(t_line(1),'color','k','LineWidth',1)
xlabel('Contra ramp peak (z)');
ylabel('Ipsi ramp peak (z)');
axis square


[chap_peak_r,chap_peak_p] = corr(contra_peaks(monkey2keep == 0), ipsi_peaks(monkey2keep == 0),'Rows','complete');

subplot(3,2,6);
hold on
plot(contra_peaks(monkey2keep == 1), ipsi_peaks(monkey2keep == 1),'.','color',CT(2,:),'MarkerSize',15);
xlim([1,5]);
ylim([1,5]);
xticks([1,3,5]);
yticks([1,3,5]);
t_line = lsline;
set(t_line(1),'color','k','LineWidth',1)
xlabel('Contra ramp peak (z)');
ylabel('Ipsi ramp peak (z)');
axis square


[george_peak_r, george_peak_p] = corr(contra_peaks(monkey2keep == 1), ipsi_peaks(monkey2keep == 1),'Rows','complete');



end % of function