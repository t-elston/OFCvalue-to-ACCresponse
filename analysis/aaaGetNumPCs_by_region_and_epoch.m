
% aaaGetNumPCs_by_region_and_epoch

bhv_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\bhv\';
rec_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\recdata\';
alignment = 'choice';
brain_area = 'ACC';



% find names of bhv files and rec files
bhv_folder_names = get_subfolder_names(bhv_dir);

n_files = numel(bhv_folder_names);

ACC_PC95=[];
OFC_PC95=[];

% now loop over each file
for f_ix = 1:n_files
    
    current_file = bhv_folder_names(f_ix);
    monkey_id = contains(current_file,'George'); % 0 = Cap; 1 = George
    
    fprintf('\n')
    disp(['***' current_file{1}]);
    
    % pull out this file's behavior
    [trialinfo] = load_behavior(bhv_dir, current_file);
    
    % pull out this file's neurons
    [ACC_choice_units, ~, choice_ts] = load_neurons(rec_dir, current_file, 'choice', 'ACC', [-500, 0]);
    [OFC_choice_units, ~, ~] = load_neurons(rec_dir, current_file, 'choice', 'OFC', [-500, 0]);
    
    [ACC_pics_units, ~, choice_ts] = load_neurons(rec_dir, current_file, 'pics', 'ACC', [0, 500]);
    [OFC_pics_units, ~, ~] = load_neurons(rec_dir, current_file, 'pics', 'OFC', [0, 500]);
    
    
    % now get the mean firing rates for the PC analysis
    f_ACC_PC95=[];
    f_OFC_PC95=[];

    
    if ~isempty(ACC_choice_units)
        
        n_acc_units = numel(ACC_choice_units(1,1,:));
        
        ACC_choice_units = squeeze(nanmean(ACC_choice_units,1));
        ACC_pics_units = squeeze(nanmean(ACC_pics_units,1));
        
        [~, ~, ~, ~, choice_expvar, ~] = pca(ACC_choice_units);
        [~, ~, ~, ~, pics_expvar, ~] = pca(ACC_pics_units);
        f_ACC_PC95(1,1) = find(cumsum(pics_expvar)>95,1);
        f_ACC_PC95(1,2) = find(cumsum(choice_expvar)>95,1);
        f_ACC_PC95(1,3) = n_acc_units;
        f_ACC_PC95(1,4) = monkey_id;

        ACC_PC95 = [ACC_PC95 ; f_ACC_PC95];

    end
    
    if ~isempty(OFC_choice_units)
        
        n_ofc_units = numel(OFC_choice_units(1,1,:));
        
        OFC_choice_units = squeeze(nanmean(OFC_choice_units,1));
        OFC_pics_units = squeeze(nanmean(OFC_pics_units,1));
        
        [~, ~, ~, ~, choice_expvar, ~] = pca(OFC_choice_units);
        [~, ~, ~, ~, pics_expvar, ~] = pca(OFC_pics_units);
        f_OFC_PC95(1,1) = find(cumsum(pics_expvar)>95,1);
        f_OFC_PC95(1,2) = find(cumsum(choice_expvar)>95,1);
        f_OFC_PC95(1,3) = n_ofc_units;
        f_OFC_PC95(1,4) = monkey_id;

        OFC_PC95 = [OFC_PC95 ; f_OFC_PC95];

    end
 
end % of looping over files
%------------------------------------------------------------------------------

% get means and SEMs for n_PCs by brain region and task epoch for each animal
[chap_OFC_mean, chap_OFC_sem] = GetMeanCI(OFC_PC95(OFC_PC95(:,4)==0,1:2),'sem');
[chap_ACC_mean, chap_ACC_sem] = GetMeanCI(ACC_PC95(ACC_PC95(:,4)==0,1:2),'sem');

[george_OFC_mean, george_OFC_sem] = GetMeanCI(OFC_PC95(OFC_PC95(:,4)==1,1:2),'sem');
[george_ACC_mean, george_ACC_sem] = GetMeanCI(ACC_PC95(ACC_PC95(:,4)==1,1:2),'sem');

% let's do an anova for each animal
% create the factors
nPCs = [ OFC_PC95(:,1) ; ACC_PC95(:,1) ; OFC_PC95(:,2) ; ACC_PC95(:,2)];
brain_area = [ones(size(OFC_PC95(:,1))) ; ones(size(ACC_PC95(:,1)))*-1 ; 
              ones(size(OFC_PC95(:,1))) ; ones(size(ACC_PC95(:,1)))*-1];
          
alignment = [ones(size(OFC_PC95(:,1))) ; ones(size(ACC_PC95(:,1))) ;
             ones(size(OFC_PC95(:,1)))*-1 ; ones(size(ACC_PC95(:,1)))*-1];
         
animal_ix = [ OFC_PC95(:,4) ; ACC_PC95(:,4) ; OFC_PC95(:,4) ; ACC_PC95(:,4)];

[~,george_anova] = anovan(nPCs(animal_ix==1),{brain_area(animal_ix==1), alignment(animal_ix==1)},...
    'model','interaction','Display','off','varnames',{'brain_area','alignment'});

[~,chap_anova] = anovan(nPCs(animal_ix==0),{brain_area(animal_ix==0), alignment(animal_ix==0)},...
    'model','interaction','Display','off','varnames',{'brain_area','alignment'});


figure; 
subplot(3,2,1);
hold on
errorbar(chap_OFC_mean, chap_OFC_sem,'LineWidth',2,'CapSize',0);
errorbar(chap_ACC_mean, chap_ACC_sem,'LineWidth',2,'CapSize',0);
ylim([1,4.5]);
yticks([2,3,4]);
xticks([1,2]);
xlim([.8,2.2]);
xticklabels({'Pics','Choice'});
ylabel('#PCs to Explain 95% Var.');
xlabel('Alignment');
legend('OFC','ACC');
title('Animal C');
axis square

subplot(3,2,2);
hold on
errorbar(george_OFC_mean, george_OFC_sem,'LineWidth',2,'CapSize',0);
errorbar(george_ACC_mean, george_ACC_sem,'LineWidth',2,'CapSize',0);
xticks([1,2]);
xlim([.8,2.2]);
ylim([1,4.5]);
yticks([2,3,4]);
xticklabels({'Pics','Choice'});
title('Animal G');
axis square



subplot(3,2,3);
hold on
plot(OFC_PC95(OFC_PC95(:,4)==0,3), OFC_PC95(OFC_PC95(:,4)==0,1),'.','Markersize',20);
plot(ACC_PC95(ACC_PC95(:,4)==0,3), ACC_PC95(ACC_PC95(:,4)==0,1),'.','Markersize',20);
yticks([1:4])
ylim([1.5,4.5])
lsline
xlabel('Ensemble Size (num. units)');
ylabel('#PCs to Explain 95% Var.');
axis square

[chap_pics_OFC_r, chap_pics_OFC_p] = corr(OFC_PC95(OFC_PC95(:,4)==0,3), OFC_PC95(OFC_PC95(:,4)==0,1));
[chap_pics_ACC_r, chap_pics_ACC_p] = corr(ACC_PC95(OFC_PC95(:,4)==0,3), ACC_PC95(OFC_PC95(:,4)==0,1));



subplot(3,2,4);
hold on
plot(OFC_PC95(OFC_PC95(:,4)==1,3), OFC_PC95(OFC_PC95(:,4)==1,1),'.','Markersize',20);
plot(ACC_PC95(ACC_PC95(:,4)==1,3), ACC_PC95(ACC_PC95(:,4)==1,1),'.','Markersize',20);
yticks([1:4])
ylim([1.5,4.5])
lsline
axis square

[george_pics_OFC_r, george_pics_OFC_p] = corr(OFC_PC95(OFC_PC95(:,4)==1,3), OFC_PC95(OFC_PC95(:,4)==1,1));
[george_pics_ACC_r, george_pics_ACC_p] = corr(ACC_PC95(ACC_PC95(:,4)==1,3), ACC_PC95(ACC_PC95(:,4)==1,1));



subplot(3,2,5);
hold on
plot(OFC_PC95(OFC_PC95(:,4)==0,3), OFC_PC95(OFC_PC95(:,4)==0,2),'.','Markersize',20);
plot(ACC_PC95(ACC_PC95(:,4)==0,3), ACC_PC95(ACC_PC95(:,4)==0,2),'.','Markersize',20);
yticks([1:4])
ylim([1.5,4.5])
lsline
xlabel('Ensemble Size (num. units)');
ylabel('#PCs to Explain 95% Var.');
axis square

[chap_choice_OFC_r, chap_choice_OFC_p] = corr(OFC_PC95(OFC_PC95(:,4)==0,3), OFC_PC95(OFC_PC95(:,4)==0,2));
[chap_choice_ACC_r, chap_choice_ACC_p] = corr(ACC_PC95(OFC_PC95(:,4)==0,3), ACC_PC95(OFC_PC95(:,4)==0,2));


subplot(3,2,6);
hold on
plot(OFC_PC95(OFC_PC95(:,4)==1,3), OFC_PC95(OFC_PC95(:,4)==1,2),'.','Markersize',20);
plot(ACC_PC95(ACC_PC95(:,4)==1,3), ACC_PC95(ACC_PC95(:,4)==1,2),'.','Markersize',20);
yticks([1:4])
ylim([1.5,4.5])
lsline
axis square

[george_choice_OFC_r, george_choice_OFC_p] = corr(OFC_PC95(OFC_PC95(:,4)==1,3), OFC_PC95(OFC_PC95(:,4)==1,2))
[george_choice_ACC_r, george_choice_ACC_p] = corr(ACC_PC95(OFC_PC95(:,4)==1,3), ACC_PC95(OFC_PC95(:,4)==1,2))




















