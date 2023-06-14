function [all_bhv] = load_all_bhv_v01(bhv_dir)

% find names of bhv files and rec files
bhv_folder_names = get_subfolder_names(bhv_dir);

n_files = numel(bhv_folder_names);

all_bhv = table;
% now loop over each file
for f_ix = 1:n_files
    
    current_file = bhv_folder_names(f_ix);
    monkey_id = contains(current_file,'George'); % 0 = Cap; 1 = George
    
    fprintf('\n') 
    disp(['***' current_file{1}]);
    
    % pull out this file's behavior
    [trialinfo] = load_behavior(bhv_dir, current_file);
    
    % find out whether he picked optimally with respect to expected value
    l_subjval = trialinfo.subjval_expval(:,1); l_subjval(isnan(l_subjval))=0;
    r_subjval = trialinfo.subjval_expval(:,2); r_subjval(isnan(r_subjval))=0;

    
    best_side = 2*((r_subjval - l_subjval)>0)-1;
    
    
    % unpack these data into columns of a table
    f_bhv = table;
    f_bhv.file_num = ones(size(trialinfo.trialtype))*f_ix;
    f_bhv.monkey = ones(size(trialinfo.trialtype))*monkey_id;
    f_bhv.trialtype = trialinfo.trialtype; % 2 = free choices
    f_bhv.picked_best = double(best_side == trialinfo.lever);
    f_bhv.best_side = best_side;
    f_bhv.lever = trialinfo.lever;
    f_bhv.rt = trialinfo.rt;
    f_bhv.left_expval = trialinfo.subjval_expval(:,1);
    f_bhv.right_expval = trialinfo.subjval_expval(:,2);
    f_bhv.left_amnt = trialinfo.amnt(:,1);
    f_bhv.right_amnt = trialinfo.amnt(:,2);
    f_bhv.left_prob = trialinfo.prob(:,1);
    f_bhv.right_prob = trialinfo.prob(:,2);
    f_bhv.saccN = trialinfo.saccN;
    f_bhv.first_sacc_loc = trialinfo.saccloc(:,1);
    f_bhv.sacc1_rt = trialinfo.MLsacc(:,1);
    f_bhv.sacc2_rt = trialinfo.MLsacc(:,2);

    
    % now concatenate these data into a bigger table
    all_bhv = [all_bhv ; f_bhv];

end % of loading files












end % of function 