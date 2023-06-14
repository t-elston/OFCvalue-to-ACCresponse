function [units, u_hemi_ix, ts] = load_neurons(rec_dir, current_file, alignment, brain_area, offsets)

units = []; % unit data: n_trials x n_times x n_units
u_hemi_ix = []; % which hemi are the units in? n_units x 1; -1 = left hemi; 1 = right hemi
ts = []; % timesteps - when, relative to the aligment

folder_contents = dir([rec_dir, current_file{1}]);
% look for the non-directory contents
sub_directory_ix = [folder_contents(:).isdir]; % returns logical vector
file_names = {folder_contents(~sub_directory_ix).name}';

% now find the files with the proper alignment
aligned_ix = contains(file_names,alignment);

% now find files with the proper brain area
brain_area_ix = contains(file_names,brain_area);

% use these indices to get the names of the relevant file(s)
relevant_files = file_names(brain_area_ix & aligned_ix);

% initialize some empty arrays for the data from each hemisphere (if
% applicable)
left_hemi_data = [];
left_hemi_ix = [];
right_hemi_data = [];
right_hemi_ix = [];

% now load the neurons

% % do we have left hemisphere data?
if any(contains(relevant_files,'L'))
    left_hemi_data = cell2mat(struct2cell(load([rec_dir, current_file{1},'/',relevant_files{contains(relevant_files,'L')}], 'SPKfrnorm_units')));
    ts = cell2mat(struct2cell(load([rec_dir, current_file{1},'/',relevant_files{contains(relevant_files,'L')}], 't_mids')));
    left_hemi_ix = ones(numel(left_hemi_data(1,1,:)),1)*-1;
    
end

%%% do we have any right hemisphere data?
if any(contains(relevant_files,'R'))
    right_hemi_data = cell2mat(struct2cell(load([rec_dir, current_file{1},'/',relevant_files{contains(relevant_files,'R')}], 'SPKfrnorm_units')));
    right_hemi_ix = ones(numel(right_hemi_data(1,1,:)),1);
end

% do we have left hemisphere data?
% if any(contains(relevant_files,'L'))
%     left_hemi_data = cell2mat(struct2cell(load([rec_dir, current_file{1},'/',relevant_files{contains(relevant_files,'L')}], 'SPKfr_units')));
%     ts = cell2mat(struct2cell(load([rec_dir, current_file{1},'/',relevant_files{contains(relevant_files,'L')}], 't_mids')));
%     left_hemi_ix = ones(numel(left_hemi_data(1,1,:)),1)*-1;
%     
% end
% 
% % do we have any right hemisphere data?
% if any(contains(relevant_files,'R'))
%     right_hemi_data = cell2mat(struct2cell(load([rec_dir, current_file{1},'/',relevant_files{contains(relevant_files,'R')}], 'SPKfr_units')));
%     right_hemi_ix = ones(numel(right_hemi_data(1,1,:)),1);
% end

ts = ts+1;

units = cat(3, left_hemi_data, right_hemi_data);
u_hemi_ix = [left_hemi_ix ; right_hemi_ix];

% only save the firing rates relative to the offsets specified in the offsets input
if ~isnan(offsets)
    [~,save_start] = min(abs(ts - offsets(1)));
    [~,save_end] = min(abs(ts - offsets(2)));
else
    save_start = 1;
    save_end = numel(ts);
end

ts = ts(save_start:save_end);
units=units(:,save_start:save_end,:);


end % of function 