function [outdata] = load_decoder_output_v01(decoder_dir, brain_area, alignment)


units = []; % unit data: n_trials x n_times x n_units
u_hemi_ix = []; % which hemi are the units in? n_units x 1; -1 = left hemi; 1 = right hemi
ts = []; % timesteps - when, relative to the aligment

folder_contents = dir([decoder_dir, current_file{1}]);
% look for the non-directory contents
sub_directory_ix = [folder_contents(:).isdir]; % returns logical vector
file_names = {folder_contents(~sub_directory_ix).name}';

% now find the files with the proper alignment
aligned_ix = contains(file_names,alignment);

% now find files with the proper brain area
brain_area_ix = contains(file_names,brain_area);

% use these indices to get the names of the relevant file(s)
relevant_files = file_names(brain_area_ix & aligned_ix);


end % of function