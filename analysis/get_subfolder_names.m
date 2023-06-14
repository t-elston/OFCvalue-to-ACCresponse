function [folder_names] = get_subfolder_names(root_dir)

d = dir(root_dir);
isub = [d(:).isdir]; %# returns logical vector
folder_names = {d(isub).name}';

folder_names(ismember(folder_names,{'.','..'})) = [];

end