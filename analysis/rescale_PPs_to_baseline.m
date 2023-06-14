function [rescaled_PPs] = rescale_PPs_to_baseline(PPs, PP_baselines)

[n_classes, n_times] = size(PP_baselines);

rescaled_PPs = NaN(size(PPs));

for c = 1:n_classes
    
    for t = 1:n_times
        
        rescaled_PPs(:,t,c) = (PPs(:,t,c) / PP_baselines(c,t))*100;
        
    end % of looping over times
    
end % of looping over classes




end % of function