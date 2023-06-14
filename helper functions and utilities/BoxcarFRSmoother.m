function [FRs,winTS] = BoxcarFRSmoother(rasters,window,step)

% rasters should be an array of ones and zeros with a sampling rate of
% milliseconds
winTS=[];
FRs=[];

ctr=0;
for i = step:step:numel(rasters)
ctr = ctr+1;
    
    if i-window/2 < 0
        wStart = i;
    else
        wStart = i - ceil(window/2);
    end
    
    if i+window/2 > numel(rasters)
        wEnd = numel(rasters);
    else
        wEnd = i + floor(window/2);
    end
    
FRs(ctr) = sum(rasters(wStart:wEnd))*ceil(1000/window);
winTS(ctr) = i;    
    
    
end % of looping through steps

end % of function