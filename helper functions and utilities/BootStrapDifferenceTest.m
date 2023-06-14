function [pval] = BootStrapDifferenceTest(data1, data2, nboots)

for b = 1:nboots
    
    randD1_ix = randi(numel(data1),[numel(data1),1]);
    randD2_ix = randi(numel(data2),[numel(data2),1]);
    
    D1boot = mean(data1(randD1_ix(1:round(.85*numel(randD1_ix)))));
    D2boot = mean(data2(randD2_ix(1:round(.85*numel(randD2_ix)))));

   bootdiffs(b) = D1boot - D2boot;
    
    
end % of cycling through bootstraps

pval = sum(bootdiffs>=0) / nboots;


end % of function