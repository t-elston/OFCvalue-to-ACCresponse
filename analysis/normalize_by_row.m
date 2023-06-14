function [ipsi_outdata, contra_outdata] = normalize_by_row(ipsi_indata, contra_indata)

ipsi_outdata = NaN(size(ipsi_indata));
contra_outdata = NaN(size(contra_indata));

for r = 1:numel(ipsi_outdata(:,1))
    
%     r_min = min([ipsi_indata(r,:), contra_indata(r,:)]);
%     r_max = max([ipsi_indata(r,:), contra_indata(r,:)]);
%     
%     ipsi_outdata(r,:) = (ipsi_indata(r,:) - r_min) / (r_max - r_min);
%     contra_outdata(r,:) = (contra_indata(r,:) - r_min) / (r_max - r_min);
% 
    r_mean = nanmean([ipsi_indata(r,:), contra_indata(r,:)]);
    r_std = nanstd([ipsi_indata(r,:), contra_indata(r,:)]);
    
    ipsi_outdata(r,:) = (ipsi_indata(r,:) - r_mean) / (r_std);
    contra_outdata(r,:) = (contra_indata(r,:) - r_mean) / (r_std);
        
%     i_r_mean = nanmean(ipsi_indata(r,:));
%     i_r_std = nanstd(ipsi_indata(r,:));
%     
%     c_r_mean = nanmean(contra_indata(r,:));
%     c_r_std = nanstd(contra_indata(r,:));
%     
%     ipsi_outdata(r,:) = (ipsi_indata(r,:) - i_r_mean) / (i_r_std);
%     contra_outdata(r,:) = (contra_indata(r,:) - c_r_mean) / (c_r_std);
%       
      
%       
end % of looping over rows

end % of function