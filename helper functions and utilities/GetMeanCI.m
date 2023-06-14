function [Xmean,XCI] = GetMeanCI(X,method)
% gets column means and 95% CIs
Xmean=[];
XCI=[];

[nVals,nCols] = size(X);
     
Xmean = nanmean(X); 


if contains(method,'tdist')
% t dist method
SEM = nanstd(X)/sqrt(nVals);          % Standard Error
ts = tinv([0.025  0.975],nVals-1);    % T-Score
tCI = abs(Xmean + ts(2)*SEM);         % Confidence Interval
XCI = tCI; 
end

if contains(method,'percentile')
% percentile method
CIFcn = prctile(X,abs([0,100]-(100-95)/2));
pCI = abs(Xmean - CIFcn(2,:));
XCI = pCI; % use percent
end

if contains(method,'bootstrap')
% bootstrap method
try
bootCIs = bootci(1000,@nanmean,X);
bCI = nanmean([abs(Xmean - bootCIs(1,:)); abs(Xmean - bootCIs(2,:))],1);   
XCI = bCI;
catch
XCI = 0;
end
end

if contains(method,'sem')
% calculate sem 
XCI = (nanstd(X) / sqrt(nVals));
end


end % of function
