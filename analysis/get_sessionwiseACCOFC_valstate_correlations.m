function [ch_rval, unch_rval, corr_animal_id, t_mids] = get_sessionwiseACCOFC_valstate_correlations


load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\ACC_valdec_output.mat');

ACCval = valdata_pics;
t_mids = valdata_pics.t_mids;

load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\OFC_valdec_output.mat');
OFCval = valdata_pics;



[ACC_mean, ACC_sem] = GetMeanCI(ACCval.ch_ppd,'sem');
[OFC_mean, OFC_sem] = GetMeanCI(OFCval.ch_ppd,'sem');

OFC_s_ids = unique(OFCval.session);
ACC_s_ids = unique(ACCval.session);

OFC_val_pps = [];
ACC_val_pps = [];


simultaneous_recs = intersect(OFC_s_ids, ACC_s_ids);
for s = 1:numel(simultaneous_recs)
    
    s_OFC_chval_pps = OFCval.ch_ppd(contains(OFCval.session,simultaneous_recs{s}),:);
    s_ACC_chval_pps = ACCval.ch_ppd(contains(ACCval.session,simultaneous_recs{s}),:);
    s_OFC_unchval_pps = OFCval.unch_ppd(contains(OFCval.session,simultaneous_recs{s}),:);
    s_ACC_unchval_pps = ACCval.unch_ppd(contains(ACCval.session,simultaneous_recs{s}),:);
    
    corr_animal_id(s) = contains(simultaneous_recs(s),'George');
    

    
    
    for t = 1:numel(ACC_mean)
        
        [ch_rval(s,t), ch_pval(s,t)] = corr(s_OFC_chval_pps(:,t), s_ACC_chval_pps(:,t),'rows','complete');
        [unch_rval(s,t), unch_pval(s,t)] = corr(s_OFC_unchval_pps(:,t), s_ACC_unchval_pps(:,t),'rows','complete');

        
    end % of looping over times
    
end % of looping over sessions

% Fisher-transform the correlation coefficients
ch_rval = atanh(ch_rval);
unch_rval = atanh(unch_rval);

end % of function
