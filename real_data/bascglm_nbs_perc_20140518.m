clear all

path = '/Users/pyeror/Work/transfert/Bascglm/results/nbs/';


study = {'cobre','blind','alloego'};
scale = {'sci10_scg10_scf10','sci100_scg100_scf100'};
contrast = {'szVScont_age_sex_FD','CBvsSC_conf','rest2VS1_avg'};
level = {'0.01','0.2'};


for i = 1:length(study)
    for j = 1:length(scale)
        for k = 1:length(level)
        
        [hdr,pos] = niak_read_vol(strcat(path,study{i},'_',contrast{i},'_',scale{j},'_perc_disc_nbs_',level{k},'_pos.nii.gz'));
        [hdr,neg] = niak_read_vol(strcat(path,study{i},'_',contrast{i},'_',scale{j},'_perc_disc_nbs_',level{k},'_neg.nii.gz'));
        new = zeros(size(pos));
        new = pos + neg;
        hdr.file_name = strcat(path,study{i},'_',contrast{i},'_',scale{j},'_perc_disc_nbs_',level{k},'.nii.gz');
        niak_write_vol(hdr,new);
        
        end
    end
end

