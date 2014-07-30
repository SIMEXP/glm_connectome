
%% syneeg extraction findings (2014/02/27)

clear all

%% input

path_results =  '/Users/pyeror/Work/transfert/Bascglm/results/cobre_glm_bascglm_20140405/';
path_scale =    {'sci100_scg100_scf100'};
%path_scale =    {'sci10_scg10_scf10'};
path_contrast = {'szVScont_age_sex_FD','sz_age_sex_FD','cont_age_sex_FD'};
path_overlap = {'szVScont_age_sex_FD','sz_age_sex_FD','cont_age_sex_FD'};
data_contrast = {'szVScont_age_sex_FD','sz_age_sex_FD','cont_age_sex_FD'};
data_overlap = {'szVScont_age_sex_FD','sz_age_sex_FD','cont_age_sex_FD'};

data_seed     = {'29_thal','73_pos_php','81_ant_php_amyg','38_temp_pole'};
data_cluster =   [29 73 81 38];
% data_seed     = {'1_bg','3_Temp','6_pfc_hp','10_vis_hp'};
% data_cluster =   [1 3 6 10];

data_type = {'fdr','ttest','effect'};



%% map de clusters séparés avec valeur à 1 et d'un clusters compound avec valeurs 1, 2, 3, etc (à partir de networks)

for i = 1:length(path_scale)
    for k = 1:length(data_cluster)
        
[hdr,mask] = niak_read_vol(strcat(path_results,'/',path_scale{i},'/networks/networks_',path_scale{i},'.nii.gz'));
submask = zeros(size(mask));
cluster = data_cluster(k);
submask(mask==cluster) = 1;
hdr.file_name = strcat(path_results,'/',path_scale{i},'/networks/cluster_',data_seed{k},'_networks_',path_scale{i},'.nii.gz');
niak_write_vol(hdr,submask);

    end
end

[hdr,mask] = niak_read_vol(strcat(path_results,'/',path_scale{1},'/networks/networks_',path_scale{1},'.nii.gz'));
submask = zeros(size(mask));
for k = 1:length(data_cluster)
    cluster = data_cluster(k);
    submask(mask==cluster) = k;
end
hdr.file_name = strcat(path_results,'/',path_scale{1},'/networks/cluster_composite_networks_',path_scale{1},'.nii.gz');
niak_write_vol(hdr,submask);




%% fdr threshold - write .csv with minimal t values surviving FDR correction

for i = 1:length(path_scale)
    for j = 1:length(path_contrast)
        for k = 1:length(data_cluster)

[hdr,fdr] = niak_read_vol(strcat(path_results,'/',path_scale{i},'/',path_contrast{j},'/fdr_',path_contrast{j},'_',path_scale{i},'.nii.gz'));
cluster = data_cluster(k);
fdr_mask = abs(fdr(:,:,:,cluster));
fdr_ok_min = fdr_mask;
fdr_ok_max = fdr_mask;
fdr_ok_min(fdr_mask==0) = 10;
fdr_ok_max(fdr_mask==0) = 0;
output_fdr_min = min(fdr_ok_min(:));    % probleme, ne fonctionne pas pour min!!!
output_fdr_max = max(fdr_ok_max(:));
data_values_min(k,j) = output_fdr_min;
data_values_max(k,j) = output_fdr_max;

        end
    end
end

opt.labels_y = data_contrast;
opt.labels_x = data_seed;
opt.precision = 2;
niak_write_csv(strcat(path_results,'/',path_scale{i},'/fdr_values_min.csv'),data_values_min,opt)
niak_write_csv(strcat(path_results,'/',path_scale{i},'/fdr_values_max.csv'),data_values_max,opt)



%% single FDR, TTEST or EFFECT maps (for use with mixed colormap)

for t = 1:length(data_type)

for i = 1:length(path_scale)
    for j = length(path_scale)
        for k = 1:length(data_cluster)
        
[hdr,con] = niak_read_vol(strcat(path_results,'/',path_scale{i},'/',path_contrast{j},'/',data_type{t},'_',path_overlap{j},'_',path_scale{i},'.nii.gz'));
cluster = data_cluster(k);
con_mask = con(:,:,:,cluster);
con_new = zeros(size(con_mask));

if t < 3
con_mask(con_mask<=-4) = -4;
con_mask(con_mask>=4) = 4;
con_new(con_mask<-1) = abs(con_mask(con_mask<-1));
con_new(con_mask>1) = con_mask(con_mask>1) + 4.001;
else
con_mask(con_mask<=-0.2) = -0.2;
con_mask(con_mask>=0.2) = 0.2;
con_new(con_mask<0) = abs(con_mask(con_mask<0));
con_new(con_mask>0) = con_mask(con_mask>0) + 0.201;     
end

hdr.file_name = strcat(path_results,'/',path_scale{i},'/mixed/mixed_',data_type{t},'_',path_contrast{j},'_',data_seed{k},'_',path_scale{i},'.nii.gz');
niak_write_vol(hdr,con_new);

        end
    end
end

end
