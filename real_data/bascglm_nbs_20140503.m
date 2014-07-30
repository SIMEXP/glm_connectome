clear all
addpath(genpath('/usr/local/quarantine/niak-boss-0.12.4/'));
addpath(genpath('/usr/local/svn/projects/'));
addpath(genpath('/usr/local/svn/psom/'));


path = '/media/database3/bascglm/results/';
path_output = '/media/database3/bascglm/results/nbs/new/';

study = {'cobre','blind','alloego'};
scale = {'sci10_scg10_scf10','sci100_scg100_scf100'};
contrast = {'szVScont_age_sex_FD','CBvsSC_conf','rest2VS1_avg'};


for i = 1:length(study)
    for j = 1:length(scale)
        
        in.glm = strcat(path,study{i},'_glm_bascglm_20140405/',scale{j},'/',contrast{i},'/glm_',contrast{i},'_',scale{j},'.mat');
        in.network = strcat(path,study{i},'_glm_bascglm_20140405/',scale{j},'/networks_',scale{j},'.nii.gz');
        out.perc_discovery = strcat(path_output,study{i},'_',contrast{i},'_',scale{j},'_perc_disc_nbs_0.2_pos.nii.gz');
        out.nbs = strcat(path_output,study{i},'_',contrast{i},'_',scale{j},'_nbs_0.2_pos.nii.gz'); % the equivalent of fdr_xxx but with NBS
        
        opt.direction = 'positive';
        opt.alpha = 0.2; % the cluster-level significance
%         opt.thresh = 2; % the uncorrected threshold on connexel-wise t-test

        glmc_brick_glm2nbs(in,out,opt);
    end
end

for i = 1:length(study)
    for j = 1:length(scale)
        
        in.glm = strcat(path,study{i},'_glm_bascglm_20140405/',scale{j},'/',contrast{i},'/glm_',contrast{i},'_',scale{j},'.mat');
        in.network = strcat(path,study{i},'_glm_bascglm_20140405/',scale{j},'/networks_',scale{j},'.nii.gz');
        out.perc_discovery = strcat(path_output,study{i},'_',contrast{i},'_',scale{j},'_perc_disc_nbs_0.2_neg.nii.gz');
        out.nbs = strcat(path_output,study{i},'_',contrast{i},'_',scale{j},'_nbs_0.2_neg.nii.gz'); % the equivalent of fdr_xxx but with NBS
        
        opt.direction = 'negative';
        opt.alpha = 0.2; % the cluster-level significance
%         opt.thresh = 2; % the uncorrected threshold on connexel-wise t-test

        glmc_brick_glm2nbs(in,out,opt);
    end
end


% for i = 1:length(study)
%     for j = 1:length(scale)
%         
%         in.glm = strcat(path,study{i},'_glm_bascglm_20140405/',scale{j},'/',contrast{i},'/glm_',contrast{i},'_',scale{j},'.mat');
%         in.network = strcat(path,study{i},'_glm_bascglm_20140405/',scale{j},'/networks_',scale{j},'.nii.gz');
%         out.perc_discovery = strcat(path_output,study{i},'_',contrast{i},'_',scale{j},'_perc_disc_nbs_0.1_pos.nii.gz');
%         out.nbs = strcat(path_output,study{i},'_',contrast{i},'_',scale{j},'_nbs_0.1_pos.nii.gz'); % the equivalent of fdr_xxx but with NBS
%         
%         opt.direction = 'positive';
%         opt.alpha = 0.1; % the cluster-level significance
% %         opt.thresh = 2; % the uncorrected threshold on connexel-wise t-test
% 
%         glmc_brick_glm2nbs(in,out,opt);
%     end
% end
% 
% for i = 1:length(study)
%     for j = 1:length(scale)
%         
%         in.glm = strcat(path,study{i},'_glm_bascglm_20140405/',scale{j},'/',contrast{i},'/glm_',contrast{i},'_',scale{j},'.mat');
%         in.network = strcat(path,study{i},'_glm_bascglm_20140405/',scale{j},'/networks_',scale{j},'.nii.gz');
%         out.perc_discovery = strcat(path_output,study{i},'_',contrast{i},'_',scale{j},'_perc_disc_nbs_0.1_neg.nii.gz');
%         out.nbs = strcat(path_output,study{i},'_',contrast{i},'_',scale{j},'_nbs_0.1_neg.nii.gz'); % the equivalent of fdr_xxx but with NBS
%         
%         opt.direction = 'negative';
%         opt.alpha = 0.1; % the cluster-level significance
% %         opt.thresh = 2; % the uncorrected threshold on connexel-wise t-test
% 
%         glmc_brick_glm2nbs(in,out,opt);
%     end
% end