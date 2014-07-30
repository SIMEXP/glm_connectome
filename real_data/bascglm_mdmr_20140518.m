clear all
addpath(genpath('/usr/local/quarantine/niak-boss-0.12.4/'));
addpath(genpath('/usr/local/svn/projects/'));
addpath(genpath('/usr/local/svn/psom/'));


path = '/media/database3/bascglm/';
path_output = '/media/database3/bascglm/results_peu/mdmr/';

study = {'cobre','blind'}; % pas d'intrasujet, 'alloego'
contrast = {'szVScont_age_sex_FD','CBvsSC_conf'}; % pas d'intrasujet 'rest2VS1_avg'
scale = {'sci10_scg10_scf10','sci100_scg100_scf100'};


opt_g.min_nb_vol = 60;     % The minimum number of volumes for an fMRI dataset to be included. This option is useful when scrubbing is used, and the resulting time series may be too short.
opt_g.min_xcorr_func = 0.5; % The minimum xcorr score for an fMRI dataset to be included. This metric is a tool for quality control which assess the quality of non-linear coregistration of functional images in stereotaxic space. Manual inspection of the values during QC is necessary to properly set this threshold.
opt_g.min_xcorr_anat = 0.5; % The minimum xcorr score for an fMRI dataset to be included. This metric is a tool for quality control which assess the quality of non-linear coregistration of the anatomical image in stereotaxic space. Manual inspection of the values during QC is necessary to properly set this threshold.
opt_g.type_files = 'glm_connectome'; % Specify to the grabber to prepare the files for the glm_connectome pipeline
        

for i = 1:length(study) % i is also valid for contrast
    
files_in.fmri = niak_grab_fmri_preprocess([path study{i} '_fmri_preprocess_bascglm_20140405/'],opt_g).fmri; % Replace the folder by the path where the results of the fMRI preprocessing pipeline were stored. 
    
    for j = 1:length(scale)
        
        files_in.glm = [path 'results_peu/' study{i} '_glm_bascglm_20140405/' scale{j} '/' contrast{i},'/glm_',contrast{i},'_',scale{j},'.mat']; 
        files_in.network = [path 'results_peu/' study{i} '_glm_bascglm_20140405/' scale{j} '/networks_' scale{j} '.nii.gz'];
        files_out.perc_discovery = [path_output study{i} '_' contrast{i} '_' scale{j} '_perc_disc_mdmr_p0.05.nii.gz'];
        files_out.mdmr = [path_output study{i} '_' contrast{i} '_' scale{j} '_mdmr_p0.05.nii.gz'];
        
        opt.nb_permutation = 10000;
        opt.p_seed = 0.05;
        opt.fdr = 0.05;
        glmc_brick_glm2mdmr(files_in,files_out,opt);
        
    end
end

