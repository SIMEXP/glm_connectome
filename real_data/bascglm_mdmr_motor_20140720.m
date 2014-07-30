clear all
addpath(genpath('/sb/project/gsf-624-aa/quarantaine/niak-boss-0.12.14/'));
addpath(genpath('/home/pbellec/svn/projects/glm_connectome/'));

path_data      = '/gs/scratch/porban/bascglm/';
path_glm       = '/gs/scratch/porban/bascglm/results/';
path_confounds = '/sb/project/gsf-624-aa/database/bascglm/models/';
path_output    = '/gs/scratch/pbellec/bascglm/results/mdmr/';
path_model     = '/sb/project/gsf-624-aa/database/bascglm/models/';

study = {'alloego'}; 
contrast = {'rest2VS1_avg'};
scale = {'sci10_scg10_scf10','sci100_scg100_scf100'};


opt_g.min_nb_vol = 60;     % The minimum number of volumes for an fMRI dataset to be included. This option is useful when scrubbing is used, and the resulting time series may be too short.
opt_g.min_xcorr_func = 0.5; % The minimum xcorr score for an fMRI dataset to be included. This metric is a tool for quality control which assess the quality of non-linear coregistration of functional images in stereotaxic space. Manual inspection of the values during QC is necessary to properly set this threshold.
opt_g.min_xcorr_anat = 0.5; % The minimum xcorr score for an fMRI dataset to be included. This metric is a tool for quality control which assess the quality of non-linear coregistration of the anatomical image in stereotaxic space. Manual inspection of the values during QC is necessary to properly set this threshold.
opt_g.type_files = 'glm_connectome'; % Specify to the grabber to prepare the files for the glm_connectome pipeline
        
pipe = struct();
list_p = [0.05, 0.2];
for pp = 1:length(list_p)
label_p = sprintf('%i',100*list_p(pp));
for i = 1:length(study) % i is also valid for contrast    
    files_in.fmri = niak_grab_fmri_preprocess([path_data study{i} '_fmri_preprocess_bascglm_20140405/'],opt_g).fmri; % Replace the folder by the path where the results of the fMRI preprocessing pipeline were stored. 
    for j = 1:length(scale)
        
        files_in.glm = [path_glm study{i} '_glm_bascglm_20140405/' scale{j} '/' contrast{i},'/glm_',contrast{i},'_',scale{j},'.mat']; 
        files_in.network = [path_glm study{i} '_glm_bascglm_20140405/' scale{j} '/networks_' scale{j} '.mnc.gz'];
        files_in.confounds = [path_model 'alloego_mdmr_confounds.csv'];
        files_out.perc_discovery = [path_output study{i} '_' contrast{i} '_' scale{j} '_perc_disc_mdmr_p' label_p '.mnc.gz'];
        files_out.mdmr = [path_output study{i} '_' contrast{i} '_' scale{j} '_mdmr_p' label_p '.mnc.gz'];
        
        opt.nb_permutation = 10000;
        opt.p_seed = list_p(pp);
        opt.fdr = 0.05;
        opt.session1 = 'session1';
        opt.session2 = 'session1';
        opt.run1 = 'run1';
        opt.run2 = 'run2';
        pipe = psom_add_job(pipe,['mdmr_' study{i} '_' scale{j} '_' label_p],'glmc_brick_glm2mdmr_diff',files_in,files_out,opt);
    end
end
end
opt_p.path_logs = [path_output 'logs'];
psom_run_pipeline(pipe,opt_p);
