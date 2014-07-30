clear all

addpath(genpath('/home/porban/quarantaine/niak-boss-0.12.4/'))

%%%%%%%%%%%%
%% Grabbing the results from BASC
%%%%%%%%%%%%
opt_g_basc.level = 'group';
opt_g_basc.flag_tseries = false;
files_in = niak_grab_stability_rest('/gs/scratch/porban/bascglm/alloego_basc_bascglm_20140405/',opt_g_basc);

%%%%%%%%%%%%
%% Grabbing the results from the NIAK fMRI preprocessing pipeline
%%%%%%%%%%%%
opt_g.min_nb_vol = 60;     % The minimum number of volumes for an fMRI dataset to be included. This option is useful when scrubbing is used, and the resulting time series may be too short.
opt_g.min_xcorr_func = 0.5; % The minimum xcorr score for an fMRI dataset to be included. This metric is a tool for quality control which assess the quality of non-linear coregistration of functional images in stereotaxic space. Manual inspection of the values during QC is necessary to properly set this threshold.
opt_g.min_xcorr_anat = 0.5; % The minimum xcorr score for an fMRI dataset to be included. This metric is a tool for quality control which assess the quality of non-linear coregistration of the anatomical image in stereotaxic space. Manual inspection of the values during QC is necessary to properly set this threshold.
%opt_g.exclude_subject = {'subject1','subject2'}; % If for whatever reason some subjects have to be excluded that were not caught by the quality control metrics, it is possible to manually specify their IDs here.
opt_g.type_files = 'glm_connectome'; % Specify to the grabber to prepare the files for the glm_connectome pipeline
%opt_g.filter.session = {'session1'}; % Just grab session 1
%opt_g.filter.run = {'rest'}; % Just grab the "rest" run
files_in.fmri = niak_grab_fmri_preprocess('/gs/scratch/porban/bascglm/alloego_fmri_preprocess_bascglm_20140405',opt_g).fmri; % Replace the folder by the path where the results of the fMRI preprocessing pipeline were stored. 

%%%%%%%%%%%%
%% Set the model
%%%%%%%%%%%%

%% Group
files_in.model.group = '/home/porban/database/bascglm/models/alloego_model_group.csv';

%% Subjects

data.group           = {'alloego'};
data.subs_alloego    = {'01','03','09','11','15','17','20','21','24','25','28','30','31','02','04','06','07','10','12','13','18','19','22','23','26','27',...
                        '33','35','36','39','43','44','47','48','51','52','55','56','59','61','32','34','37','38','41','42','49','50','53','54','57','58','60','62'}; 
data.runs            = {'session1_run1','session1_run2'};


for n1 = 1:length(data.subs_alloego)
    data.covariates_group_subs{n1} = strcat(data.group{1},data.subs_alloego{n1});
end




for n = 1:length(data.covariates_group_subs)
    
    files_in.model.individual.(data.covariates_group_subs{n}).inter_run = '/home/porban/database/bascglm/models/alloego_model_inter_run.csv';
    
end


%%%%%%%%%%%%
%% Options 
%%%%%%%%%%%%
opt.folder_out = '/gs/scratch/porban/bascglm/results/alloego_glm_global_bascglm_20140518'; % Where to store the results
opt.fdr = 0.05; % The maximal false-discovery rate that is tolerated both for individual (single-seed) maps and whole-connectome discoveries, at each particular scale (multiple comparisons across scales are addressed via permutation testing)
opt.type_fdr = 'global';
opt.fwe = 0.05; % The overall family-wise error, i.e. the probablity to have the observed number of discoveries, agregated across all scales, under the global null hypothesis of no association.
% opt.nb_samps = 1000; % The number of samples in the permutation test. This number has to be multiplied by OPT.NB_BATCH below to get the effective number of samples
% opt.nb_batch = 10; % The permutation tests are separated into NB_BATCH independent batches, which can run on parallel if sufficient computational resources are available
opt.nb_samps = 1000;
opt.nb_batch = 10;
opt.flag_rand = false; % if the flag is false, the pipeline is deterministic. Otherwise, the random number generator is initialized based on the clock for each job.


%%%%%%%%%%%%
%% Tests
%%%%%%%%%%%%



%% AVG RESTS 1 & 2

opt.test.rest1_avg.group.contrast.intercept      = 1;
opt.test.rest1_avg.group.contrast.age           = 0;
opt.test.rest1_avg.group.contrast.sex           = 0;
opt.test.rest1_avg.group.contrast.FDrun1           = 0;
opt.test.rest1_avg.inter_run.select(1).label     = 'run1';
opt.test.rest1_avg.inter_run.select(1).values    = 1;

opt.test.rest2_avg.group.contrast.intercept      = 1;
opt.test.rest2_avg.group.contrast.age           = 0;
opt.test.rest2_avg.group.contrast.sex           = 0;
opt.test.rest2_avg.group.contrast.FDrun2           = 0;
opt.test.rest2_avg.inter_run.select(1).label     = 'run2';
opt.test.rest2_avg.inter_run.select(1).values    = 1;




%% REST 2-1

opt.test.rest2VS1_avg.group.contrast.intercept      = 1;
opt.test.rest2VS1_avg.group.contrast.age           = 0;
opt.test.rest2VS1_avg.group.contrast.sex           = 0;
opt.test.rest2VS1_avg.group.contrast.FDrun1           = 0;
opt.test.rest2VS1_avg.group.contrast.FDrun2           = 0;
opt.test.rest2VS1_avg.inter_run.select(1).label     = 'run2vs1';
opt.test.rest2VS1_avg.inter_run.select(1).values    = [-1 1];
opt.test.rest2VS1_avg.inter_run.contrast.run2vs1      = 1;







%%%%%%%%%%%%
%% Run the pipeline
%%%%%%%%%%%%
opt.flag_test = false; % Put this flag to true to just generate the pipeline without running it. Otherwise the region growing will start. 
%opt.psom.max_queued = 24; % Uncomment and change this parameter to set the number of parallel threads used to run the pipeline
[pipeline,opt] = niak_pipeline_glm_connectome(files_in,opt);

