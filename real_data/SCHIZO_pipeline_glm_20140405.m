% Run the multiscale GLM-connectome pipeline (with group FDR) on the SCHIZO dataset.

clear all

addpath(genpath('/home/porban/quarantaine/niak-boss-0.12.2/'))

%%%%%%%%%%%%
%% Grabbing the results from BASC
%%%%%%%%%%%%
opt_g_basc.level = 'group';
opt_g_basc.flag_tseries = false;
files_in = niak_grab_stability_rest('/gs/scratch/porban/bascglm/cobre_basc_bascglm_20140405',opt_g_basc);

%%%%%%%%%%%%
%% Grabbing the results from the NIAK fMRI preprocessing pipeline
%%%%%%%%%%%%
opt_g.min_nb_vol = 60;      % irrelevant (see exclude list) % The minimum number of volumes for an fMRI dataset to be included. This option is useful when scrubbing is used, and the resulting time series may be too short.
opt_g.min_xcorr_func = 0.5; % The minimum xcorr score for an fMRI dataset to be included. This metric is a tool for quality control which assess the quality of non-linear coregistration of functional images in stereotaxic space. Manual inspection of the values during QC is necessary to properly set this threshold.
opt_g.min_xcorr_anat = 0.5; % The minimum xcorr score for an fMRI dataset to be included. This metric is a tool for quality control which assess the quality of non-linear coregistration of the anatomical image in stereotaxic space. Manual inspection of the values during QC is necessary to properly set this threshold.
%opt_g.exclude_subject = {'subject1','subject2'}; % If for whatever reason some subjects have to be excluded that were not caught by the quality control metrics, it is possible to manually specify their IDs here.
opt_g.type_files = 'glm_connectome'; % Specify to the grabber to prepare the files for the glm_connectome pipeline
%opt_g.filter.session = {'session1'}; % Just grab session 1
%opt_g.filter.run = {'rest'}; % Just grab the "rest" run
files_in.fmri = niak_grab_fmri_preprocess('/gs/scratch/porban/bascglm/cobre_fmri_preprocess_bascglm_20140405',opt_g).fmri; % Replace the folder by the path where the results of the fMRI preprocessing pipeline were stored. 

%%%%%%%%%%%%
%% Set the model
%%%%%%%%%%%%

%% Group
files_in.model.group = '/home/porban/database/bascglm/models/cobre_model_group.csv';


%%%%%%%%%%%%
%% Options 
%%%%%%%%%%%%
opt.folder_out = '/gs/scratch/porban/bascglm/results/cobre_glm_bascglm_20140405'; % Where to store the results
opt.fdr = 0.05; % The maximal false-discovery rate that is tolerated both for individual (single-seed) maps and whole-connectome discoveries, at each particular scale (multiple comparisons across scales are addressed via permutation testing)
%opt.type_fdr = 'group+family';
opt.fwe = 0.05; % The overall family-wise error, i.e. the probablity to have the observed number of discoveries, agregated across all scales, under the global null hypothesis of no association.
% opt.nb_samps = 1000; % The number of samples in the permutation test. This number has to be multiplied by OPT.NB_BATCH below to get the effective number of samples
% opt.nb_batch = 10; % The permutation tests are separated into NB_BATCH independent batches, which can run on parallel if sufficient computational resources are available
opt.nb_samps = 1000;
opt.nb_batch = 10;
opt.flag_rand = false; % if the flag is false, the pipeline is deterministic. Otherwise, the random number generator is initialized based on the clock for each job.


%%%%%%%%%%%%
%% Tests
%%%%%%%%%%%%

%% Group

% szVScont

opt.test.szVScont.group.contrast.sz                 = 1;

opt.test.szVScont_age.group.contrast.sz           = 1;
opt.test.szVScont_age.group.contrast.age          = 0;

opt.test.szVScont_sex.group.contrast.sz           = 1;
opt.test.szVScont_sex.group.contrast.sex          = 0;

opt.test.szVScont_FD.group.contrast.sz            = 1;
opt.test.szVScont_FD.group.contrast.FD            = 0;

opt.test.szVScont_age_sex_FD.group.contrast.sz    = 1;
opt.test.szVScont_age_sex_FD.group.contrast.age   = 0;
opt.test.szVScont_age_sex_FD.group.contrast.sex   = 0;
opt.test.szVScont_age_sex_FD.group.contrast.FD    = 0;




%% Avg

% sz

opt.test.sz.group.select(1).label                   = 'sz';
opt.test.sz.group.select(1).values                  = 1; 
opt.test.sz.group.contrast.intercept                = 1;

opt.test.sz_age.group.select(1).label              = 'sz';
opt.test.sz_age.group.select(1).values             = 1; 
opt.test.sz_age.group.contrast.intercept           = 1;
opt.test.sz_age.group.contrast.age                 = 0;

opt.test.sz_sex.group.select(1).label                    = 'sz';
opt.test.sz_sex.group.select(1).values                   = 1; 
opt.test.sz_sex.group.contrast.intercept           = 1;
opt.test.sz_sex.group.contrast.sex                 = 0;

opt.test.sz_FD.group.select(1).label                     = 'sz';
opt.test.sz_FD.group.select(1).values                    = 1; 
opt.test.sz_FD.group.contrast.intercept            = 1;
opt.test.sz_FD.group.contrast.FD                   = 0;

opt.test.sz_age_sex_FD.group.select(1).label             = 'sz';
opt.test.sz_age_sex_FD.group.select(1).values            = 1; 
opt.test.sz_age_sex_FD.group.contrast.intercept    = 1;
opt.test.sz_age_sex_FD.group.contrast.age          = 0;
opt.test.sz_age_sex_FD.group.contrast.sex          = 0;
opt.test.sz_age_sex_FD.group.contrast.FD           = 0;


% cont

opt.test.cont.group.select(1).label                   = 'sz';
opt.test.cont.group.select(1).values                  = 0; 
opt.test.cont.group.contrast.intercept          = 1;

opt.test.cont_age.group.select(1).label              = 'sz';
opt.test.cont_age.group.select(1).values             = 0; 
opt.test.cont_age.group.contrast.intercept           = 1;
opt.test.cont_age.group.contrast.age                 = 0;

opt.test.cont_sex.group.select(1).label                    = 'sz';
opt.test.cont_sex.group.select(1).values                   = 0; 
opt.test.cont_sex.group.contrast.intercept           = 1;
opt.test.cont_sex.group.contrast.sex                 = 0;

opt.test.cont_FD.group.select(1).label                     = 'sz';
opt.test.cont_FD.group.select(1).values                    = 0; 
opt.test.cont_FD.group.contrast.intercept            = 1;
opt.test.cont_FD.group.contrast.FD                   = 0;

opt.test.cont_age_sex_FD.group.select(1).label             = 'sz';
opt.test.cont_age_sex_FD.group.select(1).values            = 0; 
opt.test.cont_age_sex_FD.group.contrast.intercept    = 1;
opt.test.cont_age_sex_FD.group.contrast.age          = 0;
opt.test.cont_age_sex_FD.group.contrast.sex          = 0;
opt.test.cont_age_sex_FD.group.contrast.FD           = 0;


%%%%%%%%%%%%
%% Run the pipeline
%%%%%%%%%%%%
opt.flag_test = false; % Put this flag to true to just generate the pipeline without running it. Otherwise the region growing will start. 
%opt.psom.max_queued = 24; % Uncomment and change this parameter to set the number of parallel threads used to run the pipeline
[pipeline,opt] = niak_pipeline_glm_connectome(files_in,opt);

