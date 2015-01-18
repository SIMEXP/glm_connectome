%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setting input/output files %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

addpath(genpath('/home/porban/quarantaine/niak-boss-0.12.2/'))

opt_g.min_nb_vol = 60; % irrelevant, as subjects are in the 'exclude' list if less than 60 vols are left after scrubbing...
root_path = '/gs/scratch/porban/bascglm/';

opt_g.min_xcorr_func = 0.5;
opt_g.min_xcorr_anat = 0.5;
files_in = niak_grab_fmri_preprocess([root_path 'cobre_fmri_preprocess_bascglm_20140405/'],opt_g);

%%%%%%%%%%%%%
%% Options %%
%%%%%%%%%%%%%
opt.folder_out = [root_path 'cobre_region_growing_bascglm_20140405/']; % Where to store the results
opt.flag_roi = true; % Only generate the ROI parcelation
opt.region_growing.thre_size = 1000; % The critical size for regions

%% Run the pipeline
opt.flag_test = false;
[pipeline_rg,opt] = niak_pipeline_stability_rest(files_in,opt);
