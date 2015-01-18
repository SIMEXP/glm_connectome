clear all

addpath(genpath('/home/porban/quarantaine/niak-boss-0.12.18/'))

root_path = '/gs/scratch/porban/bascglm/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setting input/output files %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

files_in = niak_grab_region_growing([root_path 'blind_region_growing_bascglm_20140405/rois/']);

%%%%%%%%%%%%%
%% Options %%
%%%%%%%%%%%%%

opt.folder_out = [root_path 'blind_basc_msteps_bascglm_20141101/']; % Where to store the results
%opt.grid_scales = [2:49 50:5:200 210:10:500 500:50:950]; % Search in the range 2-950 clusters
opt.grid_scales = [10:10:300];
opt.scales_maps = [10 7 7;...
20 18 16;...
50 40 40;...
90 81 77;...
180 162 165;...
280 308 313]; % The scales that will be used to generate the maps of brain clusters and stability
%opt.scales_maps = repmat ([10:10:300]',1,3);
opt.stability_tseries.nb_samps = 100; % Number of bootstrap samples at the individual level. 100: the CI on indidividual stability is +/-0.1
opt.stability_group.nb_samps = 500; % Number of bootstrap samples at the group level. 500: the CI on group stability is +/-0.05

opt.flag_mixed = false;
opt.flag_ind = false;


%%%%%%%%%%%%%%%%%%%%%%
%% Run the pipeline %%
%%%%%%%%%%%%%%%%%%%%%%
%opt.psom.max_queued = 10;
opt.psom.qsub_options = '-q sw -l nodes=1:ppn=4,walltime=03:00:00';
pipeline = niak_pipeline_stability_rest(files_in,opt);