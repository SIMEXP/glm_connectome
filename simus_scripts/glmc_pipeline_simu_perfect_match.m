% A script to generate a pipeline to run multiscale simulations for the 
% GLM_connectome project

clear

%% Inputs: time series from Cambridge - 1000 functional connectomes
path_roi = '/gs/scratch/pbellec/glm_connectome/cambridge_region_growing_05scrubb/rois';
in.tseries = dir([path_roi filesep 'tseries_rois_*_session1_rest.mat']);
in.tseries = {in.tseries.name};
for ff = 1:length(in.tseries)
    in.tseries{ff} = [path_roi filesep in.tseries{ff}];
end
in.hier = [path_roi filesep 'hier_avg_connectome.mat'];

%% Output
folder_out = '/gs/scratch/pbellec/glm_connectome/simu_multiscale_xp7b/';
type_fdr = 'global';
nb_samps = 100;
%param.sc = {0,10,30,100};

%% At scale 4, cluster 4 is 38.7% of the grey matter, which makes for a percentage of true difference close to 15% (14.98%)
%% At scale 7, cluster 7 is 20.5% of the grey matter, which makes for a percentage of true difference close to 4% (4.2%)
param.sc      = {0,4,7};
param.cluster = {0,4,7};

%param.nsub = {20,90};
param.nsub   = { 40, 100};
param.alpha2 = {0.1, 0.2};


%% Build jobs
num_seed = 0;
pipe = struct;
for num_s = 1:nb_samps
    opt.list_scales = [2:10 15:5:100 110:10:200];
    %opt.list_scales = [2 5 10 50];
    opt.nb_samps = 100;    
    opt.nb_replication = 10;    
    opt.fdr = 0.05;
    opt.type_fdr = type_fdr;
    opt.perc_rand = 0;
    for num_sc = 1:length(param.sc)
    for num_nsub = 1:length(param.nsub)
    for num_alpha2 = 1:length(param.alpha2)
        opt.alpha2      = param.alpha2{num_alpha2};
        opt.nb_subject  = param.nsub{num_nsub};
        opt.scale_ref   = param.sc{num_sc};
        opt.cluster_ref = param.cluster{num_sc};
        name_job = sprintf('simu_alpha2x%i_nb_subject%i_sc%i_samp%i',ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref,num_s);
        out = sprintf('%s/%s.mat',folder_out,name_job);
        opt.rand_seed = num_seed;
        num_seed = num_seed+1;        
        pipe = psom_add_job(pipe,name_job,'glmc_brick_simu_multiscale',in,out,opt);
    end
    end
    end
end
opt_p.path_logs = [folder_out 'logs'];
opt_p.flag_pause = false;
opt_p.qsub_options = '-q sw -l walltime=05:00:00';
