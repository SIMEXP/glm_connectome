% A script to generate a pipeline to run multiscale simulations for the 
% GLM_connectome project - case of perfect match between the simulation and test clusters
clear

warning('This script will generate a lot of results in the current folder. Press CTRL-C now to interrupt !')
pause

%% Set up paths
path_curr = pwd;
path_roi  = [path_curr filesep 'rois']; % Where to save the real regional time series
path_out  = [path_curr filesep 'simu_perfect_match']; % Where to store the results of the simulation
path_logs = [path_out filesep 'logs']; % Where to save the logs of the pipeline
psom_mkdir(path_out);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Download the Cambridge - 1000 functional connectomes - time series %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~psom_exist(path_roi)
    mkdir(path_roi)
    cd(path_roi)
    fprintf('Could not find the Cambridge time series. Downloading from figshare ...\n')
    instr_dwnld = 'wget http://downloads.figshare.com/article/public/1159331';
    [status,msg] = system(instr_dwnld);
    if status~=0
        psom_clean(path_rois)
        error('Could not download the necessary data from figshare. The command was: %s. The error message was: %s',instr_dwnld,msg);
    end
    instr_unzip = 'unzip 1159331';
    [status,msg] = system(instr_unzip);
    if status~=0
        psom_clean(path_roi)
        error('Could not unzip the necessary data. The command was: %s. The error message was: %s',instr_unzip,msg);
    end
    psom_clean('1159331');
    cd(path_curr)
end
    
in.tseries = dir([path_roi filesep 'tseries_rois_*_session1_rest.mat']);
in.tseries = {in.tseries.name};
for ff = 1:length(in.tseries)
    in.tseries{ff} = [path_roi filesep in.tseries{ff}];
end
in.hier = [path_roi filesep 'hier_avg_connectome.mat'];

%%%%%%%%%%%%%%%%
%% PARAMETERS %% 
%%%%%%%%%%%%%%%%

% The number of simulation samples. 
% Note that each simulation job is performing itself several simulations. 
% The final number of simulations is NB_SAMPS*NB_REPLICATION
param.nb_samps       = 20; 
param.nb_replication = 50;

% The list of scales to be tested 
param.list_scales = [7 16 25 55 114 199 328];

% The number of permutation samples for the omnibus test
param.nb_perm = 100;

% The FDR thresholds
param.list_fdr = [0.01 0.05 0.1 0.2];

% The percentage of mismatch between the reference and test clusters
param.perc_rand = [0 0.3];

% Selection of the cluster of reference
% At scale 0, there is not simulated signal (global null hypothesis)
% At scale 8, cluster 3 is 10.4% of the grey matter, which makes for a percentage of true difference close to 1% (1.0%)
% At scale 7, cluster 3 is 13.7% of the grey matter, which makes for a percentage of true difference close to 2% (1.9%)
% At scale 5, cluster 5 is 22.6% of the grey matter, which makes for a percentage of true difference close to 5% (5.1%)
% At scale 4, cluster 3 is 30.8% of the grey matter, which makes for a percentage of true difference close to 10% (9.6%)
param.sc      = { 0 , 8 , 7 , 5 , 4 }; % The scale of reference that use used to define the ground truth
param.cluster = { 0 , 3 , 3 , 5 , 3 }; % The number of the cluster that is used at the scale of reference. The fact that SC and CLUSTER are identical is a coincidence.

% Number of subjects per group
param.nsub           = {  40 , 100 }; 

% The alpha2 parameters, which sets the effect size 
% the actual effect size is scale dependent, and will be evaluated directly from the simulation
param.alpha2         = { 0.1 , 0.2 }; 

% save the paramaters of the simulation
save([path_out filesep 'simu_param.mat'],'-struct','param')

%%%%%%%%%%%%%%%%%%%%%%%%
%% Build the pipeline %%
%%%%%%%%%%%%%%%%%%%%%%%%
num_seed = 0;  % num_seed will be used the seed the random number generate of matlab, and ensure that all simulations are 100% reproducible
pipe = struct; % Initialization of the pipeline. All simulations will be added as "jobs" in the structure PIPE, using the PSOM framework http://psom.simexp-lab.org

%% set the option parameters common to all simulations
opt.list_scales = param.list_scales;
opt.nb_samps = param.nb_perm; % The number of permutation sample for the omnibus test
opt.nb_replication = param.nb_replication; % The number of simulations inside the job. The final number of simulations is NB_SAMPS*NB_REPLICATION
opt.list_fdr = param.list_fdr; % The FDR thresholds
    
for num_s = 1:param.nb_samps % Loop over all simulation samples. Note that each simulation job is performing itself several simulations. The final number of simulations is NB_SAMPS*NB_REPLICATION 
    for num_sc = 1:length(param.sc) % Loop over the scale of reference
        for num_nsub = 1:length(param.nsub) % Loop over the number of subjects per group
            for num_alpha2 = 1:length(param.alpha2) % Loop over the effect size
                for num_p = 1:length(param.perc_rand); % Loop over the percentage of mismatch between the true and test clusters
                    % Job-specific parameters
                    opt.alpha2      = param.alpha2{num_alpha2}; % Effect size
                    opt.nb_subject  = param.nsub{num_nsub};     % # subject per group
                    opt.scale_ref   = param.sc{num_sc};         % Scale of reference (sets the percentage of true positives)
                    opt.cluster_ref = param.cluster{num_sc};    % Scale of reference, continued
                    opt.perc_rand   = param.perc_rand(num_p);
                    % The name of the simulation
                    name_job = sprintf('simu_a2%i_nsub%i_sc%i_perc%i_samp%i',ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref,ceil(100*opt.perc_rand),num_s);
                    % The output file
                    out = sprintf('%s/%s.mat',path_out,name_job);
                    opt.rand_seed = num_seed;
                    num_seed = num_seed+1;        
                    pipe = psom_add_job(pipe,name_job,'glmc_brick_simu_multiscale',in,out,opt);
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%
%% Run the pipeline %%
%%%%%%%%%%%%%%%%%%%%%%
opt_p.path_logs = path_logs;
opt_p.flag_pause = false;
opt_p.qsub_options = '-q sw -l walltime=05:00:00';
%psom_run_pipeline(pipe,opt_p);