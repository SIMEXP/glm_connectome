% A script to generate a pipeline to run multiscale simulations for the 
% GLM_connectome project

clear

warning('This script will generate a lot of results in the current folder. Press CTRL-C now to interrupt !')
pause
path_curr = pwd;

%% Make sure time series from Cambridge - 1000 functional connectomes - are available
path_roi = [path_curr filesep 'rois'];
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

%% Output
folder_out = [path_curr filesep 'simu_perfect_match']; % Where to store the results of the simulation
nb_samps = 100; % The number of simulation samples. Note that each simulation job is performing itself several simulations. The final number of simulations is NB_SAMPS*NB_REPLICATION

%% At scale 0, there is not simulated signal (global null hypothesis)
%% At scale 4, cluster 4 is 38.7% of the grey matter, which makes for a percentage of true difference close to 15% (14.98%)
%% At scale 7, cluster 7 is 20.5% of the grey matter, which makes for a percentage of true difference close to 4% (4.2%)
param.sc      = { 0 , 4 , 7 }; % The scale of reference that use used to define the ground truth
param.cluster = { 0 , 4 , 7 }; % The number of the cluster that is used at the scale of reference. The fact that SC and CLUSTER are identical is a coincidence.

param.nsub   = {  40 , 100 }; % Number of subjects per group
param.alpha2 = { 0.1 , 0.2 }; % The alpha2 parameters, which sets the effect size (although the actual effect size is scale dependent, and will be evaluated directly from the simulation


%% Build jobs
num_seed = 0;  % num_seed will be used the seed the random number generate of matlab, and ensure that all simulations are 100% reproducible
pipe = struct; % Initialization of the pipeline. All simulations will be added as "jobs" in the structure PIPE, using the PSOM framework http://psom.simexp-lab.org
for num_s = 1:nb_samps % Loop over all simulation samples. Note that each simulation job is performing itself several simulations. The final number of simulations is NB_SAMPS*NB_REPLICATION
    %opt.list_scales = [10:10:300]; % The list of scales to be tested 
    opt.list_scales = [10:10:20];
    opt.nb_samps = 5;            % The number of permutation sample for the omnibus test
    opt.nb_replication = 1;       % The number of simulations inside the job. The final number of simulations is NB_SAMPS*NB_REPLICATION
    opt.list_fdr = [0.01 0.05 0.1 0.2]; % The FDR thresholds
    opt.perc_rand = 0;
    for num_sc = 1:length(param.sc) % Loop over the scale of reference
        for num_nsub = 1:length(param.nsub) % Loop over the number of subjects per group
            for num_alpha2 = 1:length(param.alpha2) % Loop over the effect size
                opt.alpha2      = param.alpha2{num_alpha2}; % Effect size
                opt.nb_subject  = param.nsub{num_nsub};     % # subject per group
                opt.scale_ref   = param.sc{num_sc};         % Scale of reference (sets the percentage of true positives)
                opt.cluster_ref = param.cluster{num_sc};    % Scale of reference, continued
                % The name of the simulation
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
