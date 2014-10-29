%% Prepare the figures for the GLM-connectome simulation
%% histogram of p-values

clear

%% Set up paths
path_results = pwd;
if ~psom_exist([path_results filesep 'simu_param.mat'])
    error('Could not find the results of the simulations')
end
path_fig = [path_results filesep 'fig_simus_sens'];
psom_mkdir(path_fig);

%% Load parameters
thre_omnibus = 0.05;
param = load([path_results filesep 'simu_param.mat']);
nb_fdr = length(param.list_fdr);
nb_scale = length(param.list_scales);
list_color = {'b','g','r','k'};

%% Build figures
for num_p = 1:length(param.perc_rand)
    for num_sc = 2:length(param.sc) % Loop over reference clusters
        for num_nsub = 1:length(param.nsub) % Loop over the number of subjects
            for num_alpha2 = 1:length(param.alpha2) % Loop over effect size
                hf = figure;    
                % Load data
                name_fig = sprintf('simu_perc%i_sc%i_a2%i_nsub%i',ceil(100*param.perc_rand(num_p)),param.sc{num_sc},ceil(100*param.alpha2{num_alpha2}),param.nsub{num_nsub});                            
                fprintf('%s\n',name_fig)
                nb_samps = 0;
                for num_s = 1:param.nb_samps % Loop over simulation jobs
                    niak_progress(num_s,param.nb_samps);
                    name_job = sprintf('simu_a2%i_nsub%i_sc%i_perc%i_samp%i',ceil(100*param.alpha2{num_alpha2}),param.nsub{num_nsub},param.sc{num_sc},ceil(100*param.perc_rand(num_p)),num_s);                    
                    file_data = [path_results filesep name_job '.mat'];
                    if num_s == 1
                        sens = zeros(nb_scale,nb_fdr); % The within-scale sensitivity, as a function of the FDR threshold
                    end                            
                    if ~psom_exist(file_data)
                        warning('I could not find the simulation results for one sample (%i). I am going to skip it.',num_s)
                        continue % Be robust to missing files, to check intermediate results of the pipeline, as it's running
                    end
                    data = load(file_data);
                    nb_samps = nb_samps+param.nb_replication;
                    for rr = 1:param.nb_replication % Loop over replications inside the job
                        all_true = zeros(nb_scale,nb_fdr);
                        pce = data.results(rr).p_nb_disc; 
                        all_true = 0;
                        for ss = 1:nb_scale % Loop over tested scales
                            truth = data.results(rr).mask_truth{ss};
                            truth = truth*(truth')>0;
                            nb_true = sum(truth(:));
                            all_true = all_true + nb_true;
                            for ff = 1:nb_fdr % Loop over tested FDR thresholds
                                if (pce(ff)<=thre_omnibus)
                                    true_disc(ss,ff) = sum(data.results(rr).test_q{ss,ff}(truth));
                                    sens(ss,ff) = sens(ss,ff) + (true_disc(ss,ff)/nb_true);
                                end
                            end
                        end
                    end
                end
                sens = sens / nb_samps;
                
                %% Make the figure
                hold on    
                for ff = 1:length(param.list_fdr)
                    plot(param.list_scales(:),sens(:,ff),list_color{ff});
                end
                legend(cellstr(num2str(param.list_fdr'))')
                axis([0 max(param.list_scales) 0 1]);
                xlabel('Scale')
                ylabel('Sensitivity')
                title(strrep(name_fig,'_',' '))
                print(hf,[path_fig filesep name_fig '.pdf'],'-dpdf')
                close(hf)
            end
        end
    end
end