%% Prepare the figures for the GLM-connectome simulation
%% histogram of p-values

clear

%% Set threshold for omnibus test 
thre_omnibus = Inf;

%% Set up paths
path_results = pwd;
if ~psom_exist([path_results filesep 'simu_param.mat'])
    error('Could not find the results of the simulations')
end
path_fig = [path_results filesep sprintf('fig_simus_fdr_p%i',ceil(100*thre_omnibus))];
psom_mkdir(path_fig);

%% Load parameters
param = load([path_results filesep 'simu_param.mat']);
nb_fdr = length(param.list_fdr);
nb_scale = length(param.list_scales);

%% Build figures
num_p = 1;
for num_sc = 1:length(param.sc) % Loop over reference clusters
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
                    fdr = zeros(nb_scale,nb_fdr); % The within-scale FDR
                    fdr_all = zeros(1,nb_fdr); % The between-scales FDR
                end                            
                if ~psom_exist(file_data)
                    warning('I could not find the simulation results for one sample (%i). I am going to skip it.',num_s)
                    continue % Be robust to missing files, to check intermediate results of the pipeline, as it's running
                end
                data = load(file_data);
                nb_samps = nb_samps+param.nb_replication;
                for rr = 1:param.nb_replication % Loop over replications inside the job
                    all_false = zeros(nb_scale,nb_fdr);
                    all_disc = zeros(nb_scale,nb_fdr);
                    pce = data.results(rr).p_nb_disc; 
                    for ss = 1:nb_scale % Loop over tested scales
                        truth = data.results(rr).mask_truth{ss};
                        truth = truth*(truth')>0;
                        for ff = 1:nb_fdr % Loop over tested FDR thresholds
                            all_false(ss,ff) = sum(data.results(rr).test_q{ss,ff}(~truth));
                            all_disc(ss,ff) = sum(data.results(rr).test_q{ss,ff}(:));
                            if (all_disc(ss,ff)>0)
                                fdr(ss,ff) = fdr(ss,ff) + (all_false(ss,ff)/all_disc(ss,ff));
                            end 
                        end
                    end
                    for ff = 1:nb_fdr
                        if any(all_disc(:,ff))&&(pce(ff)<=thre_omnibus)
                            fdr_all(ff) = fdr_all(ff) + sum(all_false(:,ff))/sum(all_disc(:,ff),1);            
                        end
                    end
                end
            end
            fdr = fdr / nb_samps;
            fdr_all = fdr_all / nb_samps;
            
            hold on
            plot([0; max(param.list_fdr)],[0; max(param.list_fdr)],'k')            
            plot(param.list_fdr(:),fdr','b')
            plot(param.list_fdr(:),fdr_all','r')            
            axis([0 max(param.list_fdr) 0 0.2]);
            xlabel('Nominal FDR')
            ylabel('Effective FDR')
            title(strrep(name_fig,'_',' '))
            print(hf,[path_fig filesep name_fig '.pdf'],'-dpdf')
            close(hf)
        end
    end
end
