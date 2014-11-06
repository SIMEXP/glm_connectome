%% Prepare the figures for the GLM-connectome simulation
%% histogram of p-values

clear

%% Set threshold for omnibus test 
list_pce = [0.05,Inf];

%% Set up paths
path_results = pwd;
if ~psom_exist([path_results filesep 'simu_param.mat'])
    error('Could not find the results of the simulations')
end
path_fig = [path_results filesep 'figures'];
psom_mkdir(path_fig);

%% Load parameters
param = load([path_results filesep 'simu_param.mat']);
nb_fdr = length(param.list_fdr);
nb_scale = length(param.list_scales);

%% Create figure handles and plot indices
hfdr = figure; % Create a figure for the FDR plots
num_plot_fdr = 1; % keep track of the number of the current FDR plot

%% Build figures
num_p = 1;
for num_alpha2 = 1:length(param.alpha2) % Loop over effect size
    for num_nsub = 1:length(param.nsub) % Loop over the number of subjects
        for num_sc = [1 3 2] % Loop over reference clusters
            % Load data
            name_xp = sprintf('simu_perc%i_sc%i_a2%i_nsub%i',ceil(100*param.perc_rand(num_p)),param.sc{num_sc},ceil(100*param.alpha2{num_alpha2}),param.nsub{num_nsub});                            
            fprintf('%s\n',name_xp)
            nb_samps = 0;
            for num_s = 1:param.nb_samps % Loop over simulation jobs
%            for num_s = 1 % Loop over simulation jobs
                niak_progress(num_s,param.nb_samps);
                name_job = sprintf('simu_a2%i_nsub%i_sc%i_perc%i_samp%i',ceil(100*param.alpha2{num_alpha2}),param.nsub{num_nsub},param.sc{num_sc},ceil(100*param.perc_rand(num_p)),num_s);                    
                file_data = [path_results filesep name_job '.mat'];
                if num_s == 1
                    fdr = zeros(nb_scale,nb_fdr); % The within-scale FDR
                    fdr_all = zeros(nb_fdr,length(list_pce)); % The between-scales FDR
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
                    for pp = 1:length(list_pce)
                        for ff = 1:nb_fdr
                            if any(all_disc(:,ff))&&(pce(ff)<=list_pce(pp))
                                fdr_all(ff,pp) = fdr_all(ff,pp) + sum(all_false(:,ff))/sum(all_disc(:,ff),1);            
                            end
                        end
                    end
                end
            end
            fdr = fdr / nb_samps;
            fdr_all = fdr_all / nb_samps;
            switch num_sc
                case 1
                    name_fig = 'True Positive Rate 0%';
                case 2
                    name_fig = 'True Positive Rate \sim15%';
                case 3
                    name_fig = 'True Positive Rate \sim5%';
            end
            hp = subplot(length(param.nsub)*length(param.alpha2),length(param.sc),num_plot_fdr);
            FN = findall(hp,'-property','FontName');
            set(FN,'FontName','/usr/share/fonts/truetype/dejavu/DejaVuSerifCondensed.ttf');
            FS = findall(hp,'-property','FontSize');
            set(FS,'FontSize',8);
            hold on
            plot([0; max(param.list_fdr)],[0; max(param.list_fdr)],'k','LineWidth',4)            
            plot(param.list_fdr(:),fdr','b')
            plot(param.list_fdr(:),fdr_all(:,1),'r','LineWidth',4)     
            plot(param.list_fdr(:),fdr_all(:,2),'g','LineWidth',4)                      
            axis([0 max(param.list_fdr) 0 0.2]);
            if (num_nsub == length(param.nsub))&&(num_alpha2==length(param.alpha2))
                xlabel('Nominal FDR');
            end
            if (num_nsub == 1)&&(num_alpha2==1)
                title(name_fig)
            end
            if num_sc==1
                ylabel(sprintf('Effective FDR\n(N=%i,eff=%1.2f)',param.nsub{num_nsub},param.alpha2{num_alpha2}))
            end
            if num_plot_fdr==(length(param.nsub)*length(param.alpha2)*length(param.sc))
                print(hfdr,[path_fig filesep 'fig_fdr.png'],'-dpng')
            end
            num_plot_fdr = num_plot_fdr+1;
        end
    end
end
%close(hfdr)