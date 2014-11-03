%% Prepare the figures for the GLM-connectome simulation
%% histogram of p-values

clear

name_xp = 'sens';

%% Set up paths
path_results = pwd;
if ~psom_exist([path_results filesep 'simu_param.mat'])
    error('Could not find the results of the simulations')
end
path_fig = [path_results filesep 'figures'];
psom_mkdir(path_fig);

%% Load parameters
thre_omnibus = 0.05;
param = load([path_results filesep 'simu_param.mat']);
nb_fdr = length(param.list_fdr);
nb_scale = length(param.list_scales);
list_color = {'b','g','r','k'};

%% Create figure handles and plot indices
hf = figure; % Create a figure for the FDR plots
num_plot = 1; % keep track of the number of the current sensitivity plot

%% Build figures
for num_alpha2 = 1:length(param.alpha2) % Loop over effect size
    for num_nsub = 1:length(param.nsub) % Loop over the number of subjects
        for num_p = 1:length(param.perc_rand) % Loop over the degree of matching between test and ground truth clusters
            for num_sc = [3 2] % Loop over reference clusters
        
                % Load data
                fprintf('simu_perc%i_sc%i_a2%i_nsub%i\n',ceil(100*param.perc_rand(num_p)),param.sc{num_sc},ceil(100*param.alpha2{num_alpha2}),param.nsub{num_nsub});                            
                nb_samps = 0;
                for num_s = 1:param.nb_samps % Loop over simulation jobs
                    niak_progress(num_s,param.nb_samps);
                    name_job = sprintf('simu_a2%i_nsub%i_sc%i_perc%i_samp%i',ceil(100*param.alpha2{num_alpha2}),param.nsub{num_nsub},param.sc{num_sc},ceil(100*param.perc_rand(num_p)),num_s);                    
                    file_data = [path_results filesep name_job '.mat'];
                    if num_s == 1
                        sens = zeros(nb_scale,nb_fdr); % The within-scale sensitivity, as a function of the FDR threshold
                        perc_disc = zeros(nb_scale,nb_fdr); % The percentage of discovery
                        tpr = zeros(nb_scale);
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
                                    true_disc = sum(data.results(rr).test_q{ss,ff}(truth));
                                    sens(ss,ff) = sens(ss,ff) + (true_disc/nb_true);
                                    perc_disc(ss,ff) = perc_disc(ss,ff) + sum(data.results(rr).test_q{ss,ff}(:))/length(data.results(rr).test_q{ss,ff}(:));
                                    if (ff == 1)&&(num_s==1)&&(rr==1)
                                        tpr(ss) = sum(truth(:))/length(truth(:));
                                    end
                                end
                            end
                        end
                    end
                end
                sens = sens / nb_samps;
                perc_disc = perc_disc / nb_samps;
                
                %% Make the figure
                figure(hf);
                switch num_sc
                case 2
                    name_fig = 'TPR \sim15%';
                case 3
                    name_fig = 'TPR \sim5%';
                end
                if param.perc_rand(num_p)==0
                    name_fig = sprintf('%s, true/test match',name_fig);
                else 
                    name_fig = sprintf('%s, true/test mismatch',name_fig);
                end
                hp = subplot(length(param.nsub)*length(param.alpha2),4,num_plot);
                hold on
                FN = findall(hp,'-property','FontName');
                set(FN,'FontName','/usr/share/fonts/truetype/dejavu/DejaVuSerifCondensed.ttf');
                FS = findall(hp,'-property','FontSize');
                set(FS,'FontSize',8);
                
                switch name_xp
                case 'sens'
                    for ff = 1:length(param.list_fdr)
                        plot(param.list_scales(:),sens(:,ff),list_color{ff},'LineWidth',4);
                    end
                    axis([0 max(param.list_scales) 0 1]);
                case 'perc_disc'
                    for ff = 1:length(param.list_fdr)
                        plot(param.list_scales(:),perc_disc(:,ff),list_color{ff},'LineWidth',4);
                    end
                    plot(param.list_scales(:),tpr,'y','LineWidth',4)
                otherwise
                    error('%s is an unknown experiment',name_xp)
                end
                if (num_nsub == length(param.nsub))&&(num_alpha2==length(param.alpha2))
                    xlabel('Scale')
                end
                if (num_nsub == 1)&&(num_alpha2==1)
                    title(strrep(name_fig,'_',' '))
                end
                if (num_sc==3)&&(num_p==1)
                    ylabel(sprintf('Sensitivity\n(N=%i,eff=%1.2f)',param.nsub{num_nsub},param.alpha2{num_alpha2}))
                end
                if num_plot==(length(param.nsub)*length(param.alpha2)*2*length(param.perc_rand))
                    print(hf,[path_fig filesep 'fig_' name_xp '.png'],'-dpng')
                end
                num_plot = num_plot+1;
                
            end
        end
    end
end