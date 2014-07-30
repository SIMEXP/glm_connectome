%% Prepare the figures for the GLM-connectome simulation
%% histogram of p-values

%  xp7a simu_alpha2x10_nb_subject40_sc4
%     sensitivity max: 0.43176
%     sensitivity tot: 0.06250
%     sensitivity NBS: 0.33220
%     sensitivity MDMR: 0.10170
%  xp7b simu_alpha2x10_nb_subject40_sc4
%     sensitivity max: 0.64080
%     sensitivity tot: 0.02955
%     sensitivity NBS: 0.33220
%     sensitivity MDMR: 0.10170
%  xp7a simu_alpha2x20_nb_subject40_sc4
%     sensitivity max: 0.92344
%     sensitivity tot: 0.76666
%     sensitivity NBS: 0.90259
%     sensitivity MDMR: 0.82325
%  xp7b simu_alpha2x20_nb_subject40_sc4
%     sensitivity max: 0.99867
%     sensitivity tot: 0.67620
%     sensitivity NBS: 0.90259
%     sensitivity MDMR: 0.82325
%  xp7a simu_alpha2x10_nb_subject100_sc4
%     sensitivity max: 0.86020
%     sensitivity tot: 0.48091
%     sensitivity NBS: 0.74970
%     sensitivity MDMR: 0.57784
%  xp7b simu_alpha2x10_nb_subject100_sc4
%     sensitivity max: 0.97978
%     sensitivity tot: 0.37144
%     sensitivity NBS: 0.74970
%     sensitivity MDMR: 0.57784
%  xp7a simu_alpha2x20_nb_subject100_sc4
%     sensitivity max: 0.95533
%     sensitivity tot: 0.98782
%     sensitivity NBS: 0.97314
%     sensitivity MDMR: 0.99466
%  xp7b simu_alpha2x20_nb_subject100_sc4
%     sensitivity max: 1.00000
%     sensitivity tot: 0.98161
%     sensitivity NBS: 0.97314
%     sensitivity MDMR: 0.99466
%  xp7a simu_alpha2x10_nb_subject40_sc7
%     sensitivity max: 0.12409
%     sensitivity tot: 0.00750
%     sensitivity NBS: 0.05924
%     sensitivity MDMR: 0.02628
%  xp7b simu_alpha2x10_nb_subject40_sc7
%     sensitivity max: 0.24728
%     sensitivity tot: 0.00263
%     sensitivity NBS: 0.05924
%     sensitivity MDMR: 0.02628
%  xp7a simu_alpha2x20_nb_subject40_sc7
%     sensitivity max: 0.80404
%     sensitivity tot: 0.34151
%     sensitivity NBS: 0.63577
%     sensitivity MDMR: 0.44688
%  xp7b simu_alpha2x20_nb_subject40_sc7
%     sensitivity max: 0.95902
%     sensitivity tot: 0.20363
%     sensitivity NBS: 0.63577
%     sensitivity MDMR: 0.44688
%  xp7a simu_alpha2x10_nb_subject100_sc7
%     sensitivity max: 0.55628
%     sensitivity tot: 0.11928
%     sensitivity NBS: 0.31981
%     sensitivity MDMR: 0.19514
%  xp7b simu_alpha2x10_nb_subject100_sc7
%     sensitivity max: 0.80048
%     sensitivity tot: 0.05674
%     sensitivity NBS: 0.31981
%     sensitivity MDMR: 0.19514
%  xp7a simu_alpha2x20_nb_subject100_sc7                                                                                   
%     sensitivity max: 0.95301
%     sensitivity tot: 0.90550
%     sensitivity NBS: 0.94711
%     sensitivity MDMR: 0.94675
%  xp7b simu_alpha2x20_nb_subject100_sc7
%     sensitivity max: 0.99286
%     sensitivity tot: 0.84871
%     sensitivity NBS: 0.94711
%     sensitivity MDMR: 0.94675


clear

%% Folders
list_xp = {'xp7a','xp7b'};
path_results = '/home/pbellec/database/glm_connectome/simu_multiscale_';
path_fig = '/home/pbellec/database/glm_connectome/figures/fig_sens_xp7ab/';
path_xp6 = '/home/pbellec/database/glm_connectome/figures/fig_sens_xp6a/';
path_xp8 = '/home/pbellec/database/glm_connectome/figures/fig_sens_xp8a/';
psom_mkdir(path_fig);

%% Parameters
nb_samps = 100;
%param.sc = {0,10,30,100};

%% At scale 3, cluster 3 is 69.5% of the grey matter, which makes for a percentage of true difference close to 50% (48%)
%% At scale 7, cluster 7 is 20.5% of the grey matter, which makes for a percentage of true difference close to 4% (4.2%)S
param.sc      = {4,7};
param.cluster = {4,7};

param.nsub = {40,100};
param.alpha2 = {0.1,0.2};
opt.nb_replication = 10;


%% Build figures
for num_sc = 1:length(param.sc)
    for num_nsub = 1:length(param.nsub)
        for num_alpha2 = 1:length(param.alpha2)
            hf = figure;    
            % Load data
            opt.alpha2     = param.alpha2{num_alpha2};        
            opt.nb_subject = param.nsub{num_nsub};        
            opt.scale_ref  = param.sc{num_sc};
            tot_sens = 0;
            for num_xp = 1:length(list_xp)
                name_xp = list_xp{num_xp};
                fprintf('%s simu_alpha2x%i_nb_subject%i_sc%i\n',name_xp,ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref)
                for num_s = 1:nb_samps        
                    niak_progress(num_s,nb_samps);
                    name_job = sprintf('simu_alpha2x%i_nb_subject%i_sc%i_samp%i',ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref,num_s);
                    file_data = [path_results name_xp filesep name_job '.mat'];
                    data = load(file_data);
                    if num_s == 1
                        nb_scale = length(data.results(1).test_q);
                        list_scale = zeros(nb_scale,1);
                        for ss = 1:nb_scale
                            list_scale(ss) = size(data.results(1).test_q{ss},1);
                        end
                        sens = zeros(nb_scale,1);
                        sens_max = 0;
                    end        
                    for rr = 1:length(data.results)
                        perc_discovery = zeros(nb_scale,1);
                        tot_discovery = 0;
                        tot_true_discovery = 0;
                        tot_true = 0;
                        pce = data.results(rr).p_nb_disc; 
                        for ss = 1:nb_scale
                            truth = data.results(rr).mask_truth{ss};
                            truth = truth*(truth')>0;
                            tot_discovery = tot_discovery + sum(data.results(rr).test_q{ss}(:));
                            tot_true_discovery = tot_true_discovery + sum(data.results(rr).test_q{ss}(:)&truth(:));
                            tot_true = tot_true + sum(truth(:));
                            perc_discovery(ss) = sum(data.results(rr).test_q{ss}(:))/length(data.results(rr).test_q{ss}(:));
                            if any(truth(:))&&(pce<=0.05)        
                                sens(ss) = sens(ss) + (sum(data.results(rr).test_q{ss}(truth))/sum(truth(:)));
                            end
                        end
                        if (tot_true>0)&&(pce<=0.05)
                            tot_sens = tot_sens + (tot_true_discovery/tot_true);
                        end
                        [val,ind] = max(perc_discovery);
                        truth = data.results(rr).mask_truth{ind};
                        truth = truth*(truth')>0;
                        if any(truth(:))&&(pce<=0.05)
                            sens_max = sens_max + (sum(data.results(rr).test_q{ind}(truth))/sum(truth(:)));
                        end
                    end
                end
                name_fig = sprintf('simu_alpha2x%i_nb_subject%i_sc%i',ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref);
                try
                    file_xp6 = [path_xp6 filesep name_fig '.mat'];
                    data_xp6 = load(file_xp6);
                    sens_nbs = data_xp6.sens(end);
                catch
                    sens_nbs = NaN;
                end
                
                try
                    file_xp8 = [path_xp8 filesep name_fig '.mat'];
                    data_xp8 = load(file_xp8);
                    sens_mdmr = data_xp8.sens(end);
                catch
                    sens_mdmr = NaN;
                end


                sens = sens / (opt.nb_replication * nb_samps);
                tot_sens = tot_sens / (opt.nb_replication * nb_samps);
                sens_max = sens_max / (opt.nb_replication * nb_samps);
                
                fprintf('   sensitivity max: %1.5f\n',sens_max);
                fprintf('   sensitivity tot: %1.5f\n',tot_sens);
                fprintf('   sensitivity NBS: %1.5f\n',sens_nbs);
                fprintf('   sensitivity MDMR: %1.5f\n',sens_mdmr);
            end
        end
    end
end        
