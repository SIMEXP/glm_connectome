%% Prepare the figures for the GLM-connectome simulation
%% histogram of p-values

clear

%% Folders
list_xp = {'xp7a','xp7b'};
path_results = '/home/pbellec/database/glm_connectome/simu_multiscale_';
path_fig = '/home/pbellec/database/glm_connectome/figures/fig_eff_size_xp7/';
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

bins = 0:0.05:1;
delta = bins(2)-bins(1);

%% Build figures
for num_sc = 1:length(param.sc)
    for num_nsub = 1:length(param.nsub)
        for num_alpha2 = 1:length(param.alpha2)
            hf = figure;    
            clf
            % Load data    
            nb_disc = zeros(nb_samps*opt.nb_replication,1);
            npce = 1;
            opt.alpha2     = param.alpha2{num_alpha2};        
            opt.nb_subject = param.nsub{num_nsub};        
            opt.scale_ref  = param.sc{num_sc};
            for num_xp = 1:length(list_xp)
                name_xp = list_xp{num_xp};
                fprintf('simu_alpha2x%i_nb_subject%i_sc%i\n',ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref)
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
                        eff_size = zeros(nb_samps*opt.nb_replication,length(list_scale));
                    end
                    for rr = 1:length(data.results)                        
                        for ss = 1:length(data.results(rr).test_q)
                            if ~isempty(data.results(rr).mean_diff{ss})
                                mat_tmp = data.results(rr).mean_diff{ss}./data.results(rr).std_iir{ss};
                                if any((data.results(rr).std_iir{ss}(:)==0)&(data.results(rr).mean_diff{ss}(:)~=0))
                                    error('Huho')
                                end
                                mat_tmp(data.results(rr).std_iir{ss}(:)==0) = 0;
                                eff_size(npce,ss) = mean(mean(mat_tmp));  
                            end
                        end
                        npce = npce+1;
                    end
                end
            
                m_eff_size = mean(eff_size,1);
                s_eff_size = std(eff_size,[],1)/sqrt(size(eff_size,1));
                hold on
                if num_xp == 1
                    plot(list_scale',m_eff_size,'b')
                else
                    plot(list_scale',m_eff_size,'r')
                    axis([0 max(list_scale) 0 3]);
                    ha = gca;
                    set(ha,'ygrid','on')
                    name_fig = sprintf('simu_alpha2x%i_nb_subject%i_sc%i',ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref);
                    title(strrep(name_fig,'_',' '))
                    print(hf,[path_fig name_fig '.pdf'],'-dpdf')
                    close(hf)
                end
            end
        end
    end
end

        
