%% Prepare the figures for the GLM-connectome simulation
%% histogram of p-values

clear

%% Folders
path_results = '/home/pbellec/database/glm_connectome/simu_multiscale_xp6a/';
path_fig = '/home/pbellec/database/glm_connectome/figures/fig_fdr_xp6a/';
psom_mkdir(path_fig);

%% Parameters
nb_samps = 10;
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
            for num_s = 1:nb_samps        
                name_job = sprintf('simu_alpha2x%i_nb_subject%i_sc%i_samp%i',ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref,num_s);
                file_data = [path_results name_job '.mat'];
                data = load(file_data);
                if num_s == 1
                    nb_scale = length(data.results(1).test_q);
                    list_scale = zeros(nb_scale,1);
                    for ss = 1:nb_scale
                        list_scale(ss) = size(data.results(1).test_q{ss},1);
                    end
                    fdr = zeros(nb_scale,1);
                    fdr_max = 0;
                    fdr_all = 0;
                end        
                for rr = 1:length(data.results)
                    perc_discovery = zeros(nb_scale,1);
                    all_false = zeros(nb_scale,1);
                    all_disc = zeros(nb_scale,1);
                    pce = 0; 
                    for ss = 1:nb_scale
                        truth = data.results(rr).mask_truth{ss};
                        truth = truth*(truth')>0;
                        perc_discovery(ss) = sum(data.results(rr).test_q{ss}(:))/length(data.results(rr).test_q{ss}(:));
                        if ~any(truth(:))
                            if (pce<=0.05)
                                fdr(ss) = fdr(ss) + any(data.results(rr).test_q{ss}(:));
                            end
                            all_false(ss) = sum(data.results(rr).test_q{ss}(:));
                            all_disc(ss) = sum(data.results(rr).test_q{ss}(:));
                        elseif any(data.results(rr).test_q{ss}(:))
                            if (pce<=0.05)
                                fdr(ss) = fdr(ss) + (sum(data.results(rr).test_q{ss}(~truth))/sum(data.results(rr).test_q{ss}(:)));
                            end
                            all_false(ss) = sum(data.results(rr).test_q{ss}(~truth));
                            all_disc(ss) = sum(data.results(rr).test_q{ss}(:));
                        end                
                    end
                    if (sum(all_disc)>0)&&(pce<=0.05)
                        fdr_all = fdr_all + sum(all_false)/sum(all_disc);            
                    end
                    [val,ind] = max(perc_discovery);
                    truth = data.results(rr).mask_truth{ind};
                    truth = truth*(truth')>0;
                    if ~any(truth(:))
                        if (pce<=0.05)
                            fdr_max = fdr_max + any(data.results(rr).test_q{ind}(:));
                        end
                    elseif any(data.results(rr).test_q{ind}(:))
                        if (pce<=0.05)
                            fdr_max = fdr_max + (sum(data.results(rr).test_q{ind}(~truth))/sum(data.results(rr).test_q{ind}(:)));
                        end
                    end
                end
            end
            fdr = fdr / (opt.nb_replication * nb_samps);
            fdr_max = fdr_max / (opt.nb_replication * nb_samps);
            fdr_all = fdr_all / (opt.nb_replication * nb_samps);
            hold on
            plot(list_scale,repmat(fdr_all,size(list_scale)),'g')
            plot(list_scale,repmat(fdr_max,size(list_scale)),'r')    
            plot(list_scale,fdr,'b')
            name_fig = sprintf('simu_alpha2x%i_nb_subject%i_sc%i',ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref);
            title(strrep(name_fig,'_',' '))
            axis([-10 max(list_scale) -0.01 0.2]);
            print(hf,[path_fig name_fig '.pdf'],'-dpdf')    
            close(hf)
            name_res = sprintf('simu_alpha2x%i_nb_subject%i_sc%i.mat',ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref);
            save([path_fig name_res],'list_scale','fdr')
        end
    end
end
        
