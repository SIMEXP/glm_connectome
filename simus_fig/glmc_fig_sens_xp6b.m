%% Prepare the figures for the GLM-connectome simulation
%% histogram of p-values

clear

%% Folders
path_results = '/home/pbellec/database/glm_connectome/simu_multiscale_xp6b/';
path_fig = '/home/pbellec/database/glm_connectome/figures/fig_sens_xp6b/';
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
    tot_sens = 0;
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
            sens = zeros(nb_scale,1);
            sens_max = 0;
        end        
        for rr = 1:length(data.results)
            perc_discovery = zeros(nb_scale,1);
            tot_discovery = 0;
            tot_true = 0;
            pce = 0;
            for ss = 1:nb_scale
                truth = data.results(rr).mask_truth{ss};
                truth = truth*(truth')>0;
                tot_discovery = tot_discovery + sum(data.results(rr).test_q{ss}(:));
                tot_true = tot_true + sum(truth(:));
                perc_discovery(ss) = sum(data.results(rr).test_q{ss}(:))/length(data.results(rr).test_q{ss}(:));
                if any(truth(:))&&(pce<=0.05)        
                    sens(ss) = sens(ss) + (sum(data.results(rr).test_q{ss}(truth))/sum(truth(:)));
                end
            end
            if (tot_true>0)&&(pce<=0.05)
                tot_sens = tot_sens + (tot_discovery/tot_true);
            end
            [val,ind] = max(perc_discovery);
            truth = data.results(rr).mask_truth{ind};
            truth = truth*(truth')>0;
            if any(truth(:))&&(pce<=0.05)
                sens_max = sens_max + (sum(data.results(rr).test_q{ind}(truth))/sum(truth(:)));
            end
        end
    end
    sens = sens / (opt.nb_replication * nb_samps);
    tot_sens = tot_sens / (opt.nb_replication * nb_samps);
    sens_max = sens_max / (opt.nb_replication * nb_samps);
    hold on
    plot(list_scale,repmat(sens_max,size(list_scale)),'r')
    plot(list_scale,repmat(tot_sens,size(list_scale)),'g')
    plot(list_scale,sens,'b')    
    axis([0 max(list_scale)+1 -0.025 1.025])
    name_fig = sprintf('simu_alpha2x%i_nb_subject%i_sc%i',ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref);
    title(strrep(name_fig,'_',' '))
    print(hf,[path_fig name_fig '.pdf'],'-dpdf')
    close(hf)
    name_res = sprintf('simu_alpha2x%i_nb_subject%i_sc%i.mat',ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref);
    save([path_fig name_res],'list_scale','sens')
end
end
end
        
