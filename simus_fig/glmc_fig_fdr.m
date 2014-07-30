%% Prepare the figures for the GLM-connectome simulation
%% histogram of p-values

clear

%% Folders
path_results = '/home/pbellec/database/glm_connectome/simu_multiscale/';
path_fig = '/home/pbellec/database/glm_connectome/fig_fdr/';
psom_mkdir(path_fig);

%% Parameters
param.sc = {0,5,10};
%param.sc = {0};
param.nsub = {20,90};
param.alpha2 = {0.1,0.2};
nb_samps = 10;
opt.nb_replication = 100;

%% Build figures
for num_sc = 1:length(param.sc)
for num_nsub = 1:length(param.nsub)
for num_alpha2 = 1:length(param.alpha2)
    hf = figure;    
    % Load data
    pce = zeros(nb_samps*opt.nb_replication,1);
    npce = 1;
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
            for ss = 1:nb_scale
                truth = data.results(rr).mask_truth{ss};
                truth = truth*(truth')>0;
                perc_discovery(ss) = sum(data.results(rr).test_q{ss}(:))/length(data.results(rr).test_q{ss}(:));
                if ~any(truth(:))
                    fdr(ss) = fdr(ss) + any(data.results(rr).test_q{ss}(:));
                    all_false(ss) = sum(data.results(rr).test_q{ss}(:));
                    all_disc(ss) = sum(data.results(rr).test_q{ss}(:));
                elseif any(data.results(rr).test_q{ss}(:))
                    fdr(ss) = fdr(ss) + (sum(data.results(rr).test_q{ss}(~truth))/sum(data.results(rr).test_q{ss}(:)));
                    all_false(ss) = sum(data.results(rr).test_q{ss}(~truth));
                    all_disc(ss) = sum(data.results(rr).test_q{ss}(:));
                end                
            end
            if sum(all_disc)>0
                fdr_all = fdr_all + sum(all_false)/sum(all_disc);            
            end
            [val,ind] = max(perc_discovery);
            truth = data.results(rr).mask_truth{ind};
            truth = truth*(truth')>0;
            if ~any(truth(:))
                fdr_max = fdr_max + any(data.results(rr).test_q{ind}(:));
            elseif any(data.results(rr).test_q{ind}(:))
                fdr_max = fdr_max + (sum(data.results(rr).test_q{ind}(~truth))/sum(data.results(rr).test_q{ind}(:)));
            end
        end
    end
    fdr = fdr / (opt.nb_replication * nb_samps);
    fdr_max = fdr_max / (opt.nb_replication * nb_samps);
    fdr_all = fdr_all / (opt.nb_replication * nb_samps);
    plot([-1 ; 0 ; list_scale],[fdr_all ; fdr_max ; fdr],'*')
    name_fig = sprintf('simu_alpha2x%i_nb_subject%i_sc%i',ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref);
    title(strrep(name_fig,'_',' '))
    axis([-10 max(list_scale) 0 0.2]);
    print(hf,[path_fig name_fig '.png'],'-dpng')    
    close(hf)
end
end
end
        