%% Prepare the figures for the GLM-connectome simulation
%% histogram of p-values

clear

%% Folders
path_results = '/home/pbellec/database/glm_connectome/simu_multiscale_xp7a/';
path_fig = '/home/pbellec/database/glm_connectome/figures/fig_perc_disc_xp7/';
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
for num_nsub = 1
for num_alpha2 = 1
    hf = figure;    
    % Load data
    opt.alpha2     = param.alpha2{num_alpha2};        
    opt.nb_subject = param.nsub{num_nsub};        
    opt.scale_ref  = param.sc{num_sc};
    num_e = 1;
    fprintf('simu_alpha2x%i_nb_subject%i_sc%i\n',ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref)
    for num_s = 1:nb_samps        
        niak_progress(num_s,nb_samps);
        name_job = sprintf('simu_alpha2x%i_nb_subject%i_sc%i_samp%i',ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref,num_s);
        file_data = [path_results name_job '.mat'];
        data = load(file_data);
        if num_s == 1
            nb_scale = length(data.results(1).test_q);
            list_scale = zeros(nb_scale,1);
            for ss = 1:nb_scale
                list_scale(ss) = size(data.results(1).test_q{ss},1);
            end
            perc_disc = zeros(nb_scale,nb_samps*length(data.results));
            truth = zeros(nb_scale,nb_samps*length(data.results));
        end        
        for rr = 1:length(data.results)           
            for ss = 1:nb_scale
                truth(ss,num_e) = sum(sum(((data.results(rr).mask_truth{ss})*(data.results(rr).mask_truth{ss})')>0))/length(data.results(rr).test_q{ss}(:));
                perc_disc(ss,num_e) = sum(data.results(rr).test_q{ss}(:))/length(data.results(rr).test_q{ss}(:));                
            end            
            num_e = num_e+1;
        end
    end    
    perc_disc = sort(perc_disc,2);    
    hold on    
    %jbfill(list_scale(:)',perc_disc(:,floor(size(perc_disc,2)*0.1)),perc_disc(:,floor(size(perc_disc,2)*0.9)),[1 0 0],[1 0 0],0,0.25);
    %jbfill(list_scale(:)',perc_disc(:,floor(size(perc_disc,2)*0.9))',perc_disc(:,floor(size(perc_disc,2)*0.1))',[0 0 1],[0 0 0.25],0,0.25);
    %plot(list_scale,perc_disc(:,floor(size(perc_disc,2)*0.1)),'r');
    %plot(list_scale,perc_disc(:,floor(size(perc_disc,2)*0.9)),'r');
    %plot(list_scale,perc_disc(:,floor(size(perc_disc,2)/2)),'-');
    plot(list_scale,truth(:,1),'b');    
    axis([0 max(list_scale)+1 0 0.3])    
    name_fig = sprintf('simu_alpha2x%i_nb_subject%i_sc%i',ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref);
    title(strrep(name_fig,'_',' '))
    print(hf,[path_fig name_fig '.pdf'],'-dpdf')
    close(hf)
end
end
end
        
