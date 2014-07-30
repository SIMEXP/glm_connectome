%% Prepare the figures for the GLM-connectome simulation
%% histogram of p-values

clear

%% Folders
path_results = '/home/pbellec/database/glm_connectome/simu_multiscale/';
path_fig = '/home/pbellec/database/glm_connectome/fig_pce/';
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
    nb_disc = zeros(nb_samps*opt.nb_replication,1);
    npce = 1;
    opt.alpha2     = param.alpha2{num_alpha2};        
    opt.nb_subject = param.nsub{num_nsub};        
    opt.scale_ref  = param.sc{num_sc};
    for num_s = 1:nb_samps        
        name_job = sprintf('simu_alpha2x%i_nb_subject%i_sc%i_samp%i',ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref,num_s);
        file_data = [path_results name_job '.mat'];
        data = load(file_data);                
        for rr = 1:length(data.results)
            pce(npce) = data.results(rr).p_nb_disc;  
            for ss = 1:length(data.results(rr).test_q)
                nb_disc(npce) = nb_disc(ss) + sum(data.results(rr).test_q{ss}(:));
            end
            npce = npce+1;
        end
    end
    hist(pce,100);        
    name_fig = sprintf('simu_alpha2x%i_nb_subject%i_sc%i',ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref);
    title(strrep(name_fig,'_',' '))
    print(hf,[path_fig name_fig '.png'],'-dpng')
    close(hf)
    if num_sc == 1
        fprintf('%s\n pr(p<0.1)  = %1.4f\n pr(p<0.05) = %1.4f\n pr(p<0.01) = %1.4f\n\n',name_fig,sum(pce<0.1)/length(pce),sum(pce<0.05)/length(pce),sum(pce<0.01)/length(pce));
    end
end
end
end
        