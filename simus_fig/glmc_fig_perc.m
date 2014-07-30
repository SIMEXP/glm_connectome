%% Prepare the figures for the GLM-connectome simulation
%% percentage of discovery

clear

%% Folders
path_results = '/home/pbellec/database/glm_connectome/simu_multiscale/';
path_fig = '/home/pbellec/database/glm_connectome/fig_perc/';
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
    opt.alpha2     = param.alpha2{num_alpha2};        
    opt.nb_subject = param.nsub{num_nsub};        
    opt.scale_ref  = param.sc{num_sc};
    num_exp = 1;
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
            perc = zeros(nb_scale,opt.nb_replication * nb_samps);
        end        
        for rr = 1:length(data.results)
            for ss = 1:nb_scale
                perc(ss,num_exp) = sum(data.results(rr).test_q{ss}(:))/length(data.results(rr).test_q{ss}(:));
                num_exp = num_exp+1;
            end
        end        
    end
    hold on
    perc = sort(perc,2);
    plot(list_scale,perc(:,round(length(perc)/2)),'-')
    plot(list_scale,perc(:,round(length(perc)*0.05)),'--')
    plot(list_scale,perc(:,round(length(perc)*0.95)),'--')
    %for ss = 1:nb_scale
    %    plot(0.5*randn(size(perc(3,:)))+list_scale(ss),perc(ss,:),'.')
    %end
    %errorbar(list_scale,mean(perc,2),std(perc,[],2))
    name_fig = sprintf('simu_alpha2x%i_nb_subject%i_sc%i',ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref);
    title(strrep(name_fig,'_',' '))
    print(hf,[path_fig name_fig '.png'],'-dpng')
    close(hf)
end
end
end
        