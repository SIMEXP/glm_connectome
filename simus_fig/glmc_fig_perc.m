%% Prepare the figures for the GLM-connectome simulation
%% percentage of discovery

clear

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
hf = figure; % Create a figure for the plots
num_plot = 1; % keep track of the number of the current plot

%% Build plots
for num_alpha2 = 1:length(param.alpha2) % Loop over effect size
    for num_nsub = 1:length(param.nsub) % Loop over the number of subjects
        for num_sc = 1:length(param.sc) % Loop over reference clusters
            % Load data
            name_xp = sprintf('simu_perc%i_sc%i_a2%i_nsub%i',ceil(100*param.perc_rand(num_p)),param.sc{num_sc},ceil(100*param.alpha2{num_alpha2}),param.nsub{num_nsub});                            
            fprintf('%s\n',name_xp)
            nb_samps = 0;
%            for num_s = 1:param.nb_samps % Loop over simulation jobs
            for num_s = 1 % Loop over simulation jobs
                niak_progress(num_s,param.nb_samps);
                name_job = sprintf('simu_a2%i_nsub%i_sc%i_perc%i_samp%i',ceil(100*param.alpha2{num_alpha2}),param.nsub{num_nsub},param.sc{num_sc},ceil(100*param.perc_rand(num_p)),num_s);                    
                file_data = [path_results filesep name_job '.mat'];
                data = load(file_data);
                perc = zeros(nb_scale,opt.nb_replication * nb_samps);
                for rr = 1:length(data.results)
                    for ss = 1:nb_scale
                        perc(ss,num_exp) = sum(data.results(rr).test_q{ss}(:))/length(data.results(rr).test_q{ss}(:));
                        num_exp = num_exp+1;
                    end
                end        
            end
            
            %% figure
            hold on
            perc = sort(perc,2);
            plot(list_scale,perc(:,round(length(perc)/2)),'-')
            plot(list_scale,perc(:,round(length(perc)*0.05)),'--')
            plot(list_scale,perc(:,round(length(perc)*0.95)),'--')
        end
    end
end
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
        