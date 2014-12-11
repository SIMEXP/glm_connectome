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
param.nb_samps = 1;
%% Create figure handles and plot indices
hf = figure; % Create a figure for the plots
num_plot = 1; % keep track of the number of the current plot

%% Build plots
num_p = 1;
for num_sc = 2:length(param.sc) % Loop over reference clusters
    for num_alpha2 = 1 % Loop over effect size
        for num_nsub = 1 % Do not loop over the number of subjects, the effect size does not depend on the number of subjects
            % Load data
            name_xp = sprintf('simu_perc%i_sc%i_a2%i_nsub%i',ceil(100*param.perc_rand(num_p)),param.sc{num_sc},ceil(100*param.alpha2{num_alpha2}),param.nsub{num_nsub});                            
            fprintf('%s\n',name_xp)
            nb_samps = 0;
            pi1 = zeros(nb_scale,1);
            num_s = 1;
            rr = 1;
            name_job = sprintf('simu_a2%i_nsub%i_sc%i_perc%i_samp%i',ceil(100*param.alpha2{num_alpha2}),param.nsub{num_nsub},param.sc{num_sc},ceil(100*param.perc_rand(num_p)),num_s);                    
            file_data = [path_results filesep name_job '.mat'];
            data = load(file_data);
            for ss = 1:nb_scale
                truth = data.results(rr).mask_truth{ss};
                truth = truth*(truth')>0; 
                pi1(ss) = mean(truth(:));                
            end
            
            %% figure
            hp = subplot(length(param.sc)-1,1,num_plot);               
            FN = findall(hp,'-property','FontName');
            set(FN,'FontName','/usr/share/fonts/truetype/dejavu/DejaVuSerifCondensed.ttf');
            FS = findall(hp,'-property','FontSize');
            set(FS,'FontSize',8);
            hold on
            switch num_sc
            case 1
                name_sc = 'True Positive Rate 0%';
            case 2
                name_sc = 'True Positive Rate \sim1%';
            case 3
                name_sc = 'True Positive Rate \sim2%';
            case 4
                name_sc = 'True Positive Rate \sim5%';
            case 5
                name_sc = 'True Positive Rate \sim10%';
            end
            if num_sc == 2
                title('Pi1')
            end
            if num_alpha2==1
                ylabel(name_sc)
            end
            if (num_nsub == length(param.nsub))&&(num_alpha2==length(param.alpha2))
                xlabel('Scale')
            end
            
            plot(param.list_scales',pi1,'bx-')
            axis([0 max(param.list_scales) 0.007 0.15]);
            ha = gca;
            set(ha,'ytick',[0.01 0.02 0.05 0.1])
            set(ha,'yticklabel',{'0.01','0.02','0.05','0.1','0.15'})
            set(ha,'ygrid','on')                     
            set(ha,'yscale','log')
            num_plot = num_plot+1;  
        end
    end
end
print(hf,[path_fig filesep 'fig_p1.png'],'-dpng')