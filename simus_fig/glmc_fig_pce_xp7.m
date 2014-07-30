%% Prepare the figures for the GLM-connectome simulation
%% histogram of p-values

% SCALE 0 

%  xp7a_simu_alpha2x10_nb_subject40_sc0
%   pr(p<0.1)  = 0.0860
%   pr(p<0.05) = 0.0360
%   pr(p<0.01) = 0.0090
%  
%   xp7b_simu_alpha2x10_nb_subject40_sc0
%   pr(p<0.1)  = 0.0840
%   pr(p<0.05) = 0.0385
%   pr(p<0.01) = 0.0090
%  
%   xp7a_simu_alpha2x20_nb_subject40_sc0
%   pr(p<0.1)  = 0.0990
%   pr(p<0.05) = 0.0470
%   pr(p<0.01) = 0.0130
%  
%   xp7b_simu_alpha2x20_nb_subject40_sc0
%   pr(p<0.1)  = 0.0980
%   pr(p<0.05) = 0.0465
%   pr(p<0.01) = 0.0125
%  
%   xp7a_simu_alpha2x10_nb_subject100_sc0
%   pr(p<0.1)  = 0.0960
%   pr(p<0.05) = 0.0530
%   pr(p<0.01) = 0.0140
%  
%   xp7b_simu_alpha2x10_nb_subject100_sc0
%   pr(p<0.1)  = 0.0960
%   pr(p<0.05) = 0.0525
%   pr(p<0.01) = 0.0130
%  
%   xp7a_simu_alpha2x20_nb_subject100_sc0
%   pr(p<0.1)  = 0.0890
%   pr(p<0.05) = 0.0500
%   pr(p<0.01) = 0.0080
%  
%   xp7b_simu_alpha2x20_nb_subject100_sc0
%   pr(p<0.1)  = 0.0860
%   pr(p<0.05) = 0.0445
%   pr(p<0.01) = 0.0060

% SCALE 4

%  xp7a_simu_alpha2x10_nb_subject40_sc4
%   pr(p<0.1)  = 0.8180
%   pr(p<0.05) = 0.7420
%   pr(p<0.01) = 0.5670
%  
%   xp7b_simu_alpha2x10_nb_subject40_sc4
%   pr(p<0.1)  = 0.8280
%   pr(p<0.05) = 0.7455
%   pr(p<0.01) = 0.5530
%  
%   xp7a_simu_alpha2x20_nb_subject40_sc4
%   pr(p<0.1)  = 1.0000
%   pr(p<0.05) = 1.0000
%   pr(p<0.01) = 1.0000
%  
%   xp7b_simu_alpha2x20_nb_subject40_sc4
%   pr(p<0.1)  = 1.0000
%   pr(p<0.05) = 1.0000
%   pr(p<0.01) = 1.0000
%  
%   xp7a_simu_alpha2x10_nb_subject100_sc4
%   pr(p<0.1)  = 1.0000
%   pr(p<0.05) = 0.9990
%   pr(p<0.01) = 0.9960
%  
%   xp7b_simu_alpha2x10_nb_subject100_sc4
%   pr(p<0.1)  = 1.0000
%   pr(p<0.05) = 0.9995
%   pr(p<0.01) = 0.9930
%  
%   xp7a_simu_alpha2x20_nb_subject100_sc4
%   pr(p<0.1)  = 1.0000
%   pr(p<0.05) = 1.0000
%   pr(p<0.01) = 1.0000
%  
%   xp7b_simu_alpha2x20_nb_subject100_sc4
%   pr(p<0.1)  = 1.0000
%   pr(p<0.05) = 1.0000
%   pr(p<0.01) = 1.0000

 % SCALE 7 
 
%   xp7a_simu_alpha2x10_nb_subject40_sc7
%   pr(p<0.1)  = 0.3800
%   pr(p<0.05) = 0.2560
%   pr(p<0.01) = 0.1140
%  
%   xp7b_simu_alpha2x10_nb_subject40_sc7
%   pr(p<0.1)  = 0.3955
%   pr(p<0.05) = 0.2695
%   pr(p<0.01) = 0.1200
%  
%   xp7a_simu_alpha2x20_nb_subject40_sc7
%   pr(p<0.1)  = 0.9860
%   pr(p<0.05) = 0.9720
%   pr(p<0.01) = 0.9410
%  
%   xp7b_simu_alpha2x20_nb_subject40_sc7
%   pr(p<0.1)  = 0.9895
%   pr(p<0.05) = 0.9770
%   pr(p<0.01) = 0.9310
%  
%   xp7a_simu_alpha2x10_nb_subject100_sc7
%   pr(p<0.1)  = 0.8770
%   pr(p<0.05) = 0.8140
%   pr(p<0.01) = 0.6450
%  
%   xp7b_simu_alpha2x10_nb_subject100_sc7
%   pr(p<0.1)  = 0.8915
%   pr(p<0.05) = 0.8275
%   pr(p<0.01) = 0.6285
%  
%   xp7a_simu_alpha2x20_nb_subject100_sc7
%   pr(p<0.1)  = 1.0000
%   pr(p<0.05) = 1.0000
%   pr(p<0.01) = 1.0000
%  
%   xp7b_simu_alpha2x20_nb_subject100_sc7
%   pr(p<0.1)  = 1.0000
%   pr(p<0.05) = 1.0000
%   pr(p<0.01) = 1.0000


clear

%% Folders
list_xp = {'xp7a','xp7b'};
path_results = '/home/pbellec/database/glm_connectome/simu_multiscale_';
path_fig = '/home/pbellec/database/glm_connectome/figures/fig_pce_xp7/';
psom_mkdir(path_fig);

%% Parameters
nb_samps = 100;
%param.sc = {0,10,30,100};

%% At scale 3, cluster 3 is 69.5% of the grey matter, which makes for a percentage of true difference close to 50% (48%)
%% At scale 7, cluster 7 is 20.5% of the grey matter, which makes for a percentage of true difference close to 4% (4.2%)S
param.sc      = {0,4,7};
param.cluster = {0,4,7};

param.nsub = {40,100};
param.alpha2 = {0.1,0.2};
opt.nb_replication = 10;

bins = 0.025:0.05:1;
delta = bins(2)-bins(1);
%% Build figures
for num_sc = 1:length(param.sc)
    for num_nsub = 1:length(param.nsub)
        for num_alpha2 = 1:length(param.alpha2)
            clf
            % Load data
            pce = zeros(nb_samps*opt.nb_replication,1);
            nb_disc = zeros(nb_samps*opt.nb_replication,1);
            npce = 1;
            opt.alpha2     = param.alpha2{num_alpha2};        
            opt.nb_subject = param.nsub{num_nsub};        
            opt.scale_ref  = param.sc{num_sc};
            for num_xp = 1:length(list_xp)
                hf = figure;
                name_xp = list_xp{num_xp};
                fprintf('simu_alpha2x%i_nb_subject%i_sc%i\n',ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref)
                for num_s = 1:nb_samps        
                    niak_progress(num_s,nb_samps);            
                    name_job = sprintf('simu_alpha2x%i_nb_subject%i_sc%i_samp%i',ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref,num_s);
                    file_data = [path_results name_xp filesep name_job '.mat'];
                    data = load(file_data);                
                    for rr = 1:length(data.results)
                        pce(npce) = data.results(rr).p_nb_disc;  
                        for ss = 1:length(data.results(rr).test_q)
                            nb_disc(npce) = nb_disc(ss) + sum(data.results(rr).test_q{ss}(:));
                        end
                        npce = npce+1;
                    end
                end
                [Y,X] = hist(pce,bins);
                Y = Y/(delta*sum(Y));
                Y(Y>5) = 5;
                hb = bar(X,Y);
                set(hb,'facecolor','b')
                set(hb,'edgecolor','k')
                axis([-0.05 1.05 0 5]);
                ha = gca;   
                name_fig = sprintf('%s_simu_alpha2x%i_nb_subject%i_sc%i',name_xp,ceil(100*opt.alpha2),opt.nb_subject,opt.scale_ref);
                title(strrep(name_fig,'_',' '))
                print(hf,[path_fig name_fig '.pdf'],'-dpdf')
                close(hf)                
                fprintf('%s\n pr(p<0.1)  = %1.4f\n pr(p<0.05) = %1.4f\n pr(p<0.01) = %1.4f\n\n',name_fig,sum(pce<0.1)/length(pce),sum(pce<0.05)/length(pce),sum(pce<0.01)/length(pce));                
            end
        end
    end
end
        
