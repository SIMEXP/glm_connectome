% Simulations of the scalability of FDR across multiple families for independent tests
% Warning: this script will clear the workspace
% and generate several figures in the current directory
% Copyright Pierre Bellec
% Centre de recherche de l'institut de gériatrie de Montréal, 
% Department of Computer Science and Operations Research
% University of Montreal, Québec, Canada, 2010-2014
% Maintainer : pierre.bellec@criugm.qc.ca
% See LICENSE.md for licensing information.

clear all
close all
path_res = pwd;

%% Make the simulations perfectly reproducible
psom_set_rand_seed(0);

%% Simulation parameters
list_pce = [0.05 Inf]; % The significance level for the omnibus test
list_fdr = [0.01 0.05 0.1 0.2]; % The FDR thresholds
list_M = [0 0.01 0.02 0.05 0.1]; % Number of true positive per family
list_theta = [2 3 5]; % The effect size
nb_samps = 100; % Numer of samples
nb_perm = nb_samps; % Number of permutation samples for the omnibus tests
%thre_perc = 0.005; % The threshold on acceptable percentage of discovery (0.5%)
thre_perc = 0.0005; % The threshold on acceptable percentage of discovery (0.5%)

%  % A flat number of tests per family (1000)
%  % To test scalability re number of families K (here K=2)
%  name_xp = 'variablek_2';
%  list_sc = [1000 1000];
%  list_N = list_sc;

% A flat number of tests per family (1000)
% To test scalability re number of families K (here K=5)
%  name_xp = 'variablek_5';
%  list_sc = [1000 1000 1000 1000 1000];
%  list_N = list_sc;

%  % A flat number of tests per family (1000)
%  % To test scalability re number of families K (here K=10)
%  name_xp = 'variablek_10';
%  list_sc = [1000 1000 1000 1000 1000 1000 1000 1000 1000 1000];
%  list_N = list_sc;

%  % A flat number of tests per family (1000)
%  % To test scalability re number of families K (here K=10)
%  name_xp = 'variablek_30';
%  list_sc = repmat(1000,[1 30]);
%  list_N = list_sc;

%  % A fixed number of families (K=5)
%  % To test scalability re number of tests per family 
%  lk = 10000;
%  name_xp = sprintf('variablel_%i',lk);
%  list_sc = repmat(lk,[1 5]);
%  list_N = list_sc;

%  % A fixed number of families (K=5)
%  % To test scalability re number of tests per family (here L=100 for all families)
%  name_xp = 'variablel_100';
%  list_sc = [100 100 100 100 100];
%  list_N = list_sc;

%  % A fixed number of families (K=5)
%  % To test scalability re number of tests per family (here L=10000 for all families)
%  name_xp = 'variablel_10000';
%  list_sc = [10000 10000 10000 10000 10000];
%  list_N = list_sc;

%  % Multiscale analysis with regular grid
%  name_xp = 'multiscale_50_step_10';
%  list_sc = [10:10:50];
%  list_N = list_sc.*(list_sc+1)/2;

% Multiscale analysis with regular grid
name_xp = 'multiscale_100_step_10';
list_sc = [10:10:100];
list_N = list_sc.*(list_sc+1)/2;

%  % Multiscale analysis with regular grid
%  name_xp = 'multiscale_200_step_10';
%  list_sc = [10:10:200];
%  list_N = list_sc.*(list_sc+1)/2;

%  % Multiscale analysis with regular grid
%  name_xp = 'multiscale_300_step_10';
%  list_sc = [10:10:300];
%  list_N = list_sc.*(list_sc+1)/2;

%  % Same parameters as the simulation with dependent tests
%  % i.e. replicate the number of tests in multiscale connectomes
%  % using a subset of the scales from the real SCHIZO dataset
%  name_xp = 'multiscale_msteps_schizo';
%  list_sc = [7 16 25 55 114 199 328];
%  list_N = list_sc.*(list_sc+1)/2;

%% Colors for the sensitivity plot
list_color = {'b','g','r','k'};

%% Get an estimate of excepted volume of discoveries under the null
fprintf('Estimating the volume of discoveries under the null...')
disc_null = zeros(length(list_fdr),nb_perm);
for ff = 1:length(list_fdr)
    niak_progress(ff,length(list_fdr));
    q = list_fdr(ff);
    for nn = 1:length(list_N)
        N = list_N(nn);
        z_null = randn(N,nb_perm);
        pce_null = normcdf(z_null);    
        [fdr_null,test_null] = niak_fdr(pce_null,'BH',q);        
        disc_null(ff,:) = disc_null(ff,:) + (squeeze(sum(test_null,1)));
    end
end
disc_null = sort(disc_null,2);

%% Create figure handles and plot indices
hfdr = figure; % Create a figure for the FDR plots
for pp = 1:length(list_pce)
    hsens(pp) = figure; % Create figures for the sensitivity plots (one per level of omnibus test)
end
num_plot_fdr = 1; % keep track of the number of the current FDR plot
num_plot_sens = 1; % keep track of the number of the current sensitivity plot

%% Loop over parameters
for tt = 1:length(list_theta)
    for mm = 1:length(list_M)        
        fprintf('Simulation perc of true pos %1.3f perc, effect size %1.f. ',list_M(mm),list_theta(tt));
        %% Run simulations
        pce = cell(length(list_N),1);
        mask_true = cell(length(list_N),1);
        all_true = zeros(length(list_N),nb_samps);
        for nn = 1:length(list_N)
            N = list_N(nn);
            theta = zeros(N,nb_samps);
            for ss = 1:nb_samps
                M0 = floor(list_M(mm)*N);
                R = list_M(mm)*N - M0;                
                M = M0 + (rand(1,nb_samps)<=R);
                theta(:,ss) = [ list_theta(tt)*ones(M(ss),1) ; zeros(N-M(ss),1) ];
            end
            mask_true{nn} = theta>0;
            pce{nn} = normcdf(-theta + randn(size(theta)));
            all_true(nn,:) = sum(mask_true{nn},1);
        end
        
        fdr_global  = zeros(length(list_fdr),length(list_pce));         % FDR pooled across families per FDR threshold/omnibus threshold
        q_global  = zeros(length(list_fdr),nb_samps);         % FDR pooled across families per FDR threshold/omnibus threshold      
        sens = zeros(length(list_fdr),length(list_N),length(list_pce)); % Sensitivity per FDR threshold/family
        fdr  = zeros(length(list_fdr),length(list_N));                  % FDR per FDR threshold/family
        fp   = zeros(length(list_fdr),nb_samps,length(list_N)); % # of false positives per FDR threshold/sample/family
        disc = zeros(length(list_fdr),nb_samps,length(list_N)); % # of discoveries per FDR threshold/sample/family
        tp   = zeros(length(list_fdr),nb_samps,length(list_N)); % # of true positives per FDR threshold/sample/family
        perc = zeros(length(list_fdr),nb_samps,length(list_N)); % percentage of discovery per FDR threshold/sample/family
        
        for ff = 1:length(list_fdr) % Loop over FDR thresholds
            niak_progress(ff,length(list_fdr));
            q = list_fdr(ff); % Select the FDR threshold            
            
            %% Compute the global FDR tests
            ind_n = [];
            pce_gb = [];
            for nn = 1:length(list_N)
               ind_n = [ind_n ; repmat(nn,[size(pce{nn},1) 1])];
               pce_gb = [pce_gb ; pce{nn}];              
            end
            [fdr_gb,test_gb] = niak_fdr(pce_gb,'BH',q);
            for nn = 1:length(list_N)
                test_gb_nn{nn} = test_gb(ind_n==nn,:);
            end
            
            %% Compute the FDR tests, and derive the number of false/true positives, along with the number and percentage of discoveries
            for nn = 1:length(list_N) % Loop over families
                [fdr_tmp,test] = niak_fdr(pce{nn},'BH',q); 
                test = test & test_gb_nn{nn};
                fp(ff,:,nn)   = sum(test&~mask_true{nn},1); % sum the false positives 
                tp(ff,:,nn)   = sum(test&mask_true{nn},1);  % sum the true positives 
                disc(ff,:,nn) = fp(ff,:,nn) + tp(ff,:,nn);  % # of discoveries is the # false positives + # true positives
            end
            
%              %% Compute the global FDR estimate
%              disc_tot = sum(squeeze(disc(ff,:,:)),2);
%              mask_disc = disc_tot>0;
%              df = repmat(disc_null(ff,:),[nb_samps 1]);
%              dt = repmat(disc_tot,[1 nb_perm]);
%              tmp = df./dt;
%              tmp(~mask_disc) = 0;
%              q_global(ff,:) = mean(tmp,2)';
%              
            %% Compute the familywise FDR, global FDR, and familywise sensitivity
            %% as a function of the threshold on the omnibus test
            for pp = 1:length(list_pce) % Loop over omnibus thresholds   
                if list_pce(pp)==0
                    thre_pce = list_fdr(ff); % Select a threshold on the omnibus test
                else
                    thre_pce = list_pce(pp); 
                end
                
                %% Estimate the global FDR
                tmp = sum(fp(ff,:,:),3)./sum(disc(ff,:,:),3); % Divide the sum of false positives across families by the sum of discoveries across families
                tmp(isnan(tmp))=0; % Take care of the NaN, that appear when there is no discovery (the false discovery proportion is by convention 0).
                tmp(q_global(ff,:)>thre_pce) = 0; % For samples where the omnibus test is non-significant, there is no discovery. The false discovery proportion is by convention 0.
                fdr_global(ff,pp) = mean(tmp,2); % Simply average the false discovery proportions across samples to get an estimate of the effective FDR.
                
                %% Estimate the family-wise FDR
                tmp = fp(ff,:,:)./disc(ff,:,:); % Divide the false positives per family by the sum of discoveries per family
                tmp(isnan(tmp))=0; % Take care of the NaN, that appear when there is no discovery (the false discovery proportion is by convention 0).
                fdr(ff,:) = squeeze(mean(tmp,2)); % Simply average the false discovery proportions across samples to get an estimate of the effective FDR.
                
                %% Estimate the sensitivity per family
                for nn = 1:length(list_N)
                    tmp = tp(ff,:,nn)./all_true(nn,:);
                    tmp(q_global(ff,:)>thre_pce) = 0;
                    tmp = tmp(all_true(nn,:)>0); % Cannot compute sensitivity on simulations with no true positives, which happen at low scales
                    if ~isempty(tmp)
                        sens(ff,nn,pp) = mean(tmp);
                    else
                        sens(ff,nn,pp) = NaN;
                    end
                end
            end
        end
        
        %% Now make figures
        
        % The FDR figure
        name_fig = sprintf('\\pi_1=%i%s',ceil(100*list_M(mm)),'%');
        
        figure(hfdr);
        hp = subplot(length(list_theta),length(list_M),num_plot_fdr);
        FN = findall(hp,'-property','FontName');
        set(FN,'FontName','/usr/share/fonts/truetype/dejavu/DejaVuSerifCondensed.ttf');
        FS = findall(hp,'-property','FontSize');
        set(FS,'FontSize',8);
        hold on
        plot([0; max(list_fdr)],[0; max(list_fdr)],'k','LineWidth',4)            
        plot(list_fdr,fdr,'b')
        plot(list_fdr,fdr_global(:,1),'r','LineWidth',4)
        plot(list_fdr,fdr_global(:,2),'g','LineWidth',4)
        axis([0 max(list_fdr) 0 0.2]);
        if tt == length(list_theta)
            xlabel('Nominal FDR');
        end
        if tt==1
            title(name_fig)
        end
        if mm==1
            ylabel(sprintf('Effective FDR (eff=%i)',list_theta(tt)))
        end        
        if (mm==length(list_M))&&(tt==length(list_theta))
            print(hfdr,[path_res filesep 'fig_fdr_' name_xp '.png'],'-dpng')
        end        
            
        % the sensitivity figure
        for pp = 1:length(list_pce)
            figure(hsens(pp));
            name_fig = sprintf('\\pi_1=%i%s',ceil(100*list_M(mm)),'%');
            if list_M(mm)~=0
                hp = subplot(length(list_theta),length(list_M)-1,num_plot_sens);                
                hold on    
                FN = findall(hp,'-property','FontName');
                set(FN,'FontName','/usr/share/fonts/truetype/dejavu/DejaVuSerifCondensed.ttf');
                FS = findall(hp,'-property','FontSize');
                set(FS,'FontSize',8);
        
                for ff = 1:length(list_fdr)
                    plot(list_sc,sens(ff,:,pp)',list_color{ff},'LineWidth',4)
                end
                %legend(cellstr(num2str(list_fdr'))')
                axis([0 max(list_sc) 0 1]);
                if tt == length(list_theta)
                    xlabel('Scale')
                end
                if mm==2
                    ylabel(sprintf('Sensitivity (eff=%i)',list_theta(tt)))
                end
                if tt==1
                    title(name_fig)
                end
            end
            if (mm==length(list_M))&&(tt==length(list_theta))
                print(hsens(pp),[path_res filesep sprintf('sens_omnibus%i',ceil(100*list_pce(pp))) '_' name_xp '.png'],'-dpng')
            end
        end
        
        %% Increase plot numbers
        num_plot_fdr = num_plot_fdr+1;
        if list_M(mm)~=0
            num_plot_sens = num_plot_sens+1;
        end
    end
end
%close(hfdr)
for pp = 1:length(list_pce)
   %close(hsens(pp))
end