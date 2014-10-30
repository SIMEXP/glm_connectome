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
list_M = [0 0.01 0.05 0.1]; % Number of true positive per family
list_theta = [2 3 5]; % The effect size
nb_samps = 10000; % Numer of samples
nb_perm = 10000; % Number of permutation samples for the omnibus tests
list_N = [7 16 25 55 114 199 328]; % Number of tests per family

%% Colors for the sensitivity plot
list_color = {'b','g','r','k'};

%% Get an estimate of excepted volume of discoveries under the null
fprintf('Estimating the volume of discoveries under the null...')
perc_null = zeros(length(list_fdr),nb_perm);
for ff = 1:length(list_fdr)
    niak_progress(ff,length(list_fdr));
    q = list_fdr(ff);
    for nn = 1:length(list_N)
        N = list_N(nn);
        z_null = randn(N,nb_perm);
        pce_null = normcdf(z_null);    
        [fdr_null,test_null] = niak_fdr(pce_null,'BH',q);        
        perc_null(ff,:) = perc_null(ff,:) + (squeeze(sum(test_null,1)/N));
    end
    perc_null(ff,:) = perc_null(ff,:)/length(list_N);
end
perc_null = sort(perc_null,2);

%% Loop over parameters
for mm = 1:length(list_M)
    for tt = 1:length(list_theta)
        
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
        perc_global = zeros(length(list_fdr),nb_samps);                 % percentage of discovery averaged across families per FDR threshold/sample
        sens = zeros(length(list_fdr),length(list_N),length(list_pce)); % Sensitivity per FDR threshold/family
        fdr  = zeros(length(list_fdr),length(list_N));                  % FDR per FDR threshold/family
        fp   = zeros(length(list_fdr),nb_samps,length(list_N)); % # of false positives per FDR threshold/sample/family
        disc = zeros(length(list_fdr),nb_samps,length(list_N)); % # of discoveries per FDR threshold/sample/family
        tp   = zeros(length(list_fdr),nb_samps,length(list_N)); % # of true positives per FDR threshold/sample/family
        perc = zeros(length(list_fdr),nb_samps,length(list_N)); % percentage of discovery per FDR threshold/sample/family
        
        for ff = 1:length(list_fdr) % Loop over FDR thresholds
            niak_progress(ff,length(list_fdr));
            q = list_fdr(ff); % Select the FDR threshold
            
            %% Compute the FDR tests, and derive the number of false/true positives, along with the number and percentage of discoveries
            for nn = 1:length(list_N) % Loop over families
                [fdr_tmp,test] = niak_fdr(pce{nn},'BH',q);
                fp(ff,:,nn)   = sum(test&~mask_true{nn},1); % sum the false positives 
                tp(ff,:,nn)   = sum(test&mask_true{nn},1);  % sum the true positives 
                disc(ff,:,nn) = fp(ff,:,nn) + tp(ff,:,nn);  % # of discoveries is the # false positives + # true positives
                perc(ff,:,nn) = disc(ff,:,nn)/list_N(nn);  % percentage of discoveries is the number of discoveries divided by the number of tests
            end
            
            %% Compute the omnibus test
            % i.e. counting the percentage of times a larger percentage of discoveries was observed under the null than in the actual sample. 
            % That process is repeated for each sample.
            perc_global(ff,:) = mean(perc(ff,:,:),3); % average the percentage of discoveries across families
            p_omn = zeros(1,nb_samps);
            for ss = 1:nb_samps
                p_omn(ss) = sum(perc_null(ff,:)>=perc_global(ff,ss))/nb_perm;
            end
            
            %% Compute the familywise FDR, global FDR, and familywise sensitivity
            %% as a function of the threshold on the omnibus test
            for pp = 1:length(list_pce) % Loop over omnibus thresholds          
                thre_pce = list_pce(pp); % Select a threshold on the omnibus test
                
                %% Estimate the global FDR
                tmp = sum(fp(ff,:,:),3)./sum(disc(ff,:,:),3); % Divide the sum of false positives across families by the sum of discoveries across families
                tmp(isnan(tmp))=0; % Take care of the NaN, that appear when there is no discovery (the false discovery proportion is by convention 0).
                tmp(p_omn>thre_pce) = 0; % For samples where the omnibus test is non-significant, there is no discovery. The false discovery proportion is by convention 0.
                fdr_global(ff,pp) = mean(tmp,2); % Simply average the false discovery proportions across samples to get an estimate of the effective FDR.
                
                %% Estimate the family-wise FDR
                tmp = fp(ff,:,:)./disc(ff,:,:); % Divide the false positives per family by the sum of discoveries per family
                tmp(isnan(tmp))=0; % Take care of the NaN, that appear when there is no discovery (the false discovery proportion is by convention 0).
                fdr(ff,:) = squeeze(mean(tmp,2)); % Simply average the false discovery proportions across samples to get an estimate of the effective FDR.
                
                %% Estimate the sensitivity per family
                for nn = 1:length(list_N)
                    tmp = tp(ff,:,nn)./all_true(nn,:);
                    tmp(p_omn>thre_pce) = 0;
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
        for pp = 1:length(list_pce)
            % The FDR figure
            name_fig = sprintf('omnibus%i_sc%i_eff%i',ceil(100*list_pce(pp)),ceil(100*list_M(mm)),list_theta(tt));
            hf = figure;
            hold on
            plot([0; max(list_fdr)],[0; max(list_fdr)],'k')            
            plot(list_fdr,fdr,'b')
            plot(list_fdr,fdr_global(:,pp),'r')
            axis([0 max(list_fdr) 0 0.2]);
            xlabel('Nominal FDR')
            ylabel('Effective FDR')
            title(strrep(name_fig,'_',' '))
            print(hf,[path_res filesep 'fdr_' name_fig '.pdf'],'-dpdf')
            close(hf)
            
            % the sensitivity figure
            if list_M(mm)~=0
                hf = figure;
                hold on    
                for ff = 1:length(list_fdr)
                    plot(list_N,sens(ff,:,pp)',list_color{ff})
                end
                legend(cellstr(num2str(list_fdr'))')
                axis([0 max(list_N) 0 1]);
                xlabel('Scale')
                ylabel('Sensitivity')
                title(strrep(name_fig,'_',' '))
                print(hf,[path_res filesep 'sens_' name_fig '.pdf'],'-dpdf')
                close(hf)
            end
        end
    end
end