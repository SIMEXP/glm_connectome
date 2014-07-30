function results_simu = glmc_simu_multiscale(tseries_all,hier,opt,flag_init)
% Simulate multiscale changes in connectivity between two populations
%
% SYNTAX:
% RESULTS = GLMC_SIMU_MULTISCALE(TSERIES,HIER,OPT)
%
% _________________________________________________________________________
% INPUTS
%
% TSERIES_ALL
%   (cell of strings) each entry is a variable TSERIES of size T x N. 
%
% HIER
%   (string) a .mat file with a hierarchy on the N regions of the brain.
%
% OPT           
%   (structure) with the following fields.  
%
%   TYPE_DATA
%      (string, default 'real') the type of simulation. Available options:
%      'real': use real datasets for the null hypothesis. 
%      'homogeneous': use an homogeneous correlation matrix under the null. 
%         in this case OPT.{THETA,STD_THETA,T,N} below have to be specified.
%
%   ALPHA2
%      (scalar, default 0.2) the variance of the signal added to
%      simulate an effect.
%
%   THETA2
%      (scalar, default 0.1) the variance of the signal added to create
%      a background 'homogeneous' correlation.
%
%   LIST_SCALES
%      (vector of integers) the list of scales used for testing. 
% 
%   SCALE_REF
%      (integer) the scale that is used as a reference for the simulation.
%      if left empty, or equal 0, the simu is done under the null of no effect 
%      at any scale.
%
%   CLUSTER_REF
%      (integer) the cluster that is used to simulate the effect at the scale 
%      of reference.
%
%   PERC_RAND
%      (scalar, default 0) a percentage PERC_RAND of the cluster of reference is 
%      randomly assigned to other clusters.
%
%   TYPE_FDR
%      (string, default 'LSL_sym') how the FDR is controled. 
%      See the TYPE argument of NIAK_GLM_FDR.
%
%   FDR
%      (scalar, default 0.05) the level of acceptable false-discovery rate 
%      for the t-maps.
%
%   NB_SUBJECT
%      (integer) the number of subjects to use.
%
%   NB_SAMPS
%      (integer, default 0) the number of samples for the permutation-based
%      omnibus test. If the number of samples is 0, the test is not 
%      generated. 
%
%   FLAG_VERBOSE 
%      (boolean, default 1) if the flag is 1, then the function 
%      prints some infos during the processing.
%
% _________________________________________________________________________
% OUTPUTS:
%
% RESULTS
%   (structure) with the following fields:
%   TEST_Q,P_NB_DISC,MASK_TRUTH,MEAN_DIFF,STD_IIR
%         
% OPT
%   Same as input OPT, but updated.
%
% _________________________________________________________________________
% SEE ALSO:
%
% _________________________________________________________________________
% COMMENTS
%
% Copyright (c) Pierre Bellec, 
% Centre de recherche de l'institut de gériatrie de Montréal, 2013.
% Maintainer : pierre.bellec@criugm.qc.ca
% See licensing information in the code.
% Keywords : GLM-connectome, simulation

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.

%% Options
if (nargin < 4)||flag_init
    list_fields   = { 'perc_rand' , 'type_data' , 'theta2' , 'nb_samps' , 'fdr' , 'type_fdr' , 'nb_subject' , 'alpha2' , 'list_scales' , 'scale_ref' , 'cluster_ref' , 'flag_verbose' };
    list_defaults = { 0           , 'real'      , 0.1      , 0          , 0.05  , 'LSL_sym'  , NaN          , 0.2      , NaN           , NaN         , NaN           , true           };
    opt = psom_struct_defaults(opt,list_fields,list_defaults);
end

if opt.scale_ref == 0
    opt.scale_ref = [];
end
        
%% Read partitions
part = niak_threshold_hierarchy(hier,struct('thresh',opt.list_scales));
if ~isempty(opt.scale_ref)
    mask_ref = opt.list_scales == opt.scale_ref;
    if ~any(mask_ref)&&~isempty(opt.scale_ref)
        error('I could not find the scale of reference in the list of scales !!')
    end
    cluster_ref = opt.cluster_ref;
    part_ref = part(:,mask_ref) == cluster_ref;
    part_ref_raw = part_ref;
    if opt.perc_rand>0
        part_ref_r = zeros(size(part_ref));
        ind_c = find(part_ref);
        %ind_c = ind_c(randperm(length(ind_c)));        % Use a deterministic discrepancy
        ind_nc = find(~part_ref);
        %ind_nc = ind_nc(randperm(length(ind_nc)));     % Use a deterministic discrepancy   
        nb_r = floor(length(ind_c)*opt.perc_rand);
        part_ref_r(ind_c(1:(end-nb_r))) = 1;
        part_ref_r(ind_nc(1:nb_r)) = 1;
        part_ref = part_ref_r>0;        
    end
    part_ref_raw = part_ref & part_ref_raw;
else
    warning('No scale of reference was specified, simulations are generated under the global null');
    part_ref = [];
    cluster_ref = 0;
    mask_ref = false(size(opt.list_scales));
end

%% Generate a mask of truth
for ss = 1:length(opt.list_scales)
    mask_truth{ss} = false(opt.list_scales(ss),1);
    if ~isempty(opt.scale_ref)
        for cc = 1:opt.list_scales(ss)
            mask_truth{ss}(cc) = any(part_ref_raw(part(:,ss)==cc));
        end
    end
end

%% Keep a random subset of subjects
ind_s = randperm(length(tseries_all));
ind_s = ind_s(1:opt.nb_subject);
tseries_all = tseries_all(ind_s);

%% Define a random group difference
ind1 = randperm(length(tseries_all));
ind1 = ind1(1:floor(length(tseries_all)/2));

%% Read time series and generate multiscale connectomes
if opt.flag_verbose
    fprintf('Generating (multiscale) connectomes ...\n')
end
conn = cell(length(opt.list_scales),1);
for tt = 1:length(tseries_all)
    niak_progress(tt,length(tseries_all));
    tseries = tseries_all{tt};
    
    switch opt.type_data        
        case 'real'
            %% run a circular block bootstrap on the time series inside 
            %% each network at the scale of reference
            %% to kill inter-network correlations at scales lower than the ref. 
            %% The objective is to be able to tell with confidence which connections
            %% are impacted by the addition of a signal in the simulations that follow
            if ~isempty(opt.scale_ref)
                for cc = 1:opt.scale_ref
                    tseries(:,part(:,mask_ref)==cc) = niak_bootstrap_tseries(tseries(:,part(:,mask_ref)==cc));
                end
            end
        case 'homogeneous'
            [t,n] = size(tseries);
            tseries = randn([t n]);
            w = opt.theta;
            tseries = sqrt(1-w)*tseries + sqrt(w)*repmat(randn([t 1]),[1 n]);        
        otherwise
             error('%s is an unkown type of simulations',opt.type_data);
    end
    tseries = niak_normalize_tseries(tseries);
    
    if any(ind1==tt)&&~isempty(opt.scale_ref)
        % simulate an effect !
        tseries_raw = tseries;
        tseries(:,part_ref) = sqrt(1-opt.alpha2)*tseries(:,part_ref) + sqrt(opt.alpha2)*repmat(randn([size(tseries,1) 1]),[1 sum(part_ref)]);        
    end        
        
    %% Measure the empirical effect size    
    for ss = 1:length(opt.list_scales)
        iir = niak_build_intra_inter_r(tseries,part(:,ss),true,true);
        if (tt ==1)
            conn{ss} = zeros(length(tseries_all),length(iir));
            mean_diff{ss} = zeros(sum(mask_truth{ss}),sum(mask_truth{ss}));
            sum_iir{ss} = zeros(sum(mask_truth{ss}),sum(mask_truth{ss}));
            sum_iir2{ss} = zeros(sum(mask_truth{ss}),sum(mask_truth{ss}));
        end
        conn{ss}(tt,:) = iir';
        if any(ind1==tt)&&~isempty(opt.scale_ref)
            iir_raw = niak_build_intra_inter_r(tseries_raw,part(:,ss),true,true);
            tmp = niak_lvec2mat(iir) - niak_lvec2mat(iir_raw);
            mean_diff{ss} = mean_diff{ss} + tmp(mask_truth{ss},mask_truth{ss});
            tmp_iir = niak_lvec2mat(iir_raw);
            sum_iir{ss} = sum_iir{ss} + tmp_iir(mask_truth{ss},mask_truth{ss});
            sum_iir2{ss} = sum_iir2{ss} + tmp_iir(mask_truth{ss},mask_truth{ss}).^2;
        else
            tmp_iir = niak_lvec2mat(iir);
            sum_iir{ss} = sum_iir{ss} + tmp_iir(mask_truth{ss},mask_truth{ss});
            sum_iir2{ss} = sum_iir2{ss} + tmp_iir(mask_truth{ss},mask_truth{ss}).^2;
        end
    end
end
if ~isempty(opt.scale_ref)
    for ss = 1:length(opt.list_scales)
        mean_diff{ss} = mean_diff{ss}/length(ind1);        
    end    
end
for ss = 1:length(opt.list_scales)
    std_iir{ss} = sqrt(length(tseries_all)/(length(tseries_all)-1))*sqrt(sum_iir2{ss}/length(tseries_all) - (sum_iir{ss}/length(tseries_all)).^2);
end

%% Generate stats
if opt.flag_verbose
    fprintf('Generate stats ...\n    scale')
end
mask1 = false(length(tseries_all),1);
mask1(ind1) = true;
opt_glm.test  = 'ttest' ;
opt_glm.flag_beta = true ; 
opt_glm.flag_residuals = true ;
q = opt.fdr;
nb_discovery = zeros(length(opt.list_scales),1);
perc_discovery = zeros(length(opt.list_scales),1);
test_q = cell(length(opt.list_scales),1);
vol_disc = 0;
for ss = 1:length(opt.list_scales)
    if opt.flag_verbose
        fprintf(' - %i',opt.list_scales(ss))
    end
    glm(ss).x = [ones(length(tseries_all),1) double(mask1)];
    glm(ss).c = [0 ; 1];
    glm(ss).y = conn{ss};
    results = niak_glm( glm(ss) , opt_glm );
    [fdr,test_q{ss}] = niak_glm_fdr(results.pce,opt.type_fdr,q,'correlation');
    nb_discovery(ss) = sum(sum(test_q{ss},1));
    perc_discovery(ss) = mean(nb_discovery(ss)/size(fdr,1));
    if any(niak_mat2lvec(test_q{ss}));
        vol_disc = vol_disc + sum(results.ttest(niak_mat2lvec(test_q{ss})).^2);
    else
        vol_disc = vol_disc + max(abs(results.ttest)).^2;
    end
end
nb_disc = sum(nb_discovery);

if opt.flag_verbose
    fprintf('\n')
end

%% Generate samples under the null 
if opt.nb_samps>0
    if opt.flag_verbose
        fprintf('Estimate the significance of the number of findings ...\n')
    end
    p_nb_disc = 0;
    nb_disc_null = zeros([opt.nb_samps 1]);
    vol_disc_null = zeros([opt.nb_samps 1]);    
    for num_s = 1:opt.nb_samps
        if opt.flag_verbose
            niak_progress(num_s,opt.nb_samps);
        end
        glm_null = niak_permutation_glm(glm);        
        for num_e = 1:length(glm_null)
            res_null = niak_glm(glm_null(num_e),opt_glm);
            [fdr_null,test_null] = niak_glm_fdr(res_null.pce,opt.type_fdr,opt.fdr,'correlation');
            nb_disc_null(num_s) = nb_disc_null(num_s) + double(sum(test_null(:)));
            if any(niak_mat2lvec(test_null));
                vol_disc_null(num_s) = vol_disc_null(num_s) + sum(res_null.ttest(niak_mat2lvec(test_null)).^2);
            else
                vol_disc_null(num_s) = vol_disc_null(num_s) + max(abs(res_null.ttest)).^2;
            end
        end
        %p_nb_disc = p_nb_disc + double(nb_disc_null(num_s)>=nb_disc);
        p_nb_disc = p_nb_disc + double(vol_disc_null(num_s)>=vol_disc);
    end
    p_nb_disc = p_nb_disc / opt.nb_samps;
else
    p_nb_disc = NaN;
    nb_disc_null = NaN;    
end

%% Save the results
results_simu.test_q = test_q;
results_simu.p_nb_disc = p_nb_disc;
results_simu.mask_truth = mask_truth;
results_simu.mean_diff = mean_diff;
results_simu.std_iir = std_iir;