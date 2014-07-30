function [in,out,opt] = glmc_brick_glm2mdmr(in,out,opt)
% Generate connectome-wide association study using the MDMR method
%
% SYNTAX: [IN,OUT,OPT] = GLMC_BRICK_GLM2MDMR(IN,OUT,OPT)
%
% IN.FMRI.(SUBJECT).(SESSION).(RUN) (string) a 3D+t fMRI dataset.
% IN.GLM (string) a XXX_glm.mat structure. See NIAK_PIPELINE_GLM_CONNECTOME pipeline.
% IN.NETWORK (string) the name of a 3D volume with the networks used for the GLM analysis.
%
% OUT.MDMR (string) the file name of a 4D dataset. MDMR(:,:,:,n) is the t-stat map 
%   corresponding to network n. The significance of connections is derived using 
%   the MDMR method.
% OUT.PERC_DISCOVERY (string) the file name of a 3D volume, with the map of the 
%   number of discovery associated with each network, expressed as a percentage of the 
%   number of networks - 1 (i.e. the max possible number of discoveries associated 
%   with a network).
%
% OPT.FDR (scalar, default 0.05) the FDR within each connectivity map.
% OPT.P_SEED (scalar, default 0.05) the alpha level for making omnibus test for the 
%   presence of an effect in a seed-based connectivity map.
% OPT.NB_PERMUTATION (integer, default 9999) the number of permutations.
% OPT.FLAG_TEST (boolean, default false) if the flag is true, the brick does 
%   nothing but update IN, OUT, OPT.
%
% COMMENTS:
%    MDMR does not support intra-subject contrast. It will generate the average 
%    connectivity based on concatenated data for each subject.
%
%  opt_g.min_nb_vol = 50;     % The minimum number of volumes for an fMRI dataset to be included. This option is useful when scrubbing is used, and the resulting time series may be too short.
%  opt_g.min_xcorr_func = 0; % The minimum xcorr score for an fMRI dataset to be included. This metric is a tool for quality control which assess the quality of non-linear coregistration of functional images in stereotaxic space. Manual inspection of the values during QC is necessary to properly set this threshold.
%  opt_g.min_xcorr_anat = 0; % The minimum xcorr score for an fMRI dataset to be included. This metric is a tool for quality control which assess the quality of non-linear coregistration of the anatomical image in stereotaxic space. Manual inspection of the values during QC is necessary to properly set this threshold.
%  opt_g.type_files = 'glm_connectome'; % Specify to the grabber to prepare the files for the glm_connectome pipeline
%  files_in.fmri = niak_grab_fmri_preprocess('/media/database5/pyeror/gui/cobre/fmri_preprocess_260213',opt_g).fmri; % Replace the folder by the path where the results of the fMRI preprocessing pipeline were stored. 
%  files_in.glm = '/media/database5/pyeror/gui/cobre/results/glm_connectome_1_130313/sci10_scg10_scf10/szVScont_c_age_sex_FD/glm_szVScont_c_age_sex_FD_sci10_scg10_scf10.mat';
%  files_in.network = '/media/database5/pyeror/gui/cobre/results/glm_connectome_1_130313/sci10_scg10_scf10/networks_sci10_scg10_scf10.mnc.gz';
%  files_out.mdmr = '/home/pbellec/tmp/test_glm2mdmr.mnc.gz';
%  files_out.perc_discovery = '/home/pbellec/tmp/test_glm2mdmr_perc.mnc.gz';
%  opt.nb_permutation = 10;
%  opt.p_seed = 0.05;
%  opt.fdr = 0.05;
%  glmc_brick_glm2mdmr(files_in,files_out,opt);
%
% Copyright (c) Pierre Bellec, 
% Centre de recherche de l'institut de gériatrie de Montréal, 2014.
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

%% Set default
in = psom_struct_defaults(in, ...
    { 'fmri' , 'glm' , 'network' }, ...
    { NaN    , NaN   , NaN       });
   
out = psom_struct_defaults(out, ...
    { 'mdmr' , 'perc_discovery' }, ...
    { NaN    , NaN              });

if nargin < 3
    opt = struct();
end
opt = psom_struct_defaults(opt, ...
    { 'fdr' , 'p_seed' , 'nb_permutation' , 'flag_test'}, ...
    { 0.05  , 0.05     , 9999             , false      });

if opt.flag_test
    return
end

%% Read the input model
glm_niak = load(in.glm,'model_group','pce');
pce = glm_niak.pce;
glm_niak = glm_niak.model_group;
ttest = load(in.glm,'ttest');
ttest = ttest.ttest;

%% Read the input partition
[hdr,mask] = niak_read_vol(in.network);
nb_net = max(mask(:));

%% Read time series
list_subject = glm_niak.labels_x;
[list_files,labels] = niak_fmri2cell(in.fmri);
fprintf('Loading individual subject data ...');
for ss = 1:length(list_subject)
    niak_progress(ss,length(list_subject));
    subject = list_subject{ss};
    mask_ss = ismember({labels.subject},{subject});
    if ~any(mask_ss)
        error('Could not find the fMRI data for subject %s',subject)
    end
    files_ss = list_files(mask_ss);
    for ff = 1:length(files_ss)
        [hdr_f,vol] = niak_read_vol(files_ss{ff});
        if ff == 1
            tseries = niak_build_tseries(vol,mask,struct('correction','mean_var'));
        else
            tseries = [tseries ; niak_build_tseries(vol,mask,struct('correction','mean_var'))];
        end
    end 
    glm_mdmr.tseries{ss} = tseries;
end

%% Build MDMR model
mask_keep = ~ismember(glm_niak.labels_y,'intercept');
glm_mdmr.x = glm_niak.x(:,mask_keep);
glm_mdmr.labels_x = glm_niak.labels_x;
glm_mdmr.labels_y = glm_niak.labels_y(mask_keep);
for yy = 1:length(glm_mdmr.labels_y)
    if yy == 1
        opt_mdmr.formula = glm_mdmr.labels_y{1};
    else
        opt_mdmr.formula = [opt_mdmr.formula ' + ' glm_mdmr.labels_y{yy}];
    end
end
opt_mdmr.nb_permutation = opt.nb_permutation;
opt_mdmr.factors2perm = glm_mdmr.labels_y{glm_niak.c(mask_keep)>0};

%% Run MDMR
pce_seed = glmc_cwas(glm_mdmr,opt_mdmr);
test_seed = pce_seed <= opt.p_seed;

%% Now build the results of the test
[fdr,test_q] = niak_glm_fdr(pce,'BH-local',opt.fdr);
test_q(:,~test_seed) = false;

%% Save the results: discovery maps
nb_discovery = sum(test_q,1);
perc_discovery = nb_discovery/size(test_q,1);
discovery_maps = niak_part2vol(perc_discovery,mask);
hdr.file_name = out.perc_discovery;
niak_write_vol(hdr,discovery_maps);

%% Save the results: NBS maps
mdmr_maps   = zeros([size(mask) nb_net]);
ttest_mat = niak_lvec2mat (ttest);
for num_net = 1:nb_net    
    ttest_thre = ttest_mat(:,num_net);
    ttest_thre( ~test_q(:,num_net) ) = 0;
    mdmr_maps(:,:,:,num_net) = niak_part2vol(ttest_thre,mask);
end
hdr.file_name = out.mdmr;
niak_write_vol(hdr,mdmr_maps);
