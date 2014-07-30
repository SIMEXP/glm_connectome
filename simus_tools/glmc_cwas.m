function pce_seed = glmc_cwas(glm,opt)
% Connectome-wide association study
%
% SYNTAX:
%   RESULTS = GLMC_CWAS(GLM,OPT)
%
% INPUTS:
%   GLM.TSERIES (cell of arrays) TSERIES{s}(t,n) is the time series at time t, region
%      n and subject s. The number of regions has to be identical for all subjects.
%   GLM.X (array SxK) X(s,k) it the value of the k-th explanatory variable 
%      for subject s. 
%   GLM.LABELS_X (cell of strings, length S, default {'S1',...}) the subject labels.
%   GLM.LABELS_Y (cell of strings, length K) the labels of the covariates.
%
%   OPT.FORMULA (string) columns in your model file that you want to account/measure 
%      distances for (e.g. 'group + meanFD')
%   OPT.FACTORS2PERM (string) the factors given in your --formula that you want to 
%      permute to get significance of the association between the factor and distances
%   OPT.NB_PERMUTATION (integer, default 9999) the number of permutations.
%   OPT.PATH_WORK (string, default a temporary folder) where to write the temporary 
%      files. If specified, it will be kept after processing. Otherwise, a temp 
%      folder is automatically generated, and then deleted.
%   OPT.FORKS (integer, default 1) how many parallel jobs will be used to run.
%   OPT.THREADS (integer, default 1)  how many cores do you want to use in parallel.
%   OPT.MEMLIMITS (scalar, default 1) the maximum amount of gb of RAM to use
% 
% OUTPUTS:
%   PCE_SEED (vector length N) PCE_SEED(n) is the probability that the seed n
%      is associated with signal. 
%
% COMMENTS:
%   A simple wrapper around the R CWAS tools of Zarrar Shehzad.
%
% EXAMPLE:
%   opt_s.type = 'onescale';
%   opt_s.nb_clusters = 4;
%   opt_s.n = 100;
%   opt_s.t = 200;
%   opt_s.variance = [0.2 0.2 0.2 0.2];
%   for ss = 1:40
%       if ss == 21
%           opt_s.variance = 0.2:0.15:0.65;
%       end
%       glm.tseries{ss} = niak_simus_scenario(opt_s);
%   end
%   glm.x = [ones(40,1) [ zeros(20,1) ; ones(20,1) ]];
%   glm.labels_y = {'intercept','group'};
%   opt.formula = 'group';
%   opt.factors2perm = 'group';
%   results = glmc_cwas(glm,opt);
%   
% Copyright (c) Pierre Bellec, 
% Centre de recherche de l'institut de gériatrie de Montréal, 2014.
% Maintainer : pierre.bellec@criugm.qc.ca
% See licensing information in the code.
% Keywords : connectome-wide association study. fMRI.

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

%% Get inputs, and set defaults
glm = psom_struct_defaults(glm, ...
      { 'tseries' , 'x' , 'labels_x' , 'labels_y'}, ...
      { NaN       , NaN , {}         , NaN       });
opt = psom_struct_defaults(opt, ...
      { 'nb_permutation' , 'path_work' , 'formula' , 'factors2perm' , 'forks' , 'threads' , 'memlimits' }, ...
      { 9999             , ''          , NaN       , NaN            , 1       , 1         , 1           });

s = length(glm.tseries);
[s2,k]  = size(glm.x);
if s2~=s
    error('The second dimension of TSERIES should have the same length as the first dimension of X, i.e. the number of subjects');
end
if isempty(glm.labels_x)
    for ss = 1:s
        glm.labels_x{ss} = sprintf('S%i',ss);
    end
end
s3 = length(glm.labels_x);
if s3~=s
    error('The length of LABELS_X should be equal to the length of the first dimension of X, i.e. the number of subjects');
end

%% Deal with temporary folder
if isempty(opt.path_work)
    opt.path_work = psom_path_tmp('_cwas');
    flag_tmp = true;
else
    flag_tmp = false;
    opt.path_work = niak_full_path(opt.path_work);
end

%% Check the dimensions of the arrays
for ss = 1:s
    [t,n] = size(glm.tseries{ss});
    if ss == 1
        n1 = n;
    end
    if n~=n1
         error('The number of regions in the time series should be identical for all subjects');
    end
end

%% Fake nifti header
hdr.type = 'nii';
hdr.info.machine = 'native';
hdr.info.file_parent = '';
hdr.info.precision = 'float32';
hdr.info.voxel_size = [1 1 1];
hdr.info.mat = [eye(3) zeros(3,1) ; zeros(1,3) 1];
hdr.info.dimension_order = '';
hdr.info.history = '';

%% write out individual time series
file_all_tseries = [opt.path_work 'list_tseries.txt'];
hf_all = fopen(file_all_tseries,'w');
mask = true(n,1,1,1);
hdr.info.tr = 3;
for ss = 1:s    
    file_name = [opt.path_work glm.labels_x{ss} '.txt'];
    %tseries_ss = glm.tseries{ss};
    %save(file_name,'-ascii','tseries_ss');    
    sub_write_ascii(file_name,glm.tseries{ss});
    fprintf(hf_all,'%s\n',file_name);
end
fclose(hf_all);

%% Write out a fake nifti mask
file_mask = [opt.path_work 'group_mask.nii.gz'];
hdr.info.tr = 0;
hdr.file_name = file_mask;
niak_write_vol(hdr,mask);

%% Write covariates
file_model = [opt.path_work 'covariates_glm.csv'];
cell_model = [glm.labels_y(:)' ; num2cell(glm.x)];
niak_write_csv_cell(file_model,cell_model,',')

%% Call CWAS: first compute distance
folder_subdist = [opt.path_work 'subdist' filesep];
call_cwas = [ 'connectir_subdist.R ' ...              
              ' --infuncs1 ' file_all_tseries ...
              ' --brainmask1 ' file_mask ...
              ' --ztransform' ...
              ' --bg ' file_mask ...
              ' --forks ' sprintf('%i',opt.forks) ...
              ' --memlimit ' sprintf('%i',opt.memlimits) ...
              ' ' folder_subdist];
[status,msg] = system(call_cwas);
if status~=0
    error('There was a problem calling connectir_subdist.R: %s',msg);
end

%% Second call CWAS: now run MDMR
call_cwas = [ 'connectir_mdmr.R ' ...
              ' --indir ' folder_subdist ...
              ' --formula "' opt.formula '"' ...
              ' --model ' file_model ...  
              ' --permutations ' sprintf('%i',opt.nb_permutation) ...
              ' --factors2perm "' sprintf('%s',opt.factors2perm) '"'...              
              ' --forks ' sprintf('%i',opt.forks) ...
              ' --threads ' sprintf('%i',opt.threads) ...
              ' --memlimit ' sprintf('%i',opt.memlimits) ...
              ' --save-perms' ...
              ' --ignoreprocerror' ...
              ' mdmr'];
[status,msg] = system(call_cwas);
if status~=0
    error('There was a problem calling connectir_mdmr.R: %s',msg);
end

%% Read output
file_pce = [folder_subdist 'mdmr' filesep 'pvals.bin'];
hf = fopen(file_pce);
pce_seed = fread(hf,Inf,"double");
fclose(hf);
psom_clean(opt.path_work,struct('flag_verbose',false));

function [] = sub_write_ascii(filename,data)

[nx,ny] = size(data);
hf = fopen(filename,'w');
for xx = 1:nx
    for yy = 1:ny
        if yy == ny
            fprintf(hf,'%1.15f\n',data(xx,yy))
        else
            fprintf(hf,'%1.15f ',data(xx,yy))
        end
    end
end
fclose(hf);