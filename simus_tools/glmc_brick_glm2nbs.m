function [in,out,opt] = glmc_brick_glm2nbs(in,out,opt)
% Generate network-based statistics from a NIAK GLM_CONNECTOME structure
%
% SYNTAX: [IN,OUT,OPT] = GLMC_BRICK_GLM2NBS(IN,OUT,OPT)
%
% IN.GLM (string) a XXX_glm.mat structure. See NIAK_PIPELINE_GLM_CONNECTOME pipeline.
% IN.NETWORK (string) the name of a 3D volume with the networks used for the GLM analysis.
%
% OUT.NBS (string) the file name of a 4D dataset. NBS(:,:,:,n) is the t-stat map 
%   corresponding to network n. The significance of connections is derived using 
%   the NBS method.
% OUT.PERC_DISCOVERY (string) the file name of a 3D volume, with the map of the 
%   number of discovery associated with each network, expressed as a percentage of the 
%   number of networks - 1 (i.e. the max possible number of discoveries associated 
%   with a network).
%
% OPT.THRESH (scalar, default 2) the uncorrected peak threshold on t-stats.
% OPT.ALPHA (scalar, default 0.01) the test on cluster-level p-values.
% OPT.NB_SAMPS (integer, default 1000) the number of permutation tests.
% OPT.DIRECTION (string, default 'positive') the direction of the t-test.
%   Available options: 'positive' or 'negative'.
% OPT.SIZE (string, default 'extent') how to measure the size of networks. 
%   Available options: 'extent' and 'intensity'.
% OPT.FLAG_TEST (boolean, default false) if the flag is true, the brick does 
%   nothing but update IN, OUT, OPT.
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
    { 'glm' , 'network' }, ...
    { NaN   , NaN       });
   
out = psom_struct_defaults(out, ...
    { 'nbs' , 'perc_discovery' }, ...
    { NaN   , NaN              });

if nargin < 3
    opt = struct();
end
opt = psom_struct_defaults(opt, ...
    { 'direction' , 'size'   , 'thresh' , 'alpha' , 'flag_test' , 'rand_seed' , 'nb_samps' }, ...
    { 'positive'  , 'extent' , 2        , 0.01    , false       , []          , 1000       });

if opt.flag_test
    return
end

%% Control random variations
if ~isempty(opt.rand_seed)
    psom_set_rand_seed(opt.rand_seed);
end

%% Read the input model
glm_niak = load(in.glm,'model_group');
glm_niak = glm_niak.model_group;
ttest = load(in.glm,'ttest');
ttest = ttest.ttest;

%% Read the input partition
[hdr,mask] = niak_read_vol(in.network);

%% Now run NBS
stats.thresh = opt.thresh;
stats.alpha = opt.alpha;
stats.size = opt.size;
stats.test_stat = [];
nb_net = max(mask(:));
stats.N = nb_net;
ind_upper = find(triu(ones(stats.N,stats.N),1));
glm = struct();
for num_s = 1:size(glm_niak.y,1)
    if num_s == 1
        glm.y = zeros([size(glm_niak.y,1) length(ind_upper)]);
    end
    tmp = niak_lvec2mat(glm_niak.y(num_s,:));        
    glm.y(num_s,:) = tmp(ind_upper);
end    
    
glm.X = glm_niak.x;
switch opt.direction
    case 'positive'
         glm.contrast = glm_niak.c(:)';
    case 'negative'
         glm.contrast = -glm_niak.c(:)';
    otherwise
         error('%s is an unkown option for OPT.DIRECTION',opt.direction)
end
glm.test = 'ttest';
glm.perms = opt.nb_samps;    
[ncon,con_mat] = NBSstats(stats,NaN,glm);
test_q = false(size(niak_lvec2mat(glm_niak.y(1,:))));
for num_e = 1:length(con_mat)
    test_q = test_q | full(con_mat{num_e}) | full(con_mat{num_e}');
end    

%% Save the results: discovery maps
nb_discovery = sum(test_q,1);
perc_discovery = nb_discovery/size(test_q,1);
discovery_maps = niak_part2vol(perc_discovery,mask);
hdr.file_name = out.perc_discovery;
niak_write_vol(hdr,discovery_maps);

%% Save the results: NBS maps
nbs_maps   = zeros([size(mask) nb_net]);
ttest_mat = niak_lvec2mat (ttest);
for num_net = 1:nb_net    
    ttest_thre = ttest_mat(:,num_net);
    ttest_thre( ~test_q(:,num_net) ) = 0;
    nbs_maps(:,:,:,num_net) = niak_part2vol(ttest_thre,mask);
end
hdr.file_name = out.nbs;
niak_write_vol(hdr,nbs_maps);
