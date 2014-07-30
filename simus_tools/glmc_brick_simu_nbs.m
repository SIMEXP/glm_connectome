function [files_in,files_out,opt] = glmc_brick_simu_nbs(files_in,files_out,opt)
% Simulate multiscale changes in connectivity between two populations - 
% test with network-based statistics.
%
% SYNTAX:
% [FILES_IN,FILES_OUT,OPT] = GLMC_BRICK_SIMU_NBS(FILES_IN,FILES_OUT,OPT)
%
% _________________________________________________________________________
% INPUTS
%
% FILES_IN 
%   (structure) with the following fields :
%
%   TSERIES
%      (cell of strings) each entry is the name of 
%      a .mat file with one variable TSERIES of size T x N. That entry is
%      only expected for simulations based on real data.
%
%   HIER
%      (string) a .mat file with a hierarchy on the N regions of the brain.
%
% FILES_OUT       
%   (string) the results of the simulation
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
%   PERC_RAND
%      (scalar, default 0) a percentage PERC_RAND of the cluster of reference is 
%      randomly assigned to other clusters.
%
%   NB_REPLICATION
%      (integer, default 1) the number of replications of the simulation.
%
%   RAND_SEED
%      (scalar, default []) The specified value is used to seed the random
%      number generator with PSOM_SET_RAND_SEED. If left empty, no action
%      is taken.
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
%   FLAG_TEST 
%      (boolean, default 0) if FLAG_TEST equals 1, the brick does not 
%      do anything but update the default values in FILES_IN, 
%      FILES_OUT and OPT.
%        
% _________________________________________________________________________
% OUTPUTS
%
% The structures FILES_IN, FILES_OUT and OPT are updated with default
% valued. If OPT.FLAG_TEST == 0, the specified outputs are written.
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

%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('niak:brick','syntax: [FILES_IN,FILES_OUT,OPT] = GLMC_BRICK_SIMU_MULTISCALE(FILES_IN,FILES_OUT,OPT).\n Type ''help glmc_brick_simu_multiscale'' for more info.')
end

%% Files in
list_fields   = { 'tseries' , 'hier' };
list_defaults = { NaN       , NaN    };
files_in = psom_struct_defaults(files_in,list_fields,list_defaults);

%% Options
list_fields   = { 'perc_rand' , 'nb_replication' , 'type_data' , 'theta2' , 'nb_samps' , 'fdr' , 'type_fdr' , 'nb_subject' , 'rand_seed' , 'alpha2' , 'list_scales' , 'scale_ref' , 'cluster_ref' , 'flag_verbose' , 'flag_test'  };
list_defaults = { 0           , 1                , 'real'      , 0.1      , 0          , 0.05  , 'LSL_sym'  , NaN          , []          , 0.2      , NaN           , NaN         , NaN           , true           , false        };
opt = psom_struct_defaults(opt,list_fields,list_defaults);

if opt.scale_ref == 0
    opt.scale_ref = [];
end
        
%% If the test flag is true, stop here !
if opt.flag_test == 1
    return
end

%% Random number generator
if ~isempty(opt.rand_seed)
    psom_set_rand_seed(opt.rand_seed);
end

%% Read hierarchy
hier = load(files_in.hier);
hier = hier.hier;

%% Read time series
tseries_all = cell(length(opt.list_scales),1);
for tt = 1:length(files_in.tseries)
    tseries = load(files_in.tseries{tt});
    tseries_all{tt} = tseries.tseries;    
end

opt_simu = rmfield(opt,{'flag_test','rand_seed','nb_replication'});

for rr = 1:opt.nb_replication
    results(rr) = glmc_simu_nbs(tseries_all,hier,opt_simu);
end

%% Save the results
save(files_out,'results')