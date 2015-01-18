 
% Create the inputs of and launch the NIAK_PIPELINE_FMRI_PREPROCESS on the
% specified dataset (cobre dataset)
clear all

addpath(genpath('/home/porban/quarantaine/niak-boss-0.12.2/'))

path_raw_fmri       = '/home/porban/database/blind/mnc_data_select/';
path_preprocess     = '/gs/scratch/porban/bascglm/blind_fmri_preprocess_bascglm_20140405/';



%%


groups_list = {'CB', 'SC'};
    
    for group_n = 1:size(groups_list,2)
        group = groups_list{group_n};
        path_group = [path_raw_fmri,filesep,group,filesep];
        subjects_list = dir(path_group);
        subjects_list = subjects_list(3:end);
        
        for num_s = 1:size(subjects_list,1)

            %% Subject file names
            subject = subjects_list(num_s).name

            %anat
            fmrirun = dir([path_group filesep subject filesep 'anat' filesep '*.mnc.gz']);
            anat = [path_group filesep subject filesep 'anat' filesep fmrirun.name];
            
            %func
            fmrirun = dir([path_group filesep subject filesep 'fmri' filesep '*.mnc.gz']);
            fmri.session1.run1=[path_group filesep subject filesep 'fmri' filesep fmrirun.name];
                       
            files_in.([group 'xxx' strrep(subject,'_','')]).fmri = fmri;
            files_in.([group 'xxx' strrep(subject,'_','')]).anat = anat;
            
        end
        
    end
 

    
%% Building the optional inputs
opt.folder_out = path_preprocess; 
opt.size_output = 'all';

%%%%%%%%%%%%%%%%%%%%
%% Bricks options %%
%%%%%%%%%%%%%%%%%%%%

%% Slice timing
opt.slice_timing.type_acquisition = 'interleaved ascending'; % Slice timing order (available options : 'sequential ascending', 'sequential descending', 'interleaved ascending', 'interleaved descending')
opt.slice_timing.type_scanner     = 'Siemens';               % Scanner manufacturer. Only the value 'Siemens' will actually have an impact
opt.slice_timing.delay_in_tr      = 0;                       % The delay in TR ("blank" time between two volumes)
% opt.slice_timing.flag_center      = 1;
%% Motion correction (niak_brick_motion_correction)
%opt.motion_correction.suppress_vol = 0;             % Remove the first three dummy scans

%% Linear and non-linear fit of the anatomical image in the stereotaxic
%% space 
opt.t1_preprocess.nu_correct.arg = '-distance 50'; % Parameter for non-uniformity correction. 200 is a suggested value for 1.5T images, 25 for 3T images. If you find that this stage did not work well, this parameter is usually critical to improve the results.

% T1-T2 coregistration (niak_brick_anat2func)
opt.anat2func.init = 'identity'; % An initial guess of the transform. Possible values 'identity', 'center'. 'identity' is self-explanatory. The 'center' option usually does more harm than good. Use it only if you have very big misrealignement between the two images (say, 2 cm).

%% Temporal filetring (niak_brick_time_filter)
opt.time_filter.hp = 0.01; % Apply a high-pass filter at cut-off frequency 0.01Hz (slow time drifts)
opt.time_filter.lp = Inf; % Do not apply low-pass filter. Low-pass filter induce a big loss in degrees of freedom without sgnificantly improving the SNR.

%% Correction of physiological noise (niak_pipeline_corsica)
opt.corsica.sica.nb_comp = 60;
opt.corsica.component_supp.threshold = 0.15;

%% Resampling in the stereotaxic space (niak_brick_resample_vol)
%opt.resample_vol.interpolation       = 'tricubic'; % The resampling scheme. The most accurate is 'sinc' but it is awfully slow
opt.resample_vol.voxel_size          = [3 3 3];    % The voxel size to use in the stereotaxic space

%% Spatial smoothing (niak_brick_smooth_vol)
opt.bricks.smooth_vol.fwhm = 6; % Apply an isotropic 6 mm gaussin smoothing.

%% Region growing
opt.region_growing.flag_skip = 1; % Turn on/off the region growing
%opt.template_fmri = '/home/cdansereau/svn/niak/trunk/template/roi_aal.mnc.gz';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generation of the pipeline %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt.flag_test = 0;
%opt.psom.max_queued = 24; % Please try to use the two processors of my laptop, thanks !
%opt.granularity = 'subject';

[pipeline,opt] = niak_pipeline_fmri_preprocess(files_in,opt);

