glm_connectome_paper
==============

This github repository contains the scripts used to generate the results of the following publication: 

P. Bellec, Y. Benhajali, F. Carbonell, C. Dansereau, G. Albouy, M. Pelland, C. Craddock, O. Collignon, J. Doyon, E. Stip, P. Orban _Multiscale statistical testing for connectome-wide association studies in fMRI_ [arXiv:1409.2080](http://arxiv.org/abs/1409.2080) [q-bio.QM]

This code is not meant as a toolbox for public use. The GLM-connectome method is publicly available as part of the [NIAK](https://github.com/SIMEXP/niak) package for the "boss" release 0.12.18 and latter. To run the code, please download [NIAK 0.12.18](http://www.nitrc.org/frs/download.php/7163/niak-boss-0.12.18.zip) as well as the [latest release](https://github.com/SIMEXP/glm_connectome/releases) of this repository. For any question, please contact pierre.bellec (at) criugm.qc.ca. 

Data
----

Due to restrictions imposed by the consent forms signed by the subjects, some of the data used in the paper cannot be shared. The following resources were made available:
  * The [preprocessed time series](http://figshare.com/articles/Cambridge_resting_state_fMRI_time_series_preprocessed_with_NIAK_0_12_4/1159331) (at the region level) for the Cambridge sample. These were used as a basis to simulate group changes in connectivity.
  * The [preprocessed 3D+t fMRI time series](http://figshare.com/articles/COBRE_preprocessed_with_NIAK_0_12_4/1160600) of the COBRE sample, used in the SCHIZO analysis of real data. 

Simulations
-----------
 * The [simulation script](https://github.com/SIMEXP/glm_connectome/blob/master/simus_scripts/glmc_script_simu_ind.m) for the simulations with independent tests.
 * The [group hierarchical clustering](https://github.com/SIMEXP/glm_connectome/blob/master/simus_scripts/glmc_hier_clustering_cambridge.m) of the Cambridge sample.
 * The [simulation pipeline](https://github.com/SIMEXP/glm_connectome/blob/master/simus_scripts/glmc_pipeline_simu_dep_regular_grid.m) for the simulations with dependent tests, using a regular grid of scales.
 * The [simulation pipeline](https://github.com/SIMEXP/glm_connectome/blob/master/simus_scripts/glmc_pipeline_simu_dep_MSTEPS.m) for the simulations with dependent tests, using scales selected by MSTEPS.

GLM-connectome analysis of the MOTOR dataset
-----------------------------
 * The [preprocessing](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/MOTOR_pipeline_preprocess.m) pipeline.
 * The [region growing](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/MOTOR_pipeline_region_growing.m) pipeline.
 * The [Boostrap Analysis of Stable Clusters](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/MOTOR_pipeline_basc_regular_grid.m) pipeline, with regular grid of scales.
 * The [Boostrap Analysis of Stable Clusters](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/MOTOR_pipeline_basc_MSTEPS.m) pipeline, with scales selected by MSTEPS.
 * The [Multiscale Statistical Parametric Connectome](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/MOTOR_pipeline_MSPC_regular_grid.m) pipeline.
 * The [Multiscale Statistical Parametric Connectome](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/MOTOR_pipeline_MSPC_MSTEPS.m) pipeline, with scales selected by MSTEPS.

GLM-connectome analysis of the BLIND dataset
-----------------------------
 * The [preprocessing](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/BLIND_pipeline_preprocess.m) pipeline.
 * The [region growing](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/BLIND_pipeline_region_growing.m) pipeline.
 * The [Boostrap Analysis of Stable Clusters](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/BLIND_pipeline_basc_regular_grid.m) pipeline, with regular grid of scales.
 * The [Boostrap Analysis of Stable Clusters](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/BLIND_pipeline_basc_MSTEPS.m) pipeline, with scales selected by MSTEPS.
 * The [Multiscale Statistical Parametric Connectome](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/BLIND_pipeline_MSPC_regular_grid.m) pipeline.
 * The [Multiscale Statistical Parametric Connectome](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/BLIND_pipeline_MSPC_MSTEPS.m) pipeline, with scales selected by MSTEPS.

GLM-connectome analysis of the SCHIZO dataset
-----------------------------
 * The [preprocessing](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/SCHIZO_pipeline_preprocess.m) pipeline.
 * The [region growing](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/SCHIZO_pipeline_region_growing.m) pipeline.
 * The [Boostrap Analysis of Stable Clusters](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/SCHIZO_pipeline_basc_regular_grid.m) pipeline, with regular grid of scales.
 * The [Boostrap Analysis of Stable Clusters](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/SCHIZO_pipeline_basc_MSTEPS.m) pipeline, with scales selected by MSTEPS.
 * The [Multiscale Statistical Parametric Connectome](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/SCHIZO_pipeline_MSPC_regular_grid.m) pipeline.
 * The [Multiscale Statistical Parametric Connectome](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/SCHIZO_pipeline_MSPC_MSTEPS.m) pipeline, with scales selected by MSTEPS.
