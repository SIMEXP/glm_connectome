glm_connectome_paper
==============

This github repository is a companion for the following publication:

_Multiscale statistical testing for connectome-wide association studies in fMRI_

P. Bellec<sup>1,2</sup>, Y. Benhajali<sup>1,3</sup>, F. Carbonell<sup>4</sup>, C. Dansereau<sup>1,2</sup>, Z. Shehzad<sup>5,6,7</sup>, G. Albouy<sup>1</sup>, M. Pelland<sup>1,8</sup>l, C. Craddock<sup>6,7</sup>, O. Collignon<sup>1,8</sup>, J. Doyon<sup>1,8</sup>, E. Stip<sup>1,9</sup>, P. Orban<sup>1,9</sup>

  <sup>1</sup>Functional Neuroimaging Unit, Centre de Recherche de l'Institut Universitaire de Gériatrie de Montréal  
  <sup>2</sup>Department of Computer Science and Operations Research, University of Montreal, Montreal, Quebec, Canada
  <sup>3</sup>Department of Anthropology, University of Montreal, Montreal, Quebec, Canada
  <sup>4</sup>Biospective Incorporated, Montreal, Quebec, Canada
  <sup>5</sup>Department of Psychology, Yale University, New Haven, CT, United States of America
  <sup>6</sup>Nathan Kline Institute for Psychiatric Research, Orangeburg, NY, United States of America
  <sup>7</sup>Center for the Developing Brain, Child Mind Institute, New York, NY, United States of America
  <sup>8</sup>Department of Psychology, University of Montreal, Montreal, Quebec, Canada
  <sup>9</sup>Department of Psychiatry, University of Montreal, Montreal, Quebec, Canada  
  
For all question regarding the paper, please address correspondence to Pierre Bellec, CRIUGM, 4545 Queen Mary, Montreal, QC, H3W 1W5, Canada. Email: pierre.bellec (at) criugm.qc.ca. 

Disclaimer
----------

This is a snapshot of the scripts used to generate the analysis for the above-mentioned publication at the time of production (Spring 2014) of the results reported in the paper. This code is not meant as a toolbox for public use. It is not documented, and many of the options are non-functional or have not been tested. This is a research-grade code, meant to clarify the details of implementation of the algorithms described in the aforementioned paper. The GLM-connectome method is publicly available as part of the [NIAK](https://github.com/SIMEXP/niak) package for the "boss" release 0.12.14 and latter. To run the code, please download NIAK 0.12.14 from the [NITRC site](http://www.nitrc.org/frs/?group_id=411).  Some functions also depend on the "Multivariate Distance Matrix Regression" package version 512745280a5350c527201c487a7e337041d6cfc6 (see the [tutorial](https://github.com/SIMEXP/glm_connectome/wiki/MDMR) to replicate the version used in the paper), as well as the "Network-Based Statistics" toolbox version 1.2 (see this [website](https://sites.google.com/site/bctnet/comparison/nbs) for download).

Data
----

Note that due to restrictions imposed by the consent forms signed by the subjects, some of the data used in the paper cannot be shared. The following resources have been made available:


Simulations
-----------
 * The [group hierarchicl clustering](https://github.com/SIMEXP/glm_connectome/blob/master/simus_scripts/glmc_hier_clustering_cambridge.m) of the Cambridge sample.

Analysis of the MOTOR dataset
-----------------------------
 * The [preprocessing](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/MOTOR_pipeline_preprocessed_20140405.m) pipeline.
 * The [region growing](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/MOTOR_pipeline_region_growing_20140405.m) pipeline.
 * The [Boostrap Analysis of Stable Clusters](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/MOTOR_pipeline_basc_20140405.m) pipeline.
 * The [GLM connectome (group FDR)](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/MOTOR_pipeline_glm_20140405.m) pipeline.
 * The [GLM connectome (BH FDR)](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/MOTOR_pipeline_glm_BH_FDR_20140518.m) pipeline.

Analysis of the BLIND dataset
-----------------------------
 * The [preprocessing](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/BLIND_pipeline_preprocess_20140405.m) pipeline.
 * The [region growing](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/BLIND_pipeline_region_growing_20140405.m) pipeline.
 * The [Boostrap Analysis of Stable Clusters](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/BLIND_pipeline_basc_20140405.m) pipeline.
 * The [GLM connectome (group FDR)](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/BLIND_pipeline_glm_20140405.m) pipeline.
 * The [GLM connectome (BH FDR)](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/BLIND_pipeline_glm_global_BH_FDR_20140518.m) pipeline.

Analysis of the SCHIZO dataset
-----------------------------
 * The [preprocessing](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/SCHIZO_pipeline_preprocess_20140405.m) pipeline.
 * The [region growing](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/SCHIZO_pipeline_region_growing_20140405.m) pipeline.
 * The [Boostrap Analysis of Stable Clusters](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/SCHIZO_pipeline_basc_20140405.m) pipeline.
 * The [GLM connectome (group FDR)](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/SCHIZO_pipeline_glm_20140405.m) pipeline.
 * The [GLM connectome (BH FDR)](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/SCHIZO_pipeline_glm_BH_FDR_20140518.m) pipeline.

MDMR and NBS (real data)
------------------------
 * Prepare [confounds](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/MOTOR_confounds_mdmr.m) for the MDMR analysis of the MOTOR dataset.
 * Run [MDMR](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/MDMR_MOTOR_20140720.m) for the MOTOR dataset.
 * Run [MDMR](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/MDMR_BLIND_SCHIZO_20140518.m) for the BLIND and SCHIZO datasets.
 * Run [NBS](https://github.com/SIMEXP/glm_connectome/blob/master/real_data/NBS_MOTOR_BLIND_SCHIZO_503.m) on the MOTOR, BLIND and SCHIZO datasets.
