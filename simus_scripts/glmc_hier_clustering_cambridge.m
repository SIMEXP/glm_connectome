% A script to generate a hierarchical clustering on the average 
% connectome of the Cambridge database 

clear
path_roi = '/media/database6/glm_connectome/cambridge_region_growing_05scrubb/rois';
list_files = dir([path_roi filesep 'tseries_rois_*_session1_rest.mat']);
list_files = {list_files.name};

for ff = 1:length(list_files)
     data = load([path_roi filesep list_files{ff}]);
     if ff == 1
         conn = niak_fisher(niak_build_correlation(data.tseries));
     end
     conn = conn + niak_fisher(niak_build_correlation(data.tseries));     
end
conn = conn/length(list_files);
hier = niak_hierarchical_clustering(conn);
file_conn = [path_roi filesep 'hier_avg_connectome.mat'];
save(file_conn,'hier','conn');

%% BONUS: visualization
order = niak_hier2order(hier);
opt_v.limits = [-0.5 0.5];
niak_visu_matrix(conn(order,order),opt_v);