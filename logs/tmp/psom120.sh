"octave" --silent --eval "load('/home/pbellec/database/database2/cobre/cobre_qc_report_20140405/logs/PIPE.mat','path_work'), if ~ismember(path_work,{'gb_niak_omitted','gb_psom_omitted'}), path(path_work), end,fprintf('Octave version %s\n',OCTAVE_VERSION); [status,msg] = system('echo $MINC_TOOLKIT_VERSION'); fprintf('Minc-tool-kit version %s',msg); msg = which('niak_gb_vars'); fprintf('NIAK quarantine %s\n',msg); clear msg status; psom_worker('/home/pbellec/database/database2/cobre/cobre_qc_report_20140405/logs/worker/psom120/','/home/pbellec/database/database2/cobre/cobre_qc_report_20140405/logs/',120,'30-Mar-2016 10:04:19');,exit" >"/home/pbellec/database/database2/cobre/cobre_qc_report_20140405/logs/worker/psom120/worker.log" 2>&1
touch "/home/pbellec/database/database2/cobre/cobre_qc_report_20140405/logs/worker/psom120/worker.exit"