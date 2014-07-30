%% 


path_raw  = '/home/porban/database/alloego/mnc_data/';
path_data = '/home/porban/database/alloego/mnc_data_select/';
groups_list = {'allonap','allononap','egonap','egononap'};

type = {'anat','sess1/rest1','sess1/rest2'};
newtype = {'anat','fmri/rest1','fmri/rest2'};


for g = 1:size(groups_list,2)
      
    group = groups_list{g};
    path_group = [path_raw,filesep,group,filesep];
    subjects_list = dir(path_group);
    subjects_list = subjects_list(3:end);
    
    for s = 1:size(subjects_list,1)
        
        subject = subjects_list(s).name;

        mkdir([path_data subject]);
        mkdir([path_data subject '/anat/']);
        mkdir([path_data subject '/fmri/']);
        mkdir([path_data subject '/fmri/rest1']);
        mkdir([path_data subject '/fmri/rest2']);
        
        for t = 1:length(type)
            
            cmd = ['rsync -av ' path_raw groups_list{g} '/' subject '/' type{t} '/*mnc.gz ' path_data subject '/' newtype{t} '/'];
            system(cmd);
        end
    end
end

            
            


