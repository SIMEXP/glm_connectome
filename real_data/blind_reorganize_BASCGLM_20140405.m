%% 



group = {'CB','SC'};
subj_CB = {'VD_AlCh','VD_AnBe','VD_Be_Me','VD_DiCe','VD_FrCo','VD_LL','VD_Ma_Du','VD_MaLa','VD_MoBe','VD_Na_Te','VD_SePo','VD_SoSa','VD_YP','VD_Yv_La'};
subj_SC = {'VD_ChJa','VD_CJ','VD_ClDe','VD_GeAl','VD_JeRe','VD_JM','VD_JoFr','VD_KaFo','VD_LA_LH','VD_MaSa','VD_NiLe','VD_NiMi','VD_OL','VD_PG','VD_SC','VD_SG','VD_TJ'};
type = {'anat','fmri'};


for g = 1:length(group)
    
    if g ==1
        
        for s = 1:length(subj_CB)
        
            mkdir(['/home/porban/database/blind/mnc_data_select/' group{g} '/' subj_CB{s}]);
            mkdir(['/home/porban/database/blind/mnc_data_select/' group{g} '/' subj_CB{s} '/anat/']);
            mkdir(['/home/porban/database/blind/mnc_data_select/' group{g} '/' subj_CB{s} '/fmri/']);
            
            for t = 1:length(type)
            
                cmd = ['rsync -av /home/porban/database/blind/mnc_data/' type{t} '/' group{g} '/' subj_CB{s} '/*mnc.gz ' '/home/porban/database/blind/mnc_data_select/' group{g} '/' subj_CB{s} '/' type{t} '/'];
                system(cmd);
            end
        end
    else
        
        for ss = 1:length(subj_SC)
            
            mkdir(['/home/porban/database/blind/mnc_data_select/' group{g} '/' subj_SC{ss}]);
            mkdir(['/home/porban/database/blind/mnc_data_select/' group{g} '/' subj_SC{ss} '/anat/']);
            mkdir(['/home/porban/database/blind/mnc_data_select/' group{g} '/' subj_SC{ss} '/fmri/']);
            
            for tt = 1:length(type)
            
                cmd = ['rsync -av /home/porban/database/blind/mnc_data/' type{tt} '/' group{g} '/' subj_SC{ss} '/*mnc.gz ' '/home/porban/database/blind/mnc_data_select/' group{g} '/' subj_SC{ss} '/' type{tt} '/'];
                system(cmd);
            end
        end
    end
end


