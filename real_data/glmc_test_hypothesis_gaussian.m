%% Script to test the hypothesis of Gaussian residual distribution
%% for the different datasets in the GLM connectome paper

clear all
path_data = '/home/pbellec/database/glm_connectome/test_hypotheses/';

%% GLM options
opt_glm.test = 'ttest';

%% MOTOR
list_files = { [path_data 'alloego' filesep 'glm_rest2VS1_avg_sci270_scg297_scf308.mat'] ; ...
               [path_data 'blind' filesep 'glm_CBvsSC_conf_sci280_scg308_scf313.mat'] };
name_xp = { 'motor' , 'blind' } ;
for num_f = 1:length(list_files)
    fprintf('Experiment %s. ',name_xp{num_f});
    data = load(list_files{num_f});
    res = niak_glm(data.model_group,opt_glm);
    pce_w = zeros(size(res.e,2),1);
    for num_c = 1:size(res.e,2)
        niak_progress(num_c,size(res.e,2))
        [h,pce_w(num_c)] = swtest(res.e(:,num_c));
    end
    hp = subplot(1,length(list_files),num_f);
    hist(pce_w,100);
    FN = findall(hp,'-property','FontName');
    set(FN,'FontName','/usr/share/fonts/truetype/dejavu/DejaVuSerifCondensed.ttf');
    FS = findall(hp,'-property','FontSize');
    set(FS,'FontSize',8);
    title(name_xp);
    ylabel('histogram')
    xlabel('p-value (Shapiro-Wilk)');
    %axis([0 max(param.list_scales) 0.007 0.15]);
    fprintf('Percentage of pce below 0.5: %1.3f\n',sum(pce_w <0.05)/length(pce_w));
end