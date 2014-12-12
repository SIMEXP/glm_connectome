%% Script to test the hypothesis of Gaussian residual distribution
%% for the different datasets in the GLM connectome paper

clear all
path_data = '/home/pbellec/database/glm_connectome/test_hypotheses/';

%% GLM options
opt_glm.test = 'ttest';
opt_glm.flag_residuals = true;

%% MOTOR
list_files = { [path_data 'alloego' filesep 'glm_rest2VS1_avg_sci270_scg297_scf308.mat']      ; ...
               [path_data 'cambridge' filesep 'glm_subject40_sc328']                         ; ...
               [path_data 'blind' filesep 'glm_CBvsSC_conf_sci280_scg308_scf313.mat']         ; ...
               [path_data 'cambridge' filesep 'glm_subject40_sc328']                         ; ...
               [path_data 'cobre' filesep 'glm_szVScont_age_sex_FD_sci270_scg324_scf328.mat'] ; ...
               [path_data 'cambridge' filesep 'glm_subject180_sc328'] };
name_xp = { 'motor' , 'cambridge40' , 'blind' , 'cambridge100' , 'cobre' , 'cambridge180'} ;
for num_f = 1:length(list_files)
%for num_f = 5
    fprintf('Experiment %s. ',name_xp{num_f});
    data = load(list_files{num_f});
    pce_h = niak_white_test_hetero(data.model_group);
    pce_h = pce_h(~isnan(pce_h));
    hp = subplot(3,2,num_f);
    X = 0.01:0.02:0.99; 
    Y = hist(pce_h,X);
    Y = Y/(length(pce_h)*(X(2)-X(1)));
    axis([0,1,0,3])
    FN = findall(hp,'-property','FontName');
    set(FN,'FontName','/usr/share/fonts/truetype/dejavu/DejaVuSerifCondensed.ttf');
    FS = findall(hp,'-property','FontSize');
    set(FS,'FontSize',8);
    title([name_xp{num_f} ' scale ' num2str(size(niak_lvec2mat(data.model_group.y(1,:)),2))]);
    ylabel('normalized histogram')
    if (num_f==5)||(num_f==6)
        xlabel('p-value (White)');
    end
    set(gca,'Color',[193/255 226/255 247/255]);
    
    gridxy(0:0.1:1,'color',[1 1 1],'LineStyle','-','LineWidth',3) 
    gridxy([],0:0.5:3,'color',[1 1 1],'LineStyle','-','LineWidth',3) 
    hold on
    bar(X,Y)
    %axis([0 max(param.list_scales) 0.007 0.15]);
    
    [fdr,test] = niak_fdr(pce_h (:),'BH',0.05);
    fprintf('Percentage of pce below 0.5: %1.3f. Percentage of significant tests: %1.2f\n',sum(pce_h <0.05)/length(pce_h),sum(test)/length(test));
end
print([path_data 'fig_test_hetero.pdf'],'-dpdf')

return

%%% Small script to generate random models for the Cambridge datasets
path_cambridge = '/home/pbellec/database/glm_connectome/test_hypotheses/cambridge/';
list_files = dir([path_cambridge 'tseries*.mat']);
list_files = {list_files.name};
hier = load([path_cambridge 'hier_avg_connectome.mat']);
part = niak_threshold_hierarchy(hier.hier,struct('thresh',328));
for num_file = 1:length(list_files)
    data = load(list_files{num_file});
    conn_subject = niak_build_intra_inter_r(data.tseries,part,true);
    if num_file == 1
        conn = zeros(length(list_files),length(conn_subject));        
    end
    conn(num_file,:) = conn_subject;
end

model_group.y = conn(1:180,:);
model_group.x = [ones(size(model_group.y,1),1) [zeros(size(model_group.y,1)/2,1) ; ones(size(model_group.y,1)/2,1)]];
model_group.c = [0;1];
save([path_cambridge 'glm_subject180_sc328.mat'],'model_group');

model_group.y = conn(1:100,:);
model_group.x = [ones(size(model_group.y,1),1) [zeros(size(model_group.y,1)/2,1) ; ones(size(model_group.y,1)/2,1)]];
model_group.c = [0;1];
save([path_cambridge 'glm_subject100_sc328.mat'],'model_group');

model_group.y = conn(1:40,:);
model_group.x = [ones(size(model_group.y,1),1) [zeros(size(model_group.y,1)/2,1) ; ones(size(model_group.y,1)/2,1)]];
model_group.c = [0;1];
save([path_cambridge 'glm_subject40_sc328.mat'],'model_group');
