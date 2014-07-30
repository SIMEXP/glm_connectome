
% write confounds.csv file for use in the mdmr interun analysis on the
% motor dataset

clear all

path = '/Users/pyeror/Work/transfert/Bascglm/models/';
model_name = 'alloego_model_group.csv';
output_name = 'alloego_mdmr_confounds.csv';


[tab,lx,ly] = niak_read_csv([path model_name]);


% (SUBJECT)_(SESSION)_(RUN)
for i = 1:length(lx)
    new_lx{i} = [lx{i} '_session1_run1'];  
end
for ii = length(lx)+1:length(lx)+length(lx)
    iii = ii - length(lx);
    new_lx{ii} = [lx{iii} '_session1_run2'];   
end

% VARIABLES
for j = 1:2
new_ly{j} = [ly{j+1}];
end
new_ly{3} = 'FD';

% DATA
for k = 1:length(tab) %age
    new_tab(k,1) = tab(k,2);
end
for kk = length(tab)+1:length(tab)+length(tab)
    kkk = kk - length(tab);
    new_tab(kk,1) = tab(kkk,2); 
end

for k = 1:length(tab) %sex
    new_tab(k,2) = tab(k,3);
end
for kk = length(tab)+1:length(tab)+length(tab)
    kkk = kk - length(tab);
    new_tab(kk,2) = tab(kkk,3); 
end

for k = 1:length(tab) %fd
    new_tab(k,3) = tab(k,4);
end
for kk = length(tab)+1:length(tab)+length(tab)
    kkk = kk - length(tab);
    new_tab(kk,3) = tab(kkk,5); 
end



% WRITE CSV
opt.labels_y = new_ly;
opt.labels_x = new_lx;
opt.precision = 2;
niak_write_csv(strcat(path,output_name),new_tab,opt)




