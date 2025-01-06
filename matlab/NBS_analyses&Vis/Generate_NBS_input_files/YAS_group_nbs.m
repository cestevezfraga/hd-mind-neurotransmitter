
%Peter McColgan UCL
%last edit 27.05.24
%Script for MIND paper
clear all
close all

%Define variables
m.phd = 57;
m.cont = 60;

%Load demographics - rank demographics gene carriers & then controls
demo= xlsread('yas_compare_clinical'); %covariates - age, gender, site, TIV
covars = demo(:,2:4);

% Load matrices
[ndata, s_ind, alldata] = xlsread('yas_compare_clinical','A:A');
s_ind(1,:)=[];
%s_ind = s_ind';
cd yas/
for s0 = 1 : length(s_ind)
    str     = sprintf(['./',s_ind{s0},'/mind.csv']); % loading connectivity matrix
    con_mat_lab=importdata(str); % load mat file into workspace
    con_mat=con_mat_lab.data;
    con_tot(:,:,s0) = con_mat;
end

%Create Design Matrices with group with covariates -> your contrast
%here would be [-1,1,0,0,0,0]
phd_cont_dm = [cat(1,ones(m.phd,1),zeros(m.cont,1)) cat(1,zeros(m.phd,1),ones(m.cont,1)) covars];
cont_phd_dm = [cat(1,ones(m.cont,1),zeros(m.phd,1)) cat(1,zeros(m.cont,1),ones(m.phd,1)) phd_cont_dm([m.phd + 1:end 1:m.phd],3:5)];

%Group Matrices 
pcw = con_tot;
cpw = pcw(:,:,[m.phd + 1:end 1:m.phd]);

% save connectivity matrcies
nbs_files = {'pcw' 'cpw'};
for i = 1:length(nbs_files)
str1 = [nbs_files{i}];
str  = sprintf(str1);
save(fullfile(str),str1);
end

% save design matrices
dm = {'phd_cont_dm' 'cont_phd_dm'};
for i = 1:length(dm)
str2 = [dm{i}];
str  = sprintf(str2);
save(fullfile(str),str2);
end
    

