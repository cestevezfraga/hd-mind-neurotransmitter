%pmccolgan
%last edit 30.05.24
%Script for calculating functional graph metrics
clear all
close all

% Define and add paths
root = pwd;
addpath(genpath('~/Documents/MATLAB/NBS1.2/'))

%Define variables
m.phd = 57;
m.cont = 60;
m.tot = 57+60;

%Load demographics - rank demographics gene carriers & then controls
demo= xlsread('yas_compare_clinical'); %covariates - age, gender, TIV
cag = demo(:,5);
p_nfl = demo(:,10);
log_p_nfl = demo(:,11);
c_nfl = demo(:,12);
log_c_nfl = demo(:,13);
covars = demo(:,2:4);

% Create design matrix - NfL all comers
cd yas/
phdcont_p_nfl_dm = [ones(m.tot,1) [covars p_nfl]];
save('phdcont_p_nfl_dm','phdcont_p_nfl_dm')
% Create design matrix - log_p_nfl
phdcont_log_p_nfl_dm = [ones(m.tot,1) [covars log_p_nfl]];
save('phdcont_log_p_nfl_dm','phdcont_log_p_nfl_dm')
% Create design matrix - log_c_nfl
phdcont_c_nfl_dm = [ones(m.tot,1) [covars c_nfl]];
save('phdcont_c_nfl_dm','phdcont_c_nfl_dm')
% Create design matrix - log_c_nfl
phdcont_log_c_nfl_dm = [ones(m.tot,1) [covars log_c_nfl]];
save('phdcont_log_c_nfl_dm','phdcont_log_c_nfl_dm')

% Create design matrix - NfL all comers
phd_p_nfl_dm = [ones(m.phd,1) [covars(1:m.phd,:) p_nfl(1:m.phd,:)]];
save('phd_p_nfl_dm','phd_p_nfl_dm')
% Create design matrix - log_p_nfl
phd_log_p_nfl_dm = [ones(m.phd,1) [covars(1:m.phd,:) log_p_nfl(1:m.phd,:)]];
save('phd_log_p_nfl_dm','phd_log_p_nfl_dm')
% Create design matrix - log_c_nfl
phd_c_nfl_dm = [ones(m.phd,1) [covars(1:m.phd,:) c_nfl(1:m.phd,:)]];
save('phd_c_nfl_dm','phd_c_nfl_dm')
% Create design matrix - log_c_nfl
phd_log_c_nfl_dm = [ones(m.phd,1) [covars(1:m.phd,:) log_c_nfl(1:m.phd,:)]];
save('phd_log_c_nfl_dm','phd_log_c_nfl_dm')

% Load matrices
[ndata, s_ind, alldata] = xlsread('../yas_compare_clinical','A:A');
s_ind(1,:)=[];
%s_ind = s_ind';
for s0 = 1 : length(s_ind)
    str     = sprintf(['./',s_ind{s0},'/mind.csv']); % loading connectivity matrix
    con_mat_lab=importdata(str); % load mat file into workspace
    con_mat=con_mat_lab.data;
    con_tot(:,:,s0) = con_mat;
end

%Group Matrices 
pcw_nfl = con_tot;
pw_nfl = con_tot(:,:,1:m.phd);

% save connectivity matrcies
nbs_files = {'pcw_nfl' 'pw_nfl'};
for i = 1:length(nbs_files)
str1 = [nbs_files{i}];
str  = sprintf(str1);
save(fullfile(str),str1);
end