%pmccolgan
%last edit 30.05.24
%Script for calculating functional graph metrics
clear all
close all

% Define and add paths
root = pwd;
addpath(genpath('~/Documents/MATLAB/NBS1.2/'))

%Define variables
m.phd = 103;
m.cont = 111-23;
m.tot = 103+111-23;

%Load demographics - rank demographics gene carriers & then controls
demo= xlsread('track_fsv7_prehd.xls'); %covariates - age, gender, TIV
nan = find(any(isnan(demo),2))';
demo(nan,:) = [];
covars = demo(:,2:5);
p_nfl = demo(:,6);
log_p_nfl = demo(:,7);

% Create design matrix - NfL all comers
cd prehd_track/
phdcont_p_nfl_dm = [ones(m.tot,1) [covars p_nfl]];
save('phdcont_p_nfl_dm','phdcont_p_nfl_dm')
% Create design matrix - log_p_nfl
phdcont_log_p_nfl_dm = [ones(m.tot,1) [covars log_p_nfl]];
save('phdcont_log_p_nfl_dm','phdcont_log_p_nfl_dm')


% Create design matrix - NfL all comers
phd_p_nfl_dm = [ones(m.phd,1) [covars(1:m.phd,:) p_nfl(1:m.phd,:)]];
save('phd_p_nfl_dm','phd_p_nfl_dm')
% Create design matrix - log_p_nfl
phd_log_p_nfl_dm = [ones(m.phd,1) [covars(1:m.phd,:) log_p_nfl(1:m.phd,:)]];
save('phd_log_p_nfl_dm','phd_log_p_nfl_dm')

% Load matrices
[ndata, s_ind, alldata] = xlsread('../track_fsv7_prehd','A:A');
s_ind(1,:)=[];
%s_ind = s_ind';
for s0 = 1 : length(s_ind)
    str     = sprintf(['./',s_ind{s0},'/mind.csv']); % loading connectivity matrix
    con_mat_lab=importdata(str); % load mat file into workspace
    con_mat=con_mat_lab.data;
    con_tot(:,:,s0) = con_mat;
end

%Group Matrices 
con_tot(:,:,nan) = [];
pcw_nfl = con_tot;
pw_nfl = con_tot(:,:,1:m.phd);

% save connectivity matrcies
nbs_files = {'pcw_nfl' 'pw_nfl'};
for i = 1:length(nbs_files)
str1 = [nbs_files{i}];
str  = sprintf(str1);
save(fullfile(str),str1);
end