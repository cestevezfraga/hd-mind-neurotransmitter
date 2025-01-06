%pmccolgan
%last edit 30.05.24
%Script for calculating functional graph metrics
clear all
close all

% Define and add paths
root = pwd;
addpath(genpath('~/Documents/MATLAB/NBS1.2/'))

%Define variables
m.hd = 110-26;
m.cont = 111;
m.tot = 110+111-49;

%Load demographics - rank demographics gene carriers & then controls
demo= xlsread('track_fsv7'); %covariates - age, gender, TIV
nan = find(any(isnan(demo),2))';
hd_nan = nan(1:26);
demo(nan,:) = [];
covars = demo(:,2:5);
p_nfl = demo(:,6);
log_p_nfl = demo(:,7);
c_nfl = demo(:,8);
log_c_nfl = demo(:,9);


% Create design matrix - NfL all comers
cd track/
hdcont_p_nfl_dm = [ones(m.tot,1) [covars p_nfl]];
save('hdcont_p_nfl_dm','hdcont_p_nfl_dm')
% Create design matrix - log_p_nfl
hdcont_log_p_nfl_dm = [ones(m.tot,1) [covars log_p_nfl]];
save('hdcont_log_p_nfl_dm','hdcont_log_p_nfl_dm')
% Create design matrix - log_c_nfl
hdcont_c_nfl_dm = [ones(m.tot,1) [covars c_nfl]];
save('hdcont_c_nfl_dm','hdcont_c_nfl_dm')
% Create design matrix - log_c_nfl
hdcont_log_c_nfl_dm = [ones(m.tot,1) [covars log_c_nfl]];
save('hdcont_log_c_nfl_dm','hdcont_log_c_nfl_dm')

% Create design matrix - NfL all comers
hd_p_nfl_dm = [ones(m.hd,1) [covars(1:m.hd,:) p_nfl(1:m.hd,:)]];
save('hd_p_nfl_dm','hd_p_nfl_dm')
% Create design matrix - log_p_nfl
hd_log_p_nfl_dm = [ones(m.hd,1) [covars(1:m.hd,:) log_p_nfl(1:m.hd,:)]];
save('hd_log_p_nfl_dm','hd_log_p_nfl_dm')
% Create design matrix - log_c_nfl
hd_c_nfl_dm = [ones(m.hd,1) [covars(1:m.hd,:) c_nfl(1:m.hd,:)]];
save('hd_c_nfl_dm','hd_c_nfl_dm')
% Create design matrix - log_c_nfl
hd_log_c_nfl_dm = [ones(m.hd,1) [covars(1:m.hd,:) log_c_nfl(1:m.hd,:)]];
save('hd_log_c_nfl_dm','hd_log_c_nfl_dm')

% Load matrices
[ndata, s_ind, alldata] = xlsread('../track_fsv7','A:A');
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
pw_nfl = con_tot(:,:,1:m.hd);

% save connectivity matrcies
nbs_files = {'pcw_nfl' 'pw_nfl'};
for i = 1:length(nbs_files)
str1 = [nbs_files{i}];
str  = sprintf(str1);
save(fullfile(str),str1);
end