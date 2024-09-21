
%Peter McColgan UCL
%last edit 27.05.24
%Script for MIND paper
clear all
close all

%Define variables
m.phd = 110;
m.cont = 111;

%Load demographics - rank demographics gene carriers & then controls
demo= xlsread('track_fsv7'); %covariates - age, gender, site, TIV
covars = demo(:,2:5);

% Load matrices
[ndata, s_ind, alldata] = xlsread('track_fsv7','A:A');
s_ind(1,:)=[];
%s_ind = s_ind';
cd track/
for s0 = 1 : length(s_ind)
    str     = sprintf(['./',s_ind{s0},'/mind.csv']); % loading connectivity matrix
    con_mat_lab=importdata(str); % load mat file into workspace
    con_mat=con_mat_lab.data;
    con_tot(:,:,s0) = con_mat;
end

%Graph theory metrics
addpath('~/Documents/MATLAB/BCT/2019_03_03_BCT/')
addpath('~/Documents/MATLAB/Graph_analysis_01.03.14/')
for n = 1:length(s_ind)
gma.stren(:,n)=strengths_und(con_tot(:,:,n)); 
end

% Metrics
metrics  = {'stren'};

nperm = 10000; % number of permutations
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pre-HD VS Control %%%%%%%%%%%%%%%%
design=[ones(1,m.phd) 2*ones(1,m.cont)];
N=length(design);
for i = 1:length(metrics)
    str10 = ['pc_' metrics{i}];
    str11 = [metrics{i}];
    str12 = ['c_results_' metrics{i}];
    str13 = ['c_fdr_' metrics{i}];
    gma.(str10) = gma.(str11);
    for n = 1:size(gma.(str10),1)
        y = gma.(str10)(n,:);
        [B,BINT,R] = regress(y',[covars ones(size(covars,1),1)]);
        gma.(str10)(n,:) = R;
    end
    perm = perm_test4('init',m.phd,m.cont,nperm);
    gma.(str12) = perm_test4('test',gma.(str10)(:,1:m.phd)',gma.(str10)(:,m.phd + 1:end)',perm);
    gma.(str13) = gma.(str12).p_two;
end

%Write table for results
data = table;
data.Regions = con_mat_lab.textdata(2:69,1);
data.FDR = gma.c_fdr_stren; 
data.p_two = gma.c_results_stren.p_two; 
data.p_left = gma.c_results_stren.p_left;
data.p_right = gma.c_results_stren.p_right;
data.diff = gma.c_results_stren.diff';
writetable(data,'/Users/petermccolgan/Desktop/mind/track_node_results.csv')

%Write table for strength for participants for all regions
colnames = con_mat_lab.textdata(1,2:69);
ID = s_ind(1:end,1);
Group = cat(1,ones(m.phd,1),zeros(m.cont,1));
T1 = cell2table(num2cell(gma.stren'), 'VariableNames', colnames);
T2 = addvars(T1,ID,'Before',con_mat_lab.textdata(1,2:2)); 
T3 = addvars(T2,Group,'Before',con_mat_lab.textdata(1,2:2)); 
writetable(T3,'/Users/petermccolgan/Desktop/mind/track_node_strength.csv')