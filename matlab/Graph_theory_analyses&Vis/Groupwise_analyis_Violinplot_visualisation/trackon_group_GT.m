
%Peter McColgan UCL
%last edit 27.05.24
%Script for MIND paper
clear all
close all

% Permutation testing requires perm_test4.m in the path

%Define variables
m.phd = 85;
m.cont = 89;
cd ..
%Load demographics - rank demographics gene carriers & then controls
demo= xlsread('demographics'); %covariates - age, gender, site, TIV
covars = demo(:,2:5);

% Load matrices
[ndata, s_ind, alldata] = xlsread('demographics','A:A');
s_ind(1,:)=[];
%s_ind = s_ind';
cd trackon/
for s0 = 1 : length(s_ind)
    str     = sprintf(['./',s_ind{s0},'/mind.csv']); % loading connectivity matrix
    con_mat_lab=importdata(str); % load mat file into workspace
    con_mat=con_mat_lab.data;
    con_tot(:,:,s0) = con_mat;
end

%Graph theory metrics
addpath('~/Documents/MATLAB/BCT/2019_03_03_BCT/')

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
    gma.(str13) = mafdr(gma.(str12).p_two,'BHFDR','true');
end

%Write table for results
data = table;
regions = strrep(con_mat_lab.textdata(2:69,1),'_',' ');
data.Regions = regions;
data.FDR = gma.c_fdr_stren; 
data.p_two = gma.c_results_stren.p_two; 
data.p_left = gma.c_results_stren.p_left;
data.p_right = gma.c_results_stren.p_right;
data.diff = gma.c_results_stren.diff';
writetable(data,'/Users/petermccolgan/Desktop/mind/trackon_node_results.csv')

%Write table for strength for participants for all regions
colnames = regions;
ID = s_ind(1:end,1);
Group = cat(1,ones(m.phd,1),zeros(m.cont,1));
%T1 = cell2table(num2cell(gma.stren'), 'VariableNames', colnames); %raw
%strength values
T1 = cell2table(num2cell(gma.pc_stren'), 'VariableNames', colnames); %residuals
T2 = addvars(T1,ID,'Before',regions(1)); 
T3 = addvars(T2,Group,'Before',regions(1)); 
writetable(T3,'/Users/petermccolgan/Desktop/mind/trackon_node_strength.csv')

% FDR corrected
data.Regions(find(data.FDR >0.05),:); % non-significant regions
sig_regions = data.Regions(find(data.p_two <0.05),:); % significant regions (uncorrected)

%Convert group numbers to names
G = num2cell(Group);
preHD = Group == 1;
Control = Group == 0;
G(preHD) = {'HDGEC'};
G(Control) = {'Control'};
T3.Group = G;
T3.ID = []; %remove ID

%Create violin plots
T4 = stack(T3,colnames,'NewDataVariableName','Strength');
T5 = renamevars(T4,'Strength_Indicator','Brain Region');

% Create table of only significant regions
sig_inx = find(ismember(T5.("Brain Region"),sig_regions));
T6 = T5(sig_inx,:);


%Violin Plot

grpandplot(T6,"Strength",yTitle='Strength (residuals)',xFactor="Brain Region",cFactor="Group",xOrder=sig_regions,...
showXLine=true,showVln=true,showBox=false);
title('Late PreHD: HDGEC vs. Controls')
ylim([min(T6.Strength) max(T6.Strength)])
fontsize(12,"points")
% using exportgraphics()
exportgraphics(gca, '/Users/petermccolgan/Desktop/mind/trackon_highres.png','Resolution','600')