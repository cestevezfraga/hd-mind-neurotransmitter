%Create files for BNV and Circular Graph - run after running NBS analysis

% Written by Peter McColgan UCL
% last edit 04.01.25

%Add paths 
addpath(genpath('~/Documents/MATLAB/NBS1.2/'))
addpath('~/Documents/MATLAB/BCT/2019_03_03_BCT')
addpath('~/Documents/MATLAB/Graph_analysis_01.03.14/')
addpath(genpath('~/Documents/MATLAB/BrainNetViewer_20191031/'))
addpath(genpath('~/Documents/MATLAB/circularGraph-master/'))

%saveas(gcf,'NBS_yas_cont_preHD','jpg');
data = table;
% Display NBS Connections
global nbs; [i,j]=find(nbs.NBS.con_mat{1});

%% Print results and write to a file
fileID = fopen('mHD_track_nfl_p_results.csv','w');

for n=1:length(i)
i_lab=nbs.NBS.node_label{i(n)};
j_lab=nbs.NBS.node_label{j(n)};
stat=nbs.NBS.test_stat(i(n),j(n));
fprintf(fileID,'%s %s %0.2f\n',i_lab,j_lab,stat);
end

fclose(fileID);

%import coords, labels & results
coords = importdata('desikan_68_coords.txt','\t');
labels = importdata('desikan_68_labels.txt','\t');
regions = erase(labels,'ctx-');
regions = strrep(regions,'-','.');

% Create adjacency matrix
global nbs;
matrix = full(nbs.NBS.con_mat{1});

% Save full matrix for circular graph
mHD_track_matrix_nfl_p = matrix;
save('mHD_track_matrix_nfl_p','mHD_track_matrix_nfl_p')

%Concatenate in one column and remove duplicates
node_column_with_repat = cat(1,i,j);
node_column = unique(node_column_with_repat);

%Write truncated matrix to a text file - edge file
mHD_track_edge_nfl_p = matrix(node_column,node_column);
save ('mHD_track_edge_nfl_p.txt', 'mHD_track_edge_nfl_p', '-ascii')

%Create table - node file
coords_trunc = coords(node_column,:);
data = table;
data.x = coords_trunc(:,1);
data.y = coords_trunc(:,2);
data.z = coords_trunc(:,3);
data.ncolour =  ones(length(coords_trunc),1);
data.nsize =  ones(length(coords_trunc),1);
data.labels = regions(node_column);
writetable(data,'mHD_track_node_nfl_p.txt','Delimiter',' ','WriteVariableNames',0)
