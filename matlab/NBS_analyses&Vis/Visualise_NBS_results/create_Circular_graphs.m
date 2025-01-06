%% Create Circular Graph 
clear all 
close all
% Copyright 2016 The MathWorks, Inc.

%% Adapted for Desikan Atlas by Peter McColgan UCL
%Last edit 04.01.25

%% Addpaths
addpath('~/Documents/MATLAB/circularGraph-master/')

% cd into input file folder 
cd Circular_Graph_Input_files

%% Create custom colormap & rank labels by regions
figure                 % Creates a figure
set(gca,'FontSize',4) % Creates an axes and sets its FontSize to 18
colour_labels = importdata('desikan_68_colour_scheme.txt');
colour_scheme = colour_labels.data;
[~,rank] = sort(colour_scheme,'ascend');
colour = lines(max(colour_scheme));

%% Import labels & rank by region
labels = importdata('desikan_68_labels.txt','\t');
regions = erase(labels,'ctx-');
regions = strrep(regions,'-','.');
myLabel = regions(rank);

for i = 1:max(colour_scheme)
map = reshape(repelem(colour(i,:),numel(find(colour_scheme==i))),numel(find(colour_scheme==i)),3);
C{i} = map;
end
myColorMap = vertcat(C{:});

%% Group-wise graphs -  See 'Create_BNV_&_CircularGraph_input_files folder'for scripts to generate files below

% Group-wise analyses
load('mHD_matrix_cp.mat');
mHD_matrix_cp = mHD_matrix_cp(rank,rank);
load('mHD_matrix_pc.mat');
mHD_matrix_pc = mHD_matrix_pc(rank,rank);
load('YAS_matrix_cp.mat');
YAS_matrix_cp = YAS_matrix_cp(rank,rank);

% Run circular graph function - mHD_matrix_cp
%circularGraph(mHD_matrix_cp,'Colormap',myColorMap,'Label',myLabel);
%Save image using exportgraphics()
%exportgraphics(gca, '/Users/petermccolgan/Desktop/natcomms_matlab/NBS_analyses/Visualise_NBS_results/CircularGraphFigures/mHD_matrix_cp_with_labels_highres.png','Resolution','600')
% clear gca
% 
% % Run circular graph function - mHD_matrix_pc
circularGraph(mHD_matrix_pc,'Colormap',myColorMap,'Label',myLabel);
% % Save image using exportgraphics()
exportgraphics(gca, '/Users/petermccolgan/Desktop/natcomms_matlab/NBS_analyses/Visualise_NBS_results/CircularGraphFigures/mHD_matrix_pc_with_labels_highres.png','Resolution','600')
% clear gca
% 
% % Run circular graph function - YAS_matrix_cp
%circularGraph(YAS_matrix_cp,'Colormap',myColorMap,'Label',myLabel);
% % Save image using exportgraphics()
%exportgraphics(gca, '/Users/petermccolgan/Desktop/mind/circular_graph_highres/YAS_matrix_cp_with_labels_highres.png','Resolution','600')
% clear gca
% 
% %% Nfl graphs
% 
% % Nfl analyses
load('mHD_track_matrix_nfl_p.mat');
mHD_track_matrix_nfl_p = mHD_track_matrix_nfl_p(rank,rank);
load('mHD_track_matrix_nfl_pc.mat');
mHD_track_matrix_nfl_pc = mHD_track_matrix_nfl_pc(rank,rank);
load('preHD_track_matrix_p_nfl_p.mat');
preHD_track_matrix_p_nfl_p = preHD_track_matrix_p_nfl_p(rank,rank);
load('preHD_track_matrix_p_nfl_pc.mat');
preHD_track_matrix_p_nfl_pc = preHD_track_matrix_p_nfl_pc(rank,rank);
% 
% % Run circular graph function - mHD_track_matrix_nfl_p
%circularGraph(mHD_track_matrix_nfl_p,'Colormap',myColorMap,'Label',myLabel);
% % Save image using exportgraphics()
%exportgraphics(gca, '/Users/petermccolgan/Desktop/mind/circular_graph_highres/mHD_track_matrix_nfl_p_with_labels_highres.png','Resolution','600')
% clear gca
% 
% % Run circular graph function - mHD_track_matrix_nfl_pc
%circularGraph(mHD_track_matrix_nfl_pc,'Colormap',myColorMap,'Label',myLabel);
% % Save image using exportgraphics()
%exportgraphics(gca, '/Users/petermccolgan/Desktop/mind/circular_graph_highres/mHD_track_matrix_nfl_pc_with_labels_highres.png','Resolution','600')
% clear gca
% 
% % Run circular graph function - preHD_track_matrix_p_nfl_p
%circularGraph(preHD_track_matrix_p_nfl_p,'Colormap',myColorMap,'Label',myLabel);
% % Save image using exportgraphics()
%exportgraphics(gca, '/Users/petermccolgan/Desktop/mind/circular_graph_highres/preHD_track_matrix_p_nfl_p_with_labels_highres.png','Resolution','600')
% clear gca
% 
% % Run circular graph function - preHD_track_matrix_p_nfl_pc
circularGraph(preHD_track_matrix_p_nfl_pc,'Colormap',myColorMap,'Label',myLabel);
% % Save image using exportgraphics()
exportgraphics(gca, '/Users/petermccolgan/Desktop/mind/circular_graph_highres/preHD_track_matrix_p_nfl_pc_with_labels_highres.png','Resolution','600')
% clear gca
% 
% 
% 
% 
% 
% 
% 
% 
% %delete(findobj(gca,'Type','Text')) %uncomment to delete labels
