%Create images for BNV

% Written by Peter McColgan UCL
% Last edit 04.01.25

%% Add paths 
addpath(genpath('~/Documents/MATLAB/BrainNetViewer_20191031/'))

cd ./BNV_input_files

%% Group Analyses

% YAS_cp
BrainNet_MapCfg('BrainMesh_ICBM152.nv','YAS_node_cp.node','YAS_edge_cp.edge','options_highres.mat','../BNV_figures/YAS_cp.png');

% mHD_cp
BrainNet_MapCfg('BrainMesh_ICBM152.nv','mHD_node_cp.node','mHD_edge_cp.edge','options_highres.mat','../BNV_figures/mHD_cp.png');

% mHD_pc
BrainNet_MapCfg('BrainMesh_ICBM152.nv','mHD_node_pc.node','mHD_edge_pc.edge','options_highres.mat','../BNV_figures/mHD_pc.png');

%% NfL analyses

% preHD_NfL_pc
BrainNet_MapCfg('BrainMesh_ICBM152.nv','preHD_track_node_p_nfl_pc.node','preHD_track_edge_p_nfl_pc.edge','options_highres.mat','../BNV_figures/preHD_NfL_pc.png');

% preHD_NfL_p
BrainNet_MapCfg('BrainMesh_ICBM152.nv','preHD_track_node_p_nfl_p.node','preHD_track_edge_p_nfl_p.edge','options_highres.mat','../BNV_figures/preHD_NfL_p.png');

% mHD_NfL_pc
BrainNet_MapCfg('BrainMesh_ICBM152.nv','mHD_track_node_nfl_pc.node','mHD_track_edge_nfl_pc.edge','options_highres.mat','../BNV_figures/mHD_NfL_pc.png');

% mHD_NfL_p
BrainNet_MapCfg('BrainMesh_ICBM152.nv','mHD_track_node_nfl_p.node','mHD_track_edge_nfl_p.edge','options_highres.mat','../BNV_figures/mHD_NfL_p.png');
