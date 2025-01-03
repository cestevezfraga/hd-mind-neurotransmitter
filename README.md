# hd-mind-neurotransmitter
Scripts to obtain mind connectivity and investigate the its biological underpinnings

1.	mind_github.ipynb (python)
Script to generate MIND networks
Input: T1 imaging data already parcellated  with freesurfer (examples from six particioants 000-000-001 and 000-000-002)
Output are subject-specific mind networks: mind.csv files inside each participant’s folder

2.	cohens_mind_track_node_github.R 
Script to obtain Cohens D for MIND networks in patients vs controls
Input: subject specific mind.csv files and track_fsv.csv file with the demographics
Output: cohens_d_right_transp.csv file with Cohens D across ROIs

3.	cohens_mind_track_node_correlate_github.R
Script to correlate MIND networks with disease burden score and plasma neurofilament
Input: track_fsv.csv file with the demographics, mind.csv file
Output: correlation_transp.csv file with correlation coefficients and P values an


4.	epicenter_mind_github.ipynb
Epicenter analysis from Cohens D data. Figure 4 from the paper. Using data from the mHD cohort. 
Note: In my laptop it required a different environment  (conda activate py38 with an outdated version of numpy )
Input: cohens_d_right_transp.csv file with Cohens D across ROIs
Output: 
	fc: functional connecticity. sc: structural connectivity. ctx: cortico-cortical. sctx: subcortico-cortical
	fc_ctx_track_node.png —> Figure depicting ROIs with significant associations between cortico-cortical functional connectivity in controls and MIND connectivity in HD
	fc_ctx_track_node_p_fdr.txt —> FDR corrected P values 
	fc_ctx_track_node_coefficients.txt —> Correlation coefficients
	Same with structural connectivity and subcortico-cortical connections as per the legenf

5.	contributions_mind_github.ipynb
Script to investigate the relative contribution of different organizational principles to MIND connectivity. Figure 5 from the paper.
Input: 
	Cammoun033_coords.txt —> Coordinates of the Cammoun Atlas
	data/Cammoun033/*.npy files —> Matrices for different systems (gene expression, receptor similarity etc
Output:
	Figures/Cammoun033/heatmap_disease_nncorr_dis_mind.svg —> Distance organizational principle
	Figures/Cammoun033/heatmap_disease_nncorr_scd.svg —> Structural connectivity organizational principle
	Figures/Cammoun033heatmap_disease_nncorr_NEG.svg —> Remaining organizational principles


6.	data_cortex_github.ipynb
Initial script to develop receptome gradients. 
Input: 
	PET_nifti_images
	Schaeffer 100 parcellation files
Output: 
	100Parcels7Networks_receptorprofiles.csv —> region x receptor matrix
	cort_dist_100.npy 

7.	F1_github.ipynb
Script to develop receptome gradients. Figure 6 from the paper
input:
	100Parcels7Networks_receptorprofiles.csv —> region x receptor matrix
	numpy array with receptor gradients 1, 2 and 3 (rc_g1_100.npy , rc_g2_100.npy, rc_g3_100.npy )

Output: 
	G1_receptors_green.png —> Bar plot representing the relative contribution of each neurotransmitter to each receptome gradient
	G1_on_surf_redblue.png —> Brain surface depicting the influence of each gradient across brain regions regions
	RC_scree_colour.png —> proportion of the variance explained by each receptor

8.	receptome_github.ipynb
Script to investigate the association between receptome gradients and MIND connectivity. Figure 6 from the paper.
Note: In my laptop it required a different environment  (conda activate py38 with an outdated version of numpy )
Input: 
	numpy array with receptor gradients 1, 2 and 3 (rc_g1_100.npy , rc_g2_100.npy, rc_g3_100.npy )
	cohens_d_right_transp.csv (example with the mHD cohort

Output:
	rc_g1_mind_node_green.png —> association between mind connectivity in pwHD and first receptome gradient
	rc_g2_mind_node_green.png —> association between mind connectivity in pwHD and second receptome gradient
	rc_g3_mind_node_green.png —> association between mind connectivity in pwHD and third receptome gradient


9.	Reorder_cammoun_github.ipynb
Script to reorder cohens D data from DK68 to Cammoun atlas, necessary for posterior steps
Input: cohens_d_right_transp.csv (in data —> cohens_d data from the mHD cohort
Output: cammoun_track.csv


10.	make_receptor_matrix_github.ipynb
Script to generate a region x receptor matrix with PET data
Input: 
	PET parcellated data (in data/PET_parcellated/scale068/)
Output: 
	receptor_names_pet.npy —> Numpy array with receptor matrix
	receptor_data_scale068.csv —> CSV file with receptor matrix
	zscored_receptor_data_scale068.csv —> Z scored CSV file with receptor matrix

11.	neurotransmitter_github.ipynb
Script to determine the PET neurotransmitter systems associated with MIND connectivity. Figure 7A and 7B
Input:
	modified_enigma_atrophy_node.csv —> Cohens D data across ROIs reordered with the Cammoun atlas. Each column corresponds to one cohort (early preHD, 	late preHD and mHD)
	cammoun_zscored_receptor_data_scale068.csv —> PET data from healtgy controls  
	colourmap.csv —> Colourmap for the figures
Output: 
	figures/bar_dominance_enigma_test_node.eps —> Bar plot showing the influence of neurotransmitter distribution in MIND connectivty
	figures/heatmap_dominance_enigma_test_node.eps —> Heatmap showing the relative contribution of each neurotransmitter


12.	autoradiography_mind_github.ipynb

Similar to the previous script, but with autoradiography. Figures 7C and 7D

Input: modified_enigma_atrophy_node.csv
Output: (SG = Supragranular, G = granular, IG = Infragranular)
	bar_dominance_enigma_aut_node_SG_6colors.eps —> Bar plot showing the influence of neurotransmitter distribution in MIND connectivity in the supragranular layer (G = granular, IG = Infragranular)
	heatmap_dominance_enigma_aut_node_SG_6colors.eps —> Heatmap showing the relative contribution of each neurotransmitter in the supragranular layer



13.	combat_github.ipynb

Script to perform combat harminonization, Figure S1
input: 
	mind_node_mhd.csv file —> mind node data from the mHD cohort
	demos_long_mhd —> demographic data (note: to preserve anonymization this demographic data has been generated using random numbers) 
output:
	data_combat_mhd.csv —> combat harmonized mind data
	
14.	parcellate_github.sh
Script to parcellated files already processed with recon-all

15.	mind_parcellations_github.ipynb
Script to obtain mind data with different parcellation resolutions
input: 
	999-999-999 (example subject after processing with recon-all

output: 
	mind_500.sym.aparc.csv
	mind_HCP.csv
	mind_Schaefer2018_1000Parcels_17Networks_order.csv
	mind_Schaefer2018_100Parcels_17Networks_order.csv
	mind_Schaefer2018_100Parcels_7Networks_order.csv
	mind_Schaefer2018_200Parcels_7Networks_order.csv
	mind_Schaefer2018_400Parcels_7Networks_order.csv
	mind_Schaefer2018_500Parcels_17Networks_order.csv


16.	brain_surface_github.ipynb
Script to generate brain surface fgure
Note: In my laptop it required a different environment  (conda activate py38 with an outdated version of numpy )
![image](https://github.com/user-attachments/assets/73785177-7611-4aa9-afb3-f2b4d3b5d24d)
