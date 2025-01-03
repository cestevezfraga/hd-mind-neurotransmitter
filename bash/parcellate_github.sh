# Script to parcellate files already processed with recon-all 
#
#Need to rename /mri/lh.pial.T1 as mri/lh.pial [same for rh]
#fsaverage6 folder with annot files in fsaverage space in the SUBJECTS_DIR


export FREESURFER_HOME=/Applications/freesurfer/7.1.1
export SUBJECTS_DIR=/Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos
source $FREESURFER_HOME/SetUpFreeSurfer.sh

cd $SUBJECTS_DIR
echo "Working in directory; `pwd`"
echo ""

for j in `awk '{print $1}' all4.txt` ; do
 
    echo ${j}

    
#HCP parcellation [180 ROI per side]
mri_surf2surf --srcsubject fsaverage6 --trgsubject ${j} --hemi lh --sval-annot /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/fsaverage6/label/lh.HCP.annot --tval /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/lh.HCP.annot
 
mri_surf2surf --srcsubject fsaverage6 --trgsubject ${j} --hemi rh --sval-annot /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/fsaverage6/label/rh.HCP.annot --tval /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/rh.HCP.annot
 
#The numbers after lh.annot (0) and rh.annot (1000) indicate the number added to the original #labels. Eg if insula was 23 now left insula will be 23 and right insula 1023. Only for cortical, not for #subcortical labels. 
 
mri_surf2volseg --o /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/mri/HCP.nii.gz --label-cortex --i /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/mri/aseg.mgz --threads 2 --lh-annot /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/lh.HCP.annot  0 --lh-cortex-mask /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/lh.cortex.label --lh-white /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/surf/lh.white --lh-pial /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/surf/lh.pial --rh-annot /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/rh.HCP.annot 1000 --rh-cortex-mask /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/rh.cortex.label --rh-white /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/surf/rh.white --rh-pial /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/surf/rh.pial
 
# DK318 Parcellation (160 per side)
 
mri_surf2surf --srcsubject fsaverage6 --trgsubject ${j} --hemi lh --sval-annot /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/fsaverage6/label/lh.500.sym.aparc.annot --tval /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/lh.500.sym.aparc.annot
 
mri_surf2surf --srcsubject fsaverage6 --trgsubject ${j} --hemi rh --sval-annot /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/fsaverage6/label/rh.500.sym.aparc.annot --tval /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/rh.500.sym.aparc.annot
 
mri_surf2volseg --o /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/mri/500.sym.aparc.nii.gz --label-cortex --i /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/mri/aseg.mgz --threads 2 --lh-annot /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/lh.500.sym.aparc.annot 0 --lh-cortex-mask /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/lh.cortex.label --lh-white /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/surf/lh.white --lh-pial /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/surf/lh.pial --rh-annot /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/rh.500.sym.aparc.annot 1000 --rh-cortex-mask /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/rh.cortex.label --rh-white /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/surf/rh.white --rh-pial /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/surf/rh.pial
 
 
 
# SCHAEFER 100
 
mri_surf2surf --srcsubject fsaverage6 --trgsubject ${j} --hemi lh --sval-annot /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/fsaverage6/label/lh.Schaefer2018_100Parcels_17Networks_order.annot --tval /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/lh.Schaefer2018_100Parcels_17Networks_order.annot
 
mri_surf2surf --srcsubject fsaverage6 --trgsubject ${j} --hemi rh --sval-annot /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/fsaverage6/label/rh.Schaefer2018_100Parcels_17Networks_order.annot --tval /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/rh.Schaefer2018_100Parcels_17Networks_order.annot
 
mri_surf2volseg --o /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/mri/Schaefer2018_100Parcels_17Networks_order.nii.gz --label-cortex --i /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/mri/aseg.mgz --threads 2  --lh-annot /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/lh.Schaefer2018_100Parcels_17Networks_order.annot 0 --lh-cortex-mask /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/lh.cortex.label --lh-white /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/surf/lh.white --lh-pial /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/surf/lh.pial --rh-annot /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/rh.Schaefer2018_100Parcels_17Networks_order.annot 1000 --rh-cortex-mask /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/rh.cortex.label --rh-white /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/surf/rh.white --rh-pial /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/surf/rh.pial
 
 
 
# SCHAEFER 500
 
mri_surf2surf --srcsubject fsaverage6 --trgsubject ${j} --hemi lh --sval-annot /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/fsaverage6/label/lh.Schaefer2018_500Parcels_17Networks_order.annot --tval /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/lh.Schaefer2018_500Parcels_17Networks_order.annot
 
mri_surf2surf --srcsubject fsaverage6 --trgsubject ${j} --hemi rh --sval-annot /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/fsaverage6/label/rh.Schaefer2018_500Parcels_17Networks_order.annot --tval /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/rh.Schaefer2018_500Parcels_17Networks_order.annot
 
mri_surf2volseg --o /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/mri/Schaefer2018_500Parcels_17Networks_order.nii.gz --label-cortex --i /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/mri/aseg.mgz --threads 2  --lh-annot /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/lh.Schaefer2018_500Parcels_17Networks_order.annot 0 --lh-cortex-mask /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/lh.cortex.label --lh-white /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/surf/lh.white --lh-pial /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/surf/lh.pial --rh-annot /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/rh.Schaefer2018_500Parcels_17Networks_order.annot 1000 --rh-cortex-mask /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/rh.cortex.label --rh-white /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/surf/rh.white --rh-pial /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/surf/rh.pial
 
 
# SCHAEFER 1000 
 
mri_surf2surf --srcsubject fsaverage6 --trgsubject ${j} --hemi lh --sval-annot /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/fsaverage6/label/lh.Schaefer2018_1000Parcels_17Networks_order.annot --tval /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/lh.Schaefer2018_1000Parcels_17Networks_order.annot
 
mri_surf2surf --srcsubject fsaverage6 --trgsubject ${j} --hemi rh --sval-annot /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/fsaverage6/label/rh.Schaefer2018_1000Parcels_17Networks_order.annot --tval /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/rh.Schaefer2018_1000Parcels_17Networks_order.annot
 
mri_surf2volseg --o /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/mri/Schaefer2018_1000Parcels_17Networks_order.nii.gz --label-cortex --i /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/mri/aseg.mgz --threads 2 --lh-annot /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/lh.Schaefer2018_1000Parcels_17Networks_order.annot 0 --lh-cortex-mask /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/lh.cortex.label --lh-white /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/surf/lh.white --lh-pial /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/surf/lh.pial --rh-annot /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/rh.Schaefer2018_1000Parcels_17Networks_order.annot 1000 --rh-cortex-mask /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/label/rh.cortex.label --rh-white /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/surf/rh.white --rh-pial /Users/charlie/Desktop/my_projects/neurotransmitter/imaging_metrics/thickness/track/carlos/${j}/surf/rh.pial

    
    echo ${j} 'done'
 
done
