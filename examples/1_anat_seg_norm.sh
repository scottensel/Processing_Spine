#!/bin/bash
#
# SPiCiCAP framework - Preprocessing
# June 2022
#
# This script performs automatic SC segmentation, vertebrae labeling and normalization to the PAM50 template
#
# Requirements: Spinal Cord Toolbox 5.5
#
#

# List of subjects
#declare -a sub=("sub-01" "sub-02" "sub-03" "sub-04" "sub-05" "sub-06" "sub-07" "sub-08" "sub-09" "sub-10" "sub-11" "sub-12" "sub-13" "sub-14" "sub-15" "sub-16" "sub-17" "sub-18" "sub-19" "sub-20" "sub-21" "sub-22")
declare -a sub=("SBSN_S_001")

# Path of the folder containing all data
#DIREC="/home/sce30/Documents/Research/NeuronData/SPiCiCAP_Dataset_Code/Data/"
DIREC="/mnt/d/Data/SBSN_S_001/"
#DIREC="/home/sce30/Documents/Research/Data/Spine/"

# For each subject
for s in "${sub[@]}"; do
	cd $DIREC$s"/anat/"

	tput setaf 2; echo "Segmentation started in "$s
        tput sgr0;

	# Segmentation (deep learning based method) with viewer initialization
    # label centerline of spinal cord for automatic segmentation of spine
	sct_deepseg_sc -i t2.nii.gz -c t2 -centerline viewer

	tput setaf 2; echo "Done!"
        tput sgr0;

	# Vertebral labeling (manual initialization of labeling by clicking at disc C2-C3)
    # outputs the cord segmentation labeled with vertebrae level automatically 
	# sct_label_vertebrae -i t2.nii.gz -s t2_seg.nii.gz -c t2
    # creates t2_seg_labeled, _discs, warp_curve2straight, warp_straight2curve

	# Create labels at specific vertebral levels
        # sct_label_utils -i t2_seg_labeled.nii.gz -vert-body 4,7

        # label vertebrae number 3 to 9 and create a labels.nii.gz file
        # should be at posterior edge of intervertebral discs
    sct_label_utils -i t2.nii.gz -create-viewer 3:9 -o labels.nii.gz -msg "Place labels at the posterior tip of each inter-vertebral disc. E.g. Label 3: C2/C3, Label 4: C3/C4, etc."


	#
	# NOTE: segmentation and labeling should be visually evaluated!
	#

 	tput setaf 2; echo "Normalization started in "$s
        tput sgr0;

	# Normalization to PAM50 template	
	#sct_register_to_template -i t2.nii.gz -s t2_seg.nii.gz -l labels.nii.gz -c t2 -param step=0,type=label,dof=Tx_Ty_Tz_Sz:step=1,type=seg,algo=centermassrot,metric=MeanSquares,iter=10,smooth=2,gradStep=0.5,slicewise=0,smoothWarpXY=2,pca_eigenratio_th=1.6:step=2,type=seg,algo=bsplinesyn,metric=MeanSquares,iter=3,smooth=1,slicewise=0,gradStep=0.5,smoothWarpXY=2,pca_eigenratio_th=1.6:step=3,type=im,metric=CC
    
    # register the spinal cord to the PAM50 template image using the labeled image with the t2 and segmented t2 image
    sct_register_to_template -i t2.nii.gz -s t2_seg.nii.gz -ldisc labels.nii.gz -c t2 -param step=0,type=label,dof=Tx_Ty_Tz_Sz:step=1,type=seg,algo=centermassrot,metric=MeanSquares,iter=10,smooth=2,gradStep=0.5,slicewise=0,smoothWarpXY=2,pca_eigenratio_th=1.6:step=2,type=seg,algo=bsplinesyn,metric=MeanSquares,iter=3,smooth=1,slicewise=0,gradStep=0.5,smoothWarpXY=2,pca_eigenratio_th=1.6:step=3,type=im,metric=CC

    # -ldisc Labels located at the posterior edge of the intervertebral discs

	tput setaf 2; echo "Done!"
        tput sgr0;
done
