#!/bin/bash
#
# SPiCiCAP framework - Preprocessing
# July 2020
#
# This script performs motion correction and extracts slice-wise motion parameters
#
# Requirements: Spinal Cord Toolbox 3.2.7, FSL
#
#

# List of subjects
declare -a sub=("sub-01")
#"sub-02" "sub-03" "sub-04" "sub-05" "sub-06" "sub-07" "sub-08" "sub-09" "sub-10" "sub-11" "sub-12" "sub-13" "sub-14" "sub-15" "sub-16" "sub-17" "sub-18" "sub-19" "sub-20" "sub-21" "sub-22")

# Path of the folder containing all data
#DIREC="/home/sce30/Documents/Research/NeuronData/SPiCiCAP_Dataset_Code/Data/"
DIREC="/mnt/d/NeuronData/SPiCiCAP_Dataset_Code/Data/"
#DIREC="/home/sce30/Documents/Research/Data/Spine/"

# For each subject
for s in "${sub[@]}"; do
   	cd $DIREC$s"/func/"

		# Slicewise motion correction
                tput setaf 2; echo "Moco started for "$s
                tput sgr0;
		
		# Create mask to constrain motion correction metrics		

		# Check if mask does not exist
		if [ ! -f Mask/mask_fmri.nii.gz ]; then
			mkdir -p Mask
			fslmaths fmri.nii.gz -Tmean fmri_mean.nii.gz # taking the mean of the functional images across time
                	mv fmri_mean.nii.gz Mask
                	cd Mask
                	sct_get_centerline -i fmri_mean.nii.gz -c t2 -method optic #automatically getting the centerline of the spinal cord from the mean functional images
                    # output is fmri_mean_centerline.nii.gz

                	#sct_create_mask -i fmri_mean.nii.gz -p centerline,fmri_mean_centerline_optic.nii.gz -size 30mm -o mask_fmri.nii.gz
                    sct_create_mask -i fmri_mean.nii.gz -p centerline,fmri_mean_centerline.nii.gz -size 30mm -o mask_fmri.nii.gz 
                    # creates mask of spinal cord along z direction using the centerline file

			cd ..
		fi

		mkdir -p Processing

                sct_fmri_moco -i fmri.nii.gz -m Mask/mask_fmri.nii.gz -g 1 -r 0 -param smooth=2 -ofolder Moco

                #sct_fmri_moco -i fmri.nii.gz -m Mask/mask_fmri.nii.gz -r 0 -x spline -param poly=0 -ofolder Moco
                # motion correction of the fmri data. -m is using the mask voxels to limit voxels considered
                # using spline to interpret

        # just moving and renaming some files
		mv Moco/fmri_moco.nii.gz mfmri.nii.gz
		mv Moco/fmri_moco_mean.nii.gz mfmri_mean.nii.gz
		if [ -f mfmri.nii.gz ]; then
			mv fmri.nii.gz Processing/fmri.nii.gz
		fi	

		tput setaf 2; echo "Moco done!"
                tput sgr0;

done
