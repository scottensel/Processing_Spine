#!/bin/bash

#
# SBSN framework - Preprocessing
# June 2022
#
# Preparation of functional spine data
#
# Requirements: FSL, Spinal Cord Toolbox 5.5
#
#

####################################
#set -x
# Immediately exit if error
#set -e -o pipefail

# Exit if user presses CTRL+C (Linux) or CMD+C (OSX)
trap "echo Caught Keyboard Interrupt within script. Exiting now.; exit" INT

# Save script path
#PATH_SCRIPT=$PWD

# get starting time:
start=`date +%s`
####################################

# load in function that has paths to subject
. /mnt/d/SMA/Processing_Spine/path_to_subjects.sh 

tput setaf 6; 
echo -n "Enter the index of the step to perform (1 = Prepare for FEAT, 2 = Prepare for iCAP): "
tput sgr0;
read ind

# For each subject
for s in "${sub[@]}"; do

		cd $DIREC$s"/func/"

        for d in "${myFunc[@]}"; do

            # Will print */ if no directories are available
            cd $DIREC$s"/func/func"$d"/"

            tput setaf 2; echo "Smoothing started for " $DIREC$s"/func/func"$d

            # changes the template file that it loops through
            if [ "$ind" == "1" ]; then

                echo "NOW SMOOTHING"
                mkdir smooth
                
                # split the image in t and it will add a number name to it
                sct_image -i fmri_spine_moco_denoised_plusmean.nii.gz -split t -o smooth/fmri_spine_moco_denoised_plusmean.nii.gz
    
                files_to_smooth="$PWD"/smooth/*
                
                echo $files_to_smooth
                for b in $files_to_smooth; do # Loop through all files
    
                    sct_smooth_spinalcord -i $b -s fmri_spine_moco_mean_seg_corr.nii.gz -smooth 0.849329,0.849329,2.5479870902 -o $b #smooth each file based on seg
                    
                    ## dont know which is better but they produce similar outputs Centerline goes further below and above but not needed since the manual seg goes far enough especially for analysis
                    #fmri_spine_moco_mean_seg_corr.nii.gz
                    #fmri_mean_crop_centerline.nii.gz
                done
    
                # create a list of all the files in the smooth folder
                allSmoothFiles=`${FSLDIR}/bin/imglob -extensions smooth/*`
                    
                # concatenate them all together
                sct_image -i $allSmoothFiles -concat t -o fmri_spine_moco_denoised_plusmean_smooth.nii.gz
                
                # copy header for TR
                sct_image -i fmri_spine_moco.nii.gz -copy-header fmri_spine_moco_denoised_plusmean_smooth.nii.gz -o fmri_spine_moco_denoised_plusmean_smooth.nii.gz
    
                # delete the smooth folder to save space
                rm -rf smooth

            elif [ "$ind" == "2" ]; then

                echo "NOW SMOOTHING"

                mkdir smooth2
                
                # split the image in t and it will add a number name to it
                sct_image -i fmri_spine_moco_denoised.nii.gz -split t -o smooth2/fmri_spine_moco_denoised.nii.gz
    
                files_to_smooth="$PWD"/smooth2/*
                
                echo $files_to_smooth
                for b in $files_to_smooth; do # Loop through all files
    
                    sct_smooth_spinalcord -i $b -s fmri_spine_moco_mean_seg_corr.nii.gz -smooth 0.849329,0.849329,2.5479870902 -o $b #smooth each file based on seg
                    
                    ## dont know which is better but they produce simialr outputs Centerline goes further below and above but not needed since the manual seg goes far enough especially for analysis
                    #fmri_spine_moco_mean_seg_corr.nii.gz
                    #fmri_mean_crop_centerline.nii.gz
                done
    
                # create a list of all the files in the smooth folder
                allSmoothFiles=`${FSLDIR}/bin/imglob -extensions smooth2/*`
                    
                # concatenate them all together
                sct_image -i $allSmoothFiles -concat t -o fmri_spine_moco_denoised_smooth.nii.gz
                
                # copy header for TR
                sct_image -i fmri_spine_moco.nii.gz -copy-header fmri_spine_moco_denoised_smooth.nii.gz -o fmri_spine_moco_denoised_smooth.nii.gz
    
                # delete the smooth folder to save space
                rm -rf smooth2

            fi

        done

		tput setaf 2; echo "Done!" 
        	tput sgr0;		 			
done

####################################
# Display useful info for the log
end=`date +%s`
runtime=$((end-start))
echo
echo "~~~"
echo "SCT version: `sct_version`"
echo "Ran on:      `uname -nsr`"
echo "Duration:    $(($runtime / 3600))hrs $((($runtime / 60) % 60))min $(($runtime % 60))sec"
echo "~~~"
####################################