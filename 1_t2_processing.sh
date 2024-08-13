#!/bin/bash
#
# SBSN framework - Preprocessing
# June 2022
#
# This script performs automatic SC segmentation, vertebrae labeling and normalization to the PAM50 template
#
# Requirements: Spinal Cord Toolbox 5.5
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
. /mnt/d/SBSN/Processing_Spine/path_to_subjects.sh 

# For each subject
for s in "${sub[@]}"; do
	cd $DIREC$s"/anat/"

	tput setaf 2; echo "Segmentation started in "$s
        tput sgr0;
    
    # automatically segment the spinal cord from t2 image
    sct_deepseg_sc -i t2.nii.gz -c t2 #-qc $DIREC$s"/anat/"

    # automatically segment the spinal cord from t2 image
    sct_maths -i t2_seg.nii.gz -dilate 3 -o t2_seg_dilate.nii.gz #-qc $DIREC$s"/anat/"

	tput setaf 2; echo "automatic segmentation one!"
        tput sgr0;

    # root segmentation
    sct_deepseg -i t2.nii.gz -o t2_roots.nii.gz -task seg_spinal_rootlets_t2w


    # manual create labels instead
    #sct_label_utils -i t2.nii.gz -create-viewer 1,2,3,4,5,6,7,8 -o label_discs.nii.gz -msg "Click atthe posterior tip of C1 to C7/T1 inter-vertebral disc"
    
    # automatically create vertebral labeling and disc labeling
    sct_label_vertebrae -i t2.nii.gz -s t2_seg.nii.gz -c t2 -o label_discs.nii.gz #-discfile label_discs.nii.gz #-qc $DIREC$s"/anat/"
    # now we have an automatic level of the spinal cord

	tput setaf 2; echo "vertebral labeling done!"
        tput sgr0;

    # apply registration to t2 to pam50 and back
    sct_register_to_template -i t2.nii.gz -s t2_seg.nii.gz -ldisc label_discs.nii.gz -c t2 

	tput setaf 2; echo "registered to template done!"
        tput sgr0;

    # apply tranforamtion to templates. Now they are all in subject space
    # this moves all WM, GM and CSF atlases into the subjects T2 space
    sct_warp_template -d t2.nii.gz -w warp_template2anat.nii.gz -s 1 -a 1 #-qc $DIREC$s"/anat/"


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