 #!/bin/bash
#
# SBSN framework - Preprocessing
# June 2022
#
# Preparation of functional spine data
#
# Requirements: FSL, Spinal Cord Toolbox 5.5
#
# SPINE

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

	cd $DIREC$s"/func/"

    tput setaf 2; echo "Prepare third level analysis for GLM " $s"/func/func"$d
    tput sgr0; 

    cp -r "$DIREC"reg/"" $DIREC$s"/func/level_two_force_FLOB_smooth1234.gfeat/cope1.feat/stats"
    
    # these have to be the same size across all subject because they get concatenated in the 4th dimension
    # change this to the PAM50 template

    cp level_two_force_FLOB_smooth1234.gfeat/cope1.feat/mean_func.nii.gz level_two_force_FLOB_smooth1234.gfeat/cope1.feat/stats/reg/standard.nii.gz
    cp level_two_force_FLOB_smooth1234.gfeat/cope1.feat/mean_func.nii.gz level_two_force_FLOB_smooth1234.gfeat/cope1.feat/stats/reg/example_func.nii.gz

    #cp anat2template.nii.gz level_two_force_FLOB_smooth1234.gfeat/cope1.feat/example_func.nii.gz
    #cp anat2template.nii.gz level_two_force_FLOB_smooth1234.gfeat/cope1.feat/mean_func.nii.gz

    cp -r "$DIREC"reg/"" $DIREC$s"/func/level_two_force_FLOB_smooth1234.gfeat/cope4.feat/stats"
    
    # these have to be the same size across all subject because they get concatenated in the 4th dimension
    # change this to the PAM50 template

    cp level_two_force_FLOB_smooth1234.gfeat/cope4.feat/mean_func.nii.gz level_two_force_FLOB_smooth1234.gfeat/cope4.feat/stats/reg/standard.nii.gz
    cp level_two_force_FLOB_smooth1234.gfeat/cope4.feat/mean_func.nii.gz level_two_force_FLOB_smooth1234.gfeat/cope4.feat/stats/reg/example_func.nii.gz

    #cp anat2template.nii.gz level_two_force_FLOB_smooth1234.gfeat/cope4.feat/example_func.nii.gz
    #cp anat2template.nii.gz level_two_force_FLOB_smooth1234.gfeat/cope4.feat/mean_func.nii.gz

    cp -r "$DIREC"reg/"" $DIREC$s"/func/level_two_force_FLOB_smooth1234.gfeat/cope7.feat/stats"
    
    # these have to be the same size across all subject because they get concatenated in the 4th dimension
    # change this to the PAM50 template

    cp level_two_force_FLOB_smooth1234.gfeat/cope7.feat/mean_func.nii.gz level_two_force_FLOB_smooth1234.gfeat/cope7.feat/stats/reg/standard.nii.gz
    cp level_two_force_FLOB_smooth1234.gfeat/cope7.feat/mean_func.nii.gz level_two_force_FLOB_smooth1234.gfeat/cope7.feat/stats/reg/example_func.nii.gz

    #cp anat2template.nii.gz level_two_force_FLOB_smooth1234.gfeat/cope7.feat/example_func.nii.gz
    #cp anat2template.nii.gz level_two_force_FLOB_smooth1234.gfeat/cope7.feat/mean_func.nii.gz

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
