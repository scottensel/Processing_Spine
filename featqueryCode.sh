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

#
# UPDATE WHAT THIS SCRIPT DOES TEXT
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

tput setaf 6; 
echo -n "Enter the index of the step to perform (1 = Prepare Files, 2 = Run Featquery): "
tput sgr0;
read ind

# For each subject
for s in "${sub[@]}"; do

        cd $DIREC$s"/func/"

        if [ "$ind" == "1" ]; then

            for d in "${myFunc[@]}"; do
                # Will print */ if no directories are available
                cd $DIREC$s"/func/func"$d"/"

                tput setaf 2; echo "mean in each FLOB " $DIREC$s"/func/func"$d"/"
                tput sgr0; 


                # take this from each level 1 FLOB and mean it to get a mean image
                # then-tmean them all together in the respective level 2 folder and call it mean_fun

                # first
                fslmaths level_one_FLOB.feat/filtered_func_data -Tmean level_one_FLOB.feat/subjectSpace_mean_func
                
                # now transform it to the space we want have it in
                tput setaf 1; 
                echo $DIREC$s"/func/func"$d"/level_one_FLOB.feat/mean_func"
                tput setaf 6;

                # files have already been transofrmed
                # subject space images are original images so just apply warps to them
                # no need to rename files again
                sct_apply_transfo -i "level_one_FLOB.feat/subjectSpace_mean_func.nii.gz" -d ../../../template/PAM50_t2s.nii.gz -w warp_anat2template.nii.gz -o "level_one_FLOB.feat/mean_func.nii.gz"


            done


            fslmerge -t /mnt/d/SBSN/Data/Spine/$s/func/level_two_FLOB12.gfeat/cope1.feat/mean_func /mnt/d/SBSN/Data/Spine/$s/func/func1/level_one_FLOB.feat/mean_func /mnt/d/SBSN/Data/Spine/$s/func/func2/level_one_FLOB.feat/mean_func
            fslmaths /mnt/d/SBSN/Data/Spine/$s/func/level_two_FLOB12.gfeat/cope1.feat/mean_func -Tmean /mnt/d/SBSN/Data/Spine/$s/func/level_two_FLOB12.gfeat/cope1.feat/mean_func


            fslmerge -t /mnt/d/SBSN/Data/Spine/$s/func/level_two_FLOB123.gfeat/cope1.feat/mean_func /mnt/d/SBSN/Data/Spine/$s/func/func1/level_one_FLOB.feat/mean_func /mnt/d/SBSN/Data/Spine/$s/func/func2/level_one_FLOB.feat/mean_func /mnt/d/SBSN/Data/Spine/$s/func/func3/level_one_FLOB.feat/mean_func
            fslmaths /mnt/d/SBSN/Data/Spine/$s/func/level_two_FLOB123.gfeat/cope1.feat/mean_func -Tmean /mnt/d/SBSN/Data/Spine/$s/func/level_two_FLOB123.gfeat/cope1.feat/mean_func


            fslmerge -t /mnt/d/SBSN/Data/Spine/$s/func/level_two_FLOB1234.gfeat/cope1.feat/mean_func /mnt/d/SBSN/Data/Spine/$s/func/func1/level_one_FLOB.feat/mean_func /mnt/d/SBSN/Data/Spine/$s/func/func2/level_one_FLOB.feat/mean_func /mnt/d/SBSN/Data/Spine/$s/func/func3/level_one_FLOB.feat/mean_func /mnt/d/SBSN/Data/Spine/$s/func/func4/level_one_FLOB.feat/mean_func
            fslmaths /mnt/d/SBSN/Data/Spine/$s/func/level_two_FLOB1234.gfeat/cope1.feat/mean_func -Tmean /mnt/d/SBSN/Data/Spine/$s/func/level_two_FLOB1234.gfeat/cope1.feat/mean_func


            fslmerge -t /mnt/d/SBSN/Data/Spine/$s/func/level_two_FLOB12345.gfeat/cope1.feat/mean_func /mnt/d/SBSN/Data/Spine/$s/func/func1/level_one_FLOB.feat/mean_func /mnt/d/SBSN/Data/Spine/$s/func/func2/level_one_FLOB.feat/mean_func /mnt/d/SBSN/Data/Spine/$s/func/func3/level_one_FLOB.feat/mean_func /mnt/d/SBSN/Data/Spine/$s/func/func4/level_one_FLOB.feat/mean_func /mnt/d/SBSN/Data/Spine/$s/func/func5/level_one_FLOB.feat/mean_func
            fslmaths /mnt/d/SBSN/Data/Spine/$s/func/level_two_FLOB12345.gfeat/cope1.feat/mean_func -Tmean /mnt/d/SBSN/Data/Spine/$s/func/level_two_FLOB12345.gfeat/cope1.feat/mean_func


            fslmerge -t /mnt/d/SBSN/Data/Spine/$s/func/level_two_FLOB123456.gfeat/cope1.feat/mean_func /mnt/d/SBSN/Data/Spine/$s/func/func1/level_one_FLOB.feat/mean_func /mnt/d/SBSN/Data/Spine/$s/func/func2/level_one_FLOB.feat/mean_func /mnt/d/SBSN/Data/Spine/$s/func/func3/level_one_FLOB.feat/mean_func /mnt/d/SBSN/Data/Spine/$s/func/func4/level_one_FLOB.feat/mean_func /mnt/d/SBSN/Data/Spine/$s/func/func5/level_one_FLOB.feat/mean_func /mnt/d/SBSN/Data/Spine/$s/func/func6/level_one_FLOB.feat/mean_func
            fslmaths /mnt/d/SBSN/Data/Spine/$s/func/level_two_FLOB123456.gfeat/cope1.feat/mean_func -Tmean /mnt/d/SBSN/Data/Spine/$s/func/level_two_FLOB123456.gfeat/cope1.feat/mean_func

            tput setaf 6; echo "Done!" 
            tput sgr0;
        

        elif [ "$ind" == "2" ]; then


            allFLOB=`ls -d $DIREC$s/func/level_two_FLOB*`

            echo $allFLOB
            for b in $allFLOB; do # Loop through all files

                /home/sensel/fsl/bin/featquery 1 $b/cope1.feat 1  stats/pe1 featquery -p -s -b $b/cope1.feat/cluster_mask_zstat1.nii.gz
                
            done


        fi

done

# then
#fslmerge -t mean_func /mnt/d/SBSN/Data/Brain/SBSN_H_001/func/func1/level_one_FLOB.feat/reg_standard/mean_func /mnt/d/SBSN/Data/Brain/SBSN_H_001/func/func2/level_one_FLOB.feat/reg_standard/mean_func

#/mnt/d/SBSN/Data/Spine/$s/func/func1/level_one_FLOB.feat/mean_func
#/mnt/d/SBSN/Data/Spine/$s/func/func2/level_one_FLOB.feat/mean_func
##/mnt/d/SBSN/Data/Spine/$s/func/func3/level_one_FLOB.feat/mean_func
#/mnt/d/SBSN/Data/Spine/$s/func/func4/level_one_FLOB.feat/mean_func
#/mnt/d/SBSN/Data/Spine/$s/func/func5/level_one_FLOB.feat/mean_func
#/mnt/d/SBSN/Data/Spine/$s/func/func6/level_one_FLOB.feat/mean_func

#/mnt/d/SBSN/Data/Spine/$s/func/level_two_FLOB12.gfeat
#/mnt/d/SBSN/Data/Spine/$s/func/level_two_FLOB123.gfeat
#/mnt/d/SBSN/Data/Spine/$s/func/level_two_FLOB1234.gfeat
#/mnt/d/SBSN/Data/Spine/$s/func/level_two_FLOB12345.gfeat
#/mnt/d/SBSN/Data/Spine/$s/func/level_two_FLOB123456.gfeat

#fslmaths mean_func -Tmean mean_func
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