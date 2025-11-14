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
# Detection of cardiac peaks should always be checked manually (see PNM user guide)
#
# You have to run the matlab functions to format the physio data first

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
echo -n "Enter the index of the step to perform (1 = Prepare recordings, 2 = Generate EVs): "
tput sgr0;
read ind

# For each subject
for s in "${sub[@]}"; do
        cd $DIREC$s"/physio/"
        for d in "${myFunc[@]}"; do
            # Will print */ if no directories are available
            cd $DIREC$s"/physio/physio"$d"/"
            
			# 1 - PREPARE PHYSIOLOGICAL RECORDINGS
			if [ "$ind" == "1" ]; then
				tput setaf 2; echo "Prepare physiological recordings in " $s"/func/func"$d
				tput sgr0; 
                
                # this is just setting up the text files to make sure they are good
				$FSLDIR/bin/pnm_stage1 -i physio.txt -o $s -s 100 --tr=2.2 --cardiac=1 --resp=2 --trigger=3 -v

			# 2 - GENERATE EVS
			elif [ "$ind" == "2" ]; then
				tput setaf 2; echo "EVs generation in " $s"/func/func"$d
				tput sgr0; 

                if [ ! -f $DIREC$s"/func/func"$d"/fmri_spine_moco_mean_CSF_seg_corr.nii.gz" ]; then

                    cp $DIREC$s"/func/func"$d"/fmri_spine_moco_mean_CSF_seg.nii.gz" $DIREC$s"/func/func"$d"/fmri_spine_moco_mean_CSF_seg_corr.nii.gz"

                fi

                # creating EVs (regressors) to model the physiological noise, of cardiac and respiration 
                $FSLDIR/bin/pnm_evs -i $DIREC$s"/func/func"$d"/fmri_spine_moco.nii.gz" -c $s"_card.txt" -r $s"_resp.txt" -o $s --tr=2.2 --oc=4 --or=4 --multc=2 --multr=2 --csfmask=$DIREC$s"/func/func"$d"/fmri_spine_moco_mean_CSF_seg_corr.nii.gz" --sliceorder=interleaved_up --slicetiming=slice_times.txt -v
                #$FSLDIR/bin/pnm_evs -i $DIREC$s"/func/func"$d"/fmri_spine_moco.nii.gz" -c $s"_card.txt" -r $s"_resp.txt" -o $s --tr=2.2 --oc=4 --or=4 --multc=2 --multr=2 --csfmask=$DIREC$s"/func/func"$d"/fmri_spine_moco_mean_CSF_seg_corr.nii.gz" --slicetiming=slice_matrix.txt -v

                # above uses the mean fmri files
                # -c is cardiac 
                # -r is respiratory text
                # -o is output file name for EV matrix
                # -tr is tr of fmri
                # --oc order of cardiac regressos
                # --or order of resp regressos
                # --multc is order of multiplicative cardiac terms 
                # --multr is order of multiplicative resp terms 
                # --csfmask filename of csf mask image (and generate csf regressor)
                # --sliceorder specify slice ordering (up/down/interleaved_up/interleaved_down)
                # --slicedir specify slice direction (x/y/z) - default is z
                # --slicetiming   specify slice timing via an external file
			
			else
				tput setaf 1; 
				echo "Index not valid (should be 1 or 2)"
				tput sgr0; 
			fi
		
        done

		tput setaf 2; echo "Done!" 
                tput sgr0;

done

echo
echo "${sub[@]}"
echo "${myFunc[@]}"
echo

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