#!/bin/bash
#
# SPiCiCAP framework - Preprocessing
# July 2020
#
# Preparation of physiological noise regressors
#
# Requirements: FSL
#
#

#
# Detection of cardiac peaks should always be checked manually (see PNM user guide)
#
#

# List of subjects
declare -a sub=("SBSN_S_001")

# Path of the folder containing all data
DIREC="/mnt/d/SBSN/Data/"

tput setaf 6; 
echo -n "Enter the index of the step to perform (1 = Prepare recordings, 2 = Generate EVs): "
tput sgr0;
read ind

# For each subject
for s in "${sub[@]}"; do
		cd $DIREC$s"/physio/"
			# 1 - PREPARE PHYSIOLOGICAL RECORDINGS
			if [ "$ind" == "1" ]; then
				tput setaf 2; echo "Prepare physiological recordings in " $s
				tput sgr0; 

                # this is just setting up the text files to make sure they are good
				$FSLDIR/bin/fslFixText $s".txt" $s"_input.txt"
				$FSLDIR/bin/pnm_stage1 -i $s"_input.txt" -o $s -s 100 --tr=2.5 --smoothcard=0.3 --smoothresp=0.1 --resp=1 --cardiac=2 --trigger=3
				
			# 2 - GENERATE EVS
			elif [ "$ind" == "2" ]; then
				tput setaf 2; echo "EVs generation in " $s
				tput sgr0; 
                
                # creating EVs (regressors) to model the physiological noise, of cardiac and respiration 
				$FSLDIR/bin/pnm_evs -i $DIREC$s"/func/fmri_spine_moco.nii.gz" -c $s"_card.txt" -r $s"_resp.txt" -o $s --tr=2.5 --oc=4 --or=4 --multc=2 --multr=2 --csfmask=$DIREC$s"/func/Segmentation/mask_csf.nii.gz" --sliceorder=up --slicedir=z
                

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

                    # creats an ev_text list of all the ev files
					ls -1 `${FSLDIR}/bin/imglob -extensions ${DIREC}${s}/physio/${s}ev0*` > $s"_evlist.txt"
			
			else
				tput setaf 1; 
				echo "Index not valid (should be 1 or 2)"
				tput sgr0; 
			fi
		
		tput setaf 2; echo "Done!" 
                tput sgr0;
done

