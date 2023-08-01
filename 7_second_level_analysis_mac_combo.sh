#!/bin/bash

#
# SBSN framework - Preprocessing
# June 2022
#
# Preparation of functional spine data
#
# Requirements: FSL, Spinal Cord Toolbox 5.5
#
# THIS IS SPECIFICALLY FOR MAC SINCE ITS SUCKSSSSS

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

####################################################
# replace below with your path
####################################################
input="/mnt/d/SBSN/Processing_Spine/combinations.txt"

# For each subject
for s in "${sub[@]}"; do

		cd $DIREC$s"/func/"

        while IFS=$' \t\r\n' read -r -a myFunc; do

            naming=$(printf '%s' "${myFunc[@]}")

            j=1
            
            tput setaf 2; echo "Prepare template files for analysis..."
            tput sgr0;

            if [ ! -f design_leveltwo"$naming".fsf ]; then

                for runNum in {1..6}; do # this needs to be 1..6

                    ## Generate fsf file from template
                    if [ "$j" == "1" ]; then

                        for i in "../../template/second_level_template_mac.fsf"; do

                            sed -e 's@OUTDIR@'"level_two"$naming""'@g' -e 's@PATH1@'$DIREC$s"/func/func"${myFunc[((j-1))]}"/level_one.feat"'@g' \
                                -e 's@THRESH_MASK@'$DIREC"template/PAM50_cord.nii.gz"'@g' \
                                -e 's@NSUBJECTS@'${#myFunc[@]}'@g' <$i> design_leveltwo"$j".fsf 
                                #PAM50_mask_add.nii.gz
                        done

                    elif [ "$j" == "2" ]; then

                        for i in "design_leveltwo"$((j-1))".fsf"; do

                            sed -e 's@PATH2@'$DIREC$s"/func/func"${myFunc[((j-1))]}"/level_one.feat"'@g' -e 's@OUTLPATH@''@g' <$i> design_leveltwo"$j".fsf 

                        done

                        rm design_leveltwo"$((j-1))".fsf

                    # le less than or equal to
                    # lt less than
                    elif [ "$j" -le "${#myFunc[@]}" ] ; then

                        for i in "design_leveltwo"$((j-1))".fsf"; do

                            sed -e 's@'FEAT_PATH"$j"'@'"set feat_files("$j") "$DIREC$s"/func/func"${myFunc[((j-1))]}"/level_one.feat"'@g' \
                                -e 's@'EVG_PATH"$j"'@'"set fmri(evg"$j".1) 1"'@g' \
                                -e 's@'GROUPMEM_PATH"$j"'@'"set fmri(groupmem\."$j") 1"'@g' <$i> design_leveltwo"$j".fsf 

                        done

                        rm design_leveltwo"$((j-1))".fsf

                    else

                        for i in "design_leveltwo"$((j-1))".fsf"; do

                            sed -e 's@'FEAT_PATH"$j"'@''@g' \
                                -e 's@'EVG_PATH"$j"'@''@g' \
                                -e 's@'GROUPMEM_PATH"$j"'@''@g' <$i> design_leveltwo"$j".fsf 

                        done

                        rm design_leveltwo"$((j-1))".fsf 

                    fi           

                    ((j+=1));

                done

                mv design_leveltwo6.fsf design_leveltwo"$naming".fsf 

                tput setaf 2; echo "Run second level analysis"
                tput sgr0; 

                echo "Running . . . design_leveltwo"$naming".fsf"

                # Run the analysis using the fsf file
                feat design_leveltwo"$naming".fsf

                tput setaf 2; echo "Done!" 
                    tput sgr0;
            else

                echo "design_leveltwo"$naming".fsf . . . already ran"
            fi

    done < "$input"
		 			
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