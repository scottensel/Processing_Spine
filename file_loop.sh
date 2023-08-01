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

# List of subjects
declare -a sub=("SBSN_S_001" "SBSN_S_002")

# Path of the folder containing all data
DIREC="/mnt/d/SBSN/Data/"

# For each subject
for s in "${sub[@]}"; do
	cd $DIREC$s"/func/"

    for d in */; do
        # Will print */ if no directories are available
        cd $DIREC$s"/func/"$d
        echo "$DIREC$s"/func/"$d"
    done

    tput setaf 2; echo "Done!"
    tput sgr0;

done
