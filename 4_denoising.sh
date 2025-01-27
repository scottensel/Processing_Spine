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
echo -n "Enter the index of the step to perform (1 = Prepare for GLM, 2 = Prepare for GLM forces, 3 = Prepare for GLM smooth forces, 4 = Denoise from Physio, 5 = 1st Level 2nd FEAT, 6 = 1st Level 2nd FEAT smoothed): "
tput sgr0;
read ind


# For each subject
for s in "${sub[@]}"; do

		cd $DIREC$s"/func/"

        for d in "${myFunc[@]}"; do

            # Will print */ if no directories are available
            cd $DIREC$s"/func/func"$d"/"

            tput setaf 2; echo "Noise regression started for " $DIREC$s"/func/func"$d	
            if [ "$ind" -lt "5" ]; then

                # Generate EV for outliers
                #if [ ! -f outliers.png ]; then
    
                echo "Generate motion outliers..."
                tput sgr0
    
                ## generate the EVs (regressors) for moco fmri and create outliers.txt using the spinal cord mask
                fsl_motion_outliers -i fmri_spine_moco.nii.gz -o outliers.txt —m fmri_spine_moco_mean_seg_corr.nii.gz -p outliers.png --dvars --nomoco
    
                #fi
    
                tput setaf 2; echo "Prepare nuisance regressors file..."
                tput sgr0
    
                ## creates text file will all the regressor files
                ls -1 `${FSLDIR}/bin/imglob -extensions ${DIREC}${s}/physio/physio${d}/${s}ev0*` > regressors_evlist.txt

    
                # Copy header information from moco functional to moco parameters
                sct_image -i fmri_spine_moco.nii.gz -copy-header moco_params_x.nii.gz -o moco_params_x.nii.gz
                sct_image -i fmri_spine_moco.nii.gz -copy-header moco_params_y.nii.gz -o moco_params_y.nii.gz
    
                ## adds the moco_params files to the txt file
                echo "$DIREC$s"/func/func"$d"/moco_params_x.nii.gz"" >> regressors_evlist.txt # Add motion parameters (x)
                echo "$DIREC$s"/func/func"$d"/moco_params_y.nii.gz"" >> regressors_evlist.txt # Add motion parameters (y)

            fi

            # changes the template file that it loops through
            if [ "$ind" == "1" ]; then
                templateFile="template_design.fsf"
            elif [ "$ind" == "2" ]; then
                templateFile="template_design_force.fsf"
            elif [ "$ind" == "3" ]; then
                templateFile="template_design_force.fsf"
            elif [ "$ind" == "4" ]; then
                templateFile="template_design.fsf"
            elif [ "$ind" == "5" ]; then
                templateFile="template_design_FLOB.fsf"
            elif [ "$ind" == "6" ]; then
                templateFile="template_design_force_FLOB.fsf"
            elif [ "$ind" == "7" ]; then
                templateFile="template_design_force_FLOB.fsf"
            fi

            ## Generate fsf file from template
            ## path is relevant to folder we are in which is sub-XX/fun

            ## for i in "../../../template/template_design.fsf"; do
            for i in "../../../template/"$templateFile; do
                ## Include outliers as regressors if needed

                # 1 - PREPARE .fsf files properly
                if [ "$ind" == "1" ]; then
                    tput setaf 2; echo "Prepare first level analysis for GLM " $s"/func/func"$d
                    tput sgr0; 
         
                    # this is editing the text of the files
                    sed -e 's@PNMPATH@'$DIREC$s"/func/func"$d"/regressors_evlist.txt"'@g' \
                                        -e 's@OUTDIR@'"level_one"'@g' \
                                        -e 's@DATAPATH@'$DIREC$s"/func/func"$d"/fmri_spine_moco.nii.gz"'@g' \
                                        -e 's@OUTLYN@'"1"'@g' \
                                        -e 's@NPTS@'"$(fslnvols $DIREC$s"/func/func"$d"/fmri_spine_moco.nii.gz")"'@g' \
                                        -e 's@EV_TITLE@'GripTask'@g' \
                                        -e 's@EV_FILE@'$DIREC$s"/task/task"$d"/events.txt"'@g' \
                                        -e 's@OUTLPATH@'$DIREC$s"/func/func"$d"/outliers.txt"'@g' <$i> design_levelone.fsf

                elif [ "$ind" == "2" ]; then
                    tput setaf 2; echo "Prepare first level analysis for GLM " $s"/func/func"$d
                    tput sgr0; 
         
                    # this is editing the text of the files
                    sed -e 's@PNMPATH@'$DIREC$s"/func/func"$d"/regressors_evlist.txt"'@g' \
                                        -e 's@OUTDIR@'"level_one_force"'@g' \
                                        -e 's@DATAPATH@'$DIREC$s"/func/func"$d"/fmri_spine_moco.nii.gz"'@g' \
                                        -e 's@OUTLYN@'"1"'@g' \
                                        -e 's@NPTS@'"$(fslnvols $DIREC$s"/func/func"$d"/fmri_spine_moco.nii.gz")"'@g' \
                                        -e 's@EV_TITLE1@'20'@g' \
                                        -e 's@EV_FILE1@'$DIREC$s"/task/task"$d"/force20.txt"'@g' \
                                        -e 's@EV_TITLE2@'45'@g' \
                                        -e 's@EV_FILE2@'$DIREC$s"/task/task"$d"/force45.txt"'@g' \
                                        -e 's@EV_TITLE3@'70'@g' \
                                        -e 's@EV_FILE3@'$DIREC$s"/task/task"$d"/force70.txt"'@g' \
                                        -e 's@OUTLPATH@'$DIREC$s"/func/func"$d"/outliers.txt"'@g' <$i> design_levelone_force.fsf

                elif [ "$ind" == "3" ]; then
                    tput setaf 2; echo "Prepare first level analysis for GLM " $s"/func/func"$d
                    tput sgr0; 
         
                    # this is editing the text of the files
                    sed -e 's@PNMPATH@'$DIREC$s"/func/func"$d"/regressors_evlist.txt"'@g' \
                                        -e 's@OUTDIR@'"level_one_force_smooth"'@g' \
                                        -e 's@DATAPATH@'$DIREC$s"/func/func"$d"/fmri_spine_moco_smooth.nii.gz"'@g' \
                                        -e 's@OUTLYN@'"1"'@g' \
                                        -e 's@NPTS@'"$(fslnvols $DIREC$s"/func/func"$d"/fmri_spine_moco.nii.gz")"'@g' \
                                        -e 's@EV_TITLE1@'20'@g' \
                                        -e 's@EV_FILE1@'$DIREC$s"/task/task"$d"/force20.txt"'@g' \
                                        -e 's@EV_TITLE2@'45'@g' \
                                        -e 's@EV_FILE2@'$DIREC$s"/task/task"$d"/force45.txt"'@g' \
                                        -e 's@EV_TITLE3@'70'@g' \
                                        -e 's@EV_FILE3@'$DIREC$s"/task/task"$d"/force70.txt"'@g' \
                                        -e 's@OUTLPATH@'$DIREC$s"/func/func"$d"/outliers.txt"'@g' <$i> design_levelone_force_smooth.fsf

                elif [ "$ind" == "4" ]; then
                    tput setaf 2; echo "Prepare for by denoising PHYSIO " $s"/func/func"$d
                    tput sgr0; 

                    sed -e 's@PNMPATH@'$DIREC$s"/func/func"$d"/regressors_evlist.txt"'@g' \
                                        -e 's@OUTDIR@'"physio_denoised"'@g' \
                                        -e 's@DATAPATH@'$DIREC$s"/func/func"$d"/fmri_spine_moco.nii.gz"'@g' \
                                        -e 's@OUTLYN@'"1"'@g' \
                                        -e 's@NPTS@'"$(fslnvols $DIREC$s"/func/func"$d"/fmri_spine_moco.nii.gz")"'@g' \
                                        -e 's@EV_TITLE@'""'@g' \
                                        -e 's@EV_FILE@'""'@g' \
                                        -e 's@OUTLPATH@'$DIREC$s"/func/func"$d"/outliers.txt"'@g' <$i> design_denoised.fsf

                elif [ "$ind" == "5" ]; then

                    # this is editing the text of the files
                    sed -e 's@OUTDIR@'"level_one_FLOB"'@g' \
                                        -e 's@DATAPATH@'$DIREC$s"/func/func"$d"/fmri_spine_moco_denoised_plusmean_smooth.nii.gz"'@g' \
                                        -e 's@OUTLYN@'"1"'@g' \
                                        -e 's@NPTS@'"$(fslnvols $DIREC$s"/func/func"$d"/fmri_spine_moco.nii.gz")"'@g' \
                                        -e 's@EV_TITLE@'GripTask'@g' \
                                        -e 's@EV_FILE@'$DIREC$s"/task/task"$d"/events.txt"'@g' <$i> design_levelone_FLOB.fsf

                elif [ "$ind" == "6" ]; then
                    tput setaf 2; echo "Prepare for by FEAT on denoised data " $s"/func/func"$d
                    tput sgr0;

                    sed -e 's@OUTDIR@'"level_one_force_FLOB"'@g' \
                                        -e 's@DATAPATH@'$DIREC$s"/func/func"$d"/fmri_spine_moco_denoised_plusmean.nii.gz"'@g' \
                                        -e 's@NPTS@'"$(fslnvols $DIREC$s"/func/func"$d"/fmri_spine_moco.nii.gz")"'@g' \
                                        -e 's@EV_TITLE1@'20'@g' \
                                        -e 's@EV_FILE1@'$DIREC$s"/task/task"$d"/force20.txt"'@g' \
                                        -e 's@EV_TITLE2@'45'@g' \
                                        -e 's@EV_FILE2@'$DIREC$s"/task/task"$d"/force45.txt"'@g' \
                                        -e 's@EV_TITLE3@'70'@g' \
                                        -e 's@EV_FILE3@'$DIREC$s"/task/task"$d"/force70.txt"'@g' <$i> design_levelone_force_FLOB.fsf

                elif [ "$ind" == "7" ]; then
                    tput setaf 2; echo "Prepare for by FEAT on smoothed denoised data " $s"/func/func"$d
                    tput sgr0; 

                    sed -e 's@OUTDIR@'"level_one_force_smooth_FLOB"'@g' \
                                        -e 's@DATAPATH@'$DIREC$s"/func/func"$d"/fmri_spine_moco_denoised_plusmean_smooth.nii.gz"'@g' \
                                        -e 's@NPTS@'"$(fslnvols $DIREC$s"/func/func"$d"/fmri_spine_moco.nii.gz")"'@g' \
                                        -e 's@EV_TITLE1@'20'@g' \
                                        -e 's@EV_FILE1@'$DIREC$s"/task/task"$d"/force20.txt"'@g' \
                                        -e 's@EV_TITLE2@'45'@g' \
                                        -e 's@EV_FILE2@'$DIREC$s"/task/task"$d"/force45.txt"'@g' \
                                        -e 's@EV_TITLE3@'70'@g' \
                                        -e 's@EV_FILE3@'$DIREC$s"/task/task"$d"/force70.txt"'@g' <$i> design_levelone_force_smooth_FLOB.fsf


                fi

                # created a design_levelone.fsf based on the template.fsf
            done

            if [ "$ind" == "1" ]; then
                tput setaf 2; echo "Run first level analysis for " $s"/func/func"$d
                tput sgr0; 
                
                # Run the analysis using the fsf file
                feat design_levelone.fsf
            
            elif [ "$ind" == "2" ]; then
                tput setaf 2; echo "Run first level force analysis for " $s"/func/func"$d
                tput sgr0; 
                
                # Run the analysis using the fsf file
                feat design_levelone_force.fsf

            elif [ "$ind" == "3" ]; then
                tput setaf 2; echo "Run first level smoothed force analysis for " $s"/func/func"$d
                tput sgr0; 
                
                # Run the analysis using the fsf file
                feat design_levelone_force_smooth.fsf

            elif [ "$ind" == "4" ]; then
                tput setaf 2; echo "Run first level denoising " $s"/func/func"$d
                tput sgr0; 
                
                # Run the analysis using the fsf file
                feat design_denoised.fsf

                # Copy geometry to residuals for TA
                cp physio_denoised.feat/stats/res4d.nii.gz fmri_spine_moco_denoised.nii.gz
                
                # copy header info
                sct_image -i fmri_spine_moco.nii.gz -copy-header fmri_spine_moco_denoised.nii.gz -o fmri_spine_moco_denoised.nii.gz
    
                # need to add back in the mean image to use in FEAT again but not need for iCAP so we keep separate
                fslmaths fmri_spine_moco_denoised.nii.gz -add physio_denoised.feat/mean_func fmri_spine_moco_denoised_plusmean.nii.gz
                
                # copy header info
                sct_image -i fmri_spine_moco.nii.gz -copy-header fmri_spine_moco_denoised_plusmean.nii.gz -o fmri_spine_moco_denoised_plusmean.nii.gz

            elif [ "$ind" == "5" ]; then
                tput setaf 2; echo "Run first level FLOB force analysis for " $s"/func/func"$d
                tput sgr0; 
                
                # Run the analysis using the fsf file
                feat design_levelone_FLOB.fsf

            elif [ "$ind" == "6" ]; then
                tput setaf 2; echo "Run first level FLOB force analysis for " $s"/func/func"$d
                tput sgr0; 
                
                # Run the analysis using the fsf file
                feat design_levelone_force_FLOB.fsf

            elif [ "$ind" == "7" ]; then
                tput setaf 2; echo "Run first level smoothed FLOB force analysis for " $s"/func/func"$d
                tput sgr0; 
                
                # Run the analysis using the fsf file
                feat design_levelone_force_smooth_FLOB.fsf

            fi

        done	 			
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