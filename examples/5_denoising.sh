#!/bin/bash

#
# SPiCiCAP framework - Preprocessing
# July 2020
#
# Noise regression (motion, pnm, etc.)
#
# Requirements: FSL
#
#

# List of subjects
#declare -a sub=("sub-01" "sub-02" "sub-03" "sub-04" "sub-05" "sub-06" "sub-07" "sub-08" "sub-09" "sub-10" "sub-11" "sub-12" "sub-13" "sub-14" "sub-15" "sub-16" "sub-17" "sub-18" "sub-19")
declare -a sub=("sub-01")

# Path of the folder containing all data
#DIREC="/home/sce30/Documents/Research/NeuronData/SPiCiCAP_Dataset_Code/Data/"
DIREC="/mnt/d/NeuronData/SPiCiCAP_Dataset_Code/Data"
#DIREC="/home/sce30/Documents/Research/Data/Spine/"

# For each subject
for s in "${sub[@]}"; do
		cd $DIREC$s"/func/"

		tput setaf 2; echo "Noise regression started for " $s	

                echo "Generate motion outliers..."
                tput sgr0
		
                # Generate EV for outliers
                if [ ! -f outliers.png ]; then
                         ## generate the EVs (regressors) for moco fmri and create outliers.txt using the spinal cord mask
                         fsl_motion_outliers -i mfmri.nii.gz -o outliers.txt —m ../Segmentation/mask_sc.nii.gz -p outliers.png --dvars --nomoco
                fi

		# Copy header information from moco functional to moco parameters
		fslcpgeom mfmri.nii.gz Moco/moco_params_x.nii.gz
        fslcpgeom mfmri.nii.gz Moco/moco_params_y.nii.gz
        # copying header information from mfmri to moco params

                tput setaf 2; echo "Prepare nuisance regressors file..."
                tput sgr0

        ## creates text file will all the regressor files
		ls -1 `${FSLDIR}/bin/imglob -extensions ${DIREC}${s}/physio/${s}ev0*` > regressors_evlist.txt
        
        ## adds the moco_params files to the txt file
        echo "$DIREC$s"/func/Moco/moco_params_x.nii.gz"" >> regressors_evlist.txt # Add motion parameters (x)
        echo "$DIREC$s"/func/Moco/moco_params_y.nii.gz"" >> regressors_evlist.txt # Add motion parameters (y) 
        

      
		## Generate fsf file from templat
		##
		## NOTE: adapt path for template
        ## path is relevant to folder we are in which is sub-XX/func
		for i in "../../template/template_noiseregression.fsf"; do
			## Include outliers as regressors if needed
			if [ -f $DIREC$s"/func/outliers.txt" ]; then
            
            # this section confuses me a lot
            # this is editing the text of the files
			sed -e 's@PNMPATH@'$DIREC$s"/func/regressors_evlist.txt"'@g' \
                                -e 's@OUTDIR@'"noise_regression"'@g' \
                                -e 's@DATAPATH@'$DIREC$s"/func/mfmri.nii.gz"'@g' \
                            	-e 's@FILT@'"0"'@g' \
                                -e 's@OUTLYN@'"1"'@g' \
                            	-e 's@NPTS@'"$(fslnvols $DIREC$s"/func/mfmri.nii.gz")"'@g' \
                                -e 's@OUTLPATH@'$DIREC$s"/func/outliers.txt"'@g'  <$i> design_noiseregression.fsf
			else
			sed -e 's@PNMPATH@'$DIREC$s"/func/Regressors/regressors_evlist.txt"'@g' \
                                -e 's@OUTDIR@'"noise_regression"'@g' \
                                -e 's@DATAPATH@'$DIREC$s"/func/mfmri.nii.gz"'@g' \
                                -e 's@FILT@'"0"'@g' \
                                -e 's@OUTLYN@'"0"'@g' \
                                -e 's@NPTS@'"$(fslnvols $DIREC$s"/func/mfmri.nii.gz")"'@g'  <$i> design_noiseregression.fsf
			fi

            # created a design_noiseregression.fsf based on the template.fsf
 		done

		# Run the analysis using the fsf file
 		fsl5.0-feat design_noiseregression.fsf
 	
		# Copy geometry to residuals
		cp noise_regression.feat/stats/res4d.nii.gz mfmri_denoised.nii.gz
		fslcpgeom mfmri.nii.gz mfmri_denoised.nii.gz

		tput setaf 2; echo "Done!" 
        	tput sgr0;				 			
done

