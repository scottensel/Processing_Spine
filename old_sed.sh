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


# List of subjects
declare -a sub=("SBSN_S_002")

# Path of the folder containing all data
DIREC="/mnt/d/SBSN/Data/"

tput setaf 6; 
echo -n "Enter the index of the step to perform (1 = Prepare for GLM, 2 = Prepare for iCAP): "
tput sgr0;
read ind

# For each subject
for s in "${sub[@]}"; do

		cd $DIREC$s"/func/"

        for d in {1..6}; do

            # Will print */ if no directories are available
            cd $DIREC$s"/func/func"$d"/"

            tput setaf 2; echo "Noise regression started for " $s"/func/func"$d	

                    echo "Generate motion outliers..."
                    tput sgr0

                    # Generate EV for outliers
                    if [ ! -f outliers.png ]; then
                             ## generate the EVs (regressors) for moco fmri and create outliers.txt using the spinal cord mask
                             fsl_motion_outliers -i fmri_spine_moco.nii.gz -o outliers.txt —m fmri_spine_moco_mean_seg_corr.nii.gz -p outliers.png --dvars --nomoco
                    fi

                    # Copy header information from moco functional to moco parameters
                    fslcpgeom fmri_spine_moco.nii.gz moco_params_x.nii.gz
                    fslcpgeom fmri_spine_moco.nii.gz moco_params_y.nii.gz
                    # copying header information from mfmri to moco params

                    tput setaf 2; echo "Prepare nuisance regressors file..."
                    tput sgr0

            ## creates text file will all the regressor files
            ls -1 `${FSLDIR}/bin/imglob -extensions ${DIREC}${s}/physio/physio${d}/${s}ev0*` > regressors_evlist.txt

            ## adds the moco_params files to the txt file
            echo "$DIREC$s"/func/func"$d"/moco_params_x.nii.gz"" >> regressors_evlist.txt # Add motion parameters (x)
            echo "$DIREC$s"/func/func"$d"/moco_params_y.nii.gz"" >> regressors_evlist.txt # Add motion parameters (y) 

            ## Generate fsf file from template
            ## path is relevant to folder we are in which is sub-XX/func
            for i in "../../../template/template_design.fsf"; do
                ## Include outliers as regressors if needed

                # 1 - PREPARE .fsf files properly
                if [ "$ind" == "1" ]; then
                    tput setaf 2; echo "Prepare first level analysis for GLM " $s"/func/func"$d
                    tput sgr0; 

                    #-e 's@set fmri(threshmask) ""@'"set fmri(threshmask) \"$DIREC$s"/func/func"$d"/fmri_spine_moco_mean_seg_corr.nii.gz\"'@g' \

                    if [ -f $DIREC$s"/func/func"$d"/outliers.txt" ]; then

                        # this is editing the text of the files
                        sed -e 's@PNMPATH@'$DIREC$s"/func/func"$d"/regressors_evlist.txt"'@g' \
                                            -e 's@OUTDIR@'"level_one"'@g' \
                                            -e 's@DATAPATH@'$DIREC$s"/func/func"$d"/fmri_spine_moco.nii.gz"'@g' \
                                            -e 's@FILT@'"0"'@g' \
                                            -e 's@OUTLYN@'"1"'@g' \
                                            -e 's@NPTS@'"$(fslnvols $DIREC$s"/func/func"$d"/fmri_spine_moco.nii.gz")"'@g' \
                                            -e 's@2.500000@'"2.200000"'@g' \
                                            -e 's@set fmri(te) 35@'"set fmri(te) 30"'@g' \
                                            -e 's@3.1@'"1.5"'@g' \
                                            -e 's@set fmri(overwrite_yn) 0@'"set fmri(overwrite_yn) 1"'@g' \
                                            -e 's@set fmri(evtitle1) ""@'"set fmri(evtitle1) \"GripTask\""'@g' \
                                            -e 's@set fmri(shape1) 10@'"set fmri(shape1) 3"'@g' \
                                            -e 's@set fmri(convolve1) 0@'"set fmri(convolve1) 2"'@g' \
                                            -e 's@set fmri(threshmask) ""@'"set fmri(threshmask) \"$DIREC$s"/func/func"$d"/fmri_spine_moco_mean_seg_corr.nii.gz\"'@g' \
                                            -e 's@set fmri(deriv_yn1) 0@'"set fmri(deriv_yn1) 0\n\n# Custom EV file (EV 1)\nset fmri(custom1) \""$DIREC$s"/task/task"$d"/events.txt\"\n\n# Gamma sigma (EV 1)\nset fmri(gammasigma1) 3\n\n# Gamma delay (EV 1)\nset fmri(gammadelay1) 6"'@g' \
                                            -e 's@OUTLPATH@'$DIREC$s"/func/func"$d"/outliers.txt"'@g'  <$i> design_levelone.fsf


                    else

                        sed -e 's@PNMPATH@'$DIREC$s"/func/func"$d"/regressors_evlist.txt"'@g' \
                                            -e 's@OUTDIR@'"level_one"'@g' \
                                            -e 's@DATAPATH@'$DIREC$s"/func/func"$d"/fmri_spine_moco.nii.gz"'@g' \
                                            -e 's@FILT@'"0"'@g' \
                                            -e 's@OUTLYN@'"0"'@g' \
                                            -e 's@2.500000@'"2.200000"'@g' \
                                            -e 's@set fmri(te) 35@'"set fmri(te) 30"'@g' \
                                            -e 's@3.1@'"1.5"'@g' \
                                            -e 's@set fmri(overwrite_yn) 0@'"set fmri(overwrite_yn) 1"'@g' \
                                            -e 's@set fmri(evtitle1) ""@'"set fmri(evtitle1) "GripTask""'@g' \
                                            -e 's@set fmri(shape1) 10@'"set fmri(shape1) 3"'@g' \
                                            -e 's@set fmri(convolve1) 0@'"set fmri(convolve1) 2"'@g' \
                                            -e 's@set fmri(deriv_yn1) 0@'"set fmri(deriv_yn1) 0\n\n# Custom EV file (EV 1)\nset fmri(custom1) \""$DIREC$s"/task/task"$d"/events.txt\"\n\n# Gamma sigma (EV 1)\nset fmri(gammasigma1) 3\n\n# Gamma delay (EV 1)\nset fmri(gammadelay1) 6"'@g' \
                                            -e 's@NPTS@'"$(fslnvols $DIREC$s"/func/func"$d"/fmri_spine_moco.nii.gz")"'@g'  <$i> design_levelone.fsf

                    fi

                # 2 - GENERATE EVS
                elif [ "$ind" == "2" ]; then
                    tput setaf 2; echo "Prepare for iCAP " $s"/func/func"$d
                    tput sgr0; 

                    if [ -f $DIREC$s"/func/func"$d"/outliers.txt" ]; then

                        # this section confuses me a lot
                        # this is editing the text of the files
                        sed -e 's@PNMPATH@'$DIREC$s"/func/func"$d"/regressors_evlist.txt"'@g' \
                                            -e 's@OUTDIR@'"icap_prep"'@g' \
                                            -e 's@DATAPATH@'$DIREC$s"/func/func"$d"/fmri_spine_moco.nii.gz"'@g' \
                                            -e 's@FILT@'"0"'@g' \
                                            -e 's@OUTLYN@'"1"'@g' \
                                            -e 's@NPTS@'"$(fslnvols $DIREC$s"/func/func"$d"/fmri_spine_moco.nii.gz")"'@g' \
                                            -e 's@2.500000@'"2.200000"'@g' \
                                            -e 's@set fmri(te) 35@'"set fmri(te) 30"'@g' \
                                            -e 's@3.1@'"1.5"'@g' \
                                            -e 's@OUTLPATH@'$DIREC$s"/func/func"$d"/outliers.txt"'@g'  <$i> design_icapprep.fsf


                    else
                        sed -e 's@PNMPATH@'$DIREC$s"/func/func"$d"/regressors_evlist.txt"'@g' \
                                            -e 's@OUTDIR@'"icap_prep"'@g' \
                                            -e 's@DATAPATH@'$DIREC$s"/func/func"$d"/fmri_spine_moco.nii.gz"'@g' \
                                            -e 's@FILT@'"0"'@g' \
                                            -e 's@OUTLYN@'"0"'@g' \
                                            -e 's@2.500000@'"2.200000"'@g' \
                                            -e 's@set fmri(te) 35@'"set fmri(te) 30"'@g' \
                                            -e 's@3.1@'"1.5"'@g' \
                                            -e 's@NPTS@'"$(fslnvols $DIREC$s"/func/func"$d"/fmri_spine_moco.nii.gz")"'@g'  <$i> design_icapprep.fsf

                    fi

                fi

                # created a design_levelone.fsf based on the template.fsf
            done

            if [ "$ind" == "1" ]; then
                tput setaf 2; echo "Run first level analysis for " $s"/func/func"$d
                tput sgr0; 

                # Run the analysis using the fsf file
                feat design_levelone.fsf
            
            elif [ "$ind" == "2" ]; then
                tput setaf 2; echo "Run noise regression for " $s"/func/func"$d
                tput sgr0; 

                # Run the analysis using the fsf file
                feat design_icapprep.fsf

                # Copy geometry to residuals for TA
                cp icap_prep.feat/stats/res4d.nii.gz fmri_spine_moco_denoised.nii.gz
                fslcpgeom fmri_spine_moco.nii.gz fmri_spine_moco_denoised.nii.gz

            fi

        
        done
		tput setaf 2; echo "Done!" 
        	tput sgr0;				 			
done

