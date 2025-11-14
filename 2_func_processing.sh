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
. /mnt/d/SMA/Processing_Spine/path_to_subjects.sh 


tput setaf 6; 
echo -n "Enter the index of the step to perform (1 = Crop Spine, 2 = MOCO, 3 = Initial Segmentation, 4 = Label Vertebrae , 5 = Registering to template, 6 = TSNR/SSNR): "
tput sgr0;
read ind


# For each subject
for s in "${sub[@]}"; do
		cd $DIREC$s"/func/"

        for d in "${myFunc[@]}"; do
            # Will print */ if no directories are available
            cd $DIREC$s"/func/func"$d"/"

			# 1 - CROP
			if [ "$ind" == "1" ]; then
				tput setaf 2; echo "Crop spine from brain " $DIREC$s"/func/func"$d"/"
				tput sgr0; 

                # this is for seperating the spine from the brain and then does some preprocessing steps

                # create mean image
                sct_maths -i fmri.nii.gz -mean t -o fmri_mean.nii.gz

                # do original crop
                sct_crop_image -i fmri_mean.nii.gz -g 1

                # create a binary image of the entire scan and crop the bold images
                fslmaths fmri_mean_crop.nii.gz -bin fmri_mean_crop_bin.nii.gz
                sct_crop_image -i fmri.nii.gz -m fmri_mean_crop_bin.nii.gz -o fmri_spine.nii.gz
                
                # this is for creating a slice timing file later for physio correction
                sct_crop_image -i fmri.nii.gz -m fmri_mean_crop_bin.nii.gz -o fmri_spine_slices.nii.gz -b 1000
                sct_maths -i fmri_spine_slices.nii.gz -mean t -o fmri_spine_slices_mean.nii.gz

                # get the centerline
                sct_get_centerline -i fmri_mean_crop.nii.gz -c t2

                # create a binary mask centerline
                sct_create_mask -i fmri_mean_crop.nii.gz -p centerline,fmri_mean_crop_centerline.nii.gz -size 35mm -o fmri_spine_mask.nii.gz


			# 2 - MOCO
			elif [ "$ind" == "2" ]; then
				tput setaf 2; echo "MOCO in " "$DIREC$s"/func/"$d"
				tput sgr0; 

                # motion correction flirt first
                mcflirt -in fmri_spine.nii.gz -out fmri_spine_mcflirt -refvol 0 -plots -report

                # do motion correction
                sct_fmri_moco -i fmri_spine_mcflirt.nii.gz -m fmri_spine_mask.nii.gz -g 1 -r 1
                
                cp fmri_spine_mcflirt_moco.nii.gz fmri_spine_moco.nii.gz 

                # remove extra volumes
                if [ $(fslnvols fmri_spine_moco.nii.gz) -eq 163 ]; then
                    #sct_image -i fmri_spine.nii.gz -remove-vol 0,1,2,3,4
                    sct_image -i fmri_spine_moco.nii.gz -remove-vol 0,1,2,3,4
                    sct_image -i moco_params_x.nii.gz -remove-vol 0,1,2,3,4
                    sct_image -i moco_params_y.nii.gz -remove-vol 0,1,2,3,4
                fi

                # now take the mean of the motion corrected spine image
                sct_maths -i fmri_spine_moco.nii.gz -mean t -o fmri_spine_moco_mean.nii.gz


			# 4 - INITIAL SEGMENTATION
			elif [ "$ind" == "3" ]; then
				tput setaf 2; echo "Create automatic segmentation and then manually correct " "$DIREC$s"/func/"$d"
				tput sgr0; 
            
                #sct_deepseg_sc -i fmri_spine_moco_mean.nii.gz -c t2s #-qc "$DIREC$s"/func/"$d"
                # now YOU MUST go in and manually correct this autosegmentation

                # use this as the segementation as its better
                sct_deepseg -i fmri_spine_moco_mean.nii.gz -o test_epi_seg.nii.gz -task seg_sc_epi

                # rename because getting CSF creates same file name
                #mv fmri_spine_moco_mean_seg.nii.gz fmri_spine_moco_mean_deepseg.nii.gz 

                # get CSF segmentation
                sct_propseg -i fmri_spine_moco_mean.nii.gz -c t2s -CSF
                
                #replace propseg with deepseg one
                rm fmri_spine_moco_mean_seg.nii.gz

                # just to make sure file gets overwritten properly
                sleep 20

                #mv fmri_spine_moco_mean_deepseg.nii.gz fmri_spine_moco_mean_seg.nii.gz  
                mv test_epi_seg.nii.gz fmri_spine_moco_mean_seg.nii.gz

                # dilate the mask and add one so that way we can add it to the CSF mask
                # this will allow us to bin it and remove overlapping segmentations and make it easier
                # to corrrect
                sct_maths -i fmri_spine_moco_mean_seg.nii.gz -dilate 1 -shape disk -dim 2 -o fmri_spine_moco_mean_seg_dilate.nii.gz #-qc $DIREC$s"/anat/"

                fslmaths fmri_spine_moco_mean_CSF_seg.nii.gz -sub fmri_spine_moco_mean_seg_dilate.nii.gz fmri_spine_moco_mean_CSF_seg.nii.gz

                sct_maths -i fmri_spine_moco_mean_CSF_seg.nii.gz -bin 0.5 -o fmri_spine_moco_mean_CSF_seg.nii.gz

                cp fmri_spine_moco_mean_seg.nii.gz fmri_spine_moco_mean_seg_corr.nii.gz
                cp fmri_spine_moco_mean_CSF_seg.nii.gz fmri_spine_moco_mean_CSF_seg_corr.nii.gz
                # now manually correct the CSF mask
                				
                tput setaf 6; echo "NOW YOU MUST MANUAL CORRECT THESE SEGMENTATIONS"
				tput sgr0; 

			# 5 - LABEL THE VERTEBRAE
			elif [ "$ind" == "4" ]; then
				tput setaf 2; echo "Use manual segmentation in " "$DIREC$s"/func/"$d"
				tput sgr0; 

                # this is the first steps that work for registering the epi to the pam 50
                # Label at C2/C3 and C6/C7 discs (here we specify the known Z, which can vary across subjects)
                #sct_label_utils -i fmri_spine_moco_mean.nii.gz -create-viewer 3,4,5,6,7,8 -o fmri_labels.nii.gz \
                     #-msg "Click atthe posterior tip of C2/C3 to C7/T1 inter-vertebral disc"
                sct_label_utils -i fmri_spine_moco_mean.nii.gz -create-viewer 4,8 -o fmri_labels.nii.gz \
                      -msg "Click at the posterior tip of C3/C4 & C7/T1 inter-vertebral disc"

			# 6 - REGISTER TO TEMPLATE
			elif [ "$ind" == "5" ]; then
				tput setaf 2; echo "Use manual segmentation in " "$DIREC$s"/func/"$d"
				tput sgr0; 
                
                sct_image -i fmri_spine_moco_mean.nii.gz -set-qform-to-sform
                sct_image -i fmri_spine_moco_mean_seg_corr.nii.gz -set-qform-to-sform
                sct_image -i fmri_labels.nii.gz -set-qform-to-sform

                sct_image -i fmri_spine_moco_mean.nii.gz -set-sform-to-qform
                sct_image -i fmri_spine_moco_mean_seg_corr.nii.gz -set-sform-to-qform
                sct_image -i fmri_labels.nii.gz -set-sform-to-qform


                sct_label_vertebrae -i fmri_spine_moco_mean.nii.gz -s fmri_spine_moco_mean_seg_corr.nii.gz -c t2 -discfile fmri_labels.nii.gz #-qc $DIREC$s"/anat/"
               
                #sct_register_to_template -i fmri_spine_moco_mean.nii.gz -s fmri_spine_moco_mean_seg_corr.nii.gz -l fmri_labels.nii.gz -c t2s -ref subject -param step=1,type=seg,algo=centermass:step=2,type=im,algo=syn,metric=CC,slicewise=1,smooth=0,iter=3
                #sct_register_to_template -i fmri_spine_moco_mean.nii.gz -s fmri_spine_moco_mean_seg_corr.nii.gz -ldisc fmri_labels.nii.gz -c t2 -param step=1,type=seg,algo=centermass:step=2,type=im,algo=syn,metric=CC,slicewise=1,smooth=0,iter=3 #-qc "$DIREC$s"/func/"$d"
                sct_register_to_template -i fmri_spine_moco_mean.nii.gz -s fmri_spine_moco_mean_seg_corr.nii.gz -ldisc fmri_labels.nii.gz -ref subject -c t2 -param step=1,type=seg,algo=centermass:step=2,type=im,algo=syn,metric=CC,slicewise=1,smooth=0,iter=3


            # 3 - TSNR
            elif [ "$ind" == "6" ]; then
                tput setaf 2; echo "Compute TSNR and sSNR in " "$DIREC$s"/func/"$d"
                tput sgr0; 

                # compute tsnr around just spine
                sct_fmri_compute_tsnr -i fmri_spine_moco.nii.gz -o fmri_spine_moco_mean_tsnr_subject.nii.gz

                sct_apply_transfo -i fmri_spine_moco_mean_tsnr_subject.nii.gz -d ../../../template/PAM50_t2s.nii.gz -w warp_anat2template.nii.gz -o fmri_spine_moco_mean_tsnr_PAM50.nii.gz

                sct_crop_image -i fmri_spine_moco_mean_tsnr_PAM50.nii.gz -m ../../../template/PAM50_cervical_cord_all.nii.gz -b 0 -o fmri_spine_moco_mean_tsnr_PAM50.nii.gz

                #sct_compute_snr -i fmri_spine_moco.nii.gz -m fmri_spine_moco_mean_seg_corr.nii.gz \
                #    -method mult -o fmri_spine_moco_mean_ssnr_mult -v 1

                #sct_compute_snr -i fmri_spine_moco.nii.gz -m fmri_spine_moco_mean_seg_corr.nii.gz \
                #    -method diff -vol 81,82 -o fmri_spine_moco_mean_ssnr_diff -v 1

            else

				tput setaf 1; 
                echo "Index not valid (should be 1 to 6)"
				tput sgr0; 

			fi
		
        done

		tput setaf 6; echo "Done!" 
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