 #!/bin/bash

#
# SPiCiCAP framework - Preprocessing
# July 2020
#
# Smoothing along the cord with a 2x2x6mm FWHM Gaussian kernel
#
# Requirements: Spinal Cord Toolbox 3.2.7, FSL
#
#

# List of subjects
#declare -a sub=("sub-01" "sub-02" "sub-03" "sub-04" "sub-05" "sub-06" "sub-07" "sub-08" "sub-09" "sub-10" "sub-11" "sub-12" "sub-13" "sub-14" "sub-15" "sub-16" "sub-17" "sub-18" "sub-19")
declare -a sub=("sub-01")

# Path of the folder containing all data
#DIREC="/home/sce30/Documents/Research/NeuronData/SPiCiCAP_Dataset_Code/Data/"
DIREC="/home/sce30/Documents/Research/Data/Spine/"

# For each subject
for s in "${sub[@]}"; do
	cd $DIREC$s"/func/"

		tput setaf 2; echo "Smoothing started in " $s
                tput sgr0;
                
		cd Processing
		
		mkdir -p Smoothing
		cd Smoothing

		# First, isolate slice init (get the first slice)
                if [ ! -f slice_init.nii.gz ]; then
                    fslroi ../../mfmri_denoised slice_init 0 -1 0 -1 0 1 0 -1 # extract a region of interest <input> <output> <xmin> <xsize> <ymin> <ysize> <zmin> <zsize> <tmin> <tsize> . . . -1 for size sets to full image extent
                    # if size is 1 then it pulls only that frame. EX. 10 3 would give slices 10 11 12
                fi
	
                # getting last slice index
                dim=$(fslsize ../../mfmri_denoised.nii.gz)
                arrdim=($dim)
                z=${arrdim[5]}
                maxz=$(( z-1 ))
	
                # Then, isolate slice last
                if [ ! -f slice_last.nii.gz ]; then
                        fslroi ../../mfmri_denoised slice_last 0 -1 0 -1 $maxz 1
                fi

        
		# getting first and last slice for mask as well
                if [ ! -f ../../Segmentation/mask_init.nii.gz ]; then
                        fslroi ../../Segmentation/mask_sc.nii.gz ../../Segmentation/mask_init 0 -1 0 -1 0 1 0 -1
                fi
		
                if [ ! -f ../../Segmentation/mask_last.nii.gz ]; then
                        fslroi ../../Segmentation/mask_sc.nii.gz ../../Segmentation/mask_last 0 -1 0 -1 $maxz 1
                fi

		#rm residuals_native_padded.nii.gz

        #########
        # now its creating an image of the denoised mfmri and padding it with slices in the top and bottom of z
        # doing the same thing for the mask
        #######


		if [ ! -f mfmri_denoised_padded.nii.gz ]; then
            # fslmerge <direction> <output> <file1 file2 . . . >
			fslmerge -z mfmri_denoised_padded ../../mfmri_denoised slice_last slice_last slice_last slice_last slice_last
			fslmerge -z mfmri_denoised_padded slice_init slice_init slice_init slice_init slice_init mfmri_denoised_padded
		fi

                if [ ! -f mask_padded.nii.gz ]; then
                        fslmerge -z mask_padded ../../Segmentation/mask_sc ../../Segmentation/mask_last ../../Segmentation/mask_last ../../Segmentation/mask_last ../../Segmentation/mask_last ../../Segmentation/mask_last
                        fslmerge -z mask_padded ../../Segmentation/mask_init ../../Segmentation/mask_init ../../Segmentation/mask_init ../../Segmentation/mask_init ../../Segmentation/mask_init mask_padded
                fi

                # Create temporary folder
                rm -r tmp
		mkdir -p tmp

                # Copy input data to tmp folder & split
                cp mfmri_denoised_padded.nii.gz tmp
                cd tmp
                # splitting the file along the time dimension The suffix _DIM+NUMBER will be added to the intput file name.
                sct_image -i mfmri_denoised_padded.nii.gz -split t -o mfmri_denoised_padded.nii.gz

        # remove the 4d image
		rm mfmri_denoised_padded.nii.gz

                files_to_smooth="$PWD"/*
                for b in $files_to_smooth; do # Loop through all files
                        sct_smooth_spinalcord -i $b -s ../mask_padded.nii.gz -smooth 0.849329,0.849329,2.5479870902 #smooth each file based on centerline in mask_padded
                        rm $b
                done
		
		# Concat, put back in initial folder and delete temporary folder
                fslmerge -tr s_mfmri_denoised_padded_2x2x6 mfmri_denoised_padded*_smooth.nii 2.5 # combine all the smoothed finals in a temporal resolution of 2.5 seconds
                mv s_mfmri_denoised_padded_2x2x6.nii.gz ..
                cd ..

		#sct_crop_image -i s_mfmri_denoised_padded_2x2x6.nii.gz -start 5 -end -5 -dim 2 -o s_mfmri_denoised_2x2x6.nii.gz
        # now we get rid of the padding we had in the start
        sct_crop_image -i s_mfmri_denoised_padded_2x2x6.nii.gz -zmin 5 -zmax -5 -o s_mfmri_denoised_2x2x6.nii.gz

		mv s_mfmri_denoised_2x2x6.nii.gz ../../s_mfmri_denoised_2x2x6.nii.gz
		
		cd ../..

                # Copy header info
                fslcpgeom mfmri_denoised.nii.gz s_mfmri_denoised_2x2x6.nii.gz

		tput setaf 2; echo "Done!"
                tput sgr0;

done



