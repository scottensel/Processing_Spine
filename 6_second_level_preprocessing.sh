


#!/bin/bash
#
# SBSN framework - Preprocessing
# June 2022
#
# Preparation of functional spine data
#
# Requirements: FSL, Spinal Cord Toolbox 5.5
#
# SPINE

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
echo -n "Enter the index of the step to perform (0 = level 1 zfstat, 1 = Prepare for FLOB GLM (this one), 2 = Prepare for FORCE FLOB GLM (this one), 3 = Prepare for iCAP, 4 = Prepare 3rd level flip): "
tput sgr0;
read ind


# For each subject
for s in "${sub[@]}"; do

	cd $DIREC$s"/func/"

    for d in "${myFunc[@]}"; do

        MAX_JOBS=8
        job_count=0

        if [ "$ind" == "0" ]; then # this is the one

            tput setaf 2; echo "Moving files " $s"/func/func"$d
            tput sgr0; 

            cd $DIREC$s"/func/func"$d"/"
            if [ "$d" == "0" ]; then # this is the one 

                # do two different if statements for both the cope and var cope to avoid outlier cases of overwriting
                if [ -f "level_one_FLOB.feat/stats/subjectSpace_zstat1.nii.gz" ]; then

                    tput setaf 1; 
                    echo $DIREC$s"/func/func"$d"/level_one_FLOB.feat/stats/zstat1.nii.gz"
                    tput setaf 6;
                    # files have already been transofrmed
                    # subject space images are original images so just apply warps to them
                    # no need to rename files again
                    sct_apply_transfo -i "level_one_FLOB.feat/stats/subjectSpace_zstat1.nii.gz" -d ../../../template/PAM50_t2s.nii.gz -w warp_anat2template.nii.gz -o "level_one_FLOB.feat/stats/zstat1.nii.gz"

                else 
                    echo "No subject space cope"

                    # rename file then apply a transform to it so its located in PAM50 space
                    mv "level_one_FLOB.feat/stats/zstat1.nii.gz" "level_one_FLOB.feat/stats/subjectSpace_zstat1.nii.gz"

                    tput setaf 1; 
                    echo $DIREC$s"/func/func"$d"/level_one_FLOB.feat/stats/zstat1.nii.gz"

                    tput setaf 6;
                    # files have already been transofrmed
                    # subject space images are original images so just apply warps to them
                    # no need to rename files again
                    sct_apply_transfo -i "level_one_FLOB.feat/stats/subjectSpace_zstat1.nii.gz" -d ../../../template/PAM50_t2s.nii.gz -w warp_anat2template.nii.gz -o "level_one_FLOB.feat/stats/zstat1.nii.gz"

                fi

            else 

                # do two different if statements for both the cope and var cope to avoid outlier cases of overwriting
                if [ -f "level_one_force_FLOB.feat/stats/subjectSpace_zfstat1.nii.gz" ]; then

                    tput setaf 1; 
                    echo $DIREC$s"/func/func"$d"/level_one_force_FLOB.feat/stats/zfstat1.nii.gz"
                    tput setaf 6;
                    # files have already been transofrmed
                    # subject space images are original images so just apply warps to them
                    # no need to rename files again
                    sct_apply_transfo -i "level_one_force_FLOB.feat/stats/subjectSpace_zfstat1.nii.gz" -d ../../../template/PAM50_t2s.nii.gz -w warp_anat2template.nii.gz -o "level_one_force_FLOB.feat/stats/zfstat1.nii.gz"

                else 
                    echo "No subject space cope"

                    # rename file then apply a transform to it so its located in PAM50 space
                    mv "level_one_force_FLOB.feat/stats/zfstat1.nii.gz" "level_one_force_FLOB.feat/stats/subjectSpace_zfstat1.nii.gz"

                    tput setaf 1; 
                    echo $DIREC$s"/func/func"$d"/level_one_force_FLOB.feat/stats/zfstat1.nii.gz"

                    tput setaf 6;
                    # files have already been transofrmed
                    # subject space images are original images so just apply warps to them
                    # no need to rename files again
                    sct_apply_transfo -i "level_one_force_FLOB.feat/stats/subjectSpace_zfstat1.nii.gz" -d ../../../template/PAM50_t2s.nii.gz -w warp_anat2template.nii.gz -o "level_one_force_FLOB.feat/stats/zfstat1.nii.gz"

                fi

            fi

        elif [ "$ind" == "1" ]; then # this is the one

            tput setaf 2; echo "Prepare second level analysis for GLM " $s"/func/func"$d
            tput sgr0; 

            cd $DIREC$s"/func/func"$d"/"

            tput setaf 2; echo "Moving files " $s"/func/func"$d	

            cp -r "$DIREC"reg/"" $DIREC$s"/func/func"$d"/level_one_FLOB.feat/"
            
            # these have to be the same size across all subject because they get concatenated in the 4th dimension
            # change this to the PAM50 template
            cp ../../../template/PAM50_t2s.nii.gz level_one_FLOB.feat/reg/standard.nii.gz
            cp ../../../template/PAM50_t2s.nii.gz level_one_FLOB.feat/reg/example_func.nii.gz

            cp anat2template.nii.gz level_one_FLOB.feat/example_func.nii.gz
            cp anat2template.nii.gz level_one_FLOB.feat/mean_func.nii.gz
            
            mkdir level_one_FLOB.feat/reg_standard
            mkdir level_one_FLOB.feat/reg_standard/reg
            mkdir level_one_FLOB.feat/reg_standard/stats

            totalCopes=(1 2 3)
            for copeNum in "${totalCopes[@]}"; do
                
                # do two different if statements for both the cope and var cope to avoid outlier cases of overwriting
                if [ -f "level_one_FLOB.feat/stats/subjectSpace_cope"$copeNum".nii.gz" ]; then

                    tput setaf 1; 
                    echo $DIREC$s"/func/func"$d"/level_one_FLOB.feat/stats/cope"$copeNum

                    tput setaf 6;
                    # files have already been transofrmed
                    # subject space images are original images so just apply warps to them
                    # no need to rename files again
                    sct_apply_transfo -i "level_one_FLOB.feat/stats/subjectSpace_cope"$copeNum".nii.gz" -d ../../../template/PAM50_t2s.nii.gz -w warp_anat2template.nii.gz -o "level_one_FLOB.feat/stats/cope"$copeNum".nii.gz"

                    tput setaf 1;  
                    echo "Error checking cope"$copeNum
                    tput setaf 6;

                    # this is where second level always fails
                    # so I run the command myslef to error check and reapply the transoformation if this fails
                    # it will create a new file if passes and nothing if fails
                    flirt -ref level_one_FLOB.feat/reg/standard -in "level_one_FLOB.feat/stats/cope"$copeNum".nii.gz" -out "level_one_FLOB.feat/reg_standard/stats/cope"$copeNum".nii.gz" -applyxfm -init level_one_FLOB.feat/reg/example_func2standard.mat -interp trilinear -datatype float


                else 
                    echo "No subject space cope"

                    # rename file then apply a transform to it so its located in PAM50 space
                    mv "level_one_FLOB.feat/stats/cope"$copeNum".nii.gz" "level_one_FLOB.feat/stats/subjectSpace_cope"$copeNum".nii.gz"

                    tput setaf 1; 
                    echo $DIREC$s"/func/func"$d"/level_one_FLOB.feat/stats/cope"$copeNum

                    tput setaf 6;
                    # files have already been transofrmed
                    # subject space images are original images so just apply warps to them
                    # no need to rename files again
                    sct_apply_transfo -i "level_one_FLOB.feat/stats/subjectSpace_cope"$copeNum".nii.gz" -d ../../../template/PAM50_t2s.nii.gz -w warp_anat2template.nii.gz -o "level_one_FLOB.feat/stats/cope"$copeNum".nii.gz"

                    tput setaf 1;  
                    echo "Error checking cope"$copeNum
                    tput setaf 6;
                    
                    # this is where second level always fails
                    # so I run the command myslef to error check and reapply the transoformation if this fails
                    # it will create a new file if passes and nothing if fails
                    flirt -ref level_one_FLOB.feat/reg/standard -in "level_one_FLOB.feat/stats/cope"$copeNum".nii.gz" -out "level_one_FLOB.feat/reg_standard/stats/cope"$copeNum".nii.gz" -applyxfm -init level_one_FLOB.feat/reg/example_func2standard.mat -interp trilinear -datatype float

                fi

                # another if statement that does same as above except for varcope file
                if [ -f "level_one_FLOB.feat/stats/subjectSpace_varcope"$copeNum".nii.gz" ]; then

                    tput setaf 1;  
                    echo $DIREC$s"/func/func"$d"/level_one_FLOB.feat/stats/varcope"$copeNum

                    tput setaf 6;
                    # files have already been transofrmed
                    # subject space images are original images so just apply warps to them
                    # no need to rename files again
                    sct_apply_transfo -i "level_one_FLOB.feat/stats/subjectSpace_varcope"$copeNum".nii.gz" -d ../../../template/PAM50_t2s.nii.gz -w warp_anat2template.nii.gz -o "level_one_FLOB.feat/stats/varcope"$copeNum".nii.gz"

                    tput setaf 1;  
                    echo "Error checking varcope"$copeNum
                    tput setaf 6;   
                       
                    # this is where second level always fails
                    # so I run the command myslef to error check and reapply the transoformation if this fails
                    # it will create a new file if passes and nothing if fails
                    flirt -ref level_one_FLOB.feat/reg/standard -in "level_one_FLOB.feat/stats/varcope"$copeNum".nii.gz" -out "level_one_FLOB.feat/reg_standard/stats/varcope"$copeNum".nii.gz" -applyxfm -init level_one_FLOB.feat/reg/example_func2standard.mat -interp trilinear -datatype float

                else

                    echo "No subject space varcope"

                    mv "level_one_FLOB.feat/stats/varcope"$copeNum".nii.gz" "level_one_FLOB.feat/stats/subjectSpace_varcope"$copeNum".nii.gz"

                    tput setaf 1;  
                    echo $DIREC$s"/func/func"$d"/level_one_FLOB.feat/stats/varcope"$copeNum

                    tput setaf 6;
                    # files have already been transofrmed
                    # subject space images are original images so just apply warps to them
                    # no need to rename files again
                    sct_apply_transfo -i "level_one_FLOB.feat/stats/subjectSpace_varcope"$copeNum".nii.gz" -d ../../../template/PAM50_t2s.nii.gz -w warp_anat2template.nii.gz -o "level_one_FLOB.feat/stats/varcope"$copeNum".nii.gz"

                    tput setaf 1;  
                    echo "Error checking varcope"$copeNum
                    tput setaf 6;   
             
                    # this is where second level always fails
                    # so I run the command myslef to error check and reapply the transoformation if this fails
                    # it will create a new file if passes and nothing if fails
                    flirt -ref level_one_FLOB.feat/reg/standard -in "level_one_FLOB.feat/stats/varcope"$copeNum".nii.gz" -out "level_one_FLOB.feat/reg_standard/stats/varcope"$copeNum".nii.gz" -applyxfm -init level_one_FLOB.feat/reg/example_func2standard.mat -interp trilinear -datatype float

                fi

            done

            if [ -f level_one_FLOB.feat/stats/subjectSpace_mask.nii.gz ]; then

                cp ../../../template/PAM50_cervical_cord_all.nii.gz level_one_FLOB.feat/mask.nii.gz
                # why we dont use this transform below and use the one above instead
                #sct_apply_transfo -i level_one_force_FLOB.feat/subjectSpace_mask.nii.gz -d ../../../template/PAM50_t2s.nii.gz -w warp_anat2template.nii.gz -o level_one_force_FLOB.feat/mask.nii.gz

            else

                # i think this mask may have to be replaced with the PAM50 mask that we want to use
                mv level_one_FLOB.feat/mask.nii.gz level_one_FLOB.feat/subjectSpace_mask.nii.gz
                cp ../../../template/PAM50_cervical_cord_all.nii.gz level_one_FLOB.feat/mask.nii.gz

            fi

        elif [ "$ind" == "2" ]; then

            tput setaf 2; echo "Prepare second level analysis for GLM " $s"/func/func"$d
            tput sgr0; 

            cd $DIREC$s"/func/func"$d"/"

            tput setaf 2; echo "Moving files " $s"/func/func"$d 

            cp -r "$DIREC"reg/"" $DIREC$s"/func/func"$d"/level_one_force_FLOB.feat/"
            
            # these have to be the same size across all subject because they get concatenated in the 4th dimension
            # change this to the PAM50 template
            cp ../../../template/PAM50_t2s.nii.gz level_one_force_FLOB.feat/reg/standard.nii.gz
            cp ../../../template/PAM50_t2s.nii.gz level_one_force_FLOB.feat/reg/example_func.nii.gz

            cp anat2template.nii.gz level_one_force_FLOB.feat/example_func.nii.gz
            cp anat2template.nii.gz level_one_force_FLOB.feat/mean_func.nii.gz
            
            mkdir level_one_force_FLOB.feat/reg_standard
            mkdir level_one_force_FLOB.feat/reg_standard/reg
            mkdir level_one_force_FLOB.feat/reg_standard/stats

            totalCopes=(1 2 3 4 5 6 7 8 9)
            for copeNum in "${totalCopes[@]}"; do
             

                # this is a subshell process to make this run faster in parrallel
                (


                # do two different if statements for both the cope and var cope to avoid outlier cases of overwriting
                if [ -f "level_one_force_FLOB.feat/stats/subjectSpace_cope"$copeNum".nii.gz" ]; then

                    tput setaf 1; 
                    echo $DIREC$s"/func/func"$d"/level_one_force_FLOB.feat/stats/cope"$copeNum

                    tput setaf 6;
                    # files have already been transofrmed
                    # subject space images are original images so just apply warps to them
                    # no need to rename files again
                    sct_apply_transfo -i "level_one_force_FLOB.feat/stats/subjectSpace_cope"$copeNum".nii.gz" -d ../../../template/PAM50_t2s.nii.gz -w warp_anat2template.nii.gz -o "level_one_force_FLOB.feat/stats/cope"$copeNum".nii.gz"

                    tput setaf 1;  
                    echo "Error checking cope"$copeNum
                    tput setaf 6;

                    # this is where second level always fails
                    # so I run the command myslef to error check and reapply the transoformation if this fails
                    # it will create a new file if passes and nothing if fails
                    flirt -ref level_one_force_FLOB.feat/reg/standard -in "level_one_force_FLOB.feat/stats/cope"$copeNum".nii.gz" -out "level_one_force_FLOB.feat/reg_standard/stats/cope"$copeNum".nii.gz" -applyxfm -init level_one_force_FLOB.feat/reg/example_func2standard.mat -interp trilinear -datatype float


                else 
                    echo "No subject space cope"

                    # rename file then apply a transform to it so its located in PAM50 space
                    mv "level_one_force_FLOB.feat/stats/cope"$copeNum".nii.gz" "level_one_force_FLOB.feat/stats/subjectSpace_cope"$copeNum".nii.gz"

                    tput setaf 1; 
                    echo $DIREC$s"/func/func"$d"/level_one_force_FLOB.feat/stats/cope"$copeNum

                    tput setaf 6;
                    # files have already been transofrmed
                    # subject space images are original images so just apply warps to them
                    # no need to rename files again
                    sct_apply_transfo -i "level_one_force_FLOB.feat/stats/subjectSpace_cope"$copeNum".nii.gz" -d ../../../template/PAM50_t2s.nii.gz -w warp_anat2template.nii.gz -o "level_one_force_FLOB.feat/stats/cope"$copeNum".nii.gz"

                    tput setaf 1;  
                    echo "Error checking cope"$copeNum
                    tput setaf 6;
                    
                    # this is where second level always fails
                    # so I run the command myslef to error check and reapply the transoformation if this fails
                    # it will create a new file if passes and nothing if fails
                    flirt -ref level_one_force_FLOB.feat/reg/standard -in "level_one_force_FLOB.feat/stats/cope"$copeNum".nii.gz" -out "level_one_force_FLOB.feat/reg_standard/stats/cope"$copeNum".nii.gz" -applyxfm -init level_one_force_FLOB.feat/reg/example_func2standard.mat -interp trilinear -datatype float

                fi

                # another if statement that does same as above except for varcope file
                if [ -f "level_one_force_FLOB.feat/stats/subjectSpace_varcope"$copeNum".nii.gz" ]; then

                    tput setaf 1;  
                    echo $DIREC$s"/func/func"$d"/level_one_force_FLOB.feat/stats/varcope"$copeNum

                    tput setaf 6;
                    # files have already been transofrmed
                    # subject space images are original images so just apply warps to them
                    # no need to rename files again
                    sct_apply_transfo -i "level_one_force_FLOB.feat/stats/subjectSpace_varcope"$copeNum".nii.gz" -d ../../../template/PAM50_t2s.nii.gz -w warp_anat2template.nii.gz -o "level_one_force_FLOB.feat/stats/varcope"$copeNum".nii.gz"

                    tput setaf 1;  
                    echo "Error checking varcope"$copeNum
                    tput setaf 6;   
                       
                    # this is where second level always fails
                    # so I run the command myslef to error check and reapply the transoformation if this fails
                    # it will create a new file if passes and nothing if fails
                    flirt -ref level_one_force_FLOB.feat/reg/standard -in "level_one_force_FLOB.feat/stats/varcope"$copeNum".nii.gz" -out "level_one_force_FLOB.feat/reg_standard/stats/varcope"$copeNum".nii.gz" -applyxfm -init level_one_force_FLOB.feat/reg/example_func2standard.mat -interp trilinear -datatype float

                else

                    echo "No subject space varcope"

                    mv "level_one_force_FLOB.feat/stats/varcope"$copeNum".nii.gz" "level_one_force_FLOB.feat/stats/subjectSpace_varcope"$copeNum".nii.gz"

                    tput setaf 1;  
                    echo $DIREC$s"/func/func"$d"/level_one_force_FLOB.feat/stats/varcope"$copeNum

                    tput setaf 6;
                    # files have already been transofrmed
                    # subject space images are original images so just apply warps to them
                    # no need to rename files again
                    sct_apply_transfo -i "level_one_force_FLOB.feat/stats/subjectSpace_varcope"$copeNum".nii.gz" -d ../../../template/PAM50_t2s.nii.gz -w warp_anat2template.nii.gz -o "level_one_force_FLOB.feat/stats/varcope"$copeNum".nii.gz"

                    tput setaf 1;  
                    echo "Error checking varcope"$copeNum
                    tput setaf 6;   
             
                    # this is where second level always fails
                    # so I run the command myslef to error check and reapply the transoformation if this fails
                    # it will create a new file if passes and nothing if fails
                    flirt -ref level_one_force_FLOB.feat/reg/standard -in "level_one_force_FLOB.feat/stats/varcope"$copeNum".nii.gz" -out "level_one_force_FLOB.feat/reg_standard/stats/varcope"$copeNum".nii.gz" -applyxfm -init level_one_force_FLOB.feat/reg/example_func2standard.mat -interp trilinear -datatype float

                fi

                ) &
                # this is the end of the subshell

                # this counts the jobs and lets multiple run at the same time
                ((job_count++))
                if ((job_count >= MAX_JOBS)); then
                    wait -n     # Wait for *one* job to finish before starting a new one
                    ((job_count--))
                fi

            done

            wait

            if [ -f level_one_force_FLOB.feat/stats/subjectSpace_mask.nii.gz ]; then

                cp ../../../template/PAM50_cervical_cord_all.nii.gz level_one_force_FLOB.feat/mask.nii.gz
                # why we dont use this transform below and use the one above instead
                #sct_apply_transfo -i level_one_force_FLOB.feat/subjectSpace_mask.nii.gz -d ../../../template/PAM50_t2s.nii.gz -w warp_anat2template.nii.gz -o level_one_force_FLOB.feat/mask.nii.gz

            else

                # i think this mask may have to be replaced with the PAM50 mask that we want to use
                mv level_one_force_FLOB.feat/mask.nii.gz level_one_force_FLOB.feat/subjectSpace_mask.nii.gz
                cp ../../../template/PAM50_cervical_cord_all.nii.gz level_one_force_FLOB.feat/mask.nii.gz

            fi

        elif [ "$ind" == "3" ]; then

            mkdir $DIREC$s"/func/func"$d"/iCAP/TA" -p
            cd $DIREC$s"/func/func"$d"/iCAP"

            cp $DIREC"template/PAM50_cervical_cord_all.nii.gz" PAM50_cervical_cord_all.nii.gz
            cp $DIREC"template/PAM50_cervical_cord.nii.gz" PAM50_cervical_cord.nii.gz
            cp ../fmri_spine_moco_denoised_smooth.nii.gz TA/fmri_spine_moco_denoised_smooth.nii.gz

            gunzip PAM50_cervical_cord_all.nii.gz -f
            gunzip PAM50_cervical_cord.nii.gz -f

            cd $DIREC$s"/func/func"$d"/iCAP/TA"

            # transform to space
            tput setaf 2; echo "Transform to PAM50"
                    tput sgr0;

            # apply transformation
            sct_apply_transfo -i fmri_spine_moco_denoised_smooth.nii.gz -d ../../../../../template/PAM50_t2s.nii.gz -w ../../warp_anat2template.nii.gz -o fmri_spine_moco_denoised_smooth.nii.gz

            # split the data
            tput setaf 2; echo "...Split functional data"
                    tput sgr0;

            # First, split data
            fslsplit fmri_spine_moco_denoised_smooth.nii.gz resvol -t

            for n in "$PWD"/resvol*.nii.gz; do # Loop through all files

                IFS='.' read -r volname string <<< "$n"
                gunzip "${volname##*/}".nii.gz -f

            done

            tput setaf 2; echo "...Cleaning up"
                    tput sgr0;

            rm straight_ref.nii.gz
            rm warp_curve2straight.nii.gz
            rm warp_straight2curve.nii.gz
            rm straightening.cache
            rm resvol*.nii.gz
            mv fmri_spine_moco_denoised_smooth.nii.gz ../fmri_spine_moco_denoised_smooth.nii.gz


        elif [ "$ind" == "4" ]; then

            tput setaf 2; echo "Prepare third level analysis for GLM " $s
            tput sgr0; 

            cd $DIREC$s"/func/"

            # do two different if statements for both the cope and var cope to avoid outlier cases of overwriting
            if [ -f "level_two_all_force_FLOB1234.gfeat/cope1.feat/stats/subjectSpace_cope1.nii.gz" ]; then

                tput setaf 1; 
                echo $DIREC$s"/func/"

                tput setaf 6;
                # files have already been transofrmed
                # subject space images are original images so just apply warps to them
                # no need to rename files again                
                fslswapdim "level_two_all_force_FLOB1234.gfeat/cope1.feat/stats/subjectSpace_cope1.nii.gz" -x y z "level_two_all_force_FLOB1234.gfeat/cope1.feat/stats/cope1.nii.gz"
            
            else 
                echo "No subject space cope"

                # rename file then apply a transform to it so its located in PAM50 space
                mv "level_two_all_force_FLOB1234.gfeat/cope1.feat/stats/cope1.nii.gz" "level_two_all_force_FLOB1234.gfeat/cope1.feat/stats/subjectSpace_cope1.nii.gz"

                tput setaf 1; 
                echo $DIREC$s"/func/"

                tput setaf 6;
                # files have already been transofrmed
                fslswapdim "level_two_all_force_FLOB1234.gfeat/cope1.feat/stats/subjectSpace_cope1.nii.gz" -x y z "level_two_all_force_FLOB1234.gfeat/cope1.feat/stats/cope1.nii.gz"

            fi

            # another if statement that does same as above except for varcope file
            if [ -f "level_two_all_force_FLOB1234.gfeat/cope1.feat/stats/subjectSpace_varcope1.nii.gz" ]; then

                tput setaf 1; 
                echo $DIREC$s"/func/"

                tput setaf 6;
                # files have already been transofrmed
                # subject space images are original images so just apply warps to them
                # no need to rename files again                
                fslswapdim "level_two_all_force_FLOB1234.gfeat/cope1.feat/stats/subjectSpace_varcope1.nii.gz" -x y z "level_two_all_force_FLOB1234.gfeat/cope1.feat/stats/varcope1.nii.gz"
            

            else
                echo "No subject space varcope"

                # rename file then apply a transform to it so its located in PAM50 space
                mv "level_two_all_force_FLOB1234.gfeat/cope1.feat/stats/varcope1.nii.gz" "level_two_all_force_FLOB1234.gfeat/cope1.feat/stats/subjectSpace_varcope1.nii.gz"

                tput setaf 1; 
                echo $DIREC$s"/func/"

                tput setaf 6;
                # files have already been transofrmed
                fslswapdim "level_two_all_force_FLOB1234.gfeat/cope1.feat/stats/subjectSpace_varcope1.nii.gz" -x y z "level_two_all_force_FLOB1234.gfeat/cope1.feat/stats/varcope1.nii.gz"


            fi
            
            # another if statement that does same as above except for varcope file
            if [ -f "level_two_all_force_FLOB1234.gfeat/cope1.feat/stats/subjectSpace_tdof_t1.nii.gz" ]; then

                tput setaf 1; 
                echo $DIREC$s"/func/"

                tput setaf 6;
                # files have already been transofrmed
                # subject space images are original images so just apply warps to them
                # no need to rename files again                
                fslswapdim "level_two_all_force_FLOB1234.gfeat/cope1.feat/stats/subjectSpace_tdof_t1.nii.gz" -x y z "level_two_all_force_FLOB1234.gfeat/cope1.feat/stats/tdof_t1.nii.gz"
            

            else
                echo "No subject space tdof_t"

                # rename file then apply a transform to it so its located in PAM50 space
                mv "level_two_all_force_FLOB1234.gfeat/cope1.feat/stats/tdof_t1.nii.gz" "level_two_all_force_FLOB1234.gfeat/cope1.feat/stats/subjectSpace_tdof_t1.nii.gz"

                tput setaf 1; 
                echo $DIREC$s"/func/"

                tput setaf 6;
                # files have already been transofrmed
                fslswapdim "level_two_all_force_FLOB1234.gfeat/cope1.feat/stats/subjectSpace_tdof_t1.nii.gz" -x y z "level_two_all_force_FLOB1234.gfeat/cope1.feat/stats/tdof_t1.nii.gz"


            fi
        fi
        
    done

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