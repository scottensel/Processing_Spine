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
. /mnt/d/SBSN/Processing_Spine/path_to_subjects.sh 

tput setaf 6; 
echo -n "Enter the index of the step to perform (1 = Run PPI Feat, 2 = Set up PPI to PAM50, 3 = Run Second Level PPI): "
tput sgr0;
read ind

naming=$(printf '%s' "${myFunc[@]}")

RADIUS=2

# bAreas=("L_SM" "R_SM" "R_CEB" "L_CEB")
bAreas=("Right_M1" "Left_M1" "Right_SMA" "Left_SMA" "Right_PMv" "Left_PMv")
bAreas=("Right_PMv" "Left_PMv")
# bAreas=("Left_M1" "Right_SMA")
# 
# bAreas=("Left_SMA" "Left_PMv")
# bAreas=("Left_PMv")
# bAreas=("Left_M1")
# bAreas=("Left_SMA")

#bAreas=("L_CEB" "R_CEB")
#bAreas=("R_CEB")
# For each subject


for s in "${sub[@]}"; do

    if [ "$ind" == "1" ]; then

        cd $DIREC$s"/func/"

        for d in "${myFunc[@]}"; do

            # Will print */ if no directories are available
            cd $DIREC$s"/func/func"$d"/"

            templateFile="template_design_ppi.fsf"
            ## Generate fsf file from template
            ## path is relevant to folder we are in which is sub-XX/fun

            # for ii in 0 1 2 3; do
            for ii in "${!bAreas[@]}"; do   

                # if files doesnt exist we cant run this
                if [ ! -f ${bAreas[$ii]}$RADIUS"_regressor.txt" ]; then
                    continue
                fi

                ## for i in "../../../template/template_design.fsf"; do
                for i in "../../../template/"$templateFile; do

                    tput setaf 2; echo "Prepare for by PPI on denoised data " $s"/func/func"$d"/"${bAreas[$ii]}$RADIUS"_regressor.txt"
                    tput sgr0;

                    sed -e 's@OUTDIR@'"level_one_PPI_"${bAreas[$ii]}$RADIUS'@g' \
                                        -e 's@DATAPATH@'$DIREC$s"/func/func"$d"/fmri_spine_moco_denoised_plusmean_smooth.nii.gz"'@g' \
                                        -e 's@OUTLYN@'"1"'@g' \
                                        -e 's@NPTS@'"$(fslnvols $DIREC$s"/func/func"$d"/fmri_spine_moco.nii.gz")"'@g' \
                                        -e 's@EV_TITLE1@'20'@g' \
                                        -e 's@EV_FILE1@'$DIREC$s"/task/task"$d"/force20.txt"'@g' \
                                        -e 's@EV_TITLE2@'45'@g' \
                                        -e 's@EV_FILE2@'$DIREC$s"/task/task"$d"/force45.txt"'@g' \
                                        -e 's@EV_TITLE3@'70'@g' \
                                        -e 's@EV_FILE3@'$DIREC$s"/task/task"$d"/force70.txt"'@g' \
                                        -e 's@EV_TITLE4@'${bAreas[$ii]}$RADIUS'@g' \
                                        -e 's@EV_FILE4@'$DIREC$s"/func/func"$d"/"${bAreas[$ii]}$RADIUS"_regressor.txt"'@g'<$i> design_levelone_PPI.fsf


                    # created a design_levelone.fsf based on the template.fsf
                done


                tput setaf 2; echo "Run first level PPI analysis for " $s"/func/func"$d
                tput sgr0; 
                
                # Run the analysis using the fsf file
                feat design_levelone_PPI.fsf

            done

        done            

    elif [ "$ind" == "2" ]; then

        cd $DIREC$s"/func/"

        for d in "${myFunc[@]}"; do

            MAX_JOBS=8
            job_count=0

            tput setaf 2; echo "Prepare second level analysis for GLM " $s"/func/func"$d
            tput sgr0; 

            cd $DIREC$s"/func/func"$d"/"

            # for ii in 0 1 2 3; do
            for ii in "${!bAreas[@]}"; do
                
                # this is a subshell process to make this run faster in parrallel
                (

                # if files doesnt exist we cant run this
                if [ ! -d "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat" ]; then                        
                    continue
                fi                

                tput setaf 2; echo "Moving files " $s"/func/func"$d 

                cp -r $DIREC"reg/" $DIREC$s"/func/func"$d"/level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/"
                
                # these have to be the same size across all subject because they get concatenated in the 4th dimension
                # change this to the PAM50 template
                cp ../../../template/PAM50_t2s.nii.gz "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/reg/standard.nii.gz"
                cp ../../../template/PAM50_t2s.nii.gz "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/reg/example_func.nii.gz"

                cp anat2template.nii.gz "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/example_func.nii.gz"
                cp anat2template.nii.gz "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/mean_func.nii.gz"
                
                mkdir "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/reg_standard"
                mkdir "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/reg_standard/reg"
                mkdir "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/reg_standard/stats"

                totalCopes=(11)
                for copeNum in "${totalCopes[@]}"; do
                    
                    # do two different if statements for both the cope and var cope to avoid outlier cases of overwriting
                    if [ -f "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/subjectSpace_cope"$copeNum".nii.gz" ]; then

                        tput setaf 1; 
                        echo $DIREC$s"/func/func"$d"/level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/cope"$copeNum

                        tput setaf 6;
                        # files have already been transofrmed
                        # subject space images are original images so just apply warps to them
                        # no need to rename files again
                        sct_apply_transfo -i "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/subjectSpace_cope"$copeNum".nii.gz" -d ../../../template/PAM50_t2s.nii.gz -w warp_anat2template.nii.gz -o "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/cope"$copeNum".nii.gz"

                        tput setaf 1;  
                        echo "Error checking cope"$copeNum
                        tput setaf 6;

                        # this is where second level always fails
                        # so I run the command myslef to error check and reapply the transoformation if this fails
                        # it will create a new file if passes and nothing if fails
                        flirt -ref "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/reg/standard" -in "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/cope"$copeNum".nii.gz" -out "level_one_PPI_"${bAreas[$ii]}$RADIUS"B.feat/reg_standard/stats/cope"$copeNum".nii.gz" -applyxfm -init "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/reg/example_func2standard.mat" -interp trilinear -datatype float


                    else 
                        echo "No subject space cope"

                        # rename file then apply a transform to it so its located in PAM50 space
                        mv "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/cope"$copeNum".nii.gz" "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/subjectSpace_cope"$copeNum".nii.gz"

                        tput setaf 1; 
                        echo $DIREC$s"/func/func"$d"/level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/cope"$copeNum

                        tput setaf 6;
                        # files have already been transofrmed
                        # subject space images are original images so just apply warps to them
                        # no need to rename files again
                        sct_apply_transfo -i "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/subjectSpace_cope"$copeNum".nii.gz" -d ../../../template/PAM50_t2s.nii.gz -w warp_anat2template.nii.gz -o "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/cope"$copeNum".nii.gz"

                        tput setaf 1;  
                        echo "Error checking cope"$copeNum
                        tput setaf 6;
                        
                        # this is where second level always fails
                        # so I run the command myslef to error check and reapply the transoformation if this fails
                        # it will create a new file if passes and nothing if fails
                        flirt -ref "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/reg/standard" -in "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/cope"$copeNum".nii.gz" -out "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/reg_standard/stats/cope"$copeNum".nii.gz" -applyxfm -init "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/reg/example_func2standard.mat" -interp trilinear -datatype float

                    fi

                    # another if statement that does same as above except for varcope file
                    if [ -f "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/subjectSpace_varcope"$copeNum".nii.gz" ]; then

                        tput setaf 1;  
                        echo $DIREC$s"/func/func"$d"/level_one_PPI_"${bAreas[$ii]}$RADIUS"B.feat/stats/varcope"$copeNum

                        tput setaf 6;
                        # files have already been transofrmed
                        # subject space images are original images so just apply warps to them
                        # no need to rename files again
                        sct_apply_transfo -i "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/subjectSpace_varcope"$copeNum".nii.gz" -d ../../../template/PAM50_t2s.nii.gz -w warp_anat2template.nii.gz -o "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/varcope"$copeNum".nii.gz"

                        tput setaf 1;  
                        echo "Error checking varcope"$copeNum
                        tput setaf 6;   
                           
                        # this is where second level always fails
                        # so I run the command myslef to error check and reapply the transoformation if this fails
                        # it will create a new file if passes and nothing if fails
                        flirt -ref "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/reg/standard" -in "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/varcope"$copeNum".nii.gz" -out "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/reg_standard/stats/varcope"$copeNum".nii.gz" -applyxfm -init "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/reg/example_func2standard.mat" -interp trilinear -datatype float

                    else

                        echo "No subject space varcope"

                        mv "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/varcope"$copeNum".nii.gz" "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/subjectSpace_varcope"$copeNum".nii.gz"

                        tput setaf 1;  
                        echo $DIREC$s"/func/func"$d"/level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/varcope"$copeNum

                        tput setaf 6;
                        # files have already been transofrmed
                        # subject space images are original images so just apply warps to them
                        # no need to rename files again
                        sct_apply_transfo -i "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/subjectSpace_varcope"$copeNum".nii.gz" -d ../../../template/PAM50_t2s.nii.gz -w warp_anat2template.nii.gz -o "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/varcope"$copeNum".nii.gz"

                        tput setaf 1;  
                        echo "Error checking varcope"$copeNum
                        tput setaf 6;   
                 
                        # this is where second level always fails
                        # so I run the command myslef to error check and reapply the transoformation if this fails
                        # it will create a new file if passes and nothing if fails
                        flirt -ref "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/reg/standard" -in "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/varcope"$copeNum".nii.gz" -out "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/reg_standard/stats/varcope"$copeNum".nii.gz" -applyxfm -init "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/reg/example_func2standard.mat" -interp trilinear -datatype float

                    fi

                done

                if [ -f "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/stats/subjectSpace_mask.nii.gz" ]; then

                    cp ../../../template/PAM50_cervical_cord_all.nii.gz "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/mask.nii.gz"
                    # why we dont use this transform below and use the one above instead

                else

                    # i think this mask may have to be replaced with the PAM50 mask that we want to use
                    mv "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/mask.nii.gz" "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/subjectSpace_mask.nii.gz"
                    cp ../../../template/PAM50_cervical_cord_all.nii.gz "level_one_PPI_"${bAreas[$ii]}$RADIUS".feat/mask.nii.gz"

                fi

                ) &

                ((job_count++))
                if ((job_count >= MAX_JOBS)); then
                    wait -n     # Wait for *one* job to finish before starting a new one
                    ((job_count--))
                fi

            done

            wait # DONT FORGET THIS

        done

    elif [ "$ind" == "3" ]; then

        cd $DIREC$s"/func/"
        
        # this is a counter varibale used to index each run we specificed
        # its so we can select any combination of runs to use

        tput setaf 2; echo "Prepare template files for analysis..."
        tput sgr0;

        # for ii in 0 1 2 3; do
        for ii in "${!bAreas[@]}"; do

            j=1

            # if files doesnt exist we cant run this
            if [ ! -d "func1/level_one_PPI_"${bAreas[$ii]}$RADIUS".feat" ]; then                        
                continue
            fi       

            # changes the template file that it loops through
            levelFile="level_one_PPI_"${bAreas[$ii]}$RADIUS".feat"
            outputName="level_two_PPI_"${bAreas[$ii]}$RADIUS
            templateName="second_level_PPI_template.fsf"
            threshmask='PAM50_cervical_cord_all.nii.gz'

            # here I run through the number of runs were are to combine specified when you start running this
            # always must be a minimum of two runs to do a second level analysis of per subject
            for runNum in ${myFunc[@]}; do
                
                # starts with the first run specified
                ## Generate fsf file from template
                if [ "$j" == "1" ]; then
                    
                    # here selecting the original named template file
                    for i in "../../template/"$templateName; do
                        
                        # now we set some standard paths ands masks 
                        # set path 1 as the first run
                        # and resave the file as a new template inside the folder
                        sed -e 's@OUTDIR@'"$outputName"'@g' \
                            -e 's@PATH1@'$DIREC$s"/func/func"$runNum"/"$levelFile""'@g' \
                            -e 's@THRESH_MASK@'$DIREC"template/"$threshmask""'@g' \
                            -e 's@NSUBJECTS@'${#myFunc[@]}'@g' <$i> design_leveltwo_force_FLOB"$j".fsf 

                    done
                
                # must need a second run so this is when j=2
                elif [ "$j" == "2" ]; then
                    
                    # take the template file we just named in the previous statement
                    for i in "design_leveltwo_force_FLOB"$((j-1))".fsf"; do
                        
                        # now we sepcificy the second path we chose and create a new template file
                        sed -e 's@PATH2@'$DIREC$s"/func/func"$runNum"/"$levelFile""'@g' \
                            -e 's@OUTLPATH@''@g' <$i> design_leveltwo_force_FLOB"$j".fsf 

                    done
                    
                    # remove the old template file as its not needed
                    #echo design_leveltwo_force_FLOB"$((j-1))"
                    rm design_leveltwo_force_FLOB"$((j-1))".fsf

                else
                    
                    # take the template file we just named in the previous statement                
                    for i in "design_leveltwo_force_FLOB"$((j-1))".fsf"; do
     
                        # here we have to edit the file to add in more than two runs
                        # this is so we can choose whatever amount of runs we want and add into the file instead of delete
                        sed -e 's@FEAT_PATH@'"# 4D AVW data or FEAT directory ("$j")\nset feat_files("$j") "$DIREC$s"/func/func"$runNum"/"$levelFile"\n\nFEAT_PATH"'@g' \
                            -e 's@EVG_PATH@'"# Higher-level EV value for EV 1 and input "$j"\nset fmri(evg"$j".1) 1\n\nEVG_PATH"'@g' \
                            -e 's@GROUPMEM_PATH@'"# Group membership for input "$j"\nset fmri(groupmem\."$j") 1\n\nGROUPMEM_PATH"'@g' <$i> design_leveltwo_force_FLOB"$j".fsf 

                    done

                    # remove the old template file as its not needed
                    rm design_leveltwo_force_FLOB"$((j-1))".fsf

                fi           
     
                ((j+=1));
                
                # one statement to make sure the file is named properly
                if [ "$j" -gt "${#myFunc[@]}" ]; then

                    # checking if a file exists already
                    if [ -f design_leveltwo_force_FLOB"$naming".fsf ]; then

                        # remove the old previous one because im not sure it will overwrite properly
                        rm design_leveltwo_force_FLOB"$naming".fsf
                    fi
                    
                    # set some final paths to empty
                    for i in "design_leveltwo_force_FLOB"$((j-1))".fsf"; do

                        sed -e 's@FEAT_PATH@''@g' -e 's@EVG_PATH@''@g' -e 's@GROUPMEM_PATH@''@g' <$i> design_leveltwo_force_FLOB"$naming".fsf 

                    done

                    #echo design_leveltwo"$((j-1))"
                    rm design_leveltwo_force_FLOB"$((j-1))".fsf

                fi

            done


            tput setaf 2; echo "Run second level PPI analysis on"
            tput sgr0; 

            tput setaf 2;
            echo $s
            echo ${bAreas[$ii]}$RADIUS
            tput sgr0; 

            # Run the analysis using the fsf file
            feat design_leveltwo_force_FLOB"$naming".fsf

    #######################################################################
            # move into folder so it doesnt take up space
            mv design_leveltwo_force_FLOB"$naming".fsf  "$outputName".gfeat/design_leveltwo_force_FLOB"$naming".fsf 

            tput setaf 2; echo "Done!" 
            tput sgr0;  

        done
    fi
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