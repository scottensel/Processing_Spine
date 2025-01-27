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
echo -n "Enter the index of the step to perform (1 = Prepare for GLM, 2 = Prepare for seperate force GLM, 3 = Prepare for FLOB force GLM, 4 = Prepare for FLOB smoothed force GLM): "
tput sgr0;
read ind

naming=$(printf '%s' "${myFunc[@]}")

#echo "${#myFunc[@]}"
#for d in {1..${#myFunc[@]}}; do
#${myFunc[d]}

## now I need to edit this for multiple copes instead of just one


# For each subject
for s in "${sub[@]}"; do

		cd $DIREC$s"/func/"
        
        # this is a counter varibale used to index each run we specificed
        # its so we can select any combination of runs to use
        j=1

        tput setaf 2; echo "Prepare template files for analysis..."
        tput sgr0;

	    if [ "$ind" == "1" ]; then
        
            # here I run through the number of runs were are to combine specified when you start running this
            # always must be a minimum of two runs to do a second level analysis of per subject
            for runNum in ${myFunc[@]}; do
                
                # starts with the first run specified
                ## Generate fsf file from template
                if [ "$j" == "1" ]; then
                    
                    # here selecting the original named template file
                    for i in "../../template/second_level_template.fsf"; do
                        
                        # now we set some standard paths ands masks 
                        # set path 1 as the first run
                        # and resave the file as a new template inside the folder
                        sed -e 's@OUTDIR@'"level_two"$naming""'@g' \
                            -e 's@PATH1@'$DIREC$s"/func/func"$runNum"/level_one.feat"'@g' \
                            -e 's@THRESH_MASK@'$DIREC"template/PAM50_cord.nii.gz"'@g' \
                            -e 's@NSUBJECTS@'${#myFunc[@]}'@g' <$i> design_leveltwo"$j".fsf 
    
                    done
                
                # must need a second run so this is when j=2
                elif [ "$j" == "2" ]; then
                    
                    # take the template file we just named in the previous statement
                    for i in "design_leveltwo"$((j-1))".fsf"; do
                        
                        # now we sepcificy the second path we chose and create a new template file
                        sed -e 's@PATH2@'$DIREC$s"/func/func"$runNum"/level_one.feat"'@g' \
                            -e 's@OUTLPATH@''@g' <$i> design_leveltwo"$j".fsf 
    
                    done
                    
                    # remove the old template file as its not needed
                    #echo design_leveltwo"$((j-1))"
                    rm design_leveltwo"$((j-1))".fsf
    
                else
                    
                    # take the template file we just named in the previous statement                
                    for i in "design_leveltwo"$((j-1))".fsf"; do
     
                        # here we have to edit the file to add in more than two runs
                        # this is so we can choose whatever amount of runs we want and add into the file instead of delete
                        sed -e 's@FEAT_PATH@'"# 4D AVW data or FEAT directory ("$j")\nset feat_files("$j") "$DIREC$s"/func/func"$runNum"/level_one.feat\n\nFEAT_PATH"'@g' \
                            -e 's@EVG_PATH@'"# Higher-level EV value for EV 1 and input "$j"\nset fmri(evg"$j".1) 1\n\nEVG_PATH"'@g' \
                            -e 's@GROUPMEM_PATH@'"# Group membership for input "$j"\nset fmri(groupmem\."$j") 1\n\nGROUPMEM_PATH"'@g' <$i> design_leveltwo"$j".fsf 
    
                    done
    
                    # remove the old template file as its not needed
                    rm design_leveltwo"$((j-1))".fsf
    
                fi           
     
                ((j+=1));
                
                # one statement to make sure the file is named properly
                if [ "$j" -gt "${#myFunc[@]}" ]; then
    
                    # checking if a file exists already
                    if [ -f design_leveltwo"$naming".fsf ]; then
    
                        # remove the old previous one because im not sure it will overwrite properly
                        rm design_leveltwo"$naming".fsf
                    fi
                    
                    # set some final paths to empty
                    for i in "design_leveltwo"$((j-1))".fsf"; do
    
                        sed -e 's@FEAT_PATH@''@g' -e 's@EVG_PATH@''@g' -e 's@GROUPMEM_PATH@''@g' <$i> design_leveltwo"$naming".fsf 
    
                    done
    
                    #echo design_leveltwo"$((j-1))"
                    rm design_leveltwo"$((j-1))".fsf
    
                fi
    
            done
    
    
            tput setaf 2; echo "Run second level analysis"
            tput sgr0; 
    
            # Run the analysis using the fsf file
            feat design_leveltwo"$naming".fsf
                
            # move into folder so it doesnt take up space
            mv design_leveltwo"$naming".fsf level_two"$naming".gfeat/design_leveltwo"$naming".fsf 

		    tput setaf 2; echo "Done!" 
        	    tput sgr0;	


        elif [ "$ind" == "2" ]; then

            # here I run through the number of runs were are to combine specified when you start running this
            # always must be a minimum of two runs to do a second level analysis of per subject
            for runNum in ${myFunc[@]}; do
                
                # starts with the first run specified
                ## Generate fsf file from template
                if [ "$j" == "1" ]; then
                    
                    # here selecting the original named template file
                    for i in "../../template/second_level_force_template.fsf"; do
                        
                        # now we set some standard paths ands masks 
                        # set path 1 as the first run
                        # and resave the file as a new template inside the folder
                        sed -e 's@OUTDIR@'"level_two_force"$naming""'@g' \
                            -e 's@PATH1@'$DIREC$s"/func/func"$runNum"/level_one_force.feat"'@g' \
                            -e 's@THRESH_MASK@'$DIREC"template/PAM50_cord.nii.gz"'@g' \
                            -e 's@NSUBJECTS@'${#myFunc[@]}'@g' <$i> design_leveltwo_force"$j".fsf 
    
                    done
                
                # must need a second run so this is when j=2
                elif [ "$j" == "2" ]; then
                    
                    # take the template file we just named in the previous statement
                    for i in "design_leveltwo_force"$((j-1))".fsf"; do
                        
                        # now we sepcificy the second path we chose and create a new template file
                        sed -e 's@PATH2@'$DIREC$s"/func/func"$runNum"/level_one_force.feat"'@g' \
                            -e 's@OUTLPATH@''@g' <$i> design_leveltwo_force"$j".fsf 
    
                    done
                    
                    # remove the old template file as its not needed
                    #echo design_leveltwo_force"$((j-1))"
                    rm design_leveltwo_force"$((j-1))".fsf
    
                else
                    
                    # take the template file we just named in the previous statement                
                    for i in "design_leveltwo_force"$((j-1))".fsf"; do
     
                        # here we have to edit the file to add in more than two runs
                        # this is so we can choose whatever amount of runs we want and add into the file instead of delete
                        sed -e 's@FEAT_PATH@'"# 4D AVW data or FEAT directory ("$j")\nset feat_files("$j") "$DIREC$s"/func/func"$runNum"/level_one_force.feat\n\nFEAT_PATH"'@g' \
                            -e 's@EVG_PATH@'"# Higher-level EV value for EV 1 and input "$j"\nset fmri(evg"$j".1) 1\n\nEVG_PATH"'@g' \
                            -e 's@GROUPMEM_PATH@'"# Group membership for input "$j"\nset fmri(groupmem\."$j") 1\n\nGROUPMEM_PATH"'@g' <$i> design_leveltwo_force"$j".fsf 
    
                    done
    
                    # remove the old template file as its not needed
                    rm design_leveltwo_force"$((j-1))".fsf
    
                fi           
     
                ((j+=1));
                
                # one statement to make sure the file is named properly
                if [ "$j" -gt "${#myFunc[@]}" ]; then
    
                    # checking if a file exists already
                    if [ -f design_leveltwo_force"$naming".fsf ]; then
    
                        # remove the old previous one because im not sure it will overwrite properly
                        rm design_leveltwo_force"$naming".fsf
                    fi
                    
                    # set some final paths to empty
                    for i in "design_leveltwo_force"$((j-1))".fsf"; do
    
                        sed -e 's@FEAT_PATH@''@g' -e 's@EVG_PATH@''@g' -e 's@GROUPMEM_PATH@''@g' <$i> design_leveltwo_force"$naming".fsf 
    
                    done
    
                    #echo design_leveltwo"$((j-1))"
                    rm design_leveltwo_force"$((j-1))".fsf
    
                fi
    
            done
    
    
            tput setaf 2; echo "Run second level analysis"
            tput sgr0; 
    
            # Run the analysis using the fsf file
            feat design_leveltwo_force"$naming".fsf
    
            # move into folder so it doesnt take up space
            mv design_leveltwo_force"$naming".fsf level_two_force"$naming".gfeat/design_leveltwo_force"$naming".fsf 

		    tput setaf 2; echo "Done!" 
        	    tput sgr0;	


        elif [ "$ind" -gt "2" ]; then

            # changes the template file that it loops through
            if [ "$ind" == "3" ]; then
                levelFile="level_one_force_FLOB.feat"
                outputName="level_two_force_FLOB"
                templateName="second_level_force_FLOB_template.fsf"
            elif [ "$ind" == "4" ]; then
                levelFile="level_one_force_smooth_FLOB.feat"
                outputName="level_two_force_FLOB_smooth"
                templateName="second_level_force_FLOB_template.fsf"
            elif [ "$ind" == "5" ]; then
                levelFile="level_one_FLOB.feat"
                outputName="level_two_FLOB"
                templateName="second_level_FLOB_template.fsf"
            fi

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
                        sed -e 's@OUTDIR@'"$outputName$naming"'@g' \
                            -e 's@PATH1@'$DIREC$s"/func/func"$runNum"/"$levelFile""'@g' \
                            -e 's@THRESH_MASK@'$DIREC"template/PAM50_cervical_cord.nii.gz"'@g' \
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
    
    
            tput setaf 2; echo "Run second level analysis"
            tput sgr0; 
    
            # Run the analysis using the fsf file
            feat design_leveltwo_force_FLOB"$naming".fsf

    #######################################################################
            # move into folder so it doesnt take up space
            mv design_leveltwo_force_FLOB"$naming".fsf  "$outputName$naming".gfeat/design_leveltwo_force_FLOB"$naming".fsf 

		    tput setaf 2; echo "Done!" 
        	    tput sgr0;	


        fi

			 			
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