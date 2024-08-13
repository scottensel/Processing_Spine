#!/bin/bash
#
# SBSN framework - Preprocessing
# June 2022
#
# This points to the path of your subjects

# List of subjects
#declare -a sub=("SBSN_H_001")
#declare -a sub=("SBSN_H_001" "SBSN_H_002" "SBSN_H_003" "SBSN_H_004")
#declare -a sub=("SBSN_H_003")
#declare -a sub=("SBSN_H_004")

# Path of the folder containing all data
#DIREC="/mnt/d/SBSN/Data/Spine/"
#DIREC="/mnt/c/Users/scott/Documents/"

####################################
# Here edit which paths you want to be able to choose from
# the folder after this should contain all of the SBSN_X_00X folders
myPaths=("/mnt/d/SBSN/Data/Spine/" "/mnt/c/Users/scott/Documents/Spine/")
####################################

echo "Choose a Saved Path:"
for i in ${!myPaths[@]}; do
  echo "$i -> path ${myPaths[$i]}" 
done

tput setaf 6; 
echo -n "Path Number: "
tput sgr0;
read indPath

# setting the variable to correct path
DIREC=${myPaths[$indPath]}


echo "Choose a Subject"
mySub=("SBSN_H_" "SBSN_M_" "SBSN_S_")
for i in ${!mySub[@]}; do
  echo "$i -> path ${mySub[$i]}" 
done
tput setaf 6; 
echo -n "Subject Type: "
tput sgr0;
read indSub

echo "Choose Subject Number"
tput setaf 6;
read -p "Enter Subject Numbers to Process: " -a myRuns
tput sgr0;

sub=()
for i in ${!myRuns[@]}; do

    if [ ${myRuns[$i]} -gt 9 ]; then

        if [ ${myRuns[$i]} -gt 99 ]; then

            sub+=("${mySub[$indSub]}${myRuns[$i]}")

        else

            sub+=("${mySub[$indSub]}"0"${myRuns[$i]}")

        fi
    
    else

        sub+=("${mySub[$indSub]}"00"${myRuns[$i]}")

    fi
done


echo "Choose Functional Files to run"
echo "Enter 0 for all"
tput setaf 6; 
read -p "Enter Functional Runs to Process: " -a myFunc
tput sgr0;

if [ "$myFunc" == "0" ]; then
    myFunc=(1 2 3 4 5 6)
fi

echo $DIREC
echo "${sub[@]}"
echo "${myFunc[@]}"
