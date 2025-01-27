%% template matching
clear all

%%% SPINE

% addpath('/Users/pirondinilab/spinalcordtoolbox/cbiNifti');
addpath('D:\NHP_code\cbiNifti')

% varibales to set up before
subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_007','SBSN_H_008','SBSN_H_010'}; 
zScore = 1.5;

copeFile = 'cope1.feat';
% copeFile = 'cope4.feat';
% copeFile = 'cope7.feat';

%% THINGS TO ADD

% gunzip('/Volumes/rnelshare/projects/human/brain_spine_stroke_SBSN/Data/sreya/Spine/template/PAM50_levels.nii.gz');
% [tempLevels, ~] = cbiReadNifti('/Volumes/rnelshare/projects/human/brain_spine_stroke_SBSN/Data/sreya/Spine/template/PAM50_levels.nii');
gunzip('D:\SBSN\Data\Spine\template\PAM50_levels.nii.gz');
[spinalLevels, ~] = cbiReadNifti('D:\SBSN\Data\Spine\template/PAM50_levels.nii');
gunzip('D:\SBSN\Data\Spine\template\PAM50_rl.nii.gz');
[lrLevels, ~] = cbiReadNifti('D:\SBSN\Data\Spine/template/PAM50_rl.nii');
gunzip('D:\SBSN\Data\Spine\template\PAM50_dv.nii.gz');
[dvLevels, ~] = cbiReadNifti('D:\SBSN\Data\Spine/template/PAM50_dv.nii');

allData = {};
for i = 1:length(subName)
    %direc = fullfile('/Volumes/MyPassport/Sub_Data/new_data/new_spine', subName{i}, 'func');
%     direc = fullfile('/Volumes/rnelshare/projects/human/brain_spine_stroke_SBSN/Data/sreya/Spine', subName{i}, 'func');
    direc = fullfile('D:\SBSN\Data\Spine', subName{i}, 'func');

    subjectFolder = dir(direc);

    disp(subName{i})

    allData{i, 1} = {};
    allData{i, 2} = subName{i};

    fileCounter = 1;
    for folder = 3:length(subjectFolder)

        %is dir and name contains gfeat
        if subjectFolder(folder).isdir && contains(subjectFolder(folder).name, 'level_two_force_FLOB')

            disp(subjectFolder(folder).name)

            fileName = strsplit(subjectFolder(folder).name, '.');

            if ~exist(fullfile(direc, subjectFolder(folder).name, copeFile, '/thresh_zstat1.nii'))

                gunzip(fullfile(direc, subjectFolder(folder).name,  copeFile, '/thresh_zstat1.nii.gz'));

            end
            [dataFile, ~] = cbiReadNifti(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii'));

            if contains(fileName{1}, 'smooth')
                allData{i}{fileCounter, 2} = fileName{1}(28:end);
            else
                allData{i}{fileCounter, 2} = fileName{1}(21:end);                
            end

            counter = 1;
            vertLevel2 = [];
%             for level = [4,5,6,7]
            for level = [1,2]
                % given indices of vertebral level
                vertLevel = dataFile(tempLevels == level);

                %going to combione all the logical maps to calculate dice
                %coeff
                vertLevel2 = [vertLevel2; vertLevel];

                % now must extranct number of voxels active and their
                % z-scores (vertlevel is a map of all voxels in that level)

                % number of active voxels
                allData{i}{fileCounter, 1}{1, counter} = sum(reshape(vertLevel>zScore, 1, []));

                %z-score of active voxels
                allData{i}{fileCounter, 1}{2, counter} = vertLevel(vertLevel>zScore);

                counter = counter + 1;
                allData{i}{fileCounter, 3} = vertLevel2;
                if level >= 4
                    allData{i}{fileCounter, level} = vertLevel;
                else
                    allData{i}{fileCounter, 4} = vertLevel;
                end

            end

            fileCounter = fileCounter + 1;
        end

    end
end


allDataSmooth = allData;
for i = 1:length(allData)
    allData{i,1}(6:10, :) = [];
    allDataSmooth{i,1}(1:5, :) = [];
end

% save('savedData/allDataPAMcope7', 'allData')
% save('savedData/allDataSmoothPAMcope7', 'allDataSmooth')
% 
% save('savedData/allDataRLcope7', 'allData')
% save('savedData/allDataSmoothRLcope7', 'allDataSmooth')

save('savedData/allDataDVcope7', 'allData')
save('savedData/allDataSmoothDVcope7', 'allDataSmooth')


