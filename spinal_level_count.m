 %% template matching
clear all

%%% SPINE

% addpath('/Users/pirondinilab/spinalcordtoolbox/cbiNifti');
addpath('D:\NHP_code\cbiNifti')

% varibales to set up before
subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_007','SBSN_H_008','SBSN_H_010'}; 
zScore = 2.3;

copeFile = 'cope1.feat';
% copeFile = 'cope4.feat';
% copeFile = 'cope7.feat';

%% THINGS TO ADD

% gunzip('D:\SBSN\Data\Spine/template/PAM50_cervical_cord.nii.gz');
[tempLevels, ~] = cbiReadNifti('D:\SBSN\Data\Spine/template/PAM50_cervical_cord.nii');
% gunzip('D:\SBSN\Data\Spine/template\spinal_levels/spinal_level_05.nii.gz');
[tempLevels5, ~] = cbiReadNifti('D:\SBSN\Data\Spine/template\spinal_levels/spinal_level_05.nii');
% gunzip('D:\SBSN\Data\Spine/template\spinal_levels/spinal_level_06.nii.gz');
[tempLevels6, ~] = cbiReadNifti('D:\SBSN\Data\Spine/template\spinal_levels/spinal_level_06.nii');
% gunzip('D:\SBSN\Data\Spine/template\spinal_levels/spinal_level_07.nii.gz');
[tempLevels7, ~] = cbiReadNifti('D:\SBSN\Data\Spine/template\spinal_levels/spinal_level_07.nii');
% gunzip('D:\SBSN\Data\Spine/template\spinal_levels/spinal_level_08.nii.gz');
[tempLevels8, ~] = cbiReadNifti('D:\SBSN\Data\Spine/template\spinal_levels/spinal_level_08.nii');

tempLevels5(:,:, 871:end) = 0;
tempLevels8(:,:, 1:734) = 0;

tempLevels5(tempLevels5 >= 0.0002) = 1;
tempLevels5(tempLevels5 < 0.0002) = 0; 
tempLevels6(tempLevels6 >= 0.0002) = 1;
tempLevels6(tempLevels6 < 0.0002) = 0; 
tempLevels7(tempLevels7 >= 0.0002) = 1;
tempLevels7(tempLevels7 < 0.0002) = 0; 
tempLevels8(tempLevels8 >= 0.0002) = 1;
tempLevels8(tempLevels8 < 0.0002) = 0; 

tempLevelsAll = tempLevels5 + tempLevels6 + tempLevels7 + tempLevels8;

tempLevelsAll(tempLevelsAll >= 1) = 1;
% 4, Spinal level C5, spinal_level_05.nii.gz
% 5, Spinal level C6, spinal_level_06.nii.gz
% 6, Spinal level C7, spinal_level_07.nii.gz
% 7, Spinal level C8, spinal_level_08.nii.gz
% 
% tempLevels5 = tempLevels5(tempLevels == 1);
% tempLevels8 = tempLevels8(tempLevels == 1);

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
        if subjectFolder(folder).isdir && contains(subjectFolder(folder).name, 'level_two_FLOB1234')

            disp(subjectFolder(folder).name)

            fileName = strsplit(subjectFolder(folder).name, '.');

            if ~exist(fullfile(direc, subjectFolder(folder).name, copeFile, '/thresh_zstat1.nii'))

                gunzip(fullfile(direc, subjectFolder(folder).name,  copeFile, '/thresh_zstat1.nii.gz'));

            end
            [dataFile, ~] = cbiReadNifti(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii'));

            disp(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii'))

            if contains(fileName{1}, 'smooth')
                allData{i}{fileCounter, 2} = fileName{1}(28:end);
            else
                allData{i}{fileCounter, 2} = fileName{1}(21:end);                
            end

%             dataFile5 = dataFile;
%             dataFile5(tempLevels5 >= 0.002) = 1;
%             dataFile5(dataFile5~=1) = 0;
% 
%             dataFile6 = dataFile;
%             dataFile6(tempLevels6 >= 0.002) = 1;
%             dataFile6(tempLevels6 < 0.002) = 0;            
% 
%             dataFile7 = dataFile;
%             dataFile7(tempLevels7 >= 0.002) = 1;
%             dataFile7(tempLevels7 < 0.002) = 0;
% 
%             dataFile8 = dataFile;
%             dataFile8(tempLevels8 >= 0.002) = 1;
%             dataFile8(tempLevels8 < 0.002) = 0;
%                 figure;
%                 plot(squeeze(sum(sum(dataFile5(:,:,717:888), 1), 2)))
%                 hold on
%                 plot(squeeze(sum(sum(dataFile6(:,:,717:888), 1), 2)))
%                 plot(squeeze(sum(sum(dataFile7(:,:,717:888), 1), 2)))
%                 plot(squeeze(sum(sum(dataFile8(:,:,717:888), 1), 2)))
            % number of active voxels
            cAll = dataFile.*tempLevelsAll;
            cAll(cAll>zScore) = 1;
%             c6 = dataFile.*tempLevels6;
%             c6(c6>zScore) = 1;
%             c7 = dataFile.*tempLevels7;
%             c7(c7>zScore) = 1;
%             c8 = dataFile.*tempLevels8;
%             c8(c8>zScore) = 1;

            % number of active voxels
            allData{i}{fileCounter, 1}{1, 1} = squeeze(sum(cAll,[1,2]));
%             allData{i}{fileCounter, 1}{1, 2} = squeeze(sum(c6,[1,2]));
%             allData{i}{fileCounter, 1}{1, 3} = squeeze(sum(c7,[1,2]));
%             allData{i}{fileCounter, 1}{1, 4} = squeeze(sum(c8,[1,2]));
            
%             717:888
            fileCounter = fileCounter + 1;
        end

    end
end

% 
% allDataSmooth = allData;
% for i = 1:length(allData)
%     allData{i,1}(6:10, :) = [];
%     allDataSmooth{i,1}(1:5, :) = [];
% end
% 
% 
% save('savedData/allDataSLcope7', 'allData')
% save('savedData/allDataSmoothSLcope7', 'allDataSmooth')