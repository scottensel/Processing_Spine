%% template matching
clear all

%%% SPINE
addpath('D:\NHP_code\cbiNifti')

% varibales to set up before
subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_007','SBSN_H_008','SBSN_H_010'};
% subName = {'SBSN_H_010','SBSN_H_008'};

zScore = 1.5;
%% THINGS TO ADD
% I need to remove the timings for the first x volumes
% these get removed from the original length of the images
gunzip('D:\SBSN\Data\Spine\template\PAM50_levels.nii.gz');
[pamLevels, ~] = cbiReadNifti('D:\SBSN\Data\Spine\template\PAM50_levels.nii');

allData = {};
for i = 1:length(subName)
    direc = fullfile('D:\SBSN\Data\Spine', subName{i}, 'func');
    subjectFolder = dir(direc);

    disp(subName{i})
    
    allData{i, 1} = {};
    allData{i, 2} = subName{i};
    
    fileCounter = 1;  
    for folder = 3:length(subjectFolder)
        
        %is dir and name contains gfeat
        if subjectFolder(folder).isdir && contains(subjectFolder(folder).name, 'gfeat')
            
            disp(subjectFolder(folder).name)
            
            fileName = strsplit(subjectFolder(folder).name, '.');
            
            if ~exist(fullfile(direc, subjectFolder(folder).name, '/cope1.feat/thresh_zstat1.nii'))
                
                gunzip(fullfile(direc, subjectFolder(folder).name, '/cope1.feat/thresh_zstat1.nii.gz'));
                
            end        
            [dataFile, ~] = cbiReadNifti(fullfile(direc, subjectFolder(folder).name, '/cope1.feat/thresh_zstat1.nii'));
            
            allData{i}{fileCounter, 2} = fileName{1}(10:end);
            
            counter = 1;
            vertLevel2 = [];
            for level = [4,5,6,7]
                % given indices of vertebral level
                vertLevel = dataFile(pamLevels == level);
                
                %going to combione all the logical maps to calculate dice
                %coeff
                vertLevel2 = [vertLevel2; vertLevel];
                
                % now must extranct number of voxels active and their
                % z-scores (vertlevel is a map of all voxels in that level)
                
                % number of active voxels
                allData{i}{fileCounter, 1}{1, counter} = sum(reshape(vertLevel>zScore, 1, []));
                
                %z-score of active voxels
                allData{i}{fileCounter, 1}{2, counter} = vertLevel(vertLevel>zScore);
                
                counter = counter +1;
            end
            
            allData{i}{fileCounter, 3} = vertLevel2;
            
            fileCounter = fileCounter + 1;
        end
        
    end
end


save allData.mat allData


