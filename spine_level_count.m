%%
close all

%% Frame Displacement
clear all

% addpath('/Users/pirondinilab/spinalcordtoolbox/cbiNifti');
addpath('D:\NHP_code\cbiNifti')

% varibales to set up before
subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_007','SBSN_H_008','SBSN_H_010','SBSN_H_017','SBSN_H_018','SBSN_H_019'}; 
% subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_007','SBSN_H_008'}; 

allData = {};
for i = 1:length(subName)

    allData{i, 1} = {};
    allData{i, 2} = subName{i};

    for j = 0:6
        direc = fullfile('D:\SBSN\Data\Spine', subName{i}, 'func', ['func', num2str(j)]);
        disp(fullfile(direc, 'fmri_spine_moco_mean_tsnr_PAM50.nii'))

        try
            motionData = importdata(fullfile('D:\SBSN\Data\Spine', subName{i}, 'func', ['func', num2str(j)], 'moco_params.tsv'));
    %         motionData = motionData(11:end,:);
            motionData = motionData.data;
            dataDiff = motionData(2:end,:)-motionData(1:(end-1),:);
            FD = sum(abs(dataDiff), 2);
    
            FD = FD(5:end,:);
    
            % number of active voxels
            allData{i, 1}{j+1, 1} = [mean(reshape(FD, [], 1),'omitnan'), std(reshape(FD, [], 1),'omitnan')];
            allData{i, 1}{j+1, 2} = [subName{i}, ' func', num2str(j)];
        catch
            continue
        end

    
    end
end

% this doesnt grab the rest part
for i = 1:length(allData)

    for j = 1:length(allData{i,1})

        frameD(i,j) = allData{i,1}{j,1}(1);

    end

end

% all runs but rest run
disp('task runs mean')
mean(mean(frameD(:,2:end))) 
std(mean(frameD(:,2:end)))/sqrt(length(mean(frameD(:,2:end))))

% only the rest run
disp('rest runs')
mean(frameD(:, 1))
std(frameD(:, 1))/sqrt(length(frameD(:,1)))

% means across runs
disp('all runs')
mean(frameD, 1) 
std(frameD, 1)/sqrt(length(mean(frameD)))

%% TSNR
clear all

% addpath('/Users/pirondinilab/spinalcordtoolbox/cbiNifti');
addpath('D:\NHP_code\cbiNifti')

% varibales to set up before
subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_019','SBSN_H_007','SBSN_H_008','SBSN_H_010','SBSN_H_017','SBSN_H_018'}; 
% subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_007','SBSN_H_008'}; 
subSplit = 5;

% fmri_brain_moco_mean_tsnr_MNI152.nii.gz
% gunzip('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii.gz');
% [spineLevels, ~] = cbiReadNifti('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii');
% spineLevels(spineLevels < 5) = 0;
% spineLevels(spineLevels > 8) = 0;
% spineLevels(spineLevels > 0) = spineLevels(spineLevels > 0) - 4;

% gunzip('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii.gz');
[spineLevels, ~] = cbiReadNifti('F:\SMA_HOLDER\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii');
spineLevels(spineLevels < 1) = 0;
spineLevels(spineLevels > 8) = 0;
% spineLevels(spineLevels > 0) = spineLevels(spineLevels > 0) - 1;

% gunzip('D:\SBSN\Data\Spine\template\PAM50_levels.nii.gz');
% [spineLevels, ~] = cbiReadNifti('D:\SBSN\Data\Spine\template\PAM50_levels.nii');
% spineLevels(spineLevels < 4) = 0;
% spineLevels(spineLevels > 7) = 0;
% spineLevels(spineLevels > 0) = spineLevels(spineLevels > 0) - 3;

allData = {};
for i = 1:length(subName)

    allData{i, 1} = {};
    allData{i, 2} = subName{i};

    for j = 0:6
        direc = fullfile('D:\SBSN\Data\Spine', subName{i}, 'func', ['func', num2str(j)]);
    
        subjectFolder = dir(direc);
    
        disp(subName{i})

        if ~exist(fullfile(direc, 'fmri_spine_moco_mean_tsnr_PAM50.nii'))

            gunzip(fullfile(direc, 'fmri_spine_moco_mean_tsnr_PAM50.nii.gz'));

        end

        [dataFile, ~] = cbiReadNifti(fullfile(direc, 'fmri_spine_moco_mean_tsnr_PAM50.nii'));

        disp(fullfile(direc, 'fmri_spine_moco_mean_tsnr_PAM50.nii'))

        mag = dataFile(spineLevels>=1);

        % here we will do it based on each vertebral level
        magSeperate = [];
        for k = 1:max(reshape(spineLevels, [], 1))

            magSeperate(k) = mean(dataFile(spineLevels==k));

        end

        % number of active voxels
        allData{i, 1}{j+1, 1} = [mean(reshape(mag, [], 1),'omitnan'), std(reshape(mag, [], 1),'omitnan')];
        allData{i, 1}{j+1, 2} = [subName{i}, ' func', num2str(j)];
        allData{i, 1}{j+1, 3} = [magSeperate];

    end
end

segTSNR = [];
for i = 1:length(allData)
    for j = 1:length(allData{i,1})

        tsnr(i, j) = allData{i, 1}{j, 1}(1);

        segTSNR(i, j, 1:max(reshape(spineLevels, [], 1))) = allData{i, 1}{j, 3};

    end

end

%segTSNR (x,:,:) x is across runs
%segTSNR (:,x,:) x is across subjects
%segTSNR (:,:,x) x is the vertera ( C1 - C8)
% segTSNR(:,1,8) is one run across all subjects in C8
% segTSNR(1,:,8) is on subject across all runs in C8
mean(mean(tsnr)) 


% total tSNR
plotCreator(tsnr, 1:7, subSplit);
make_pretty
xlim([0.75, 7.25])
ylabel('TSNR')
xlabel('Run Number');
title(sprintf('TSNR of Each Run'));
xticklabels({'Rest','1','2','3','4','5','6'})

saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_tsnr.png');
saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_tsnr.svg');

% all runs but rest run
mean(mean(tsnr(:, 2:end))) 
std(mean(tsnr(:, 2:end)))/sqrt(length(mean(tsnr(:,2:end))))

% only the rest run
mean(tsnr(:, 1))
std(tsnr(:, 1))/sqrt(length(tsnr(:, 1)))

% means across runs
mean(tsnr, 1) 
std(tsnr, 1)/sqrt(length(mean(tsnr)))

% Define the run combinations (1 to 6)
runCombinations = 1:6;

% Spearman's correlation for active voxels
[rho_activeVoxels, pval_activeVoxels] = corr(runCombinations', mean(tsnr(:, 2:end))', 'Type', 'Spearman');

% Display results for active voxels
disp('Spearman correlation for active voxels:');
disp('Correlation coefficients (rho):');
disp(rho_activeVoxels);
disp('P-values:');
disp(pval_activeVoxels);
% runRepeatedMeasuresANOVA(tsnr)

% Individual tSNR
%segTSNR (x,:,:) x is across runs
%segTSNR (:,x,:) x is across subjects
%segTSNR (:,:,x) x is the vertera ( C1 - C8)
% segTSNR(:,1,8) is one run across all subjects in C8
% segTSNR(1,:,8) is on subject across all runs in C8
plotCreator(mean(segTSNR(:, :, 2:4), 3), 1:7, subSplit);
make_pretty
xlim([0.75,7.25])
ylabel('TSNR')
xlabel('Run Number');
title(sprintf('TSNR of C2-C4 in Each Run'));
xticklabels({'Rest','1','2','3','4','5','6'})


plotCreator(mean(segTSNR(:, :, 7:8), 3), 1:7, subSplit);
make_pretty
xlim([0.75,7.25])
ylabel('TSNR')
xlabel('Run Number');
title(sprintf('TSNR of C7-C8 in Each Run'));
xticklabels({'Rest','1','2','3','4','5','6'})


plotCreator(mean(segTSNR(:, :, 1), 3), 1:7, subSplit);
make_pretty
xlim([0.75,7.25])
ylabel('TSNR')
xlabel('Run Number');
title(sprintf('TSNR of C1 in Each Run'));
xticklabels({'Rest','1','2','3','4','5','6'})

% plotCreator(segTSNR(:, :, 1), 1:7, subSplit);
% make_pretty
% xlim([0.75,7.25])
% ylabel('TSNR')
% xlabel('Run Number');
% title(sprintf('TSNR of C2-C3 in Each Run'));
% xticklabels({'Rest','1','2','3','4','5','6'})
% 
% 
% plotCreator(segTSNR(:, :, 7), 1:7, subSplit);
% make_pretty
% xlim([0.75,7.25])
% ylabel('TSNR')
% xlabel('Run Number');
% title(sprintf('TSNR of C7-C8 in Each Run'));
% xticklabels({'Rest','1','2','3','4','5','6'})

%% template matching
clear all

%%% SPINE

% addpath('/Users/pirondinilab/spinalcordtoolbox/cbiNifti');
addpath('D:\NHP_code\cbiNifti')

% varibales to set up before
subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_019','SBSN_H_007','SBSN_H_008','SBSN_H_010','SBSN_H_017','SBSN_H_018'}; 
% subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_008','SBSN_H_010'}; 
% subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_007','SBSN_H_008'}; 

zScore = 2.1;
subSplit = 5;

copeFile = 'cope1.feat';
analysisFile = 'level_two_all_force_FLOB';
% gunzip('D:\SBSN\Data\Spine\template\PAM50_levels.nii.gz');
% [spineLevels, ~] = cbiReadNifti('D:\SBSN\Data\Spine\template\PAM50_levels.nii');
% spineLevels(spineLevels < 4) = 0;
% spineLevels(spineLevels > 7) = 0;
% spineLevels(spineLevels > 0) = spineLevels(spineLevels > 0) - 3;

% gunzip('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii.gz');
[spineLevels, ~] = cbiReadNifti('F:\SMA_HOLDER\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii');
spineLevels(spineLevels < 5) = 0;
spineLevels(spineLevels > 8) = 0;
spineLevels(spineLevels > 0) = spineLevels(spineLevels > 0) - 4;



% gunzip('D:\SBSN\Data\Spine\template\PAM50_rl.nii.gz');
[lrLevels, ~] = cbiReadNifti('D:\SBSN\Data\Spine\template\PAM50_all_rl.nii');
% gunzip('D:\SBSN\Data\Spine\template\PAM50_dv.nii.gz');
[dvLevels, ~] = cbiReadNifti('D:\SBSN\Data\Spine\template\PAM50_all_dv.nii');

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
        if subjectFolder(folder).isdir && contains(subjectFolder(folder).name, analysisFile)

            disp(subjectFolder(folder).name)

            fileName = strsplit(subjectFolder(folder).name, '.');

            if ~exist(fullfile(direc, subjectFolder(folder).name, copeFile, 'thresh_zstat1.nii'))

                gunzip(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii.gz'));

            end

            [dataFile, ~] = cbiReadNifti(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii'));

            if i == 1
                dataFile = flip(dataFile, 1);
            end

            disp(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii'))

            numVoxels = sum(sum(sum((dataFile>=zScore).*(spineLevels>=1))));
            mag = dataFile(spineLevels>=1);

            numVoxelsSeperate = [];
            magSeperate = [];
            for j = 1:max(reshape(spineLevels, [], 1))
                numVoxelsSeperate(j) = sum(sum(sum((dataFile>=zScore).*(spineLevels==j))));
                numVoxelsSeperate(j) = numVoxelsSeperate(j)/sum(sum(sum(spineLevels==j)))*100;
                var = dataFile(spineLevels==j);
                magSeperate(j) = mean(var(var>zScore));  
            end

            for j = 1:2
                numVoxelsLR(j) = sum(sum(sum((dataFile>=zScore).*(lrLevels==j))));
                numVoxelsLR(j) = numVoxelsLR(j)/sum(sum(sum(lrLevels==j)))*100;

                numVoxelsDV(j) = sum(sum(sum((dataFile>=zScore).*(dvLevels==j))));
                numVoxelsDV(j) = numVoxelsDV(j)/sum(sum(sum(dvLevels==j)))*100;

                var = dataFile(lrLevels==j);
                magSeperateLR(j) = mean(var(var>zScore)); 

                var = dataFile(dvLevels==j);
                magSeperateDV(j) = mean(var(var>zScore)); 
                
            end

            % number of active voxels
            allData{i, 1}{fileCounter, 1} = [numVoxels/sum(sum(sum(spineLevels>=1)))*100, mean(mag(mag>zScore)), std(mag(mag>zScore))];
            allData{i, 1}{fileCounter, 2} = subjectFolder(folder).name;
            allData{i, 1}{fileCounter, 3} = [numVoxelsSeperate; magSeperate];


            allData{i, 1}{fileCounter, 4} = [numVoxelsLR; magSeperateLR];         
            allData{i, 1}{fileCounter, 5} = [numVoxelsDV; magSeperateDV];


            fileCounter = fileCounter + 1;
        end

    end
end


for i = 1:length(allData)
    for j = 1:length(allData{i,1})

        activeVoxels(i,j) = allData{i,1}{j,1}(1);
        zScores(i,j) = allData{i,1}{j,1}(2);

        actVoxSeg(i, j, 1:max(reshape(spineLevels, [], 1))) = allData{i, 1}{j, 3}(1,:);
        zSeg(i, j, 1:max(reshape(spineLevels, [], 1))) = allData{i, 1}{j, 3}(2,:);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        actVoxelsSeg2(i,:) = allData{i,1}{1,3}(1,:);
        zSeg2(i,:) = allData{i,1}{1,3}(2,:);

        actVoxelsSeg3(i,:) = allData{i,1}{2,3}(1,:);
        zSeg3(i,:) = allData{i,1}{2,3}(2,:);

        actVoxelsSeg5(i,:) = allData{i,1}{4,3}(1,:);
        zSeg5(i,:) = allData{i,1}{4,3}(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        actVoxelsSeg4(i,:) = allData{i,1}{3,3}(1,:);
        zSeg4(i,:) = allData{i,1}{3,3}(2,:);

        actVoxelsSeg6(i,:) = allData{i,1}{5,3}(1,:);
        zSeg6(i,:) = allData{i,1}{5,3}(2,:);

        %%%%% 1v 2d 1 l 2 r
        lrSeg4(i,:) = allData{i,1}{3,4}(1,:);
        zlrSeg4(i,:) = allData{i,1}{3,4}(2,:);

        dvSeg4(i,:) = allData{i,1}{3,5}(1,:);
        zdvSeg4(i,:) = allData{i,1}{3,5}(2,:);

%         lrSeg6(i,:) = allData{i,1}{5,4}(1,:);
%         zlrSeg6(i,:) = allData{i,1}{5,4}(2,:);
% 
%         dvSeg6(i,:) = allData{i,1}{5,5}(1,:);
%         zdvSeg6(i,:) = allData{i,1}{5,5}(2,:);
    end

end

% for i = 1:nRegions
% 
    [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(actVoxelsSeg2(:), actVoxelsSeg3(:), 10000, 0.05, 10);
%     brainNames{i}
    ci95, rejectNull

    [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(actVoxelsSeg2(:), actVoxelsSeg4(:), 10000, 0.05, 10);
%     brainNames{i}
    ci95, rejectNull

    [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(actVoxelsSeg2(:), actVoxelsSeg5(:), 10000, 0.05, 10);
%     brainNames{i}
    ci95, rejectNull

    [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(actVoxelsSeg2(:), actVoxelsSeg6(:), 10000, 0.05, 10);
%     brainNames{i}
    ci95, rejectNull

    [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(actVoxelsSeg3(:), actVoxelsSeg4(:), 10000, 0.05, 10);
%     brainNames{i}
    ci95, rejectNull

    [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(actVoxelsSeg3(:), actVoxelsSeg5(:), 10000, 0.05, 10);
%     brainNames{i}
    ci95, rejectNull

    [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(actVoxelsSeg3(:), actVoxelsSeg6(:), 10000, 0.05, 10);
%     brainNames{i}
    ci95, rejectNull

    [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(actVoxelsSeg4(:), actVoxelsSeg5(:), 10000, 0.05, 10);
%     brainNames{i}
    ci95, rejectNull

    [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(actVoxelsSeg4(:), actVoxelsSeg6(:), 10000, 0.05, 10);
%     brainNames{i}
    ci95, rejectNull

    [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(actVoxelsSeg5(:), actVoxelsSeg6(:), 10000, 0.05, 10);
%     brainNames{i}
    ci95, rejectNull    
% % end

figure;
hBar=barh([mean(actVoxelsSeg4)']);
X=get(hBar,'XData').'+[hBar.XOffset];
hold on  %4 runs
hEB = errorbar([mean(actVoxelsSeg4)'], X, [(std(actVoxelsSeg4)/sqrt(length(actVoxelsSeg4)))'], 'horizontal', '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4, 1), 1))/10;
scatter(actVoxelsSeg4, [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1)], 30, 'k','o','filled'); 
set (gca,'YDir','reverse')
yticks(1:length(1:4)); yticklabels({'C5','C6','C7','C8'})
ylabel('Cervical level')
xlabel('Active Voxels');
title(sprintf('Average Active Voxel 4 vs 5 Runs combined'));
% ylim([4.5,8.5])
make_pretty

% Save the plot as a PNG image
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_voxel_area.png');
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_voxel_area.svg');
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel4_area_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel4_area_z', num2str(zScore), '.svg']);


plotCreator(activeVoxels, 1:5, subSplit);
make_pretty
xlim([0.75,5.25])
xlabel('Run Combination')
ylabel('Active Voxels');
title(sprintf('Average Active Successive runs'));
xticks(1:5)
xticklabels({'1-2','1-3','1-4','1-5','1-6'})

% Save the plot as a PNG image
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_voxel_success.png');
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_voxel_success.svg');
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_success_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_success_z', num2str(zScore), '.svg']);


% Define the run combinations (1 to 6)
runCombinations = 1:5;

% Spearman's correlation for active voxels
[rho_activeVoxels, pval_activeVoxels] = corr(runCombinations', mean(activeVoxels)', 'Type', 'Spearman');

% Display results for active voxels
disp('Spearman correlation for active voxels:');
disp('Correlation coefficients (rho):');
disp(rho_activeVoxels);
disp('P-values:');
disp(pval_activeVoxels);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time = [1, 2, 3, 4, 5];
figure;
for i = 1:length(subName)
    findSlopePtsSingle(activeVoxels(i, :), time, i, subSplit)
    hold on
end
make_pretty
xlim([0.75,5.25])
xlabel('Run Combination')
ylabel('Active Voxels');
title(sprintf('Line of Best Subject'));
xticks(1:5)
xticklabels({'1-2','1-3','1-4','1-5','1-6'})
% Save the plot as a PNG image
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_SUPP_fig', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_SUPP_fig', num2str(zScore), '.svg']);



time = [1, 2, 3, 4, 5];
figure;
for i = 1:length(subName)
    findSlopePtsSingle(zScores(i, :), time, i, subSplit)
    hold on
end
make_pretty
xlim([0.75,5.25])
xlabel('Run Combination')
ylabel('zScore');
title(sprintf('Line of Best Subject'));
xticks(1:5)
xticklabels({'1-2','1-3','1-4','1-5','1-6'})
% Save the plot as a PNG image
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_SUPP_fig', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_SUPP_fig', num2str(zScore), '.svg']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotCreator(diff(activeVoxels')', 1:4, length(subName));
% % figure;
% % plot(diff(activeVoxels'), '.-r')
% % hold on
% % plot(mean(diff(activeVoxels')'), 'k')
% % errorbar(1:4,mean(diff(activeVoxels')'), std(diff(activeVoxels')')/sqrt(length(activeVoxels)), 'Color','black')
% make_pretty
% xlim([0.75,4.25])
% xlabel('Run Combination')
% ylabel('Active Voxels');
% title(sprintf('Difference Average Active Successive runs'));
% xticks(1:4)
% xticklabels({'2-3','3-4','4-5','5-6'})
% 
% % Save the plot as a PNG image
% % saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_voxel_diff.png');
% % saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_voxel_diff.svg');
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_diff_z', num2str(zScore), '.png']);
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_diff_z', num2str(zScore), '.svg']);


plotCreator(zScores, 1:5, subSplit);
make_pretty
xlim([0.75,5.25])
xlabel('Run Combination')
ylabel('Z-Score');
title(sprintf('Average Zscore Successive runs'));
xticks(1:5)
xticklabels({'1-2','1-3','1-4','1-5','1-6'})

% Save the plot as a PNG image
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_success_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_success_z', num2str(zScore), '.svg']);


% Spearman's correlation for Z-scores
[rho_ZScores, pval_ZScores] = corr(runCombinations', mean(zScores)', 'Type', 'Spearman');

% Display results for Z-scores
disp('Spearman correlation for Z-scores:');
disp('Correlation coefficients (rho):');
disp(rho_ZScores);
disp('P-values:');
disp(pval_ZScores);

% plotCreator(diff(zScores')', 1:4, subSplit);
% % figure;
% % plot(diff(zScores'), '.-r')
% % hold on
% % plot(mean(diff(zScores')'), 'k')
% % errorbar(1:4,mean(diff(zScores')'), std(diff(zScores')')/sqrt(length(zScores)), 'Color','black')
% make_pretty
% xlim([0.75,4.25])
% xlabel('Run Combination')
% ylabel('zScores');
% title(sprintf('Difference Average zScores Successive runs'));
% xticks(1:4)
% xticklabels({'2-3','3-4','4-5','5-6'})
% 
% % Save the plot as a PNG image
% % saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_zscore_diff.png');
% % saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_zscore_diff.svg');
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_diff_z', num2str(zScore), '.png']);
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_diff_z', num2str(zScore), '.svg']);

time = [1, 2, 3, 4, 5];
figure;
findSlopePts(activeVoxels, time)
make_pretty
xlim([0.75,5.25])
xlabel('Run Combination')
ylabel('Active Voxels');
title(sprintf('Line of Best Fit Active Voxels Successive runs'));
xticks(1:5)
xticklabels({'1-2','1-3','1-4','1-5','1-6'})
% Save the plot as a PNG image
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_slope_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_slope_z', num2str(zScore), '.svg']);

time = [1, 2, 3, 4, 5];
figure;
findSlopePts(zScores, time)
make_pretty
xlim([0.75,5.25])
xlabel('Run Combination')
ylabel('Z-Score');
title(sprintf('Line of Best Fit Zscore Successive runs'));
xticks(1:5)
xticklabels({'1-2','1-3','1-4','1-5','1-6'})
% Save the plot as a PNG image
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_slope_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_slope_z', num2str(zScore), '.svg']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


spineNames = {'C5','C6','C7','C8'};
for i = 1:4
    for j = 1:4

        if i ~= j && j > i
    
            [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(actVoxelsSeg4(:, i), actVoxelsSeg4(:, j), 10000, 0.05, 6);
            spineNames{i}, spineNames{j} 
            rejectNull

        end
    end

end


for i = 1:4
    for j = 1:4

        if i ~= j && j > i
    
            [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(zSeg4(:, i), zSeg4(:, j), 10000, 0.05, 6);
            spineNames{i}, spineNames{j} 
            rejectNull

        end
    end

end



for i = 1:4
    [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(actVoxelsSeg4(:, i), actVoxelsSeg6(:, i), 10000, 0.05, 8);
    spineNames{i}
    ci95, rejectNull

end
for i = 1:4
    [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(actVoxelsSeg4(:, i), actVoxelsSeg5(:, i), 10000, 0.05, 8);
    spineNames{i}
    ci95, rejectNull

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% per vertebra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spineStr = {'C1','C2','C3','C4','C5','C6','C7','C8'};
% time = [1, 2, 3, 4, 5];
% figure;
% findSlopePts(mean(actVoxSeg(:, :, 1:2), 3), time)
% hold on
% findSlopePts(mean(actVoxSeg(:, :, 7:8), 3), time)
% errorbar(time, mean(mean(actVoxSeg(:, :, 1:2), 3, 'omitnan')), std(mean(actVoxSeg(:, :, 1:2), 3), 'omitnan')/sqrt(size(actVoxSeg(:, :, 1:2), 1)), '.', 'Color', 'black', 'Marker', 'none')
% errorbar(time, mean(mean(actVoxSeg(:, :, 7:8), 3, 'omitnan')), std(mean(actVoxSeg(:, :, 7:8), 3), 'omitnan')/sqrt(size(actVoxSeg(:, :, 7:8), 1)), '.', 'Color', 'black', 'Marker', 'none')
% make_pretty
% xlim([0.75,5.25])
% xlabel('Run Combination')
% ylabel('Active Voxels');
% xticks(1:5)
% xticklabels({'1-2','1-3','1-4','1-5','1-6'})
% title(sprintf(['Line of Best Fit Active Voxels Successive runs ', spineStr{[1, 2, 7, 8]}]'));
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_vert2_slope_z', num2str(zScore), '.png']);
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_vert2_slope_z', num2str(zScore), '.svg']);
% 
% 
% % Define the run combinations (1 to 6)
% runCombinations = 1:5;
% 
% % Spearman's correlation for active voxels
% [rho_activeVoxels, pval_activeVoxels] = corr(runCombinations', mean(mean(actVoxSeg(:, :, 1:2), 3))', 'Type', 'Spearman');
% 
% % Display results for active voxels
% disp('Spearman correlation for active voxels c1-c2:');
% disp('Correlation coefficients (rho):');
% disp(rho_activeVoxels);
% disp('P-values:');
% disp(pval_activeVoxels);

% % % spineStr = {'C1','C2','C3','C4','C5','C6','C7','C8'};
% % time = [1, 2, 3, 4, 5];
% % figure;
% % findSlopePts(mean(zSeg(:, :, 1), 3, 'omitnan'), time)
% % hold on
% % findSlopePts(mean(zSeg(:, :, 7:8), 3, 'omitnan'), time)
% % make_pretty
% % xlim([0.75,5.25])
% % xlabel('Run Combination')
% % ylabel('Zscore');
% % xticks(1:5)
% % xticklabels({'1-2','1-3','1-4','1-5','1-6'})
% % title(sprintf(['Line of Best Fit Active Voxels Successive runs ', spineStr{[1, 7, 8]}]'));
% % saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_vert_slope_z', num2str(zScore), '.png']);
% % saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_vert_slope_z', num2str(zScore), '.svg']);
% % 
% % Individual tSNR
% %segTSNR (x,:,:) x is across runs
% %segTSNR (:,x,:) x is across subjects
% %segTSNR (:,:,x) x is the vertera ( C1 - C8)
% % segTSNR(:,1,8) is one run across all subjects in C8
% % segTSNR(1,:,8) is on subject across all runs in C8
% % make_pretty
% spineStr = {'C1','C2','C3','C4','C5','C6','C7','C8'};
% % UPPER LEVELS
% spineNumList = {1, 2, 3, 4, 5, 6, 7, 8, 1:2, 1:3, 1:4, 2:3, 2:4, 5:8, 6:8, 7:8};
% for i = 1:length(spineNumList)
%     spineNum = spineNumList{i};
% 
%     %ACTIVE VOXEL
%     % plotCreator(activeVoxels, 1:5, subSplit);
%     plotCreator(mean(actVoxSeg(:, :, spineNum), 3), 1:5, subSplit);
%     make_pretty
%     xlim([0.75,5.25])
%     xlabel('Run Combination')
%     ylabel('Active Voxels');
%     title(sprintf(['Average Active Successive runs Vert ', spineStr{spineNum}]));
%     xticks(1:5)
%     xticklabels({'1-2','1-3','1-4','1-5','1-6'})
%     saveas(gcf, ['D:\SBSN\Manuscript\plots\vert\Spine_voxel_success_z', num2str(zScore), '_', [spineStr{spineNum}], '.png']);
%     saveas(gcf, ['D:\SBSN\Manuscript\plots\vert\Spine_voxel_success_z', num2str(zScore), '_', [spineStr{spineNum}], '.svg']);
%     
%     plotCreator(mean(zSeg(:, :, spineNum), 3), 1:5, subSplit);
%     make_pretty
%     xlim([0.75,5.25])
%     xlabel('Run Combination')
%     ylabel('Zscore');
%     title(sprintf(['Average Zscore Successive runs Vert ', spineStr{spineNum}]));
%     xticks(1:5)
%     xticklabels({'1-2','1-3','1-4','1-5','1-6'})
%     saveas(gcf, ['D:\SBSN\Manuscript\plots\vert\Spine_zscore_success_z', num2str(zScore), '_', [spineStr{spineNum}], '.png']);
%     saveas(gcf, ['D:\SBSN\Manuscript\plots\vert\Spine_zscore_success_z', num2str(zScore), '_', [spineStr{spineNum}], '.svg']);
%     
%     time = [1, 2, 3, 4, 5];
%     figure;
%     findSlopePts(mean(actVoxSeg(:, :, spineNum), 3), time)
%     make_pretty
%     xlim([0.75,5.25])
%     xlabel('Run Combination')
%     ylabel('Active Voxels');
%     title(sprintf(['Line of Best Fit Active Voxels Successive runs ', spineStr{spineNum}]'));
%     xticks(1:5)
%     xticklabels({'1-2','1-3','1-4','1-5','1-6'})
%     saveas(gcf, ['D:\SBSN\Manuscript\plots\vert\Spine_voxel_slope_z', num2str(zScore), '_', [spineStr{spineNum}], '.png']);
%     saveas(gcf, ['D:\SBSN\Manuscript\plots\vert\Spine_voxel_slope_z', num2str(zScore), '_', [spineStr{spineNum}], '.svg']);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % LOWER LEVELS
% spineNum = 6:7;
% plotCreator(mean(actVoxSeg(:, :, spineNum), 3), 1:5, subSplit);
% make_pretty
% xlim([0.75,5.25])
% xlabel('Run Combination')
% ylabel('Active Voxels');
% title(sprintf(['Average Active Successive runs Vert ', spineStr{spineNum}]'));
% xticks(1:5)
% xticklabels({'1-2','1-3','1-4','1-5','1-6'})
% saveas(gcf, [' ', spineStr{spineNum}]Spine_voxel_success_z', num2str(zScore), '_', [spineStr{spineNum}], '.png']);
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_success_z', num2str(zScore), '_', [spineStr{spineNum}], '.svg']);
% 
% plotCreator(mean(zSeg(:, :, spineNum), 3), 1:5, subSplit);
% make_pretty
% xlim([0.75,5.25])
% xlabel('Run Combination')
% ylabel('Zscore');
% title(sprintf(['Average Zscore Successive runs Vert ', spineStr{spineNum}]));
% xticks(1:5)
% xticklabels({'1-2','1-3','1-4','1-5','1-6'})
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_success_z', num2str(zScore), '_', [spineStr{spineNum}], '.png']);
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_success_z', num2str(zScore), '_', [spineStr{spineNum}], '.svg']);
% 
% time = [1, 2, 3, 4, 5];
% findSlopePts(mean(actVoxSeg(:, :, spineNum), 3), time)
% make_pretty
% xlim([0.75,5.25])
% xlabel('Run Combination')
% ylabel('Active Voxels');
% title(sprintf(['Line of Best Fit Active Voxels Successive runs ', spineStr{spineNum}]));
% xticks(1:5)
% xticklabels({'1-2','1-3','1-4','1-5','1-6'})
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_slope_z', num2str(zScore), '_', [spineStr{spineNum}], '.png']);
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_slope_z', num2str(zScore), '_', [spineStr{spineNum}], '.svg']);
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_slope_z', num2str(zScore), '.png']);
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_slope_z', num2str(zScore), '.svg']);

% time = [1, 2, 3, 4, 5];
% findSlopePts(mean(zSeg(:, :, 2:4), 3), time)
% make_pretty
% xlim([0.75,5.25])
% xlabel('Run Combination')
% ylabel('Active Voxels');
% title(sprintf('Line of Best Fit Zscore Successive runs C2-C4'));
% xticks(1:5)
% xticklabels({'1-2','1-3','1-4','1-5','1-6'})
% Save the plot as a PNG image
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_slope_z', num2str(zScore), '.png']);
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_slope_z', num2str(zScore), '.svg']);


% time = [1, 2, 3, 4, 5];
% findSlopePts(mean(zSeg(:, :, 6:8), 3), time)
% make_pretty
% xlim([0.75,5.25])
% xlabel('Run Combination')
% ylabel('Z-Score');
% title(sprintf('Line of Best Fit Zscore Successive runs C6-C8'));
% xticks(1:5)
% xticklabels({'1-2','1-3','1-4','1-5','1-6'})
% Save the plot as a PNG image
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_slope_z', num2str(zScore), '.png']);
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_slope_z', num2str(zScore), '.svg']);
nRegions = 4;

figure;
hBar=barh([mean(actVoxelsSeg4)', mean(actVoxelsSeg6)']);
X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
hold on  %4 runs
hEB = errorbar([mean(actVoxelsSeg4)' mean(actVoxelsSeg6)'], X, [(std(actVoxelsSeg4)/sqrt(length(actVoxelsSeg4)))', (std(actVoxelsSeg6)/sqrt(length(actVoxelsSeg6)))'], 'horizontal', '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(1:subSplit,:), 1),1))/10;
Y = randVec + X(1:nRegions, 1).';              
Y2 = randVec + X(1:nRegions, 2).'; 
scatter(actVoxelsSeg4(1:subSplit,:), Y, 30, 'k','o','filled'); 
scatter(actVoxelsSeg6(1:subSplit,:), Y2, 30, 'k','o','filled'); 
% scatter(actVoxelsSeg4(1:subSplit,:), [randVec+X(1,3), randVec+X(2,3), randVec+X(3,3), randVec+X(4,3)], 30, 'k','o','filled'); 
% scatter(actVoxelsSeg6(1:subSplit,:), [randVec+X(1,5), randVec+X(2,5), randVec+X(3,5), randVec+X(4,5)], 30, 'k','o','filled'); 
randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(subSplit+1:end,:), 1),1))/10;
Y = randVec + X(1:nRegions, 1).';              
Y2 = randVec + X(1:nRegions, 2).';              

scatter(actVoxelsSeg4(subSplit+1:end,:), Y, 60, 'k','x'); 
scatter(actVoxelsSeg6(subSplit+1:end,:), Y2, 60, 'k','x'); 
% scatter(actVoxelsSeg4(subSplit+1:end,:), [randVec+X(1,3), randVec+X(2,3), randVec+X(3,3), randVec+X(4,3)], 60, 'k','x'); 
% scatter(actVoxelsSeg6(subSplit+1:end,:), [randVec+X(1,5), randVec+X(2,5), randVec+X(3,5), randVec+X(4,5)], 60, 'k','x'); 
set (gca,'YDir','reverse')
yticks(1:length(1:4)); yticklabels({'C5','C6','C7','C8'})
ylabel('Cervical level')
xlabel('Active Voxels');
title(sprintf('Average Active Voxel 4 vs 6 Runs combined'));
make_pretty

% % Save the plot as a PNG image
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_area46_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_area46_z', num2str(zScore), '.svg']);

% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_voxel_area.png');
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_voxel_area.svg');


figure;
hBar=barh([mean(zSeg4, 'omitnan')', mean(zSeg6, 'omitnan')']);
X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
hold on  %4 runs
hEB = errorbar([mean(zSeg4, 'omitnan')', mean(zSeg6, 'omitnan')'], X, [(std(zSeg4, 'omitnan')/sqrt(length(zSeg4)))',  (std(zSeg6, 'omitnan')/sqrt(length(zSeg6)))'], 'horizontal', '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
% errorbar(mean(actVoxelsSeg4), 1:subSplit, std(actVoxelsSeg4)/sqrt(length(actVoxelsSeg4)), 'horizontal', '.', 'Color','black')
randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(1:subSplit,:), 1),1))/10;
Y = randVec + X(1:nRegions, 1).';              
Y2 = randVec + X(1:nRegions, 2).';   

scatter(zSeg4(1:subSplit,:), Y, 30, 'k','o','filled'); 
scatter(zSeg6(1:subSplit,:), Y2, 30, 'k','o','filled'); 
randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(subSplit+1:end,:), 1),1))/10;
Y = randVec + X(1:nRegions, 1).';              
Y2 = randVec + X(1:nRegions, 2).'; 

scatter(zSeg4(subSplit+1:end,:), Y, 60, 'k','x'); 
scatter(zSeg6(subSplit+1:end,:), Y2, 60, 'k','x'); 

set (gca,'YDir','reverse')
yticks(1:length(1:4)); yticklabels({'C5','C6','C7','C8'})
ylabel('Cervical level')
xlabel('Z-score');
title(sprintf('Average Z-Score 4 vs 6 Runs combined'));
make_pretty
% 
% % Save the plot as a PNG image
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_zscore_area.png');
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_zscore_area.svg');
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_area46_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_area46_z', num2str(zScore), '.svg']);



% figure;
% hBar=barh([mean(actVoxelsSeg4)', mean(actVoxelsSeg6)']);
% X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
% hold on  %4 runs
% hEB = errorbar([mean(actVoxelsSeg4)', mean(actVoxelsSeg6)'], X, [(std(actVoxelsSeg4)/sqrt(length(actVoxelsSeg4)))',  (std(actVoxelsSeg6)/sqrt(length(actVoxelsSeg6)))'], 'horizontal', '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
% % randVec = (-1 + (1+1)*rand(7,1))/10;
% % scatter(actVoxelsSeg4, [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1)], 30, 'k','o','filled'); 
% % scatter(actVoxelsSeg6, [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2)], 30, 'k','o','filled');
% randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(1:subSplit,:), 1)))/10;
% scatter(actVoxelsSeg4(1:subSplit,:), [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1), randVec+X(5,1), randVec+X(6,1), randVec+X(7,1), randVec+X(8,1)], 30, 'k','o','filled'); 
% scatter(actVoxelsSeg6(1:subSplit,:), [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2), randVec+X(5,2), randVec+X(6,2), randVec+X(7,2), randVec+X(8,2)], 30, 'k','o','filled'); 
% randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(subSplit+1:end, :), 1)))/10;
% scatter(actVoxelsSeg4(subSplit+1:end,:), [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1), randVec+X(5,1), randVec+X(6,1), randVec+X(7,1), randVec+X(8,1)], 60, 'k','x'); 
% scatter(actVoxelsSeg6(subSplit+1:end,:), [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2), randVec+X(5,2), randVec+X(6,2), randVec+X(7,2), randVec+X(8,2)], 60, 'k','x'); 
% set (gca,'YDir','reverse')
% yticks(1:length(1:8)); yticklabels({'C1','C2','C3','C4','C5','C6','C7','C8'})
% ylabel('Cervical level')
% xlabel('Active Voxels');
% title(sprintf('Average Active Voxel 4 vs 6 Runs combined'));
% make_pretty
% 
% % Save the plot as a PNG image
% % saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_voxel_area.png');
% % saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_voxel_area.svg');
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_area_z', num2str(zScore), '.png']);
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_area_z', num2str(zScore), '.svg']);
% 
% 
% figure;
% hBar=barh([mean(zSeg4, 'omitnan')', mean(zSeg6, 'omitnan')']);
% X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
% hold on  %4 runs
% hEB = errorbar([mean(zSeg4, 'omitnan')', mean(zSeg6, 'omitnan')'], X, [(std(zSeg4, 'omitnan')/sqrt(length(zSeg4)))',  (std(zSeg6, 'omitnan')/sqrt(length(zSeg6)))'], 'horizontal', '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
% randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(1:subSplit,:), 1)))/10;
% scatter(zSeg4(1:subSplit,:), [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1), randVec+X(5,1), randVec+X(6,1), randVec+X(7,1), randVec+X(8,1)], 30, 'k','o','filled'); 
% scatter(zSeg6(1:subSplit,:), [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2), randVec+X(5,2), randVec+X(6,2), randVec+X(7,2), randVec+X(8,2)], 30, 'k','o','filled'); 
% randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(subSplit+1:end, :), 1)))/10;
% scatter(zSeg4(subSplit+1:end,:), [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1), randVec+X(5,1), randVec+X(6,1), randVec+X(7,1), randVec+X(8,1)], 60, 'k','x'); 
% scatter(zSeg6(subSplit+1:end,:), [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2), randVec+X(5,2), randVec+X(6,2), randVec+X(7,2), randVec+X(8,2)], 60, 'k','x'); 
% set (gca,'YDir','reverse')
% yticks(1:length(1:8)); yticklabels({'C1','C2','C3','C4','C5','C6','C7','C8'})
% ylabel('Cervical level')
% xlabel('Z-score');
% title(sprintf('Average Z-Score 4 vs 6 Runs combined'));
% make_pretty
% 
% % Save the plot as a PNG image
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_area_z', num2str(zScore), '.png']);
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_area_z', num2str(zScore), '.svg']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure;
% hBar=barh([mean(actVoxelsSeg4)', mean(actVoxelsSeg5)']);
% X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
% hold on  %4 runs
% hEB = errorbar([mean(actVoxelsSeg4)', mean(actVoxelsSeg5)'], X, [(std(actVoxelsSeg4)/sqrt(length(actVoxelsSeg4)))',  (std(actVoxelsSeg5)/sqrt(length(actVoxelsSeg5)))'], 'horizontal', '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
% % randVec = (-1 + (1+1)*rand(7,1))/10;
% % scatter(actVoxelsSeg4, [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1)], 30, 'k','o','filled'); 
% % scatter(actVoxelsSeg6, [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2)], 30, 'k','o','filled');
% randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(1:subSplit,:), 1)))/10;
% scatter(actVoxelsSeg4(1:subSplit,:), [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1), randVec+X(5,1), randVec+X(6,1), randVec+X(7,1), randVec+X(8,1)], 30, 'k','o','filled'); 
% scatter(actVoxelsSeg5(1:subSplit,:), [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2), randVec+X(5,2), randVec+X(6,2), randVec+X(7,2), randVec+X(8,2)], 30, 'k','o','filled'); 
% randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(subSplit+1:end, :), 1)))/10;
% scatter(actVoxelsSeg4(subSplit+1:end,:), [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1), randVec+X(5,1), randVec+X(6,1), randVec+X(7,1), randVec+X(8,1)], 60, 'k','x'); 
% scatter(actVoxelsSeg5(subSplit+1:end,:), [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2), randVec+X(5,2), randVec+X(6,2), randVec+X(7,2), randVec+X(8,2)], 60, 'k','x'); 
% set (gca,'YDir','reverse')
% yticks(1:length(1:8)); yticklabels({'C1','C2','C3','C4','C5','C6','C7','C8'})
% ylabel('Cervical level')
% xlabel('Active Voxels');
% title(sprintf('Average Active Voxel 4 vs 5 Runs combined'));
% ylim([4.5,8.5])
% make_pretty
% 
% % Save the plot as a PNG image
% % saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_voxel_area.png');
% % saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_voxel_area.svg');
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel45_area_z', num2str(zScore), '.png']);
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel45_area_z', num2str(zScore), '.svg']);
% 
% 
% figure;
% hBar=barh([mean(zSeg4, 'omitnan')', mean(zSeg5, 'omitnan')']);
% X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
% hold on  %4 runs
% hEB = errorbar([mean(zSeg4, 'omitnan')', mean(zSeg5, 'omitnan')'], X, [(std(zSeg4, 'omitnan')/sqrt(length(zSeg4)))',  (std(zSeg5, 'omitnan')/sqrt(length(zSeg5)))'], 'horizontal', '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
% randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(1:subSplit,:), 1)))/10;
% scatter(zSeg4(1:subSplit,:), [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1), randVec+X(5,1), randVec+X(6,1), randVec+X(7,1), randVec+X(8,1)], 30, 'k','o','filled'); 
% scatter(zSeg5(1:subSplit,:), [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2), randVec+X(5,2), randVec+X(6,2), randVec+X(7,2), randVec+X(8,2)], 30, 'k','o','filled'); 
% randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(subSplit+1:end, :), 1)))/10;
% scatter(zSeg4(subSplit+1:end,:), [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1), randVec+X(5,1), randVec+X(6,1), randVec+X(7,1), randVec+X(8,1)], 60, 'k','x'); 
% scatter(zSeg5(subSplit+1:end,:), [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2), randVec+X(5,2), randVec+X(6,2), randVec+X(7,2), randVec+X(8,2)], 60, 'k','x'); 
% set(gca,'YDir','reverse')
% yticks(1:length(1:8)); yticklabels({'C1','C2','C3','C4','C5','C6','C7','C8'})
% ylabel('Cervical level')
% xlabel('Z-score');
% ylim([4.5,8.5])
% % ylim()
% title(sprintf('Average Z-Score 4 vs 5 Runs combined'));
% make_pretty
% 
% % Save the plot as a PNG image
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore45_area_z', num2str(zScore), '.png']);
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore45_area_z', num2str(zScore), '.svg']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% 1v 2d  1 l 2 r

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% hBar = bar([mean(lrSeg4, 'omitnan'), mean(lrSeg6, 'omitnan')], 'FaceColor', 'flat');
% X = get(hBar,'XData').'+[hBar.XOffset];
% hold on  %4 runs
% hEB = errorbar(X, [mean(lrSeg4, 'omitnan'), mean(lrSeg6, 'omitnan')], [(std(lrSeg4, 'omitnan')/sqrt(length(lrSeg4))), (std(lrSeg6, 'omitnan')/sqrt(length(lrSeg6)))], 'vertical', '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
% randVec = (-1 + (1+1)*rand(2,1))/12;
% scatter([randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1)], [lrSeg4, lrSeg6],  30, 'k','o','filled'); 
% % scatter([randVec+X(2,1), randVec+X(2,2), randVec+X(2,3), randVec+X(2,4)], [dvSeg4, dvSeg6],   30, 'k','o','filled'); 
% xticks(1:subSplit); xticklabels({'L(pre)', 'R(pre)', 'L(post)', 'R(post)'})
% xlabel('Cervical Area');
% ylabel('Active Voxels');
% title(sprintf('Average Voxel Pre vs Post'));
% hBar.CData(2,:) = [0.8500 0.3250 0.0980];
% hBar.CData(4,:) = [0.8500 0.3250 0.0980];
% % set(hBar,'FaceColor','blue');
% make_pretty

figure;
hBar=bar([mean(lrSeg4, 'omitnan'); mean(dvSeg4, 'omitnan')]);
X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
hold on  %4 runs
hEB = errorbar(X, [mean(lrSeg4, 'omitnan')', mean(dvSeg4, 'omitnan')'], [(std(lrSeg4, 'omitnan')/sqrt(length(lrSeg4)))',  (std(dvSeg4, 'omitnan')/sqrt(length(dvSeg4)))'], 'vertical', '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
% randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(1:subSplit,:), 1)))/10;

randVec = (-1 + (1+1)*rand(2))/10;
% scatter([randVec+X(1,1), randVec+X(1,2)]', lrSeg4(1:subSplit,:),  30, 'k','o','filled'); 
% scatter([randVec+X(2,1), randVec+X(2,2)]', dvSeg4(1:subSplit,:),  30, 'k','o','filled'); 
% randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(subSplit+1:end, :), 1)))/10;

randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(subSplit+1:end, :), 1)))/10;
% scatter([randVec+X(1,1), randVec+X(2,1)]', lrSeg4(subSplit+1:end,:),  60, 'k','x'); 
% scatter([randVec+X(1,2), randVec+X(2,2)]', dvSeg4(subSplit+1:end,:),  60, 'k','x'); 
xticks(1:length(1:2)); xticklabels({'LR','VD'})
xlabel('Cervical Area')
ylabel('Active Voxels');
% title(sprintf('Average Active Voxel 4 vs 6 Runs combined'));
make_pretty

% figure;
% plot(lrSeg4(1:subSplit, :)', '.-b', 'Color',  [0.55, 0.55, 1], 'MarkerSize', 15)
% hold on
% plot(lrSeg4(subSplit+1:end, :)', 'x-b', 'MarkerSize', 15)
% 
% plot(mean(lrSeg4, 'omitnan'), 'k', 'LineWidth', 3)
% errorbar(1:2, mean(lrSeg4, 'omitnan'), std(lrSeg4, 'omitnan')/sqrt(length(lrSeg4)), 'Color', 'black', 'LineWidth', 2)
% xlim([0.75,2.25])
% make_pretty

% Save the plot as a PNG image
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_lrdv2_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_lrdv2_z', num2str(zScore), '.svg']);


[h,p,ci,stats] = ttest(lrSeg4(:, 1), lrSeg4(:, 2));
h,p,ci,stats
[ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(lrSeg4(:, 1), lrSeg4(:, 2), 10000, 0.05);
ci95, rejectNull
[h,p,ci,stats] = ttest(dvSeg4(:, 1), dvSeg4(:, 2));
h,p,ci,stats
[ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(dvSeg4(:, 1), dvSeg4(:, 2), 10000, 0.05);
ci95, rejectNull
% figure;
% hBar = bar([mean(dvSeg4, 'omitnan'), mean(dvSeg6, 'omitnan')], 'FaceColor', 'flat');
% X = get(hBar,'XData').'+[hBar.XOffset];
% hold on  %4 runs
% hEB = errorbar(X, [mean(dvSeg4, 'omitnan'), mean(dvSeg6, 'omitnan')], [(std(dvSeg4, 'omitnan')/sqrt(length(dvSeg4))) , (std(dvSeg6, 'omitnan')/sqrt(length(dvSeg6)))], 'vertical', '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
% randVec = (-1 + (1+1)*rand(2,1))/12;
% scatter([randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1)], [dvSeg4, dvSeg6],  30, 'k','o','filled'); 
% % scatter([randVec+X(2,1), randVec+X(2,2), randVec+X(2,3), randVec+X(2,4)], [dvSeg4, dvSeg6],   30, 'k','o','filled'); 
% xticks(1:4); xticklabels({'V(pre)', 'D(pre)', 'V(post)', 'D(post)'})
% xlabel('Cervical Area')
% ylabel('Active Voxels');
% title(sprintf('Average Voxel Pre vs Post'));
% hBar.CData(2,:) = [0.8500 0.3250 0.0980];
% hBar.CData(4,:) = [0.8500 0.3250 0.0980];
% make_pretty


figure;
hBar=bar([mean(zlrSeg4, 'omitnan')', mean(zdvSeg4, 'omitnan')']);
X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
hold on  %4 runs
hEB = errorbar(X, [mean(zlrSeg4, 'omitnan')', mean(zdvSeg4, 'omitnan')'], [(std(zlrSeg4, 'omitnan')/sqrt(length(zlrSeg4)))',  (std(zdvSeg4, 'omitnan')/sqrt(length(zdvSeg4)))'], 'vertical', '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
% randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(1:subSplit,:), 1)))/10;
% scatter([randVec+X(1,1), randVec+X(2,1)], zlrSeg4(1:subSplit,:),  30, 'k','o','filled'); 
% scatter([randVec+X(1,2), randVec+X(2,2)], zdvSeg4(1:subSplit,:),  30, 'k','o','filled'); 
randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(subSplit+1:end, :), 1)))/10;
% scatter([randVec+X(1,1), randVec+X(2,1)], zlrSeg4(subSplit+1:end,:),  60, 'k','x'); 
% scatter([randVec+X(1,2), randVec+X(2,2)], zdvSeg4(subSplit+1:end,:),  60, 'k','x');  
xticks(1:length(1:2)); xticklabels({'LR','VD'})
xlabel('Cervical Area')
ylabel('Z-Score');
% title(sprintf('Average Z-Score 4 vs 6 Runs combined'));
make_pretty

% Save the plot as a PNG image
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_zscore_lrdv.png');
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_zscore_lrdv.svg');
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_lrdv2_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_lrdv2_z', num2str(zScore), '.svg']);


% plotCreator(zdvSeg4, 1:2, subSplit);
% make_pretty
% xlim([0.75,2.25])
% ylabel('Active Voxels')
% xlabel('Group');
% title(sprintf('Mean Active Voxel'));
% % make_pretty
% xticklabels({'Task-Free','Task'})
% xticks(1:2)

% [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(zPre2, zPost2, 10000, 0.001);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is all the verts
figure;
hBar=barh([mean(actVoxelsSeg2)', mean(actVoxelsSeg3)', mean(actVoxelsSeg4)', mean(actVoxelsSeg5)', mean(actVoxelsSeg6)']);
X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
hold on  %4 runs
hEB = errorbar([mean(actVoxelsSeg2)', mean(actVoxelsSeg3)', mean(actVoxelsSeg4)', mean(actVoxelsSeg5)', mean(actVoxelsSeg6)'], X, [(std(actVoxelsSeg2)/sqrt(length(actVoxelsSeg2)))', (std(actVoxelsSeg3)/sqrt(length(actVoxelsSeg3)))', (std(actVoxelsSeg4)/sqrt(length(actVoxelsSeg4)))',  (std(actVoxelsSeg5)/sqrt(length(actVoxelsSeg5)))',  (std(actVoxelsSeg6)/sqrt(length(actVoxelsSeg6)))'], 'horizontal', '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
% randVec = (-1 + (1+1)*rand(7,1))/10;
% scatter(actVoxelsSeg4, [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1)], 30, 'k','o','filled'); 
% scatter(actVoxelsSeg6, [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2)], 30, 'k','o','filled');
randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(1:subSplit,:), 1)))/10;
% scatter(actVoxelsSeg2(1:subSplit,:), [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1)], 30, 'k','o','filled'); 
% scatter(actVoxelsSeg3(1:subSplit,:), [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2)], 30, 'k','o','filled'); 
% scatter(actVoxelsSeg4(1:subSplit,:), [randVec+X(1,3), randVec+X(2,3), randVec+X(3,3), randVec+X(4,3)], 30, 'k','o','filled'); 
% scatter(actVoxelsSeg5(1:subSplit,:), [randVec+X(1,4), randVec+X(2,4), randVec+X(3,4), randVec+X(4,4)], 30, 'k','o','filled'); 
% scatter(actVoxelsSeg6(1:subSplit,:), [randVec+X(1,5), randVec+X(2,5), randVec+X(3,5), randVec+X(4,5)], 30, 'k','o','filled'); 
randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(subSplit+1:end, :), 1)))/10;
% scatter(actVoxelsSeg2(subSplit+1:end,:), [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1)], 60, 'k','x'); 
% scatter(actVoxelsSeg3(subSplit+1:end,:), [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2)], 60, 'k','x'); 
% scatter(actVoxelsSeg4(subSplit+1:end,:), [randVec+X(1,3), randVec+X(2,3), randVec+X(3,3), randVec+X(4,3)], 60, 'k','x'); 
% scatter(actVoxelsSeg5(subSplit+1:end,:), [randVec+X(1,4), randVec+X(2,4), randVec+X(3,4), randVec+X(4,4)], 60, 'k','x'); 
% scatter(actVoxelsSeg6(subSplit+1:end,:), [randVec+X(1,5), randVec+X(2,5), randVec+X(3,5), randVec+X(4,5)], 60, 'k','x'); 
set (gca,'YDir','reverse')
yticks(1:length(1:4)); yticklabels({'C5','C6','C7','C8'})
ylabel('Cervical level')
xlabel('Active Voxels');
title(sprintf('Average Active Voxel 4 vs 6 Runs combined'));
make_pretty

% % Save the plot as a PNG image
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_all_area_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_all_area_z', num2str(zScore), '.svg']);

% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_voxel_area.png');
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_voxel_area.svg');


figure;
hBar=barh([mean(zSeg2, 'omitnan')', mean(zSeg3, 'omitnan')', mean(zSeg4, 'omitnan')', mean(zSeg5, 'omitnan')', mean(zSeg6, 'omitnan')']);
X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
hold on  %4 runs
hEB = errorbar([mean(zSeg2, 'omitnan')', mean(zSeg3, 'omitnan')', mean(zSeg4, 'omitnan')', mean(zSeg5, 'omitnan')', mean(zSeg6, 'omitnan')'], X, [(std(zSeg2, 'omitnan')/sqrt(length(zSeg2)))', (std(zSeg3, 'omitnan')/sqrt(length(zSeg3)))', (std(zSeg4, 'omitnan')/sqrt(length(zSeg4)))', (std(zSeg5, 'omitnan')/sqrt(length(zSeg5)))',  (std(zSeg6, 'omitnan')/sqrt(length(zSeg6)))'], 'horizontal', '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(1:subSplit,:), 1)))/10;
% scatter(zSeg2(1:subSplit,:), [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1)], 30, 'k','o','filled'); 
% scatter(zSeg3(1:subSplit,:), [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2)], 30, 'k','o','filled'); 
% scatter(zSeg4(1:subSplit,:), [randVec+X(1,3), randVec+X(2,3), randVec+X(3,3), randVec+X(4,3)], 30, 'k','o','filled'); 
% scatter(zSeg5(1:subSplit,:), [randVec+X(1,4), randVec+X(2,4), randVec+X(3,4), randVec+X(4,4)], 30, 'k','o','filled'); 
% scatter(zSeg6(1:subSplit,:), [randVec+X(1,5), randVec+X(2,5), randVec+X(3,5), randVec+X(4,5)], 30, 'k','o','filled'); 
randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(subSplit+1:end, :), 1)))/10;
% scatter(zSeg2(subSplit+1:end,:), [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1)], 60, 'k','x'); 
% scatter(zSeg3(subSplit+1:end,:), [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2)], 60, 'k','x'); 
% scatter(zSeg4(subSplit+1:end,:), [randVec+X(1,3), randVec+X(2,3), randVec+X(3,3), randVec+X(4,3)], 60, 'k','x'); 
% scatter(zSeg5(subSplit+1:end,:), [randVec+X(1,4), randVec+X(2,4), randVec+X(3,4), randVec+X(4,4)], 60, 'k','x'); 
% scatter(zSeg6(subSplit+1:end,:), [randVec+X(1,5), randVec+X(2,5), randVec+X(3,5), randVec+X(4,5)], 60, 'k','x'); 
set (gca,'YDir','reverse')
yticks(1:length(1:4)); yticklabels({'C5','C6','C7','C8'})
ylabel('Cervical level')
xlabel('Z-score');
title(sprintf('Average Z-Score 4 vs 6 Runs combined'));
make_pretty
% 
% % Save the plot as a PNG image
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_zscore_area.png');
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_zscore_area.svg');
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_all_area_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_all_area_z', num2str(zScore), '.svg']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % this is all the verts
% figure;
% hBar=barh([mean(actVoxelsSeg2)', mean(actVoxelsSeg3)', mean(actVoxelsSeg4)', mean(actVoxelsSeg5)', mean(actVoxelsSeg6)']);
% X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
% hold on  %4 runs
% hEB = errorbar([mean(actVoxelsSeg2)', mean(actVoxelsSeg3)', mean(actVoxelsSeg4)', mean(actVoxelsSeg5)', mean(actVoxelsSeg6)'], X, [(std(actVoxelsSeg2)/sqrt(length(actVoxelsSeg2)))', (std(actVoxelsSeg3)/sqrt(length(actVoxelsSeg3)))', (std(actVoxelsSeg4)/sqrt(length(actVoxelsSeg4)))',  (std(actVoxelsSeg5)/sqrt(length(actVoxelsSeg5)))',  (std(actVoxelsSeg6)/sqrt(length(actVoxelsSeg6)))'], 'horizontal', '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
% % randVec = (-1 + (1+1)*rand(7,1))/10;
% % scatter(actVoxelsSeg4, [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1)], 30, 'k','o','filled'); 
% % scatter(actVoxelsSeg6, [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2)], 30, 'k','o','filled');
% randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(1:subSplit,:), 1)))/10;
% scatter(actVoxelsSeg2(1:subSplit,:), [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1), randVec+X(5,1), randVec+X(6,1), randVec+X(7,1), randVec+X(8,1)], 30, 'k','o','filled'); 
% scatter(actVoxelsSeg3(1:subSplit,:), [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2), randVec+X(5,2), randVec+X(6,2), randVec+X(7,2), randVec+X(8,2)], 30, 'k','o','filled'); 
% scatter(actVoxelsSeg4(1:subSplit,:), [randVec+X(1,3), randVec+X(2,3), randVec+X(3,3), randVec+X(4,3), randVec+X(5,3), randVec+X(6,3), randVec+X(7,3), randVec+X(8,3)], 30, 'k','o','filled'); 
% scatter(actVoxelsSeg5(1:subSplit,:), [randVec+X(1,4), randVec+X(2,4), randVec+X(3,4), randVec+X(4,4), randVec+X(5,4), randVec+X(6,4), randVec+X(7,4), randVec+X(8,4)], 30, 'k','o','filled'); 
% scatter(actVoxelsSeg6(1:subSplit,:), [randVec+X(1,5), randVec+X(2,5), randVec+X(3,5), randVec+X(4,5), randVec+X(5,5), randVec+X(6,5), randVec+X(7,5), randVec+X(8,5)], 30, 'k','o','filled'); 
% randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(subSplit+1:end, :), 1)))/10;
% scatter(actVoxelsSeg2(subSplit+1:end,:), [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1), randVec+X(5,1), randVec+X(6,1), randVec+X(7,1), randVec+X(8,1)], 60, 'k','x'); 
% scatter(actVoxelsSeg3(subSplit+1:end,:), [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2), randVec+X(5,2), randVec+X(6,2), randVec+X(7,2), randVec+X(8,2)], 60, 'k','x'); 
% scatter(actVoxelsSeg4(subSplit+1:end,:), [randVec+X(1,3), randVec+X(2,3), randVec+X(3,3), randVec+X(4,3), randVec+X(5,3), randVec+X(6,3), randVec+X(7,3), randVec+X(8,3)], 60, 'k','x'); 
% scatter(actVoxelsSeg5(subSplit+1:end,:), [randVec+X(1,4), randVec+X(2,4), randVec+X(3,4), randVec+X(4,4), randVec+X(5,4), randVec+X(6,4), randVec+X(7,4), randVec+X(8,4)], 60, 'k','x'); 
% scatter(actVoxelsSeg6(subSplit+1:end,:), [randVec+X(1,5), randVec+X(2,5), randVec+X(3,5), randVec+X(4,5), randVec+X(5,5), randVec+X(6,5), randVec+X(7,5), randVec+X(8,5)], 60, 'k','x'); 
% set (gca,'YDir','reverse')
% yticks(1:length(1:8)); yticklabels({'C1','C2','C3','C4','C5','C6','C7','C8'})
% ylabel('Cervical level')
% xlabel('Active Voxels');
% title(sprintf('Average Active Voxel 4 vs 6 Runs combined'));
% make_pretty
% 
% % % Save the plot as a PNG image
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_all_area_z', num2str(zScore), '.png']);
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_voxel_all_area_z', num2str(zScore), '.svg']);
% 
% % saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_voxel_area.png');
% % saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_voxel_area.svg');
% 
% 
% figure;
% hBar=barh([mean(zSeg2, 'omitnan')', mean(zSeg3, 'omitnan')', mean(zSeg4, 'omitnan')', mean(zSeg5, 'omitnan')', mean(zSeg6, 'omitnan')']);
% X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
% hold on  %4 runs
% hEB = errorbar([mean(zSeg2, 'omitnan')', mean(zSeg3, 'omitnan')', mean(zSeg4, 'omitnan')', mean(zSeg5, 'omitnan')', mean(zSeg6, 'omitnan')'], X, [(std(zSeg2, 'omitnan')/sqrt(length(zSeg2)))', (std(zSeg3, 'omitnan')/sqrt(length(zSeg3)))', (std(zSeg4, 'omitnan')/sqrt(length(zSeg4)))', (std(zSeg5, 'omitnan')/sqrt(length(zSeg5)))',  (std(zSeg6, 'omitnan')/sqrt(length(zSeg6)))'], 'horizontal', '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
% randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(1:subSplit,:), 1)))/10;
% scatter(zSeg2(1:subSplit,:), [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1), randVec+X(5,1), randVec+X(6,1), randVec+X(7,1), randVec+X(8,1)], 30, 'k','o','filled'); 
% scatter(zSeg3(1:subSplit,:), [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2), randVec+X(5,2), randVec+X(6,2), randVec+X(7,2), randVec+X(8,2)], 30, 'k','o','filled'); 
% scatter(zSeg4(1:subSplit,:), [randVec+X(1,3), randVec+X(2,3), randVec+X(3,3), randVec+X(4,3), randVec+X(5,3), randVec+X(6,3), randVec+X(7,3), randVec+X(8,3)], 30, 'k','o','filled'); 
% scatter(zSeg5(1:subSplit,:), [randVec+X(1,4), randVec+X(2,4), randVec+X(3,4), randVec+X(4,4), randVec+X(5,4), randVec+X(6,4), randVec+X(7,4), randVec+X(8,4)], 30, 'k','o','filled'); 
% scatter(zSeg6(1:subSplit,:), [randVec+X(1,5), randVec+X(2,5), randVec+X(3,5), randVec+X(4,5), randVec+X(5,5), randVec+X(6,5), randVec+X(7,5), randVec+X(8,5)], 30, 'k','o','filled'); 
% randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(subSplit+1:end, :), 1)))/10;
% scatter(zSeg2(subSplit+1:end,:), [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1), randVec+X(5,1), randVec+X(6,1), randVec+X(7,1), randVec+X(8,1)], 60, 'k','x'); 
% scatter(zSeg3(subSplit+1:end,:), [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2), randVec+X(5,2), randVec+X(6,2), randVec+X(7,2), randVec+X(8,2)], 60, 'k','x'); 
% scatter(zSeg4(subSplit+1:end,:), [randVec+X(1,3), randVec+X(2,3), randVec+X(3,3), randVec+X(4,3), randVec+X(5,3), randVec+X(6,3), randVec+X(7,3), randVec+X(8,3)], 60, 'k','x'); 
% scatter(zSeg5(subSplit+1:end,:), [randVec+X(1,4), randVec+X(2,4), randVec+X(3,4), randVec+X(4,4), randVec+X(5,4), randVec+X(6,4), randVec+X(7,4), randVec+X(8,4)], 60, 'k','x'); 
% scatter(zSeg6(subSplit+1:end,:), [randVec+X(1,5), randVec+X(2,5), randVec+X(3,5), randVec+X(4,5), randVec+X(5,5), randVec+X(6,5), randVec+X(7,5), randVec+X(8,5)], 60, 'k','x'); 
% set (gca,'YDir','reverse')
% yticks(1:length(1:8)); yticklabels({'C1','C2','C3','C4','C5','C6','C7','C8'})
% ylabel('Cervical level')
% xlabel('Z-score');
% title(sprintf('Average Z-Score 4 vs 6 Runs combined'));
% make_pretty
% % 
% % % Save the plot as a PNG image
% % saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_zscore_area.png');
% % saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_zscore_area.svg');
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_all_area_z', num2str(zScore), '.png']);
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Spine_zscore_all_area_z', num2str(zScore), '.svg']);


%% SPINAL PROJECTION
% this is the spinal projection
% first sum the count of supra threshold voxels across all the slice we
% want into one transverse plane for each subject
% then I want to binraize it so count >= 1 is a 1
% then when I combine across all subjects I use a count to show a heatmap
% of how common a location of activation was acorss alll subjects

% SPINE

% addpath('/Users/pirondinilab/spinalcordtoolbox/cbiNifti');
addpath('D:\NHP_code\cbiNifti')

% varibales to set up before
subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_019','SBSN_H_007','SBSN_H_008','SBSN_H_010','SBSN_H_017','SBSN_H_018'}; 
% subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_008','SBSN_H_010'}; 
% subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_007','SBSN_H_008'}; 

zScore = 2.1;

copeFile = 'cope1.feat';
analysisFile = 'level_two_all_force_FLOB';

% gunzip('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii.gz');
[spineLevels, ~] = cbiReadNifti('F:\SMA_HOLDER\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii');
spineLevels(spineLevels < 5) = 0;
spineLevels(spineLevels > 8) = 0;
spineLevels(spineLevels > 0) = spineLevels(spineLevels > 0) - 4;

% gunzip('F:\SMA_HOLDER\MRI_data_upper_limb\Spine\template\PAM50_wm.nii.gz');
[wm, ~] = cbiReadNifti('F:\SMA_HOLDER\MRI_data_upper_limb\Spine\template\PAM50_wm.nii');
projectedWm = sum((wm).*(spineLevels>=1), 3);



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
        zPlane = {};

        %is dir and name contains gfeat
        if subjectFolder(folder).isdir && contains(subjectFolder(folder).name, analysisFile)

            disp(subjectFolder(folder).name)

            fileName = strsplit(subjectFolder(folder).name, '.');

            if ~exist(fullfile(direc, subjectFolder(folder).name, copeFile, 'thresh_zstat1.nii'))

                gunzip(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii.gz'));

            end

            [dataFile, ~] = cbiReadNifti(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii'));

            if i == 1
                dataFile = flip(dataFile, 1);
            end

            disp(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii'))

% this is the spinal projection
% first sum the count of supra threshold voxels across all the slice we
% want into one transverse plane for each subject
% then I want to binraize it so count >= 1 is a 1
% then when I combine across all subjects I use a count to show a heatmap
% of how common a location of activation was acorss alll subjects

            % this is where i sum the count of 1s that are in the plane
            % across these spinal level
            projectedPlane = sum((dataFile>=zScore).*(spineLevels>=1), 3);
%             projectedPlane(projectedPlane>=1) = 1;

            % NEW: z-score summary across selected slices for this run/subject
            % use max across slices so each voxel carries its strongest effect
%             x = (dataFile >= zScore) .* (spineLevels >= 1);
%             dataFile(~x) = 0;

            dataFile(spineLevels<1) = 0;
            for ii = 1:size(dataFile, 1)
                for jj = 1:size(dataFile, 1)

%                     if spineLevels(ii,jj,:)>=1

                    vals = squeeze(dataFile(ii,jj,:));   % 900x1 vector of values along the 3rd dim
                    zPlane{ii,jj} = vals(vals > zScore);  % keep only nonzero entries

%                     end
                end
            end            

%             % --- NEW: per-voxel z across selected levels (max across slices, robust masking) ---
%             inLevels = (spineLevels >= 1);
%             tmp = dataFile;
%             tmp(~inLevels) = -Inf;                  % exclude out-of-level slices
%             zPlaneMax = max(tmp, [], 3);            % [X Y]
%             zPlaneMax(~any(inLevels,3)) = 0;        % 0 where no in-level slice exis

            % number of active voxels
            allData{i, 1}{fileCounter, 1} = projectedPlane;
            allData{i, 1}{fileCounter, 2} = subjectFolder(folder).name;
%             allData{i, 1}{fileCounter, 3} = [numVoxelsSeperate; magSeperate];

            allData{i, 1}{fileCounter, 3} = zPlane;                % NEW: z map

            fileCounter = fileCounter + 1;
        end

    end
end


plane2 = zeros(141, 141);
plane3 = zeros(141, 141);
plane4 = zeros(141, 141);
plane5 = zeros(141, 141);
plane6 = zeros(141, 141);

% NEW: z-score sums and sample counts for mean z
Zsum2 = cell(141, 141); 
Zsum3 = cell(141, 141);
Zsum4 = cell(141, 141); 
Zsum5 = cell(141, 141);
Zsum6 = cell(141, 141);

% Zcnt2 = zeros(141, 141);
% Zcnt3 = zeros(141, 141);
% Zcnt4 = zeros(141, 141); 
% Zcnt5 = zeros(141, 141);
% Zcnt6 = zeros(141, 141);

for i = 1:length(allData)
% then when I combine across all subjects I use a count to show a heatmap
% of how common a location of activation was acorss alll subjects

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    plane2 = plane2 + allData{i,1}{1,1};

    plane3 = plane3 + allData{i,1}{2,1};

    plane4 = plane4 + allData{i,1}{3,1};

    plane5 = plane5 + allData{i,1}{4,1};

    plane6 = plane6 + allData{i,1}{5,1};


    for ii = 1:size(allData{i,1}{1,3}, 1)
        for jj = 1:size(allData{i,1}{1,3}, 1)
        
            Zsum2{ii, jj} = [Zsum2{ii, jj}; allData{i,1}{1,3}{ii, jj}];
        
            Zsum3{ii, jj} = [Zsum3{ii, jj}; allData{i,1}{2,3}{ii, jj}];
        
            Zsum4{ii, jj} = [Zsum4{ii, jj}; allData{i,1}{3,3}{ii, jj}];
        
            Zsum5{ii, jj} = [Zsum5{ii, jj}; allData{i,1}{4,3}{ii, jj}];
        
            Zsum6{ii, jj} = [Zsum6{ii, jj}; allData{i,1}{5,3}{ii, jj}];
        end
    end

end

% plane2 = plane2/max(plane2(:));
% plane3 = plane3/max(plane3(:));
% plane4 = plane4/max(plane4(:));
% plane5 = plane5/max(plane5(:));
% plane6 = plane6/max(plane6(:));
%%% NORMALIZE BY THE TOTAL COUNT HERE
projectedWm(projectedWm <= 20) = 0;
projectedWm(projectedWm > 20) = 1;

% projectedWm(63:79,68:77) < 1 = 100;
% Grab the submatrix
sub = projectedWm(63:79, 68:77);

% Make the change inside the submatrix
sub(sub < 1) = 100;

% Write it back (matrix size stays the same)
projectedWm(63:79, 68:77) = sub;

projectedWm(projectedWm <= 20) = 0;
projectedWm(projectedWm > 20) = 1;

% === PLOTTING: dots colored by individual z-scores (no background) ===
rows = 57:85;          % z-range
cols = 60:83;          % y-range
j    = 0.99;           % jitter in "cell units"

planes  = {plane2, plane3, plane4, plane5, plane6};
Zplanes = {Zsum2,  Zsum3,  Zsum4,  Zsum5,  Zsum6};
titles  = {'Runs 1-2','Runs 1-3','Runs 1-4','Runs 1-5','Runs 1-6'};

% dot size range (adjust to taste)
sMin = 3;      % size at z = zScore
sMax = 35;     % size at z = 5

% WM mask window (same transforms as plots)
maskSub = fliplr(rot90(projectedWm(rows, cols))) > 0;

f = figure('Units','normalized','OuterPosition',[0 0 1 1]);
t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

ax = gobjects(1,numel(planes));
for k = 1:numel(planes)
    ax(k) = nexttile(t, k); hold(ax(k), 'on');

    Msub   = fliplr(rot90(planes{k}(rows, cols)));
    Zcells = fliplr(rot90(Zplanes{k}(rows, cols)));

    [rIdx, cIdx, ~] = find(Msub);

    for m = 1:numel(rIdx)
        rr = rIdx(m);
        cc = cIdx(m);
        zvec = Zcells{rr, cc}(:);

        yJ = cc + (rand(numel(zvec),1) - 0.5)*2*j;
        zJ = rr + (rand(numel(zvec),1) - 0.5)*2*j;

        zCap = min(max(zvec, zScore), 5);
        s    = sMin + (sMax - sMin) * (zCap - zScore) / (5 - zScore);

        yCell = min(max(round(yJ), 1), size(maskSub, 2));
        zCell = min(max(round(zJ), 1), size(maskSub, 1));
        inMask = maskSub(sub2ind(size(maskSub), zCell, yCell));

        scatter(yJ(~inMask), zJ(~inMask), s(~inMask), zvec(~inMask), 'filled');
        h = scatter(yJ(inMask),  zJ(inMask),  s(inMask),  zvec(inMask),  'filled');
        set(h, 'MarkerEdgeColor','k', 'LineWidth', 1.25);
    end

    colormap(ax(k), winter);
    set(ax(k), 'YDir','reverse');
    axis(ax(k), 'equal', 'tight');
    xlim(ax(k), [0.5, size(Msub,2)+0.5]);
    ylim(ax(k), [0.5, size(Msub,1)+0.5]);
    clim(ax(k), [zScore, 5]);
    xlabel(ax(k), 'y');
    ylabel(ax(k), 'z');
    title(ax(k), titles{k});
    box(ax(k), 'on');
end

% Reserve the 6th tile for the colorbar and put the bar there
nexttile(t, 6); axis off
cb = colorbar(ax(1));              % tie colorbar to the first axes' CLim/colormap
cb.Layout.Tile = 6;                % place it in tile #6
cb.Label.String = 'Voxel z-score';

% % Optional: one shared colormap for all tiles
% colormap(t, winter);

% Set the figure units to 'normalized' and outerposition to fill the screen
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

saveas(f, ['D:\SBSN\Manuscript\plots\projection\Spine_transverse_gradient_z', num2str(zScore), '.png']);
saveas(f, ['D:\SBSN\Manuscript\plots\projection\Spine_transverse_gradient_z', num2str(zScore), '.svg']);



% ----- settings -----
rows = 57:85;          % z-range
cols = 60:83;          % y-range
j    = 0.99;           % jitter in "cell units"

planes = {plane2, plane3, plane4, plane5, plane6};
titles = {'Runs 1-2','Runs 1-3','Runs 1-4','Runs 1-5','Runs 1-6'};

% Precompute mask patch (same for all panels)
maskSub = fliplr(rot90(projectedWm(rows, cols)));
maskSub = maskSub > 0;   % logical

figure; clf;
for k = 1:numel(planes)
    subplot(2,3,k);

    % --- counts submatrix with identical transforms ---
    Msub = fliplr(rot90(planes{k}(rows, cols)));

    % --- nonzero cells and replication counts ---
    [zIdx, yIdx, vals] = find(Msub);
    rep = floor(vals);
    rep(~isfinite(rep)) = 0;
    rep = max(rep, 0);
    keep = rep > 0;
    zIdx = zIdx(keep);
    yIdx = yIdx(keep);
    rep   = rep(keep);

    % --- build repeated coordinates ---
    yAll = repelem(yIdx, rep);
    zAll = repelem(zIdx, rep);

    % --- jitter so overlapping points spread out ---
    rng(1 + k); % reproducible, different per panel
    yJ = yAll + (rand(size(yAll)) - 0.5)*2*j;
    zJ = zAll + (rand(size(zAll)) - 0.5)*2*j;

    % --- map jittered points to nearest cell and check mask membership ---
    yCell = round(yJ);  zCell = round(zJ);
    yCell = min(max(yCell, 1), size(maskSub, 2));
    zCell = min(max(zCell, 1), size(maskSub, 1));

    inMask = maskSub(sub2ind(size(maskSub), zCell, yCell));   % logical vector

    % --- plot (red = outside mask, black = inside mask) ---
    scatter(yJ(~inMask), zJ(~inMask), 10, 'filled', 'MarkerFaceColor', [.7 .7 .7]); hold on;
    scatter(yJ(inMask),  zJ(inMask),  10, 'filled', 'r');
    axis equal tight
    xlim([0.5, size(maskSub,2)+0.5]); 
    ylim([0.5, size(maskSub,1)+0.5]);
    set(gca,'YDir','reverse');  % image-like orientation
    xlabel('y'); ylabel('z');
    title(titles{k});
    box on

    % only show legend once to avoid clutter
    if k == numel(planes)
        legend({'outside mask','inside mask'}, 'Location','best');
    end
end

% optional: hide empty 6th subplot if you only have 5 panels
subplot(2,3,6); axis off

% saveas(gcf, ['D:\SBSN\Manuscript\plots\projection\Spine_transverse_z', num2str(zScore), '.png']);
% saveas(gcf, ['D:\SBSN\Manuscript\plots\projection\Spine_transverse_z', num2str(zScore), '.svg']);


%% DICE matching
clear all

%%% SPINE

% addpath('/Users/pirondinilab/spinalcordtoolbox/cbiNifti');
addpath('D:\NHP_code\cbiNifti')

% varibales to set up before
subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_019','SBSN_H_007','SBSN_H_008','SBSN_H_010','SBSN_H_017','SBSN_H_018'}; 
% subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_008','SBSN_H_010'}; 
% subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_007','SBSN_H_008'}; 
subSplit = 5;

zScore = 2.1;

copeFile = 'cope1.feat';
analysisFile = 'level_two_all_force_FLOB';
% gunzip('D:\SBSN\Data\Spine\template\PAM50_levels.nii.gz');
% [spineLevels, ~] = cbiReadNifti('D:\SBSN\Data\Spine\template\PAM50_levels.nii');
% spineLevels(spineLevels < 4) = 0;
% spineLevels(spineLevels > 7) = 0;
% spineLevels(spineLevels > 0) = spineLevels(spineLevels > 0) - 3;

% gunzip('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii.gz');
[spineLevels, ~] = cbiReadNifti('F:\SMA_HOLDER\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii');
spineLevels(spineLevels < 5) = 0;
spineLevels(spineLevels > 8) = 0;
spineLevels(spineLevels > 0) = spineLevels(spineLevels > 0) - 4;


% % gunzip('D:\SBSN\Data\Spine\template\PAM50_rl.nii.gz');
% [lrLevels, ~] = cbiReadNifti('D:\SBSN\Data\Spine\template\PAM50_all_rl.nii');
% % gunzip('D:\SBSN\Data\Spine\template\PAM50_dv.nii.gz');
% [dvLevels, ~] = cbiReadNifti('D:\SBSN\Data\Spine\template\PAM50_all_dv.nii');

allData = {};
for i = 1:length(subName)

    direc = fullfile('D:\SBSN\Data\Spine', subName{i}, 'func');

    subjectFolder = dir(direc);

    disp(subName{i})

    allData{i, 1} = {};
    allData{i, 2} = subName{i};

    fileCounter = 1;
    diceCell = {};
    for folder = 3:length(subjectFolder)

        %is dir and name contains gfeat
        if subjectFolder(folder).isdir && contains(subjectFolder(folder).name, analysisFile)

            disp(subjectFolder(folder).name)

            fileName = strsplit(subjectFolder(folder).name, '.');

            if ~exist(fullfile(direc, subjectFolder(folder).name, copeFile, 'thresh_zstat1.nii'))

                gunzip(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii.gz'));

            end

            [dataFile, ~] = cbiReadNifti(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii'));
            
            if i == 1
                dataFile = flip(dataFile, 1);
            end

            disp(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii'))


            diceCell{fileCounter} = dataFile; 

            fileCounter = fileCounter + 1;
        end

    end

    diceValues = [];
    diceValuesSeperate = [];
    for j = 1:length(diceCell)-1

        diceValues(j) = dice((diceCell{j} >= zScore).*(spineLevels>=1), (diceCell{j+1} >= zScore).*(spineLevels>=1));

        for k = 1:4
%             numVoxelsSeperate(j) = sum(sum(sum((dataFile>=zScore).*(spineLevels==j))));

            section = dice((diceCell{j} >= zScore).*(spineLevels==k), (diceCell{j+1} >= zScore).*(spineLevels==k));

            if isempty(section)

                diceValuesSeperate(j, k) = 0;

            else

                diceValuesSeperate(j, k) = section;

            end
            

        end

    end

    % number of active voxels
    allDice(i, :) = diceValues;
    allDiceSep(:, :, i) = diceValuesSeperate;

end

plotCreator(allDice, 1:4, subSplit);
make_pretty
xlim([0.75,4.25])
xlabel('Run Combination')
ylabel('DICE Score');
title(sprintf('Average Active Successive runs'));
xticks(1:4)
xticklabels({'2vs3','3vs4','4vs5','5vs6'})
saveas(gcf, ['D:\SBSN\Manuscript\plots\DICE\All_spine_area_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\DICE\All_spine_area_z', num2str(zScore), '.svg']);

spineLevels = {'C5','C6','C7','C8'};
for j = 1:4

    % run x area x subject
    allDiceSep2 = squeeze(allDiceSep(:, j, :))';
    
    
    plotCreator(allDiceSep2, 1:4, subSplit);
    make_pretty
    xlim([0.75,4.25])
    xlabel('Run Combination')
    ylabel('DICE Score');
    title(sprintf(['DICE Score Across Runs' spineLevels{j}]));
    xticks(1:4)
    xticklabels({'2vs3','3vs4','4vs5','5vs6'})
    saveas(gcf, ['D:\SBSN\Manuscript\plots\DICE\', spineLevels{j}, 'area_z', num2str(zScore), '.png']);
    saveas(gcf, ['D:\SBSN\Manuscript\plots\DICE\', spineLevels{j}, 'area_z', num2str(zScore), '.svg']);

end



%% CONTROL ANALYSIS
clear all

%%% SPINE

% addpath('/Users/pirondinilab/spinalcordtoolbox/cbiNifti');
addpath('D:\NHP_code\cbiNifti')

% varibales to set up before
subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_019','SBSN_H_007','SBSN_H_008','SBSN_H_010','SBSN_H_017','SBSN_H_018'}; 

% subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_007','SBSN_H_008'}; 
zScore = 0.0;
subSplit = 5;

%%THINGS TO ADD
% gunzip('D:\SBSN\Data\Spine\template\PAM50_levels.nii.gz');
% [spineLevels, ~] = cbiReadNifti('D:\SBSN\Data\Spine\template\PAM50_levels.nii');
% spineLevels(spineLevels < 4) = 0;
% spineLevels(spineLevels > 7) = 0;
% spineLevels(spineLevels > 0) = spineLevels(spineLevels > 0) - 3;

% gunzip('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii.gz');
[spineLevels, ~] = cbiReadNifti('F:\SMA_HOLDER\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii');
spineLevels(spineLevels < 5) = 0;
spineLevels(spineLevels > 8) = 0;
spineLevels(spineLevels > 0) = spineLevels(spineLevels > 0) - 4;

% gunzip('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii.gz');
% [spineLevels, ~] = cbiReadNifti('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii');
% % spineLevels(spineLevels < 1) = 0;
% spineLevels(spineLevels > 8) = 0;
% % spineLevels(spineLevels > 0) = spineLevels(spineLevels > 0) - 1;

allData = {};
for i = 1:length(subName)
    allData{i, 1} = {};
    allData{i, 2} = subName{i};

    for j = 0:6
        direc = fullfile('D:\SBSN\Data\Spine', subName{i}, 'func', ['func', num2str(j)]);
    
        subjectFolder = dir(direc);
    
        disp(subName{i})
    
        fileCounter = 1;
        for folder = 3:length(subjectFolder)
    
            %is dir and name contains gfeat
            if subjectFolder(folder).isdir && contains(subjectFolder(folder).name, 'level_one_force_FLOB') && (j > 0)
    
                disp(subjectFolder(folder).name)
    
                fileName = strsplit(subjectFolder(folder).name, '.');
    
                if ~exist(fullfile(direc, subjectFolder(folder).name, 'stats\zfstat1.nii')) 
    
                    gunzip(fullfile(direc, subjectFolder(folder).name, 'stats\zfstat1.nii.gz'));


                end
    
                [dataFile1, ~] = cbiReadNifti(fullfile(direc, subjectFolder(folder).name, 'stats\zfstat1.nii'));


                if i == 1
                    dataFile1 = flip(dataFile1, 1);
                end

                disp(fullfile(direc, subjectFolder(folder).name, 'stats\zfstat1.nii'))
    
                numVoxels = sum(sum(sum((dataFile1>=zScore).*(spineLevels>=1))));
                mag = dataFile1(spineLevels>=1);


                % number of active voxels
                allData{i, 1}{j+1, 1} = [numVoxels/sum(sum(sum(spineLevels>=1)))*100, mean(mag(mag>zScore)), std(mag(mag>zScore))];
                allData{i, 1}{j+1, 2} = subjectFolder(folder).name;
    
                fileCounter = fileCounter + 1;


            % this is for the rest run since the individual force levels
            % are not taken for the rest and only the overall force since
            % there should be no difference
            elseif subjectFolder(folder).isdir && contains(subjectFolder(folder).name, 'level_one_FLOB') && (j == 0)

                disp(subjectFolder(folder).name)
    
                fileName = strsplit(subjectFolder(folder).name, '.');
    
                if ~exist(fullfile(direc, subjectFolder(folder).name, 'stats\zstat1.nii')) 
    
                    gunzip(fullfile(direc, subjectFolder(folder).name, 'stats\zstat1.nii.gz'));

                end
    
                [dataFile1, ~] = cbiReadNifti(fullfile(direc, subjectFolder(folder).name, 'stats\zstat1.nii'));

                if i == 1
                    dataFile1 = flip(dataFile1, 1);

                end

                disp(fullfile(direc, subjectFolder(folder).name, 'stats\zstat1.nii'))
    
                numVoxels = sum(sum(sum((dataFile1>=zScore).*(spineLevels>=1))));
                mag = dataFile1(spineLevels>=1);

                % number of active voxels
                allData{i, 1}{j+1, 1} = [numVoxels/sum(sum(sum(spineLevels>=1)))*100, mean(mag(mag>zScore)), std(mag(mag>zScore))];
                allData{i, 1}{j+1, 2} = subjectFolder(folder).name;
    
                fileCounter = fileCounter + 1;

            end
    
        end
    end
end

for i = 1:length(allData)
    for j = 1:length(allData{i,1})

        activeVoxels(i,j) = allData{i,1}{j,1}(1);
        zScores(i,j) = allData{i,1}{j,1}(2);
         
    end
        
end


% [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(activeVoxels(:,1), mean(activeVoxels(:,2:end),2, 'omitnan'), 10000, 0.05);
% rejectNull
% [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(zScores(:,1), mean(zScores(:,2:end),2, 'omitnan'), 10000, 0.05);
% rejectNull

[h,p,ci,stats] = ttest(activeVoxels(:, 1), mean(activeVoxels(:, 2:end), 2, 'omitnan'));
h,p,ci,stats
[h,p,ci,stats] = ttest(zScores(:, 1), mean(zScores(:, 2:end), 2, 'omitnan'));
h,p,ci,stats

activeVoxels2 = [activeVoxels(:, 1), mean(activeVoxels(:, 2:end), 2, 'omitnan')];
zScores2 = [zScores(:, 1), mean(zScores(:, 2:end),2, 'omitnan')];

mean(activeVoxels2)
std(activeVoxels2)/sqrt(length(activeVoxels2))

mean(zScores2)
std(zScores2)/sqrt(length(zScores2))

plotCreator(activeVoxels2, 1:2, subSplit);
% figure;
% plot(activeVoxels2', '.-r')
% hold on
% plot(mean(activeVoxels2), 'k')
% errorbar(1:2, mean(activeVoxels2), std(activeVoxels2)/sqrt(length(activeVoxels2)), 'Color','black')
make_pretty
xlim([0.75,2.25])
ylabel('Active Voxels')
xlabel('Group');
title(sprintf('Mean Active Voxel'));
make_pretty
xticklabels({'Task-Free','Task'})
xticks(1:2)

% % Save the plot as a PNG image
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_voxel_control.png');
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_voxel_control.svg');
saveas(gcf, ['D:\SBSN\Manuscript\plots\control\Spine_voxel_control_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\control\Spine_voxel_control_z', num2str(zScore), '.svg']);


plotCreator(zScores2, 1:2, subSplit);
% figure;
% plot(zScores2', '.-r')
% hold on
% plot(mean(zScores2), 'k')
% errorbar(1:2, mean(zScores2), std(zScores2)/sqrt(length(zScores2)), 'Color','black')
make_pretty
xlim([0.75,2.25])
ylabel('Z-Score')
xlabel('Group');
title(sprintf('Mean Z-Score'));
make_pretty
xticklabels({'Task-Free','Task'})
xticks(1:2)

% Save the plot as a PNG image
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_zscore_control.png');
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_zscore_control.svg');
saveas(gcf, ['D:\SBSN\Manuscript\plots\control\Spine_zscore_control_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\control\Spine_zscore_control_z', num2str(zScore), '.svg']);

plotCreator(activeVoxels, 1:7, subSplit);
% figure;
% plot(activeVoxels', '.-r')
% hold on
% plot(mean(activeVoxels), 'k')
% errorbar(1:7, mean(activeVoxels), std(activeVoxels)/sqrt(length(activeVoxels)), 'Color','black')
make_pretty
xlim([0.75, 7.25])
ylabel('Active Voxels')
xlabel('Group');
title(sprintf('Mean Active Voxel'));
make_pretty
xticklabels({'Rest','1','2','3','4','5','6'})

% Save the plot as a PNG image
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_active_control_all.png');
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_active_control_all.svg');
saveas(gcf, ['D:\SBSN\Manuscript\plots\control\Spine_active_control_all_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\control\Spine_active_control_all_z', num2str(zScore), '.svg']);


plotCreator(zScores, 1:7, subSplit);
% figure;
% plot(zScores', '.-r')
% hold on
% plot(mean(zScores), 'k')
% errorbar(1:7, mean(zScores), std(zScores)/sqrt(length(zScores)), 'Color','black')
make_pretty
xlim([0.75, 7.25])
ylabel('Z-Score')
xlabel('Group');
title(sprintf('Mean Z-Score'));
make_pretty
xticklabels({'Rest','1','2','3','4','5','6'})

% % Save the plot as a PNG image
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_zscore_control_all.png');
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_zscore_control_all.svg');
saveas(gcf, ['D:\SBSN\Manuscript\plots\control\Spine_zscore_control_all_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\control\Spine_zscore_control_all_z', num2str(zScore), '.svg']);

%%
% %% Percent Signal Change
% clear all
% % varibales to set up before
% subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_007','SBSN_H_008','SBSN_H_010'}; 
% 
% copeFile = 'cope1.feat';
% 
% allData = {};
% for i = 1:length(subName)
% 
%     direc = fullfile('D:\SBSN\Data\Spine', subName{i}, 'func');
% 
%     subjectFolder = dir(direc);
% 
%     disp(subName{i})
% 
%     allData{i, 1} = {};
%     allData{i, 2} = subName{i};
% 
%     fileCounter = 1;
%     for folder = 3:length(subjectFolder)
% 
%         %is dir and name contains gfeat
%         if subjectFolder(folder).isdir && contains(subjectFolder(folder).name, 'level_two_FLOB')
% 
%             txtFile = importdata(fullfile(direc, subjectFolder(folder).name, copeFile, 'featquery+', 'report.txt'));
% %             fclose(fid);
% 
%             fileName = strsplit(subjectFolder(folder).name, '.');
% 
%             % number of active voxels
%             allData{i, 1}{fileCounter, 1} = abs(txtFile.data(4));
%             allData{i, 1}{fileCounter, 2} = subjectFolder(folder).name;
% 
%             fileCounter = fileCounter + 1;
%         end
% 
%     end
% end
% 
% for i = 1:length(allData)
%     for j = 1:length(allData{i,1})
% 
%         activeVoxels(i,j) = allData{i,1}{j,1};
%          
%     end
%         
% end
% 
% plotCreator(activeVoxels, 1:5);
% % figure;
% % plot(activeVoxels', '.-r')
% % hold on
% % plot(mean(activeVoxels), 'k')
% % errorbar(1:5, mean(activeVoxels), std(activeVoxels)/sqrt(length(activeVoxels)), 'Color','black')
% make_pretty
% xlim([0.75,5.25])
% ylabel('Percent Signal Change')
% xlabel('Group');
% title(sprintf('Signal Change'));
% xticks(1:5)
% xticklabels({'1-2','1-3','1-4','1-5','1-6'})
% 
% 
% % Save the plot as a PNG image
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_signal_success.png');
% saveas(gcf, 'D:\SBSN\Manuscript\plots\Spine_signal_success.svg');
% 
% % Define the run combinations (1 to 6)
% runCombinations = 1:5;
% 
% 
% % Spearman's correlation for Z-scores
% [rho_ZScores, pval_ZScores] = corr(runCombinations', mean(activeVoxels)', 'Type', 'Spearman');
% 
% % Display results for Z-scores
% disp('Spearman correlation for Z-scores:');
% disp('Correlation coefficients (rho):');
% disp(rho_ZScores);
% disp('P-values:');
% disp(pval_ZScores);

%% functions

function plotCreator(value, len, subSplit)
    % This function plots a sine wave given a frequency and duration.
    % Inputs:
    %   frequency - Frequency of the sine wave
    %   duration  - Duration (in seconds) of the sine wave to plot

    % Define the time vector [0.2 0.5 0.9 0.2]

%     Red (Full Red): [1, 0, 0]
%     Slightly Lighter Red: [1, 0.2, 0.2]
%     Lighter Red: [1, 0.4, 0.4]
%     Medium Light Red: [1, 0.6, 0.6]
%     Light Red: [1, 0.75, 0.75]
%     Very Light Red: [1, 0.85, 0.85]
%     Almost Pink (Lightest Red): [1, 0.9, 0.9]


%     Blue (Full Blue): [0, 0, 1]
%     Slightly Lighter Blue: [0.2, 0.2, 1]
%     Lighter Blue: [0.4, 0.4, 1]
%     Medium Light Blue: [0.6, 0.6, 1]
%     Light Blue: [0.75, 0.75, 1]
%     Very Light Blue: [0.85, 0.85, 1]
%     Almost White (Lightest Blue): [0.9, 0.9, 1]

%     opac = {[0.5, 0.5, 1]; [0.15, 0.15, 1]; [0.25, 0.25, 1]; [0.35, 0.35, 1]; [0.45, 0.45, 1]; [0.55, 0.55, 1]};

    figure;
    plot(value(1:subSplit, :)', '.-b', 'Color',  [0.55, 0.55, 1], 'MarkerSize', 15)
    hold on
    plot(value(subSplit+1:end, :)', 'x-b', 'MarkerSize', 15)

    plot(mean(value, 'omitnan'), 'k', 'LineWidth', 3)
    errorbar(len, mean(value, 'omitnan'), std(value, 'omitnan')/sqrt(length(value)), 'Color', 'black', 'LineWidth', 2)

end

function findSlopePts(data, time)

%     time = [1, 2, 3, 4, 5]; % Example time data
    [changepts,~] = findchangepts(mean(data), 'Statistic', 'linear', 'MaxNumChanges', 5);

    % Step 2: Plot data with detected changepoints
%     figure;
    plot(mean(data), 'o', 'DisplayName', 'Data');
    hold on;
    xline(time(changepts), '--r', 'DisplayName', 'Changepoint');
    xlabel('Time');
    ylabel('Response');
    title('Detected Changepoints');
    legend show;
    
    % Step 3: Perform piecewise linear regression for each segment
    fittedResponse = zeros(size(mean(data))); % Initialize fitted response array
    segments = [1, changepts, length(mean(data))]; % Define segments including changepoints
    
    slopes = []; % To store the slopes for each segment
    mA = mean(data);
    for i = 1:length(segments)-1
        % Get indices for the current segment
        idx = segments(i):segments(i+1);
        
        % Fit a linear model to the segment
        p = polyfit(time(idx), mA(idx), 1);
        
        % Get the fitted values for this segment
        fittedResponse(idx) = polyval(p, time(idx));
        
        % Store the slope for statistical testing later
        slopes = [slopes; p(1)]; % Store the slope (first coefficient of polyfit)
    end
    
    % Step 4: Plot the piecewise regression fit
    plot(time, fittedResponse, '-r', 'DisplayName', 'Piecewise Fit');
    legend show;
    
    % Step 5: Statistical Test for the Slopes (t-test for slopes close to zero)
    alpha = 0.05; % Significance level
    disp('Statistical Test Results for Each Segment:');
    for i = 1:length(slopes)
        fprintf('Segment %d: Slope = %.4f\n', i, slopes(i));
        
        % Perform a t-test to check if the slope is significantly different from zero
        % Null hypothesis: slope = 0 (indicating a plateau)
        [h, p_value] = ttest(slopes(i), 0, 'Alpha', alpha);
        
        if h == 0
            fprintf('Segment %d: The slope is NOT significantly different from zero (p = %.4f).\n', i, p_value);
        else
            fprintf('Segment %d: The slope is significantly different from zero (p = %.4f).\n', i, p_value);
        end
    end

end

% function findSlopePts(data, time, fitType)
%     % findSlopePts: Detects change points and fits either piecewise linear or quadratic segments
%     % Inputs:
%     %   data    - [nTrials x nTimepoints] matrix
%     %   time    - vector of time points
%     %   fitType - 'linear' or 'quadratic'
%     
%     if nargin < 3
%         fitType = 'linear'; % Default to linear if not specified
%     end
% 
%     % Step 1: Detect change points on the mean signal
%     mA = mean(data);
%     [changepts,~] = findchangepts(mA, 'Statistic', 'linear', 'MaxNumChanges', 5);
% 
%     % Step 2: Plot data with detected changepoints
%     plot(mA, 'o', 'DisplayName', 'Data');
%     hold on;
%     xline(time(changepts), '--r', 'DisplayName', 'Changepoint');
%     xlabel('Time');
%     ylabel('Response');
%     title('Detected Changepoints');
%     
%     % Step 3: Perform piecewise regression for each segment
%     fittedResponse = zeros(size(mA));
%     segments = [1, changepts, length(mA)];
%     slopes = [];
% 
%     for i = 1:length(segments)-1
%         idx = segments(i):segments(i+1);
%         
%         if strcmpi(fitType, 'linear')
%             p = polyfit(time(idx), mA(idx), 1);
%             slopes = [slopes; p(1)];
%         elseif strcmpi(fitType, 'quadratic')
%             p = polyfit(time(idx), mA(idx), 2);
%             % For simplicity, store only the linear term (slope component) from the quadratic
%             slopes = [slopes; p(2)];
%         else
%             error('fitType must be ''linear'' or ''quadratic''.');
%         end
%         
%         fittedResponse(idx) = polyval(p, time(idx));
%     end
% 
%     % Step 4: Plot the piecewise regression fit
%     plot(time, fittedResponse, '-r', 'DisplayName', [fitType ' Fit']);
%     legend show;
% 
%     % Step 5: Statistical Test for the Slopes (t-test)
%     alpha = 0.05;
%     disp('Statistical Test Results for Each Segment:');
%     for i = 1:length(slopes)
%         fprintf('Segment %d: Slope = %.4f\n', i, slopes(i));
%         [h, p_value] = ttest(slopes(i), 0, 'Alpha', alpha);
%         
%         if h == 0
%             fprintf('Segment %d: The slope is NOT significantly different from zero (p = %.4f).\n', i, p_value);
%         else
%             fprintf('Segment %d: The slope is significantly different from zero (p = %.4f).\n', i, p_value);
%         end
%     end
% end


function findSlopePtsSingle(data, time, i, subSplit)

%     time = [1, 2, 3, 4, 5]; % Example time data
    [changepts,~] = findchangepts(data, 'Statistic', 'linear', 'MaxNumChanges', 2);
    
    % Step 2: Plot data with detected changepoints
%     figure;

    if i > subSplit
        plot(data, 'xb', 'MarkerSize', 15);

    else

        plot(data, '.b', 'MarkerSize', 15);

    end
    hold on;
%     xline(time(changepts), '--r', 'DisplayName', 'Changepoint');
    xlabel('Time');
    ylabel('Response');
    title('Detected Changepoints');
%     legend show;
    
    % Step 3: Perform piecewise linear regression for each segment
    fittedResponse = zeros(size(data)); % Initialize fitted response array
    segments = [1, changepts, length(data)]; % Define segments including changepoints
    
    slopes = []; % To store the slopes for each segment
    mA = data;
    for i = 1:length(segments)-1
        % Get indices for the current segment
        idx = segments(i):segments(i+1);
        
        % Fit a linear model to the segment
        p = polyfit(time(idx), mA(idx), 1);
        
        % Get the fitted values for this segment
        fittedResponse(idx) = polyval(p, time(idx));
        
        % Store the slope for statistical testing later
        slopes = [slopes; p(1)]; % Store the slope (first coefficient of polyfit)
    end
    
    % Step 4: Plot the piecewise regression fit
    plot(time, fittedResponse, '-b');

end


function runRepeatedMeasuresANOVA(inputData)
    % This function performs a repeated measures ANOVA for 5, 6, or 7 time points
    % on the provided input data for 7 subjects. It also applies a linear contrast.
    % 
    % Inputs:
    %   inputData - A 7x5, 7x6, or 7x7 matrix where rows are subjects and columns are time points
    
    % Step 1: Validate the input data
    [numSubjects, numTimePoints] = size(inputData);
    
    if numSubjects ~= 7
        error('Input data must have 7 subjects (rows).');
    end
    
    if numTimePoints ~= 5 && numTimePoints ~= 6 && numTimePoints ~= 7
        error('Input data must have 5, 6, or 7 time points (columns).');
    end

    % Step 2: Define the subject identifiers
    subject = (1:numSubjects)';  % 7 subjects

    % Step 3: Create the data table
    % Convert the matrix to table format, with columns for each time point
    timePointNames = arrayfun(@(n) sprintf('Time%d', n), 1:numTimePoints, 'UniformOutput', false);
    dataTable = array2table(inputData, 'VariableNames', timePointNames);
    dataTable.Subject = subject;

    % Step 4: Fit the repeated measures model
    % The formula dynamically handles the number of time points
    formula = sprintf('Time1-Time%d ~ 1', numTimePoints);
    rm = fitrm(dataTable, formula, 'WithinDesign', 1:numTimePoints);

    % Step 5: Perform repeated measures ANOVA
    ranovaResults = ranova(rm);
    disp('Repeated Measures ANOVA Results:');
    disp(ranovaResults);

    % Step 6: Apply Linear Contrast
    % Define the linear contrast weights based on the number of time points
    if numTimePoints == 5
        contrastWeights = [-2 -1 0 1 2];
    elseif numTimePoints == 6
        contrastWeights = [-5 -3 -1 1 3 5];
    elseif numTimePoints == 7
        contrastWeights = [-3 -2 -1 0 1 2 3];
    end
    
    % Get the mean response for each time point
    means = mean(inputData, 'omitnan');  % Calculate the marginal means ignoring NaNs
    
    % Step 7: Calculate the linear contrast result
    % Apply the contrast weights to the means
    contrastValue = contrastWeights * means';
    
    % Step 8: Calculate standard error and test statistic
    % Standard error calculation: sqrt(sum of (contrastWeights^2 / n))
    variances = var(inputData, 'omitnan');  % Variances for each time point
    se = sqrt(sum((contrastWeights.^2 .* variances) / numSubjects));  % Standard error of contrast
    
    % Test statistic (t-value for linear contrast)
    tValue = contrastValue / se;
    
    % Convert t-value to F-value
    FValue = tValue^2;  % F-value is t-value squared for one degree of freedom

    % Degrees of freedom
    df1 = 1;  % One contrast, so numerator df is 1
    df2 = numSubjects - 1;  % Denominator df = number of subjects - 1

    % Calculate the p-value for the F-statistic
    pValue = 1 - fcdf(FValue, df1, df2);

    % Display the linear contrast result
    disp('Linear Contrast Results:');
    fprintf('Contrast Value: %.4f\n', contrastValue);
    fprintf('Standard Error: %.4f\n', se);
    fprintf('T-value: %.4f\n', tValue);
    fprintf('F-value: %.4f\n', FValue);
    fprintf('p-value: %.4f\n', pValue);
end