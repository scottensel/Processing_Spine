%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEFT HAND
% RIGHT HEMISPHERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stroke = [46, 43, 61, 60, 46, 86, 62, 62, 58, 74]
% health = [46, 46, 61, 61, 49, 71, 65, 61, 57, 70]
%% Frame Displacement
clear all

% addpath('/Users/pirondinilab/spinalcordtoolbox/cbiNifti');
addpath('D:\NHP_code\cbiNifti')

% varibales to set up before
% subName = {'SBSN_S_001','SBSN_S_002','SBSN_S_003','SBSN_S_004','SBSN_S_005', 'SBSN_S_006', 'SBSN_S_044', 'SBSN_S_055', 'SBSN_S_066','SBSN_S_077'}; 
subName = {'SBSN_H_007','SBSN_H_008','SBSN_H_010','SBSN_H_011','SBSN_H_013','SBSN_H_014','SBSN_H_015','SBSN_H_016','SBSN_H_017','SBSN_H_018',...
    'SBSN_S_001','SBSN_S_002','SBSN_S_003','SBSN_S_004','SBSN_S_005','SBSN_S_006','SBSN_S_007','SBSN_S_044','SBSN_S_008','SBSN_S_009'}; 

allData = {};
for i = 1:length(subName)

    allData{i, 1} = {};
    allData{i, 2} = subName{i};

    for j = 1:4
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
            allData{i, 1}{j, 1} = [mean(reshape(FD, [], 1),'omitnan'), std(reshape(FD, [], 1),'omitnan')];
            allData{i, 1}{j, 2} = [subName{i}, ' func', num2str(j)];
        catch
            continue
        end

    
    end
end


for i = 1:length(allData)

    frameD(i,:) = allData{i,1}{1,1}(1);

end

mean(mean(frameD)) 
std(frameD)/sqrt(length(frameD))

% means across runs
mean(frameD, 1) 
std(frameD, 1)/sqrt(length(mean(frameD)))

% %% TSNR
% clear all
% 
% % addpath('/Users/pirondinilab/spinalcordtoolbox/cbiNifti');
% addpath('D:\NHP_code\cbiNifti')
% 
% % varibales to set up before
% % subName = {'SBSN_S_001','SBSN_S_002','SBSN_S_003','SBSN_S_004','SBSN_S_005','SBSN_S_006','SBSN_S_044','SBSN_S_055','SBSN_S_066','SBSN_S_077'}; 
% subName = {'SBSN_H_007','SBSN_H_008','SBSN_H_010','SBSN_H_011','SBSN_H_013','SBSN_H_014','SBSN_H_015','SBSN_H_016','SBSN_H_017','SBSN_H_018',...
%     'SBSN_S_001','SBSN_S_002','SBSN_S_003','SBSN_S_004','SBSN_S_005','SBSN_S_006','SBSN_S_007','SBSN_S_044','SBSN_S_008','SBSN_S_009'}; 
% 
% % fmri_brain_moco_mean_tsnr_MNI152.nii.gz
% 
% % gunzip('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii.gz');
% [spineLevels, ~] = cbiReadNifti('F:\SMA_HOLDER\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii');
% spineLevels(spineLevels < 5) = 0;
% spineLevels(spineLevels > 8) = 0;
% spineLevels(spineLevels > 0) = spineLevels(spineLevels > 0) - 4;
% 
% allData = {};
% for i = 1:length(subName)
% 
%     allData{i, 1} = {};
%     allData{i, 2} = subName{i};
% 
%     for j = 1:4
%         direc = fullfile('D:\SBSN\Data\Spine', subName{i}, 'func', ['func', num2str(j)]);
%     
%         subjectFolder = dir(direc);
%     
%         disp(subName{i})
% 
%         if ~exist(fullfile(direc, 'fmri_spine_moco_mean_tsnr_PAM50.nii'))
% 
%             gunzip(fullfile(direc, 'fmri_spine_moco_mean_tsnr_PAM50.nii.gz'));
% 
%         end
% 
%         [dataFile, ~] = cbiReadNifti(fullfile(direc, 'fmri_spine_moco_mean_tsnr_PAM50.nii'));
% 
%         disp(fullfile(direc, 'fmri_spine_moco_mean_tsnr_PAM50.nii'))
% 
%         mag = dataFile(spineLevels>=1);
% 
%         % number of active voxels
%         allData{i, 1}{j, 1} = [mean(reshape(mag, [], 1),'omitnan'), std(reshape(mag, [], 1),'omitnan')];
%         allData{i, 1}{j, 2} = [subName{i}, ' func', num2str(j)];
% 
%     
%     end
% end
% 
% 
% tsnr = [];
% for i = 1:length(allData)
%     for j = 1:length(allData{i,1})
% 
%         tsnr(i,j) = allData{i,1}{j,1}(1);
% 
%     end
% 
% end
% 
% 
% plotCreator(tsnr, 1:4);
% make_pretty
% xlim([0.75,4.25])
% ylabel('TSNR')
% xlabel('Run Number')
% xticks([1:4])
% xticklabels([1:4])
% % xticklabels(strrep(subName, '_', '-'))
% title(sprintf('TSNR of Participant'));
% 
% % figure;
% % plot(tsnr', '.-r')
% % hold on
% % plot(mean(tsnr), 'k')
% % errorbar(1:7,mean(tsnr), std(tsnr)/sqrt(length(tsnr)), 'Color','black')
% % Save the plot as a PNG image
% saveas(gcf, 'D:\SBSN\stroke\Spine_tsnr.png');
% saveas(gcf, 'D:\SBSN\stroke\Spine_tsnr.svg');
% 
% % Define the run combinations (1 to 6)
% % runCombinations = 1:7;
% % 
% % % Spearman's correlation for active voxels
% % [rho_activeVoxels, pval_activeVoxels] = corr(runCombinations', mean(tsnr)', 'Type', 'Spearman');
% % 
% % % Display results for active voxels
% % disp('Spearman correlation for active voxels:');
% % disp('Correlation coefficients (rho):');
% % disp(rho_activeVoxels);
% % disp('P-values:');
% % disp(pval_activeVoxels);
% 
% 
% % runRepeatedMeasuresANOVA(tsnr)

%% template matching
% strokeAge = [46 43 61 60 46 86 62 62 58 74];
% controlAge = [46 46 61 61 49 71 65 61 57 70];
% 
% mean(strokeAge)
% std(strokeAge)
% 
% mean(controlAge)
% std(controlAge)
% close all
clear all

%%% SPINE

% addpath('/Users/pirondinilab/spinalcordtoolbox/cbiNifti');
addpath('D:\NHP_code\cbiNifti')

% varibales to set up before
subSplit = 10;
% subName = {'SBSN_S_001','SBSN_S_002','SBSN_S_003','SBSN_S_004','SBSN_S_005','SBSN_S_006','SBSN_S_044','SBSN_S_055','SBSN_S_066','SBSN_S_077'}; 
subName = {'SBSN_H_007','SBSN_H_008','SBSN_H_010','SBSN_H_011','SBSN_H_013','SBSN_H_014','SBSN_H_015','SBSN_H_016','SBSN_H_017','SBSN_H_018',...
    'SBSN_S_001','SBSN_S_002','SBSN_S_003','SBSN_S_004','SBSN_S_005','SBSN_S_006','SBSN_S_007','SBSN_S_044','SBSN_S_008','SBSN_S_009'}; 
% 8 11 13 14 17

% this is who gets flipped
% 7 8 10 11 14 15 16 17 18
% 3 5 8 9


fugl = [35, 23, 36, 29, 30, 23, 32, 32, 17, 34];
controls = ones(subSplit, 1)' * 66;

% subName = {'SBSN_S_001','SBSN_S_002','SBSN_S_003','SBSN_S_004','SBSN_S_005','SBSN_S_006','SBSN_S_044','SBSN_S_055','SBSN_S_066','SBSN_S_077'}; 
% subName = {'SBSN_S_044','SBSN_S_055','SBSN_S_066','SBSN_S_077'}; 

zScore = 2.1;

copeFile = 'cope1.feat';

% gunzip('F:\SMA_HOLDER\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii.gz');
[spineLevels, ~] = cbiReadNifti('F:\SMA_HOLDER\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii');
spineLevels(spineLevels < 5) = 0;
spineLevels(spineLevels > 8) = 0;
spineLevels(spineLevels > 0) = spineLevels(spineLevels > 0) - 4;

[spineLevelsAll, ~] = cbiReadNifti('F:\SMA_HOLDER\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii');
spineLevelsAll(spineLevelsAll > 8) = 0;
% spineLevels(spineLevels > 0) = spineLevels(spineLevels > 0);

% gunzip('D:\SBSN\Data\Spine\template\PAM50_rl.nii.gz');
[lrLevels, ~] = cbiReadNifti('D:\SBSN\Data\Spine\template\PAM50_rl.nii');
% gunzip('D:\SBSN\Data\Spine\template\PAM50_dv.nii.gz');
[dvLevels, ~] = cbiReadNifti('D:\SBSN\Data\Spine\template\PAM50_dv.nii');
% 
% % 4, WM left lateral corticospinal tract, PAM50_atlas_04.nii.gz
% % 5, WM right lateral corticospinal tract, PAM50_atlas_05.nii.gz
% % 8, WM left rubrospinal tract, PAM50_atlas_08.nii.gz
% % 9, WM right rubrospinal tract, PAM50_atlas_09.nii.gz
% % 16, WM left ventrolateral reticulospinal tract, PAM50_atlas_16.nii.gz
% % 17, WM right ventrolateral reticulospinal tract, PAM50_atlas_17.nii.gz
% % 20, WM left ventral reticulospinal tract, PAM50_atlas_20.nii.gz
% % 21, WM right ventral reticulospinal tract, PAM50_atlas_21.nii.gz
% % 22, WM left ventral corticospinal tract, PAM50_atlas_22.nii.gz
% % 23, WM right ventral corticospinal tract, PAM50_atlas_23.nii.gz
% % 26, WM left medial reticulospinal tract, PAM50_atlas_26.nii.gz
% % 27, WM right medial reticulospinal tract, PAM50_atlas_27.nii.gz
% spinalTractsNum = [4, 5, 8, 9, 16, 17, 20, 21, 22, 23, 26, 27];
% % spineCombo = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 7, 8];
% spineCombo = [1, 2, 3, 3, 4, 4, 4, 4, 5, 5, 4, 4];
% allTracts = [];
% for i = 1:length(spinalTractsNum)
% 
%     [spineTract, ~] = cbiReadNifti(['D:\SBSN\Data\Spine\template\atlas\', 'PAM50_atlas_', sprintf('%02d', spinalTractsNum(i)), '.nii']);
%     
% %     allTracts{i} = spineTract;
%     if length(allTracts) <= spineCombo(i)
% 
%         allTracts{spineCombo(i)} = spineTract;
% 
%     else
%         allTracts{spineCombo(i)} = allTracts{spineCombo(i)} + spineTract;
% 
%     end 
%     
% end


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
        if subjectFolder(folder).isdir && contains(subjectFolder(folder).name, 'level_two_all_force_FLOB1234.gfeat')

            disp(subjectFolder(folder).name)

            fileName = strsplit(subjectFolder(folder).name, '.');

            if ~exist(fullfile(direc, subjectFolder(folder).name, copeFile, 'thresh_zstat1.nii'))

                gunzip(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii.gz'));

            end

            [dataFile, ~] = cbiReadNifti(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii'));
           
            % these are the subjects who had lesions on the opposit side
%             if ismember(i, [2:7 9:subSplit subSplit+3 subSplit+5 subSplit+9 subSplit+10])
%                 dataFile = flip(dataFile, 1);
%             end

            % these are the subjects who had lesions on the opposit side
            if ismember(i, [1:4 6:subSplit subSplit+3 subSplit+5 subSplit+9 subSplit+10])
                dataFile = flip(dataFile, 1);
            end

            disp(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii'))

            numVoxels = sum(sum(sum((dataFile>=zScore).*(spineLevels>=1))));
            mag = dataFile(spineLevels>=1);

            numVoxelsSeperate = [];
            magSeperate = [];
            numVoxelsAllRegions = [];
            for j = 1:4
                numVoxelsSeperate(j) = sum(sum(sum((dataFile>=zScore).*(spineLevels==j))));
                numVoxelsSeperate(j) = numVoxelsSeperate(j)/sum(sum(sum(spineLevels==j)))*100;
                var = dataFile(spineLevels==j);
                magSeperate(j) = mean(var(var>zScore));  
            end

            for j = 1:2
                numVoxelsLR(j) = sum(sum(sum((dataFile>=zScore).*(lrLevels==j).*(spineLevelsAll>=6))));
                numVoxelsLR(j) = numVoxelsLR(j)/sum(sum(sum((lrLevels==j).*(spineLevelsAll>=6))))*100;

                numVoxelsDV(j) = sum(sum(sum((dataFile>=zScore).*(dvLevels==j).*(spineLevelsAll>=6))));
                numVoxelsDV(j) = numVoxelsDV(j)/sum(sum(sum((dvLevels==j).*(spineLevelsAll>=6))))*100;

                var = dataFile(lrLevels==j);
                magSeperateLR(j) = mean(var(var>zScore)); 

                var = dataFile(dvLevels==j);
                magSeperateDV(j) = mean(var(var>zScore)); 
                
            end
% .*(spineLevels>=3)
            for j = 1:4
                for k = 1:2
                    for k2 = 1:2
                        % L R V D
                        numVoxelsAllRegion(j, k, k2) = sum(sum(sum((dataFile>=zScore) .* (spineLevels==j) .* (lrLevels==k) .* (dvLevels==k2)  )));
                        numVoxelsAllRegion(j, k, k2) = numVoxelsAllRegion(j, k, k2)/sum(sum(sum((spineLevels==j) .* (lrLevels==k) .* (dvLevels==k2))))*100;

                    end
%                     var = dataFile(spineLevels==j);
%                     magSeperate(j) = mean(var(var>zScore));  
% 
%                     var = dataFile(lrLevels==j);
%                     magSeperateLR(k) = mean(var(var>zScore)); 
%     
%                     var = dataFile(dvLevels==j);
%                     magSeperateDV(k) = mean(var(var>zScore)); 
                end

            end

%             for j = 1:length(allTracts)
%                 index = (spineLevelsAll>=6).*(allTracts{j} > 0.5);
% 
% 
%                 numVoxelsTracts(j) = sum(sum(sum((dataFile>=zScore).*logical(index))));
%                 numVoxelsTracts(j) = numVoxelsTracts(j) / sum(sum(sum( logical(index) ))) *100;
%                 
%                 
%                 var = dataFile(logical(index));
%                 numVoxelsTractsZ(j) = mean(var(var>zScore)); 
%             end

            % number of active voxels
            allData{i, 1}{fileCounter, 1} = [numVoxels/sum(sum(sum(spineLevels>=1)))*100, mean(mag(mag>zScore)), std(mag(mag>zScore))];
            allData{i, 1}{fileCounter, 2} = subjectFolder(folder).name;
            allData{i, 1}{fileCounter, 3} = [numVoxelsSeperate; magSeperate];


            allData{i, 1}{fileCounter, 4} = [numVoxelsLR; magSeperateLR];         
            allData{i, 1}{fileCounter, 5} = [numVoxelsDV; magSeperateDV];

            allData{i, 1}{fileCounter, 6} = numVoxelsAllRegion;

%             allData{i, 1}{fileCounter, 7} = numVoxelsTracts;
%             allData{i, 1}{fileCounter, 8} = numVoxelsTractsZ;



            fileCounter = fileCounter + 1;
        end

    end
end


for i = 1:length(allData)
%     for j = 1:length(allData{i,1})

        activeVoxels(i) = allData{i,1}{1,1}(1);
        zScores(i) = allData{i,1}{1,1}(2);
        
        actVoxelsSegH(i, :) = allData{i,1}{1,3}(1,:);
        zSegH(i, :) = allData{i,1}{1,3}(2,:);

%         actVoxelsSeg6(i,:) = allData{i,1}{5,3}(1,:);
%         zSeg6(i,:) = allData{i,1}{5,3}(2,:);

        %%%%% 1v 2d 1 l 2 r
        lrSegH(i, :) = allData{i,1}{1,4}(1,:);
        zlrSegH(i, :) = allData{i,1}{1,4}(2,:);

        dvSegH(i, :) = allData{i,1}{1,5}(1,:);
        zdvSegH(i, :) = allData{i,1}{1,5}(2,:);

        quadSegH{i} = allData{i,1}{1,6};
        
%         tractsH(i, :) = allData{i,1}{1,7};
%         tractsZH(i, :) = allData{i,1}{1,8};

%         lrSeg6(i,:) = allData{i,1}{5,4}(1,:);
%         zlrSeg6(i,:) = allData{i,1}{5,4}(2,:);
% 
%         dvSeg6(i,:) = allData{i,1}{5,5}(1,:);
%         zdvSeg6(i,:) = allData{i,1}{5,5}(2,:);
%     end

end

activeVoxelsH = activeVoxels(1:subSplit);
activeVoxelsS = activeVoxels(subSplit+1:end);
zScoresH = zScores(1:subSplit);
zScoresS = zScores(subSplit+1:end);

actVoxelsSegS = actVoxelsSegH(subSplit+1:end,:);
zSegS = zSegH(subSplit+1:end,:);

actVoxelsSegH(subSplit+1:end,:) = [];
zSegH(subSplit+1:end,:) = [];

% LR
lrSegS = lrSegH(subSplit+1:end,:);
zlrSegS = zlrSegH(subSplit+1:end,:);

lrSegH(subSplit+1:end,:) = [];
zlrSegH(subSplit+1:end,:) = [];

% DV
dvSegS = dvSegH(subSplit+1:end,:);
zdvSegS = zdvSegH(subSplit+1:end,:);

dvSegH(subSplit+1:end,:) = [];
zdvSegH(subSplit+1:end,:) = [];


% QUAD
quadSegS = quadSegH(1, subSplit+1:end);
quadSegH = quadSegH(1, 1:subSplit);

% %tracts
% tractsS = tractsH(subSplit+1:end,:);
% tractsH(subSplit+1:end,:) = [];
% 
% tractsZS = tractsZH(subSplit+1:end,:);
% tractsZH(subSplit+1:end,:) = [];

figure;
plot(fugl, activeVoxelsS, '.', 'MarkerSize', 48, 'LineWidth', 1.5);
hold on
plot(controls, activeVoxelsH, '.', 'MarkerSize', 48, 'LineWidth', 1.5);

plot(mean(fugl), mean(activeVoxelsS), '.k', 'MarkerSize', 24, 'LineWidth', 1.5)
plot(mean(controls), mean(activeVoxelsH), '.k', 'MarkerSize', 24, 'LineWidth', 1.5)
errorbar( [mean(fugl), mean(controls)], [mean(activeVoxelsS, 'omitnan')', mean(activeVoxelsH, 'omitnan')'], [(std(activeVoxelsS, 'omitnan')/sqrt(length(activeVoxelsS)))',  (std(activeVoxelsH, 'omitnan')/sqrt(length(activeVoxelsH)))'], '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar

make_pretty
xlim([min(fugl)-0.25,66+0.25])
ylabel('% Active Voxels');
xlabel('Fugly-Meyer Score')
title('FM vs % Active Voxels')
% Save the plot as a PNG image
saveas(gcf, ['D:\SBSN\stroke\Spine_voxel_FM_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\stroke\Spine_voxel_FM_z', num2str(zScore), '.svg']);

figure;
plot(fugl, zScoresS, '.', 'MarkerSize', 48, 'LineWidth', 1.5);
hold on
plot(controls, zScoresH, '.', 'MarkerSize', 48, 'LineWidth', 1.5);
plot(mean(fugl), mean(zScoresS), '.k', 'MarkerSize', 24, 'LineWidth', 1.5)
plot(mean(controls), mean(zScoresH), '.k', 'MarkerSize', 24, 'LineWidth', 1.5)
errorbar( [mean(fugl), mean(controls)], [mean(zScoresS, 'omitnan')', mean(zScoresH, 'omitnan')'], [(std(zScoresS, 'omitnan')/sqrt(length(zScoresS)))',  (std(zScoresH, 'omitnan')/sqrt(length(zScoresH)))'], '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar

make_pretty
xlim([min(fugl)-0.25,66+0.25])
ylabel('Z-Score');
xlabel('Fugly-Meyer Score')
title('FM vs Z-score')
% Save the plot as a PNG image
saveas(gcf, ['D:\SBSN\stroke\Spine_zscore_FM_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\stroke\Spine_zscore_FM_z', num2str(zScore), '.svg']);

% figure;
% plot(fugl, actVoxelsSegS(:,1), '.', 'MarkerSize', 48, 'LineWidth', 1.5);
% figure;
% plot(fugl, actVoxelsSegS(:,2), '.', 'MarkerSize', 48, 'LineWidth', 1.5);
% figure;
% plot(fugl, actVoxelsSegS(:,3), '.', 'MarkerSize', 48, 'LineWidth', 1.5);
% figure;
% plot(fugl, actVoxelsSegS(:,4), '.', 'MarkerSize', 48, 'LineWidth', 1.5);



% % plotCreator(activeVoxelsS, 1:5);
% figure;
% plot(activeVoxelsS);
% hold on
% plot(activeVoxelsH);
% make_pretty
% % xlim([0.75,length(subName)+0.25])
% xlabel('Participant')
% xticklabels(strrep(subName, '_', '-'))
% ylabel('Active Voxels');
% title('Average % Active Successive runs');
% % xticks(1:5)
% % xticklabels({'1-2','1-3','1-4','1-5','1-6'})
% 
% % Save the plot as a PNG image
% saveas(gcf, ['D:\SBSN\stroke\Spine_voxel_success_z', num2str(zScore), '.png']);
% saveas(gcf, ['D:\SBSN\stroke\Spine_voxel_success_z', num2str(zScore), '.svg']);

% % plotCreator(zScores, 1:5);
% figure;
% plot(zScoresH);
% hold on
% plot(zScoresS);
% % figure;
% % plot(zScores', '.-r')
% % hold on
% % plot(mean(zScores), 'k')
% % errorbar(1:5,mean(zScores), std(zScores)/sqrt(length(zScores)), 'Color','black')
% make_pretty
% % xlim([0.75,length(subName)+0.25])
% xlabel('Participant')
% xticklabels(strrep(subName, '_', '-'))
% ylabel('Z-Score');
% title(sprintf('Average Zscore Successive runs'));
% % xticks(1:5)
% % xticklabels({'1-2','1-3','1-4','1-5','1-6'})
% 
% % Save the plot as a PNG image
% saveas(gcf, ['D:\SBSN\stroke\Spine_zscore_success_z', num2str(zScore), '.png']);
% saveas(gcf, ['D:\SBSN\stroke\Spine_zscore_success_z', num2str(zScore), '.svg']);


% Spearman's correlation for Z-scores
% [rho_ZScores, pval_ZScores] = corr(runCombinations', mean(zScores)', 'Type', 'Spearman');
% 
% % Display results for Z-scores
% disp('Spearman correlation for Z-scores:');
% disp('Correlation coefficients (rho):');
% disp(rho_ZScores);
% disp('P-values:');
% disp(pval_ZScores);
% 
% plotCreator(diff(zScores')', 1:4);
% figure;
% plot(diff(zScores'), '.-r')
% hold on
% plot(mean(diff(zScores')'), 'k')
% errorbar(1:4,mean(diff(zScores')'), std(diff(zScores')')/sqrt(length(zScores)), 'Color','black')
% make_pretty
% xlim([0.75,4.25])
% xlabel('Run Combination')
% ylabel('zScores');
% title(sprintf('Difference Average zScores Successive runs'));
% xticks(1:4)
% xticklabels({'2-3','3-4','4-5','5-6'})
% 
% % Save the plot as a PNG image
% saveas(gcf, 'D:\SBSN\stroke\Spine_zscore_diff.png');
% saveas(gcf, 'D:\SBSN\stroke\Spine_zscore_diff.svg');



% time = [1, 2, 3, 4, 5];
% findSlopePts(activeVoxels, time)
% make_pretty
% xlim([0.75,length(subName)+0.25])
% xlabel('Run Combination')
% ylabel('Active Voxels');
% title(sprintf('Difference Average Active Voxels Successive runs'));
% xticks(1:5)
% xticklabels({'1-2','1-3','1-4','1-5','1-6'})
% % Save the plot as a PNG image
% % saveas(gcf, 'D:\SBSN\stroke\Spine_voxel_slope.png');
% % saveas(gcf, 'D:\SBSN\stroke\Spine_voxel_slope.svg');
% 
% 
% 
% time = [1, 2, 3, 4, 5];
% findSlopePts(zScores, time)
% make_pretty
% xlim([0.75,length(subName)+0.25])
% xlabel('Run Combination')
% ylabel('Z-Score');
% title(sprintf('Average Zscore Successive runs'));
% xticks(1:5)
% xticklabels({'1-2','1-3','1-4','1-5','1-6'})
% Save the plot as a PNG image
% saveas(gcf, 'D:\SBSN\stroke\Spine_zscore_slope.png');
% saveas(gcf, 'D:\SBSN\stroke\Spine_zscore_slope.svg');

% 
% runRepeatedMeasuresANOVA(activeVoxels)
% runRepeatedMeasuresANOVA(zScores)
% actVoxelsSeg6 = actVoxelsSegS(3:end,:);
% zSeg6 = zSegS(3:end,:);
% 
% actVoxelsSegS(3:end,:) = [];
% zSegS(3:end,:) = [];

figure;
hBar = barh([mean(actVoxelsSegH)', mean(actVoxelsSegS)']);
X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
hold on  %4 runs
hEB = errorbar([mean(actVoxelsSegH)', mean(actVoxelsSegS)'], X, [(std(actVoxelsSegH)/sqrt(length(actVoxelsSegH)))',  (std(actVoxelsSegS)/sqrt(length(actVoxelsSegS)))'], 'horizontal','.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
randVec = (-1 + (1+1)*rand(subSplit,1))/10;
scatter(actVoxelsSegH, [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1)], 30, 'k','o','filled'); 
randVec = (-1 + (1+1)*rand(length(subName)-subSplit,1))/10;
scatter(actVoxelsSegS, [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2)], 30, 'k','o','filled'); 
% scatter(actVoxelsSeg6, [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2)], 20, 'k','o','filled'); 
set (gca,'YDir','reverse')
yticks(1:length(1:4)); yticklabels({'C5','C6','C7','C8'})
ylabel('Spinal Level')
xlabel('Active Voxels');
title(sprintf('Average Active Voxel 4 Runs Control vs Stroke'));
make_pretty

% Save the plot as a PNG image
saveas(gcf, ['D:\SBSN\stroke\Spine_voxel_area_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\stroke\Spine_voxel_area_z', num2str(zScore), '.svg']);

for idx = 1:4
    [H,P] = ttest2(actVoxelsSegH(:, idx), actVoxelsSegS(:, idx));
%     spineNames{idx}
    H, P
    
    
    [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(actVoxelsSegH(:, idx), actVoxelsSegS(:, idx), 10000, 0.05, 4);
%     spineNames{idx}
    ci95, rejectNull

end

figure;
hBar = barh([mean(zSegH, 'omitnan')', mean(zSegS, 'omitnan')']);
X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
hold on  %4 runs
hEB = errorbar([mean(zSegH, 'omitnan')', mean(zSegS, 'omitnan')'], X, [(std(zSegH, 'omitnan')/sqrt(length(zSegH)))',  (std(zSegS, 'omitnan')/sqrt(length(zSegS)))'], 'horizontal','.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
randVec = (-1 + (1+1)*rand(subSplit,1))/10;
scatter(zSegH, [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1)], 30, 'k','o','filled'); 
randVec = (-1 + (1+1)*rand(length(subName)-subSplit,1))/10;
scatter(zSegS, [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2)], 30, 'k','o','filled'); 
% scatter(zSeg6, [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2)], 20, 'k','o','filled'); 
set (gca,'YDir','reverse')
yticks(1:length(1:4)); yticklabels({'C5','C6','C7','C8'})
ylabel('Spinal Level')
xlabel('Z-score');
title(sprintf('Average Z-Score 4 Runs Control vs Stroke'));
make_pretty

% Save the plot as a PNG image
saveas(gcf, ['D:\SBSN\stroke\Spine_zscore_area_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\stroke\Spine_zscore_area_z', num2str(zScore), '.svg']);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% 1v 2d  1 l 2 r
figure;
hBar=bar([mean(lrSegH, 'omitnan'); mean(dvSegH, 'omitnan')]);
X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
hold on  %4 runs
hEB = errorbar(X, [mean(lrSegH, 'omitnan'); mean(dvSegH, 'omitnan')], [(std(lrSegH, 'omitnan')/sqrt(length(lrSegH))); (std(dvSegH, 'omitnan')/sqrt(length(dvSegH)))], 'vertical', '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
randVec = (-1 + (1+1)*rand(subSplit, 1))/10;
scatter([randVec+X(1,1), randVec+X(1,2), randVec+X(2,1), randVec+X(2,2)], [lrSegH, dvSegH],  30, 'k','o','filled'); 
% scatter([randVec+X(2,1), randVec+X(2,2), randVec+X(2,3), randVec+X(2,4)], [dvSegH, dvSeg6],   20, 'k','o','filled'); 
xticks(1:2); xticklabels({'L R','V D'})
xlabel('Cervical Area')
ylabel('Active Voxels');
title(sprintf('Controls'))
% title(sprintf('Average Active Voxel 4 vs 6 Runs combined'));
make_pretty

% Save the p[lot as a PNG image
saveas(gcf, ['D:\SBSN\stroke\Spine_voxel_lrdv_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\stroke\Spine_voxel_lrdv_z', num2str(zScore), '.svg']);



leftIdx = 1;
rightIdx = 2;
regionNames = {'Left vs Right'};


% Calculate Laterality Index for Pre and Post
LI_pre = (lrSegH(:,leftIdx) - lrSegH(:,rightIdx)) ./ ...
         (lrSegH(:,leftIdx) + lrSegH(:,rightIdx));
LI_post = (lrSegS(:,leftIdx) - lrSegS(:,rightIdx)) ./ ...
          (lrSegS(:,leftIdx) + lrSegS(:,rightIdx));

% Average and SEM
meanLI_pre = mean(LI_pre, 1);
semLI_pre = std(LI_pre, 0, 1) / sqrt(size(LI_pre, 1));
meanLI_post = mean(LI_post, 1);
semLI_post = std(LI_post, 0, 1) / sqrt(size(LI_post, 1));

% Plotting
figure;
hBar = bar([meanLI_pre; meanLI_post]', 'grouped');
hold on;

% Add error bars
ngroups = size([meanLI_pre; meanLI_post]', 1);
nbars = size([meanLI_pre; meanLI_post]', 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));


% Prepare the data for errorbar plotting
yvals = [meanLI_pre; meanLI_post];
errors = [semLI_pre; semLI_post];

% Then plot:
% for i = 1:nbars
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
% end

% X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
x=get(hBar,'XData').'+[hBar.XOffset];

randVec = (-1 + (1+1)*rand(length(subName)-subSplit, 1))/10;
errorbar(x, yvals, errors, 'k.', 'LineWidth', 1, 'Marker', 'none');

scatter([randVec+x(1), randVec+x(2)], [LI_pre, LI_post],  30, 'k','o','filled'); 


ylabel('Laterality Index (LI)');
xlabel('Brain Region');
xticks(1:length(regionNames));
xticklabels(regionNames);
legend({'Control', 'Stroke'}, 'Location', 'northeast');
title('Laterality Index');
ylim([-1, 1]);
make_pretty  % Optional custom styling


% % Save the figure
% Save the plot as a PNG image
% saveas(gcf, ['D:\SBSN\stroke\LI_area_z', num2str(zScore), '.png']);
% saveas(gcf, ['D:\SBSN\stroke\LI_area_z', num2str(zScore), '.svg']);
disp('control vs stroke')
% [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(LI_pre, LI_post, 10000, 0.05);
% 
ranksum(LI_pre, LI_post)


% Save the plot as a PNG image
saveas(gcf, ['D:\SBSN\stroke\LI_spine_area_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\stroke\LI_spine_area_z', num2str(zScore), '.svg']);



leftIdx = 1;
rightIdx = 2;
regionNames = {'ventral vs dorsal'};

% Calculate Laterality Index for Pre and Post
LI_pre = (dvSegH(:,leftIdx) - dvSegH(:,rightIdx)) ./ ...
         (dvSegH(:,leftIdx) + dvSegH(:,rightIdx));
LI_post = (dvSegS(:,leftIdx) - dvSegS(:,rightIdx)) ./ ...
          (dvSegS(:,leftIdx) + dvSegS(:,rightIdx));

% Average and SEM
meanLI_pre = mean(LI_pre, 1);
semLI_pre = std(LI_pre, 0, 1) / sqrt(size(LI_pre, 1));
meanLI_post = mean(LI_post, 1);
semLI_post = std(LI_post, 0, 1) / sqrt(size(LI_post, 1));

ranksum(LI_pre, LI_post)

% Plotting
figure;
hBar = bar([meanLI_pre; meanLI_post]', 'grouped');
hold on;

% Add error bars
ngroups = size([meanLI_pre; meanLI_post]', 1);
nbars = size([meanLI_pre; meanLI_post]', 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));


% Prepare the data for errorbar plotting
yvals = [meanLI_pre; meanLI_post];
errors = [semLI_pre; semLI_post];

% Then plot:
% for i = 1:nbars
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
% end

% X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
x=get(hBar,'XData').'+[hBar.XOffset];

randVec = (-1 + (1+1)*rand(length(subName)-subSplit, 1))/10;
errorbar(x, yvals, errors, 'k.', 'LineWidth', 1, 'Marker', 'none');

scatter([randVec+x(1), randVec+x(2)], [LI_pre, LI_post],  30, 'k','o','filled'); 


ylabel('Laterality Index (LI)');
xlabel('Brain Region');
xticks(1:length(regionNames));
xticklabels(regionNames);
legend({'Control', 'Stroke'}, 'Location', 'northeast');
title('Laterality Index Pre vs Post');
ylim([-1, 1]);
make_pretty  % Optional custom styling

% Save the plot as a PNG image
saveas(gcf, ['D:\SBSN\stroke\LI_spine_zscore_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\stroke\LI_spine_zscore_z', num2str(zScore), '.svg']);




%%%%% 1v 2d  1 l 2 r
figure;
hBar=bar([mean(zlrSegH, 'omitnan'); mean(zdvSegH, 'omitnan')]);
X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
hold on  %4 runs
hEB = errorbar(X, [mean(zlrSegH, 'omitnan'); mean(zdvSegH, 'omitnan')], [(std(zlrSegH, 'omitnan')/sqrt(length(zlrSegH))); (std(zdvSegH, 'omitnan')/sqrt(length(zdvSegH)))], 'vertical', '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
randVec = (-1 + (1+1)*rand(subSplit, 1))/10;
scatter([randVec+X(1,1), randVec+X(1,2), randVec+X(2,1), randVec+X(2,2)], [zlrSegH, zdvSegH],  30, 'k','o','filled'); 
% scatter([randVec+X(2,1), randVec+X(2,2), randVec+X(2,3), randVec+X(2,4)], [zdvSegH, zdvSeg6],   20, 'k','o','filled'); 
xticks(1:length(1:2)); xticklabels({'L R','V D'})
title(sprintf('Controls'))
xlabel('Cervical Area')
ylabel('Z-score');
% title(sprintf('Average Z-Score 4 vs 6 Runs combined'));
make_pretty

% Save the plot as a PNG image
saveas(gcf, ['D:\SBSN\stroke\Spine_zscore_lrdv_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\stroke\Spine_zscore_lrdv_z', num2str(zScore), '.svg']);


figure;
hBar=bar([mean(lrSegS, 'omitnan'); mean(dvSegS, 'omitnan')]);
X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
hold on  %4 runs
hEB = errorbar(X, [mean(lrSegS, 'omitnan'); mean(dvSegS, 'omitnan')], [(std(lrSegS, 'omitnan')/sqrt(length(lrSegS))); (std(dvSegS, 'omitnan')/sqrt(length(dvSegS)))], 'vertical', '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
randVec = (-1 + (1+1)*rand(length(subName)-subSplit, 1))/10;
scatter([randVec+X(1,1), randVec+X(1,2), randVec+X(2,1), randVec+X(2,2)], [lrSegS, dvSegS],  30, 'k','o','filled'); 
% scatter([randVec+X(2,1), randVec+X(2,2), randVec+X(2,3), randVec+X(2,4)], [dvSegS, dvSeg6],   20, 'k','o','filled'); 
xticks(1:2); xticklabels({'L R','V D'})
title(sprintf('Stroke'))
xlabel('Cervical Area')
ylabel('Active Voxels');
% title(sprintf('Average Active Voxel 4 vs 6 Runs combined'));
make_pretty

% Save the p[lot as a PNG image
saveas(gcf, ['D:\SBSN\stroke\SpineS_voxel_lrdv_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\stroke\SpineS_voxel_lrdv_z', num2str(zScore), '.svg']);

%%%%% 1v 2d  1 l 2 r
figure;
hBar=bar([mean(zlrSegS, 'omitnan'); mean(zdvSegS, 'omitnan')]);
X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
hold on  %4 runs
hEB = errorbar(X, [mean(zlrSegS, 'omitnan'); mean(zdvSegS, 'omitnan')], [(std(zlrSegS, 'omitnan')/sqrt(length(zlrSegS))); (std(zdvSegS, 'omitnan')/sqrt(length(zdvSegS)))], 'vertical', '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
randVec = (-1 + (1+1)*rand(length(subName)-subSplit, 1))/10;
scatter([randVec+X(1,1), randVec+X(1,2), randVec+X(2,1), randVec+X(2,2)], [zlrSegS, zdvSegS],  30, 'k','o','filled'); 
% scatter([randVec+X(2,1), randVec+X(2,2), randVec+X(2,3), randVec+X(2,4)], [zdvSegS, zdvSeg6],   20, 'k','o','filled'); 
xticks(1:length(1:2)); xticklabels({'L R','V D'})
title(sprintf('Stroke'))
xlabel('Cervical Area')
ylabel('Z-score');
% title(sprintf('Average Z-Score 4 vs 6 Runs combined'));
make_pretty

% Save the plot as a PNG image
saveas(gcf, ['D:\SBSN\stroke\SpineS_zscore_lrdv_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\stroke\SpineS_zscore_lrdv_z', num2str(zScore), '.svg']);

% [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(zPre2, zPost2, 10000, 0.001);



% 
% 
% Quad names and mapping to subplot positions
quadNames = {'Left Ventral', 'Left Dorsal'; 'Right Ventral', 'Right Dorsal'};
subplotPos = [1, 3; 2, 4];  % top-left=1, bottom-left=2, top-right=3, bottom-right=4

N_H = numel(quadSegH);
N_S = numel(quadSegS);
[L, R, C] = size(quadSegH{1});   % Expect 4x2x2
spinalLevels = {'C5','C6','C7','C8'};

figure;
for r = 1:R
    for c = 1:C
        % Extract values for this quadrant
        Hvals = zeros(N_H, L);
        for p = 1:N_H
            A = quadSegH{p};
            Hvals(p, :) = A(:, r, c).';
        end
        Svals = zeros(N_S, L);
        for p = 1:N_S
            A = quadSegS{p};
            Svals(p, :) = A(:, r, c).';
        end

        % Means and SEMs
        meanH = mean(Hvals, 1)';                   
        meanS = mean(Svals, 1)';                   
        semH  = std(Hvals,0,1)'/sqrt(N_H);      
        semS  = std(Svals,0,1)'/sqrt(N_S);      

        % Subplot in correct position
        subplot(2, 2, subplotPos(r, c));
        hold on;

        % Grouped horizontal bars
        hBar = barh([meanH, meanS]);
%         if numel(hBar) >= 2
%             hBar(1).FaceColor = [0.75 0.75 0.95]; % H group
%             hBar(2).FaceColor = [0.95 0.75 0.75]; % S group
%         end

        xHc = hBar(1).XEndPoints;  % centers for H bars
        xSc = hBar(2).XEndPoints;  % centers for S bars

        X = [xHc; xSc]';
        % Errorbars
        errorbar([meanH, meanS], X, [semH, semS], 'horizontal', '.', ...
                 'Color', 'k', 'Marker', 'none', 'LineWidth', 1);
        
        % Scatter points with jitter
        randVecH = (-1 + 2*rand(N_H, 1))/10;
        yH = [randVecH + X(1,1); randVecH + X(2,1); randVecH + X(3,1); randVecH + X(4,1)];
        scatter(Hvals(:), yH, 30, 'k', 'o', 'filled');

        randVecS = (-1 + 2*rand(N_S, 1))/10;
        yS = [randVecS + X(1,2); randVecS + X(2,2); randVecS + X(3,2); randVecS + X(4,2)];
        scatter(Svals(:), yS, 30, 'k', 'o', 'filled');

        % Axes, labels, title
        set(gca, 'YDir', 'reverse');
        yticks(1:L); yticklabels(spinalLevels);
        xlim([0 100])
        ylabel('Spinal Level'); xlabel('Active Voxels');
        title(quadNames{r,c});
        grid on; box off;
        hold off;

    end
end
make_pretty
sgtitle('Average Active Voxels by Level and Quadrant (H vs S) with Individual Subjects');

% Save the plot as a PNG image
saveas(gcf, ['D:\SBSN\stroke\Spine_voxel_quad_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\stroke\Spine_voxel_quad_z', num2str(zScore), '.svg']);


%% template matching
% strokeAge = [46 43 61 60 46 86 62 62 58 74];
% controlAge = [46 46 61 61 49 71 65 61 57 70];
% 
% mean(strokeAge)
% std(strokeAge)
% 
% mean(controlAge)
% std(controlAge)
% close all
clear all

%%% SPINE

% addpath('/Users/pirondinilab/spinalcordtoolbox/cbiNifti');
addpath('D:\NHP_code\cbiNifti')

% varibales to set up before
subSplit = 10;
% subName = {'SBSN_S_001','SBSN_S_002','SBSN_S_003','SBSN_S_004','SBSN_S_005','SBSN_S_006','SBSN_S_044','SBSN_S_055','SBSN_S_066','SBSN_S_077'}; 
subName = {'SBSN_H_007','SBSN_H_008','SBSN_H_010','SBSN_H_011','SBSN_H_013','SBSN_H_014','SBSN_H_015','SBSN_H_016','SBSN_H_017','SBSN_H_018',...
    'SBSN_S_001','SBSN_S_002','SBSN_S_003','SBSN_S_004','SBSN_S_005','SBSN_S_006','SBSN_S_007','SBSN_S_044','SBSN_S_008','SBSN_S_009'}; 
% 8 11 13 14 17

% this is who gets flipped
% 7 8 10 11 14 15 16 17 18
% 3 5 8 9


fugl = [35, 23, 36, 29, 30, 23, 32, 32, 17, 34];
controls = ones(subSplit, 1)' * 66;

% subName = {'SBSN_S_001','SBSN_S_002','SBSN_S_003','SBSN_S_004','SBSN_S_005','SBSN_S_006','SBSN_S_044','SBSN_S_055','SBSN_S_066','SBSN_S_077'}; 
% subName = {'SBSN_S_044','SBSN_S_055','SBSN_S_066','SBSN_S_077'}; 

zScore = 2.1;

copeFile = 'cope1.feat';

% gunzip('F:\SMA_HOLDER\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii.gz');
[spineLevels, ~] = cbiReadNifti('F:\SMA_HOLDER\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii');
spineLevels(spineLevels < 5) = 0;
spineLevels(spineLevels > 8) = 0;
spineLevels(spineLevels > 0) = spineLevels(spineLevels > 0) - 4;

[spineLevelsAll, ~] = cbiReadNifti('F:\SMA_HOLDER\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii');
spineLevelsAll(spineLevelsAll > 8) = 0;
% spineLevels(spineLevels > 0) = spineLevels(spineLevels > 0);

% gunzip('D:\SBSN\Data\Spine\template\PAM50_rl.nii.gz');
[lrLevels, ~] = cbiReadNifti('D:\SBSN\Data\Spine\template\PAM50_rl.nii');
% gunzip('D:\SBSN\Data\Spine\template\PAM50_dv.nii.gz');
[dvLevels, ~] = cbiReadNifti('D:\SBSN\Data\Spine\template\PAM50_dv.nii');

% 4, WM left lateral corticospinal tract, PAM50_atlas_04.nii.gz
% 5, WM right lateral corticospinal tract, PAM50_atlas_05.nii.gz
% 8, WM left rubrospinal tract, PAM50_atlas_08.nii.gz
% 9, WM right rubrospinal tract, PAM50_atlas_09.nii.gz
% 16, WM left ventrolateral reticulospinal tract, PAM50_atlas_16.nii.gz
% 17, WM right ventrolateral reticulospinal tract, PAM50_atlas_17.nii.gz
% 20, WM left ventral reticulospinal tract, PAM50_atlas_20.nii.gz
% 21, WM right ventral reticulospinal tract, PAM50_atlas_21.nii.gz
% 22, WM left ventral corticospinal tract, PAM50_atlas_22.nii.gz
% 23, WM right ventral corticospinal tract, PAM50_atlas_23.nii.gz
% 26, WM left medial reticulospinal tract, PAM50_atlas_26.nii.gz
% 27, WM right medial reticulospinal tract, PAM50_atlas_27.nii.gz
spinalTractsNum = [4, 5, 8, 9, 16, 17, 20, 21, 22, 23, 26, 27];
% spineCombo = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 7, 8];
spineCombo = [1, 2, 3, 3, 4, 4, 4, 4, 5, 5, 4, 4];
allTracts = [];
for i = 1:length(spinalTractsNum)

    [spineTract, ~] = cbiReadNifti(['D:\SBSN\Data\Spine\template\atlas\', 'PAM50_atlas_', sprintf('%02d', spinalTractsNum(i)), '.nii']);
    
%     allTracts{i} = spineTract;
    if length(allTracts) <= spineCombo(i)

        allTracts{spineCombo(i)} = spineTract;

    else
        allTracts{spineCombo(i)} = allTracts{spineCombo(i)} + spineTract;

    end 
    
end


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
        if subjectFolder(folder).isdir && contains(subjectFolder(folder).name, 'level_two_all_force_FLOB1234.gfeat')

            disp(subjectFolder(folder).name)

            fileName = strsplit(subjectFolder(folder).name, '.');

            if ~exist(fullfile(direc, subjectFolder(folder).name, copeFile, 'thresh_zstat1.nii'))

                gunzip(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii.gz'));

            end

            [dataFile, ~] = cbiReadNifti(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii'));
           
            % these are the subjects who had lesions on the opposit side
%             if ismember(i, [2:7 9:subSplit subSplit+3 subSplit+5 subSplit+9 subSplit+10])
%                 dataFile = flip(dataFile, 1);
%             end

            % these are the subjects who had lesions on the opposit side
            if ismember(i, [1:4 6:subSplit subSplit+3 subSplit+5 subSplit+9 subSplit+10])
                dataFile = flip(dataFile, 1);
            end

            disp(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii'))

            numVoxels = sum(sum(sum((dataFile>=zScore).*(spineLevels>=1))));
            mag = dataFile(spineLevels>=1);

            numVoxelsSeperate = [];
            magSeperate = [];
            numVoxelsAllRegions = [];
            for j = 1:4
                numVoxelsSeperate(j) = sum(sum(sum((dataFile>=zScore).*(spineLevels==j))));
                numVoxelsSeperate(j) = numVoxelsSeperate(j)/sum(sum(sum(spineLevels==j)))*100;
                var = dataFile(spineLevels==j);
                magSeperate(j) = mean(var(var>zScore));  
            end

            for j = 1:2
                numVoxelsLR(j) = sum(sum(sum((dataFile>=zScore).*(lrLevels==j).*(spineLevelsAll>=6))));
                numVoxelsLR(j) = numVoxelsLR(j)/sum(sum(sum((lrLevels==j).*(spineLevelsAll>=6))))*100;

                numVoxelsDV(j) = sum(sum(sum((dataFile>=zScore).*(dvLevels==j).*(spineLevelsAll>=6))));
                numVoxelsDV(j) = numVoxelsDV(j)/sum(sum(sum((dvLevels==j).*(spineLevelsAll>=6))))*100;

                var = dataFile(lrLevels==j);
                magSeperateLR(j) = mean(var(var>zScore)); 

                var = dataFile(dvLevels==j);
                magSeperateDV(j) = mean(var(var>zScore)); 
                
            end
% .*(spineLevels>=3)
            for j = 1:4
                for k = 1:2
                    for k2 = 1:2
                        % L R V D
                        numVoxelsAllRegion(j, k, k2) = sum(sum(sum((dataFile>=zScore) .* (spineLevels==j) .* (lrLevels==k) .* (dvLevels==k2)  )));
                        numVoxelsAllRegion(j, k, k2) = numVoxelsAllRegion(j, k, k2)/sum(sum(sum((spineLevels==j) .* (lrLevels==k) .* (dvLevels==k2))))*100;

                    end
%                     var = dataFile(spineLevels==j);
%                     magSeperate(j) = mean(var(var>zScore));  
% 
%                     var = dataFile(lrLevels==j);
%                     magSeperateLR(k) = mean(var(var>zScore)); 
%     
%                     var = dataFile(dvLevels==j);
%                     magSeperateDV(k) = mean(var(var>zScore)); 
                end

            end

            for j = 1:length(allTracts)
                index = (spineLevelsAll>=6).*(allTracts{j} > 0.5);


                numVoxelsTracts(j) = sum(sum(sum((dataFile>=zScore).*logical(index))));
                numVoxelsTracts(j) = numVoxelsTracts(j) / sum(sum(sum( logical(index) ))) *100;
                
                
                var = dataFile(logical(index));
                numVoxelsTractsZ(j) = mean(var(var>zScore)); 
            end

            % number of active voxels
            allData{i, 1}{fileCounter, 1} = [numVoxels/sum(sum(sum(spineLevels>=1)))*100, mean(mag(mag>zScore)), std(mag(mag>zScore))];
            allData{i, 1}{fileCounter, 2} = subjectFolder(folder).name;
            allData{i, 1}{fileCounter, 3} = [numVoxelsSeperate; magSeperate];


            allData{i, 1}{fileCounter, 4} = [numVoxelsLR; magSeperateLR];         
            allData{i, 1}{fileCounter, 5} = [numVoxelsDV; magSeperateDV];

            allData{i, 1}{fileCounter, 6} = numVoxelsAllRegion;

            allData{i, 1}{fileCounter, 7} = numVoxelsTracts;
            allData{i, 1}{fileCounter, 8} = numVoxelsTractsZ;



            fileCounter = fileCounter + 1;
        end

    end
end


for i = 1:length(allData)
%     for j = 1:length(allData{i,1})

        activeVoxels(i) = allData{i,1}{1,1}(1);
        zScores(i) = allData{i,1}{1,1}(2);
        
        actVoxelsSegH(i, :) = allData{i,1}{1,3}(1,:);
        zSegH(i, :) = allData{i,1}{1,3}(2,:);

%         actVoxelsSeg6(i,:) = allData{i,1}{5,3}(1,:);
%         zSeg6(i,:) = allData{i,1}{5,3}(2,:);

        %%%%% 1v 2d 1 l 2 r
        lrSegH(i, :) = allData{i,1}{1,4}(1,:);
        zlrSegH(i, :) = allData{i,1}{1,4}(2,:);

        dvSegH(i, :) = allData{i,1}{1,5}(1,:);
        zdvSegH(i, :) = allData{i,1}{1,5}(2,:);

        quadSegH{i} = allData{i,1}{1,6};
        
        tractsH(i, :) = allData{i,1}{1,7};
        tractsZH(i, :) = allData{i,1}{1,8};

%         lrSeg6(i,:) = allData{i,1}{5,4}(1,:);
%         zlrSeg6(i,:) = allData{i,1}{5,4}(2,:);
% 
%         dvSeg6(i,:) = allData{i,1}{5,5}(1,:);
%         zdvSeg6(i,:) = allData{i,1}{5,5}(2,:);
%     end

end

activeVoxelsH = activeVoxels(1:subSplit);
activeVoxelsS = activeVoxels(subSplit+1:end);
zScoresH = zScores(1:subSplit);
zScoresS = zScores(subSplit+1:end);

actVoxelsSegS = actVoxelsSegH(subSplit+1:end,:);
zSegS = zSegH(subSplit+1:end,:);

actVoxelsSegH(subSplit+1:end,:) = [];
zSegH(subSplit+1:end,:) = [];

% LR
lrSegS = lrSegH(subSplit+1:end,:);
zlrSegS = zlrSegH(subSplit+1:end,:);

lrSegH(subSplit+1:end,:) = [];
zlrSegH(subSplit+1:end,:) = [];

% DV
dvSegS = dvSegH(subSplit+1:end,:);
zdvSegS = zdvSegH(subSplit+1:end,:);

dvSegH(subSplit+1:end,:) = [];
zdvSegH(subSplit+1:end,:) = [];


% QUAD
quadSegS = quadSegH(1, subSplit+1:end);
quadSegH = quadSegH(1, 1:subSplit);

%tracts
tractsS = tractsH(subSplit+1:end,:);
tractsH(subSplit+1:end,:) = [];

tractsZS = tractsZH(subSplit+1:end,:);
tractsZH(subSplit+1:end,:) = [];


% tractNames = { ...
%     'L-lateral-corticospinal', ...   % 4
%     'R-lateral-corticospinal', ...   % 5
%     'L-rubrospinal', ...             % 8
%     'R-rubrospinal', ...             % 9
%     'L-ventrolateral-reticulospinal', ... % 16
%     'R-ventrolateral-reticulospinal', ... % 17
%     'L-ventral-reticulospinal', ...  % 20
%     'R-ventral-reticulospinal', ...  % 21
%     'L-ventral-corticospinal', ...   % 22
%     'R-ventral-corticospinal', ...   % 23
%     'L-medial-reticulospinal', ...   % 26
%     'R-medial-reticulospinal' ...    % 27
% };
% spineCombo = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 7, 8];
% MED IS BILATERAL FIBERS flexor muscles
% PONTINE IS IPSILATERAL extensor
% both act on interneurons
tractNames = { ...
    'L-lateral-corticospinal', ...   % 4
    'R-lateral-corticospinal', ...   % 5
    'L-rubrospinal', ...             % 8
    'R-rubrospinal', ...             % 9
    'L Med reticulospinal', ... % 16
    'R Med reticulospinal', ...     
    'L Pont reticulospinal', ... % 16
    'R Pont reticulospinal', ...
    'L-ventral-corticospinal', ...   % 22
    'R-ventral-corticospinal', ...   % 23
};
tractNames = { ...
    'L-lateral-corticospinal', ...   % 4
    'R-lateral-corticospinal', ...   % 5
    'rubrospinal', ...             % 8
    'reticulospinal', ...
    'L-ventral-corticospinal', ...   % 22
};
figure;
hBar = barh([mean(tractsH, 'omitnan')', mean(tractsS, 'omitnan')']);
X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
hold on  %4 runs
hEB = errorbar([mean(tractsH, 'omitnan')', mean(tractsS, 'omitnan')'], X, [(std(tractsH, 'omitnan')/sqrt(length(tractsH)))',  (std(tractsS, 'omitnan')/sqrt(length(tractsS)))'], 'horizontal','.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
randVec = (-1 + (1+1)*rand(subSplit,1))/10;
scatter(tractsH, [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1), randVec+X(5,1)], 30, 'k','o','filled'); 
randVec = (-1 + (1+1)*rand(length(subName)-subSplit,1))/10;
scatter(tractsS, [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2), randVec+X(5,2)], 30, 'k','o','filled'); 
% scatter(actVoxelsSeg6, [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2)], 20, 'k','o','filled'); 
set (gca,'YDir','reverse')
yticks(1:length(1:12)); yticklabels(tractNames)
ylabel('Tract')
xlabel('Active Voxels');
title(sprintf('Average Active Voxel 4 Runs Control vs Stroke'));
make_pretty

% Save the plot as a PNG image
saveas(gcf, ['D:\SBSN\stroke\Spine_voxel_tracts_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\stroke\Spine_voxel_tracts_z', num2str(zScore), '.svg']);

disp('control vs stroke')
for i = 1:length(tractNames)
    [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(tractsH(:, i), tractsS(:, i), 10000, 0.05);
    tractNames{i}
    ci95, rejectNull

end
disp('control vs stroke')
for i = 1:length(tractNames)
%     [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(tractsZH(:, i), tractsZS(:, i), 10000, 0.05);
%     ci95, rejectNull
    try
        [P,H] = ranksum(tractsH(:, i), tractsS(:, i),'tail', 'both') ;
    catch
        tractNames{i}
    end
    P, H
    tractNames{i}

end

% figure;
% hBar = barh([mean(tractsZH, 'omitnan')', mean(tractsZS, 'omitnan')']);
% X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
% hold on  %4 runs
% hEB = errorbar([mean(tractsZH, 'omitnan')', mean(tractsZS, 'omitnan')'], X, [(std(tractsZH, 'omitnan')/sqrt(length(tractsZH)))',  (std(tractsZS, 'omitnan')/sqrt(length(tractsZS)))'], 'horizontal','.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
% randVec = (-1 + (1+1)*rand(subSplit,1))/10;
% scatter(tractsZH, [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1), randVec+X(5,1), randVec+X(6,1), randVec+X(7,1), randVec+X(8,1), randVec+X(9,1), randVec+X(10,1)], 30, 'k','o','filled'); 
% randVec = (-1 + (1+1)*rand(length(subName)-subSplit,1))/10;
% scatter(tractsZS, [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2), randVec+X(5,2), randVec+X(6,2), randVec+X(7,2), randVec+X(8,2), randVec+X(9,2), randVec+X(10,2)], 30, 'k','o','filled'); 
% % scatter(actVoxelsSeg6, [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2)], 20, 'k','o','filled'); 
% set (gca,'YDir','reverse')
% yticks(1:length(1:12)); yticklabels(tractNames)
% ylabel('Tract')
% xlabel('Zscore');
% title(sprintf('Average Active Voxel 4 Runs Control vs Stroke'));
% make_pretty
% 
% % Save the plot as a PNG image
% saveas(gcf, ['D:\SBSN\stroke\Spine_zscore_tracts_z', num2str(zScore), '.png']);
% saveas(gcf, ['D:\SBSN\stroke\Spine_zscore_tracts_z', num2str(zScore), '.svg']);
% 
% disp('control vs stroke')
% for i = 1:length(tractNames)
% %     [ci95, rejectNull, diffSampMeans] = bootstrapCompMeans(tractsZH(:, i), tractsZS(:, i), 10000, 0.05);
% %     ci95, rejectNull
%     try
%         [P,H] = ranksum(tractsZH(:, i), tractsZS(:, i),'tail', 'both') ;
%     catch
%         tractNames{i}
%     end
%     P, H
%     tractNames{i}
% 
% end


%% SPINAL PROJECTION
% this is the spinal projection
% first sum the count of supra threshold voxels across all the slice we
% want into one transverse plane for each subject
% then I want to binraize it so count >= 1 is a 1
% then when I combine across all subjects I use a count to show a heatmap
% of how common a location of activation was acorss alll subjects

%%% SPINE

% addpath('/Users/pirondinilab/spinalcordtoolbox/cbiNifti');
addpath('D:\NHP_code\cbiNifti')

% variables to set up before
subSplit = 10;
% subName = {'SBSN_S_001','SBSN_S_002','SBSN_S_003','SBSN_S_004','SBSN_S_005','SBSN_S_006','SBSN_S_044','SBSN_S_055','SBSN_S_066','SBSN_S_077'}; 
subName = {'SBSN_H_007','SBSN_H_008','SBSN_H_010','SBSN_H_011','SBSN_H_013','SBSN_H_014','SBSN_H_015','SBSN_H_016','SBSN_H_017','SBSN_H_018',...
    'SBSN_S_001','SBSN_S_002','SBSN_S_003','SBSN_S_004','SBSN_S_005','SBSN_S_006','SBSN_S_007','SBSN_S_044','SBSN_S_008','SBSN_S_009'}; 

zScore = 3.1;

copeFile = 'cope1.feat';
% analysisFile = 'level_two_all_force_FLOB';

% --- SPINAL LEVELS MASK (SAME AS FIRST SCRIPT) ---
[spineLevels, ~] = cbiReadNifti('F:\SMA_HOLDER\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii');
spineLevels(spineLevels < 5) = 0;
spineLevels(spineLevels > 8) = 0;
spineLevels(spineLevels > 0) = spineLevels(spineLevels > 0) - 4;

% --- WHITE MATTER PROJECTION (SAME AS FIRST SCRIPT) ---
[wm, ~] = cbiReadNifti('F:\SMA_HOLDER\MRI_data_upper_limb\Spine\template\PAM50_wm.nii');
projectedWm = sum((wm).*(spineLevels>=1), 3);

projectedWm(projectedWm <= 20) = 0;
projectedWm(projectedWm > 20) = 1;

% tweak the WM mask in the same ROI as before
sub = projectedWm(63:79, 68:77);
sub(sub < 1) = 100;
projectedWm(63:79, 68:77) = sub;

projectedWm(projectedWm <= 20) = 0;
projectedWm(projectedWm > 20) = 1;

% PER-SUBJECT: projected plane + per-voxel z-lists

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

        % is dir and name contains gfeat
        if subjectFolder(folder).isdir && contains(subjectFolder(folder).name, 'level_two_all_force_FLOB1234.gfeat')

            disp(subjectFolder(folder).name)

            fileName = strsplit(subjectFolder(folder).name, '.');

            if ~exist(fullfile(direc, subjectFolder(folder).name, copeFile, 'thresh_zstat1.nii'))
                gunzip(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii.gz'));
            end

            [dataFile, ~] = cbiReadNifti(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii'));

            if ismember(i, [1:4 6:subSplit subSplit+3 subSplit+5 subSplit+9 subSplit+10])
                dataFile = flip(dataFile, 1);
            end
            disp(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii'))

            % this is the spinal projection
            % first sum the count of supra threshold voxels across all the slice we
            % want into one transverse plane for each subject
            % then I want to binraize it so count >= 1 is a 1
            % then when I combine across all subjects I use a count to show a heatmap
            % of how common a location of activation was acorss alll subjects

            % sum the count of 1s that are in the plane across these spinal levels
            projectedPlane = sum((dataFile>=zScore).*(spineLevels>=1), 3);
%           projectedPlane(projectedPlane>=1) = 1;

            % collect per-voxel suprathreshold z-scores (like in the first script)
            dataFile(spineLevels<1) = 0;
            zPlane = cell(size(dataFile,1), size(dataFile,2));
            for ii = 1:size(dataFile, 1)
                for jj = 1:size(dataFile, 2)
                    vals = squeeze(dataFile(ii,jj,:));   % values along the 3rd dim
                    zPlane{ii,jj} = vals(vals > zScore); % suprathreshold z's only
                end
            end            

            % store for this subject
            allData{i, 1}{fileCounter, 1} = projectedPlane;          % count plane
            allData{i, 1}{fileCounter, 2} = subjectFolder(folder).name;
            allData{i, 1}{fileCounter, 3} = zPlane;                  % per-voxel z

            fileCounter = fileCounter + 1;
        end

    end
end

% GROUP-LEVEL: CONTROL vs STROKE planes + pooled z-values

planeH = zeros(141, 141);   % Control
planeS = zeros(141, 141);   % Stroke

ZsumH = cell(141, 141);     % Control z-scores per voxel
ZsumS = cell(141, 141);     % Stroke z-scores per voxel

for i = 1:length(allData)
    % combine across all subjects: count how common activation is

    if isempty(allData{i,1})
        continue;   % skip subjects with no data
    end

    if i <= subSplit
        % --- CONTROL GROUP ---
        planeH = planeH + allData{i,1}{1,1};

        zPlane = allData{i,1}{1,3};
        for ii = 1:size(zPlane, 1)
            for jj = 1:size(zPlane, 2)
                ZsumH{ii, jj} = [ZsumH{ii, jj}; zPlane{ii, jj}];
            end
        end

    else
        % --- STROKE GROUP ---
        planeS = planeS + allData{i,1}{1,1};

        zPlane = allData{i,1}{1,3};
        for ii = 1:size(zPlane, 1)
            for jj = 1:size(zPlane, 2)
                ZsumS{ii, jj} = [ZsumS{ii, jj}; zPlane{ii, jj}];
            end
        end
    end

end

% FIGURE 1: z-colored dots (Control vs Stroke), same style as earlier script

rows = 57:85;          % z-range
cols = 60:83;          % y-range
j    = 0.99;           % jitter in "cell units"

planes  = {planeH, planeS};
Zplanes = {ZsumH,  ZsumS};
titles  = {'Control','Stroke'};

% dot size range
sMin = 3;      % size at z = zScore
sMax = 35;     % size at z = 5

% WM mask window (same transforms as earlier)
maskSub = fliplr(rot90(projectedWm(rows, cols))) > 0;

f = figure('Units','normalized','OuterPosition',[0 0 1 1]);
t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

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
        if isempty(zvec)
            continue;
        end

        % jittered coordinates
        yJ = cc + (rand(numel(zvec),1) - 0.5)*2*j;
        zJ = rr + (rand(numel(zvec),1) - 0.5)*2*j;

        % size by z-score (clamped)
        zCap = min(max(zvec, zScore), 5);
        s    = sMin + (sMax - sMin) * (zCap - zScore) / (5 - zScore);

        % classify by WM mask
        yCell = min(max(round(yJ), 1), size(maskSub, 2));
        zCell = min(max(round(zJ), 1), size(maskSub, 1));
        inMask = maskSub(sub2ind(size(maskSub), zCell, yCell));

        % outside WM: filled dots, no edge
        scatter(yJ(~inMask), zJ(~inMask), s(~inMask), zvec(~inMask), 'filled');

        % inside WM: black-edged dots
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

% Reserve the 3rd tile for the colorbar (like tile 6 in your original)
nexttile(t, 3); axis off
cb = colorbar(ax(1));              % tie colorbar to the first axes' CLim/colormap
cb.Layout.Tile = 3;
cb.Label.String = 'Voxel z-score';

saveas(f, ['D:\SBSN\stroke\Spine_transverse_gradient_Control_vs_Stroke_z', num2str(zScore), '.png']);
saveas(f, ['D:\SBSN\stroke\Spine_transverse_gradient_Control_vs_Stroke_z', num2str(zScore), '.svg']);

%FIGURE 2: count-based density dots (Control vs Stroke), same logic as earlier

figure; clf;
maskSub = fliplr(rot90(projectedWm(rows, cols))) > 0;

for k = 1:numel(planes)
    subplot(1,2,k);

    % counts submatrix with identical transforms
    Msub = fliplr(rot90(planes{k}(rows, cols)));

    % nonzero cells and replication counts
    [zIdx, yIdx, vals] = find(Msub);
    rep = floor(vals);                 % *** no extra normalization, same logic ***
    rep(~isfinite(rep)) = 0;
    rep = max(rep, 0);
    keep = rep > 0;
    zIdx = zIdx(keep);
    yIdx = yIdx(keep);
    rep   = rep(keep);

    % build repeated coordinates
    yAll = repelem(yIdx, rep);
    zAll = repelem(zIdx, rep);

    % jitter so overlapping points spread out
    rng(1 + k); % reproducible, different per panel
    yJ = yAll + (rand(size(yAll)) - 0.5)*2*j;
    zJ = zAll + (rand(size(zAll)) - 0.5)*2*j;

    % map jittered points to nearest cell and check mask membership
    yCell = round(yJ);  zCell = round(zJ);
    yCell = min(max(yCell, 1), size(maskSub, 2));
    zCell = min(max(zCell, 1), size(maskSub, 1));

    inMask = maskSub(sub2ind(size(maskSub), zCell, yCell));   % logical vector

    % plot (gray = outside mask, red = inside mask)
    scatter(yJ(~inMask), zJ(~inMask), 10, 'filled', 'MarkerFaceColor', [.7 .7 .7]); hold on;
    scatter(yJ(inMask),  zJ(inMask),  10, 'filled', 'r');
    axis equal tight
    xlim([0.5, size(maskSub,2)+0.5]); 
    ylim([0.5, size(maskSub,1)+0.5]);
    set(gca,'YDir','reverse');  % image-like orientation
    xlabel('y'); ylabel('z');
    title(titles{k});
    box on

    if k == numel(planes)
        legend({'outside mask','inside mask'}, 'Location','best');
    end
end

saveas(gcf, ['D:\SBSN\stroke\Spine_transverse_density_Control_vs_Stroke_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\stroke\Spine_transverse_density_Control_vs_Stroke_z', num2str(zScore), '.svg']);

% % ORIGINAL HEATMAP-STYLE FIGURE (unchanged logic: normalize here only)
% 
% planeH = planeH/max(planeH(:));
% planeS = planeS/max(planeS(:));
% 
% figure;
% subplot(1, 2, 1)
% imagesc(fliplr(rot90(planeH(57:85,60:83))))
% colormap('hot')
% title('Control')
% 
% subplot(1, 2, 2)
% imagesc(fliplr(rot90(planeS(57:85,60:83))))
% colormap('hot')
% title('Stroke')
% 
% colormap('hot')
% h = colorbar('Position', [0.93 0.11 0.02 0.815]); % [x y width height]
% 
% saveas(gcf, ['D:\SBSN\stroke\Spine_transverse_z', num2str(zScore), '.png']);
% saveas(gcf, ['D:\SBSN\stroke\Spine_transverse_z', num2str(zScore), '.svg']);

%%
function plotCreator(value,len)
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
    plot(value(1,:)', '.-b', 'MarkerSize', 15)
    hold on
    for i = 2:length(value) 

        plot(value(i,:)', '.-', 'Color', 'b', 'MarkerSize', 15)
        hold on

    end
    hold on
    plot(mean(value, 'omitnan'), 'k')
    errorbar(len, mean(value, 'omitnan'), std(value, 'omitnan')/sqrt(length(value)), 'Color','black')

end

function findSlopePts(data, time)

%     time = [1, 2, 3, 4, 5]; % Example time data
    [changepts,~] = findchangepts(mean(data), 'Statistic', 'linear', 'MaxNumChanges', 5);
    
    % Step 2: Plot data with detected changepoints
    figure;
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