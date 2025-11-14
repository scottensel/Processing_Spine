%% Brain HMAT Analysis (12 regions from single HMAT.nii)
% Adds Template Matching and DICE sections using HMAT (labels 1..12)

close all; clear;

addpath('D:\NHP_code\cbiNifti')   % cbiReadNifti

% -----------------------------
% Study/sample configuration
% -----------------------------
subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_007','SBSN_H_008','SBSN_H_010','SBSN_H_017','SBSN_H_018'};
zScore  = 3.1;
copeFile = 'cope1.feat';

% -----------------------------
% HMAT atlas (single file, 12 labels)
% -----------------------------
% Index mapping (HMAT.nii):
% 1 R_M1, 2 L_M1, 3 R_S1, 4 L_S1,
% 5 R_SMA, 6 L_SMA, 7 R_preSMA, 8 L_preSMA,
% 9 R_PMd, 10 L_PMd, 11 R_PMv, 12 L_PMv
hmatPath = 'D:\SBSN\Data\Brain\template\HMAT_website\HMAT_2mm.nii';
[HMAT, ~] = cbiReadNifti(hmatPath);

regionNames = {'Right M1','Left M1','Right S1','Left S1', ...
               'Right SMA','Left SMA','Right preSMA','Left preSMA', ...
               'Right PMd','Left PMd','Right PMv','Left PMv'};
nRegions = numel(regionNames);   % 12

% ======================================================================
% A) LEVEL-TWO (successive-run) summary inside HMAT
% ======================================================================
allData = {};
for i = 1:length(subName)
    direc = fullfile('D:\SBSN\Data\Brain', subName{i}, 'func');
    subjectFolder = dir(direc);

    disp(subName{i})

    allData{i, 1} = {};
    allData{i, 2} = subName{i};

    fileCounter = 1;
    for folder = 3:length(subjectFolder)

        if subjectFolder(folder).isdir && contains(subjectFolder(folder).name, 'level_two_force_FLOB')

            disp(subjectFolder(folder).name)

            if ~exist(fullfile(direc, subjectFolder(folder).name, copeFile, 'thresh_zstat1.nii'))
                gunzip(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii.gz'));
            end

            [dataFile, ~] = cbiReadNifti(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii'));

            if i == 1
                dataFile = flip(dataFile, 1);
            end

            numVoxelsHMAT = sum(sum(sum((dataFile >= zScore) .* (HMAT >= 1))));
            magHMAT = dataFile(HMAT >= 1);

            numVoxelsSeperate = [];
            magSeperate = [];
            for r = 1:nRegions
                msk = (HMAT == r);
                numVoxelsSeperate(r) = sum(sum(sum((dataFile(msk) >= zScore))));
                numVoxelsSeperate(r) = numVoxelsSeperate(r) / max(1, sum(msk(:))) * 100;

                vr = dataFile(msk);
                vr = vr(vr > zScore);
                if ~isempty(vr)
                    magSeperate(r) = mean(vr);
                end
            end

            allData{i,1}{fileCounter, 1} = [ numVoxelsHMAT / sum(HMAT(:) >= 1) * 100,  mean(magHMAT(magHMAT>zScore)), std(magHMAT(magHMAT>zScore)) ];
            allData{i,1}{fileCounter, 2} = subjectFolder(folder).name;
            allData{i,1}{fileCounter, 3} = [numVoxelsSeperate; magSeperate];

            fileCounter = fileCounter + 1;
        end
    end
end

for i = 1:length(allData)
    for j = 1:length(allData{i,1})
        activeVoxels(i,j) = allData{i,1}{j,1}(1);
        zScores(i,j)      = allData{i,1}{j,1}(2);

        if j == 1
            actVoxelsSeg2(i,:) = allData{i,1}{j,3}(1,:);  zSeg2(i,:) = allData{i,1}{j,3}(2,:);
        elseif j == 2
            actVoxelsSeg3(i,:) = allData{i,1}{j,3}(1,:);  zSeg3(i,:) = allData{i,1}{j,3}(2,:);
        elseif j == 3
            actVoxelsSeg4(i,:) = allData{i,1}{j,3}(1,:);  zSeg4(i,:) = allData{i,1}{j,3}(2,:);
        elseif j == 4
            actVoxelsSeg5(i,:) = allData{i,1}{j,3}(1,:);  zSeg5(i,:) = allData{i,1}{j,3}(2,:);
        elseif j == 5
            actVoxelsSeg6(i,:) = allData{i,1}{j,3}(1,:);  zSeg6(i,:) = allData{i,1}{j,3}(2,:);
        end
    end
end

% Example figure (whole HMAT)
plotCreator(activeVoxels, 1:5);
make_pretty
xlim([0.75,5.25])
xlabel('Run Combination')
ylabel('Active Voxels (%)');
title(sprintf('Average Active Successive runs (HMAT)'));
xticks(1:5)
xticklabels({'1-2','1-3','1-4','1-5','1-6'})

figure;
hBar=barh([mean(actVoxelsSeg4)', mean(actVoxelsSeg6)']);
X=cell2mat(get(hBar,'XData')).'+[hBar.XOffset];
hold on  %4 runs
hEB = errorbar([mean(actVoxelsSeg4)', mean(actVoxelsSeg6)'], X, [(std(actVoxelsSeg4)/sqrt(length(actVoxelsSeg4)))',  (std(actVoxelsSeg6)/sqrt(length(actVoxelsSeg6)))'], 'horizontal', '.', 'Color', 'black', 'Marker', 'none');  % add the errorbar
% errorbar(mean(actVoxelsSeg4), 1:4, std(actVoxelsSeg4)/sqrt(length(actVoxelsSeg4)), 'horizontal', '.', 'Color','black')
randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(1:4,:), 1),1))/10;
scatter(actVoxelsSeg4(1:4,:), [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1), randVec+X(5,1), randVec+X(6,1), randVec+X(7,1), randVec+X(8,1), randVec+X(9,1), randVec+X(10,1), randVec+X(11,1), randVec+X(12,1)], 30, 'k','o','filled'); 
scatter(actVoxelsSeg6(1:4,:), [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2), randVec+X(5,2), randVec+X(6,2), randVec+X(7,2), randVec+X(8,2), randVec+X(9,2), randVec+X(10,2), randVec+X(11,2), randVec+X(12,2)], 30, 'k','o','filled'); 
randVec = (-1 + (1+1)*rand(size(actVoxelsSeg4(5:end,:), 1),1))/10;
scatter(actVoxelsSeg4(5:end,:), [randVec+X(1,1), randVec+X(2,1), randVec+X(3,1), randVec+X(4,1), randVec+X(5,1), randVec+X(6,1), randVec+X(7,1), randVec+X(8,1), randVec+X(9,1), randVec+X(10,1), randVec+X(11,1), randVec+X(12,1)], 60, 'k','x'); 
scatter(actVoxelsSeg6(5:end,:), [randVec+X(1,2), randVec+X(2,2), randVec+X(3,2), randVec+X(4,2), randVec+X(5,2), randVec+X(6,2), randVec+X(7,2), randVec+X(8,2), randVec+X(9,2), randVec+X(10,2), randVec+X(11,2), randVec+X(12,2)], 60, 'k','x'); 
set (gca,'YDir','reverse')
yticks(1:length(1:12)); yticklabels({'M1 (R)','M1 (L)','S1 (R)','S1 (L)', 'SMA (R)','SMA (L)','preSMA (R)','preSMA (L)', 'PMd (R)','PMd (L)','PMv (R)','PMv (L)'});
ylabel('Brain Area')
xlabel('Active Voxels');
title(sprintf('Average Active Voxel 4 vs 6 Runs combined'));
make_pretty
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Brain_voxel_success_HMAT_z', num2str(zScore), '.png']);
% saveas(gcf, ['D:\SBSN\Manuscript\plots\Brain_voxel_success_HMAT_z', num2str(zScore), '.svg']);

% ======================================================================
% B) TEMPLATE MATCHING (Level-ONE per-run inside HMAT)
%     Rest j=0 uses zstat1.nii; Task j>0 uses zfstat1.nii (your style)
% ======================================================================
allData_TM = {};
for i = 1:length(subName)

    allData_TM{i, 1} = {};
    allData_TM{i, 2} = subName{i};

    for j = 0:6
        direc = fullfile('D:\SBSN\Data\Brain', subName{i}, 'func', ['func', num2str(j)]);
        subjectFolder = dir(direc);

        disp([subName{i} ' func' num2str(j)])

        for folder = 3:length(subjectFolder)
            if subjectFolder(folder).isdir && contains(subjectFolder(folder).name, 'level_one_force_FLOB') && (j > 0)

                if ~exist(fullfile(direc, subjectFolder(folder).name, 'stats\zfstat1.nii'))
                    gunzip(fullfile(direc, subjectFolder(folder).name,  'stats\zfstat1.nii.gz'));
                end

                [dataFile, ~] = cbiReadNifti(fullfile(direc, subjectFolder(folder).name,  'stats\zfstat1.nii'));
                if i == 1
                    dataFile = flip(dataFile, 1);
                end

                numVoxels = sum(sum(sum((dataFile>=zScore).*(HMAT>=1))));
                mag       = dataFile(HMAT>=1);

                allData_TM{i, 1}{j+1, 1} = [numVoxels/sum(HMAT(:)>=1)*100, mean(mag(mag>zScore)), std(mag(mag>zScore))];
                allData_TM{i, 1}{j+1, 2} = subjectFolder(folder).name;

            elseif subjectFolder(folder).isdir && contains(subjectFolder(folder).name, 'level_one_FLOB') && (j == 0)

                if ~exist(fullfile(direc, subjectFolder(folder).name, 'stats\zstat1.nii'))
                    gunzip(fullfile(direc, subjectFolder(folder).name,  'stats\zstat1.nii.gz'));
                end

                [dataFile, ~] = cbiReadNifti(fullfile(direc, subjectFolder(folder).name,  'stats\zstat1.nii'));
                if i == 1
                    dataFile = flip(dataFile, 1);
                end

                numVoxels = sum(sum(sum((dataFile>=zScore).*(HMAT>=1))));
                mag       = dataFile(HMAT>=1);

                allData_TM{i, 1}{j+1, 1} = [numVoxels/sum(HMAT(:)>=1)*100, mean(mag(mag>zScore)), std(mag(mag>zScore))];
                allData_TM{i, 1}{j+1, 2} = subjectFolder(folder).name;

            end
        end
    end
end

for i = 1:length(allData_TM)
    for j = 1:length(allData_TM{i,1})
        activeVoxels_TM(i,j) = allData_TM{i,1}{j,1}(1);
        zScores_TM(i,j)      = allData_TM{i,1}{j,1}(2);
    end
end

% Rest vs Task (mean over task runs)
activeVoxels_TM2 = [activeVoxels_TM(:,1), mean(activeVoxels_TM(:,2:end),2,'omitnan')];
zScores_TM2      = [zScores_TM(:,1),      mean(zScores_TM(:,2:end),2,'omitnan')];

plotCreator(activeVoxels_TM2, 1:2);
make_pretty
xlim([0.75,2.25])
ylabel('Active Voxels (%)')
xlabel('Group');
title(sprintf('Mean Active Voxels — Rest vs Task (HMAT)'));
xticklabels({'Rest','Task'})
xticks(1:2)
saveas(gcf, ['D:\SBSN\Manuscript\plots\Brain_voxel_control_HMAT_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\Brain_voxel_control_HMAT_z', num2str(zScore), '.svg']);

plotCreator(zScores_TM2, 1:2);
make_pretty
xlim([0.75,2.25])
ylabel('Z-Score')
xlabel('Group');
title(sprintf('Mean Z-Score — Rest vs Task (HMAT)'));
xticklabels({'Rest','Task'})
xticks(1:2)
saveas(gcf, ['D:\SBSN\Manuscript\plots\Brain_zscore_control_HMAT_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\Brain_zscore_control_HMAT_z', num2str(zScore), '.svg']);

% ======================================================================
% C) DICE (successive level-two maps) within HMAT, and per HMAT region
% ======================================================================
clear diceCell

allDice = [];
allDiceSep = [];

for i = 1:length(subName)
    direc = fullfile('D:\SBSN\Data\Brain', subName{i}, 'func');
    subjectFolder = dir(direc);

    disp(['DICE: ' subName{i}])

    fileCounter = 1;
    diceCell = {};
    for folder = 3:length(subjectFolder)
        if subjectFolder(folder).isdir && contains(subjectFolder(folder).name, 'level_two_force_FLOB')

            if ~exist(fullfile(direc, subjectFolder(folder).name, copeFile, 'thresh_zstat1.nii'))
                gunzip(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii.gz'));
            end
            [dataFile, ~] = cbiReadNifti(fullfile(direc, subjectFolder(folder).name,  copeFile, 'thresh_zstat1.nii'));
            if i == 1
                dataFile = flip(dataFile, 1);
            end

            diceCell{fileCounter} = dataFile;
            fileCounter = fileCounter + 1;
        end
    end

    diceValues = [];
    diceValuesSeperate = [];
    for j = 1:length(diceCell)-1
        A = (diceCell{j}   >= zScore);
        B = (diceCell{j+1} >= zScore);

        % whole HMAT
        diceValues(j) = dice(A & (HMAT>=1), B & (HMAT>=1));

        % per HMAT region (12)
        for k = 1:nRegions
            Ak = A & (HMAT == k);
            Bk = B & (HMAT == k);
            section = dice(Ak, Bk);
            if isempty(section)
                diceValuesSeperate(j, k) = 0;
            else
                diceValuesSeperate(j, k) = section;
            end
        end
    end

    allDice(i,:)        = diceValues;
    allDiceSep(:,:, i)  = diceValuesSeperate;
end

% Whole HMAT DICE
plotCreator(allDice, 1:4);
make_pretty
xlim([0.75,4.25])
xlabel('Run Combination')
ylabel('DICE Score');
title(sprintf('DICE (HMAT) — Successive Level-Two Maps'));
xticks(1:4)
xticklabels({'2vs3','3vs4','4vs5','5vs6'})
saveas(gcf, ['D:\SBSN\Manuscript\plots\Brain_DICE_success_HMAT_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\Brain_DICE_success_HMAT_z', num2str(zScore), '.svg']);

% Per-region DICE figures
brainNames = regionNames;
for j = 1:nRegions
    allDiceSep2 = squeeze(allDiceSep(:, j, :))';
    plotCreator(allDiceSep2, 1:4);
    make_pretty
    xlim([0.75,4.25])
    xlabel('Run Combination')
    ylabel('DICE Score');
    title(sprintf(['DICE Across Runs — ' brainNames{j}]));
    xticks(1:4)
    xticklabels({'2vs3','3vs4','4vs5','5vs6'})
    saveas(gcf, ['D:\SBSN\Manuscript\plots\Brain_DICE_' brainNames{j} '_HMAT_z' num2str(zScore) '.png']);
    saveas(gcf, ['D:\SBSN\Manuscript\plots\Brain_DICE_' brainNames{j} '_HMAT_z' num2str(zScore) '.svg']);
end

% ======================================================================
% D) Laterality Index (optional example using HMAT pairs, 4 vs 6)
% ======================================================================
leftIdx  = [2 4 6 8 10 12];
rightIdx = [1 3 5 7 9 11];
regionPairs = {'M1','S1','SMA','preSMA','PMd','PMv'};

LI_pre  = (actVoxelsSeg4(:,leftIdx) - actVoxelsSeg4(:,rightIdx)) ./ ...
          (actVoxelsSeg4(:,leftIdx) + actVoxelsSeg4(:,rightIdx));
LI_post = (actVoxelsSeg6(:,leftIdx) - actVoxelsSeg6(:,rightIdx)) ./ ...
          (actVoxelsSeg6(:,leftIdx) + actVoxelsSeg6(:,rightIdx));

meanLI_pre  = mean(LI_pre, 1, 'omitnan');  semLI_pre  = std(LI_pre, 0, 1, 'omitnan')/sqrt(size(LI_pre,1));
meanLI_post = mean(LI_post,1, 'omitnan');  semLI_post = std(LI_post,0, 1, 'omitnan')/sqrt(size(LI_post,1));

figure;
hBar = bar([meanLI_pre; meanLI_post]', 'grouped');
hold on;
ngroups = size([meanLI_pre; meanLI_post]',1);
nbars   = size([meanLI_pre; meanLI_post]',2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
yvals = [meanLI_pre; meanLI_post];
errors = [semLI_pre; semLI_post];
for b = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*b-1) * groupwidth / (2*nbars);
    errorbar(x, yvals(b,:), errors(b,:), 'k.', 'LineWidth', 1);
end
ylabel('Laterality Index (LI)'); xlabel('Region Pair');
xticks(1:length(regionPairs)); xticklabels(regionPairs);
legend({'4 runs','6 runs'}, 'Location','northeast');
title('Laterality Index (4 vs 6) — HMAT');
ylim([-1, 1]); make_pretty
saveas(gcf, ['D:\SBSN\Manuscript\plots\LI_HMAT_4v6_z', num2str(zScore), '.png']);
saveas(gcf, ['D:\SBSN\Manuscript\plots\LI_HMAT_4v6_z', num2str(zScore), '.svg']);


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

    opac = {[1, 0.5, 0.5]; [1, 0.15, 0.15]; [1, 0.25, 0.25]; [1, 0.35, 0.35]; [1, 0.45, 0.45]; [1, 0.55, 0.55]};

    figure;
    plot(value(1,:)', '.-r', 'MarkerSize', 15)
    hold on
    for i = 2:length(value) 
        
        if i > 4
            plot(value(i,:)', 'x-r', 'MarkerSize', 15)
            hold on
        else
            plot(value(i,:)', '.-r', 'MarkerSize', 15)
            hold on
        end

    end
    hold on
    plot(mean(value, 'omitnan'), 'k')
    errorbar(len, mean(value, 'omitnan'), std(value, 'omitnan')/sqrt(length(value)), 'Color', 'black')

end


function findSlopePts(data, time)

%     time = [1, 2, 3, 4, 5]; % Example time data
    [changepts,~] = findchangepts(mean(data), 'Statistic', 'linear', 'MaxNumChanges', 2);
    
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


function findSlopePtsSingle(data, time, i)

%     time = [1, 2, 3, 4, 5]; % Example time data
    [changepts,~] = findchangepts(data, 'Statistic', 'linear', 'MaxNumChanges', 2);
    
    % Step 2: Plot data with detected changepoints
%     figure;

    if i > 4
        plot(data, 'xr', 'MarkerSize', 15);

    else

        plot(data, '.r', 'MarkerSize', 15);

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
    plot(time, fittedResponse, '-r');

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