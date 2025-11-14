% SPINAL PROJECTION
% this is the spinal projection
% first sum the count of supra threshold voxels across all the slice we
% want into one transverse plane for each subject
% then I want to binraize it so count >= 1 is a 1
% then when I combine across all subjects I use a count to show a heatmap
% of how common a location of activation was acorss alll subjects

%%% SPINE

% addpath('/Users/pirondinilab/spinalcordtoolbox/cbiNifti');
addpath('D:\NHP_code\cbiNifti')

% varibales to set up before
subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_007','SBSN_H_008','SBSN_H_010','SBSN_H_017','SBSN_H_018'}; 
% subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_008','SBSN_H_010'}; 
% subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_007','SBSN_H_008'}; 

zScore = 3.1;

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

figure; clf;
for k = 1:numel(planes)
    subplot(2,3,k); hold on;

    Msub   = fliplr(rot90(planes{k}(rows, cols)));        % counts (N×N)
    Zcells = fliplr(rot90(Zplanes{k}(rows, cols)));       % cell array (N×N)

    [rIdx, cIdx, ~] = find(Msub);                         % active cells

    for m = 1:numel(rIdx)
        rr = rIdx(m);
        cc = cIdx(m);
        zvec = Zcells{rr, cc}(:);                         % one z per dot

        % jittered positions
        yJ = cc + (rand(numel(zvec),1) - 0.5)*2*j;
        zJ = rr + (rand(numel(zvec),1) - 0.5)*2*j;

        % per-dot sizes mapped to z in [zScore, 5]
        zCap = min(max(zvec, zScore), 5);                 % cap to color/size range
        s    = sMin + (sMax - sMin) * (zCap - zScore) / (5 - zScore);

        % classify jittered dots relative to WM mask (nearest cell)
        yCell = min(max(round(yJ), 1), size(maskSub, 2));
        zCell = min(max(round(zJ), 1), size(maskSub, 1));
        inMask = maskSub(sub2ind(size(maskSub), zCell, yCell));

        % plot outside-WM dots (no outline)
        scatter(yJ(~inMask), zJ(~inMask), s(~inMask), zvec(~inMask), 'filled');

        % plot inside-WM dots (black outline)
        h = scatter(yJ(inMask), zJ(inMask), s(inMask), zvec(inMask), 'filled');
        set(h, 'MarkerEdgeColor','k', 'LineWidth', 1.25);   % outline
    end


    colormap(winter);
%     colormap(turbo);
    axis equal tight
    set(gca,'YDir','reverse');
    xlim([0.5, size(Msub,2)+0.5]);
    ylim([0.5, size(Msub,1)+0.5]);
    clim([zScore, 5])                  % anything above 5 uses max color
    colorbar
    xlabel('y');
    ylabel('z');
    title(titles{k});
    box on
end

subplot(2,3,6); axis off
