% This file is to create a slice_timing text file and a physio text file
% the slice timing will be a conitnous file of the time it take each slice
% of each volume was acquired while the physio is the value of the
% respiration and heart rate at each of those slices
clear all

%%% SPINE

% varibales to set up before
subName = 'SBSN_H_018';
disp(subName)

volRemoved = 5;
TR = 2.2;

% adding paths to the subject
addpath('D:\NHP_code\cbiNifti')

direc = fullfile('D:\SBSN\Data\Spine', subName, 'physio');
physio_folders = dir(direc);

direc2 = fullfile('D:\SBSN\Data\Spine',subName,'func');
slice_number_folder = dir(direc2);

%% THINGS TO ADD
% I need to remove the timings for the first x volumes
% these get removed from the original length of the images

for folder = 3:length(physio_folders)

%     if contains(physio_folders(folder).name, 'physio0')
%         continue
%     end
    
    % here its just opening the txt file we need to parse through which has
    % all of our information
    fid = fopen(fullfile(direc, physio_folders(folder).name, 'test.txt'),'rt');
    physio = textscan(fid, '%[^\n]','HeaderLines', 3);
    fclose(fid);
    
    % these are two common strings that we will use to find our spots in
    % the files
    check = "ACQ_TIME_TICS  CHANNEL  VALUE  SIGNAL";
    check2 = "VOLUME   SLICE   ACQ_START_TICS  ACQ_FINISH_TICS  ECHO";
    
    
    % unzip the files from the .gz extention or they dont read into matlab
    gunzip(fullfile(direc2, slice_number_folder(folder).name, 'fmri_spine_slices_mean.nii.gz'))
    gunzip(fullfile(direc2, slice_number_folder(folder).name, 'fmri_spine_moco_mean.nii.gz'))
    % here i load in the two seperate files created in the beginning during segmentation
    [sliceFileWhole, ~] = cbiReadNifti(fullfile(direc2, slice_number_folder(folder).name, 'fmri_spine_slices_mean.nii')); 
    [sliceFile, ~] = cbiReadNifti(fullfile(direc2, slice_number_folder(folder).name, 'fmri_spine_moco_mean.nii')); 
    
    % now delete the .nii files because they can mess us up in later steps
    % of other scripts
    delete(fullfile(direc2, slice_number_folder(folder).name, 'fmri_spine_slices_mean.nii'))
    delete(fullfile(direc2, slice_number_folder(folder).name, 'fmri_spine_moco_mean.nii'))
    
    
    % the slice whole has values of 1000 where it was cropped from
    % so find the first one that doesnt have a mean value of 1000 and that means that row
    % was included in the crop
    checkMean = 1000;
    idx = 1;
    while mean(mean(sliceFileWhole(:, :, idx))) == 1000
        
        idx = idx + 1;
        
    end

    % here is checking for the specific strings we know will be in the file and
    % grabs the index of them
    indexes = [];
    count = 1;
    for i = 1:length(physio{1, 1})

        if strcmp(check, physio{1,1}{i, 1}) ||  strcmp(check2, physio{1,1}{i, 1}) 

            indexes(count) = i; 
            count = count + 1;

        end

        % I forget what this part does
        % doesnt seem to be attached to anything though but gonna leave in
        if contains(physio{1,1}{i, 1}, 'FirstTime')
            lastIdx = i;
        end       
        
    end

    % we then seperate out the resp and puls files
    puls_section = {physio{1,1}{indexes(1)+1:indexes(2)-1, 1}}';
    resp_section = {physio{1,1}{indexes(2)+1:indexes(3)-1, 1}}';

    % remove the bad text at the end of each file
    rem_puls = [];
    count = 1;
    for i = length(puls_section)-20:length(puls_section)
        rem_puls(count) = contains(puls_section{i,1}, 'PULS');
        count = count+1;
    end
    puls_section = {puls_section{1:end-sum(nonzeros(~rem_puls)), 1}}';

    rem_resp = [];
    count = 1;
    for i = length(resp_section)-20:length(resp_section)
        rem_resp(count) = contains(resp_section{i,1}, 'RESP');
        count = count+1;
    end
    resp_section = {resp_section{1:end-sum(nonzeros(~rem_resp)), 1}}';


    %% this is rambling below
    % we should take the timings from the bottom of the list. and then use them
    % The manually specified slicetiming file needs to be a text file that 
    % specifies the timing (in seconds) of every slice relative to the start of 
    % the volume acquisition (t=0). The file must have one number per slice and be 
    % arranged in one line, separated by spaces, and each number is the timing value 
    % for that slice (e.g. third number gives the timing of the third slice,
    % where numbering of slices is consistent with the NIFTI voxel coordinate;
    % i.e. the way that FSLView/FSLeyes shows them).

    % first time of first slice of first volume subtracted 
    % from time of start of the first slice of next volume = 880 tics
    % which sampled at 400 Hz is 2.2 seconds.
    % so I can find the points in the physio that match up with this.
    % at 100 Hz I can match the triggers for the starts of the first slices and match these up
    
    %%
    % here we have to set up a timing file per physio since the slice timings are
    % interleaved due to the multiband factor and accelration used in the sequence
    % we have to know what slice we cut off for the spine and use those timings to help with the regression
    % this next loop creates the timing file of the slices
    idxSub = 0;
    if ~contains(physio{1,1}{end}, 'LastTime')
        idxSub = 1;
    end
%     endIndex = contains(physio{1,1}{indexes(3)+30:end, 1}, 'FirstTime');
    slice_timing = {physio{1,1}{indexes(3)+1:end-2-idxSub, 1}}';
    sliceTimings = zeros(100, 163);
    mult = 1;
    for i = 1:length(slice_timing)

         slice_split = split(string(slice_timing{i,1}));


         if i == 1
             startTimer = str2double(slice_split(3));

             sliceTimings(str2double(slice_split(2))+1, mult) = (str2double(slice_split(3)) - startTimer) / 400;


             startTime = str2double(slice_split(3));
         else

             sliceTimings(str2double(slice_split(2))+1, mult) = (str2double(slice_split(3)) - startTimer) / 400;
             stopTime = str2double(slice_split(4));
         end


        if mod(i, 100) == 0
            mult = mult + 1;
        end        
    end
    
    % where i need to find the correct slices of the cropped file to grab
    % I think it needs to be 
%     sliceTimings = sliceTimings(size(sliceFileWhole, 3) - size(sliceFile, 3) + 1:size(sliceFileWhole, 3)-idx+1, :);
    sliceTimings = sliceTimings(idx:idx+size(sliceFile, 3)-1, :);

    % now its in the slices x volumes shape
    % we just get rid of the first x vloumes and subtract that last time
    % which is TR*num volumes removed which is 5
    % from all of them making the first slice time 0 instead
    sliceTimings = sliceTimings(:, volRemoved+1:end) - volRemoved*TR;
    
    % save slice timing matrix
    % The slicetiming file needs to be a text matrix that specifies the timing (in seconds) 
%     of every slice relative to the start of the volume acquisition (t=0). 
%     The matrix must have dimensions: (number of slices) by (number of volumes). 
%     The values represent the timings for one particular slice, across all volumes 
%     (e.g. timing of slice 2 in volume 1, slice 2 in volume 2, slice 2 in volume 3, etc.)
%     which allows for a variable slice timing per volume. 


    writematrix(sliceTimings', fullfile(direc, physio_folders(folder).name, 'slice_times.txt'), 'Delimiter', 'space');
    %NumSlicesx158 timings of each slice from 0 (in seconds) from t=0
    
% %     The slicetiming file needs to be a text matrix that specifies the timing (in seconds)
% %     of every slice relative to the start of the volume acquisition (t=0). 
% %     The matrix must have dimensions: (number of slices) by (number of volumes). 
% %     The values represent the timings for one particular slice, across all volumes 
% %     (e.g. timing of slice 2 in volume 1, slice 2 in volume 2, slice 2 in volume 3, etc.)
% %     which allows for a variable slice timing per volume.



%     % i dont know whats correct anymore tbh
%     sliceMatrix = sliceTimings;
%     sliceMatrix(:,2:end) = repmat(sliceMatrix(:,1), 1, size(sliceMatrix(:,2:end), 2));
%     writematrix(sliceMatrix, fullfile(direc, physio_folders(folder).name, 'slice_matrix.txt'), 'Delimiter', 'space');
%     
%     %reshape into one column format
%     sliceTimings(:,2:end) = repmat(sliceTimings(:,1), 1, size(sliceTimings(:,2:end), 2));
% 
%     sliceTimings = reshape(sliceTimings, 1, []);
%     
%     % save the vector
%     writematrix(sliceTimings, fullfile(direc, physio_folders(folder).name, 'slice.txt'), 'Delimiter', 'space');

    %%%%%%%%%%%%%%%%%%%%%
    
    % this next 2 blocks just gets the values of the puls and resp
    puls_value = [];
    count = 1;
    for i = 1:length(puls_section)

         puls_split = split(string(puls_section{i,1}));

         if str2double(puls_split(1)) >= startTime && str2double(puls_split(1)) <= stopTime
             puls_value(count, 1) = str2double(puls_split(1));

             puls_value(count, 2) = str2double(puls_split(3));

             count = count + 1;
         end
    end


    resp_value = [];
    count = 1;
    for i = 1:length(resp_section)
         resp_split = split(string(resp_section{i,1}));

         if str2double(resp_split(1)) >= startTime && str2double(resp_split(1)) <= stopTime

             resp_value(count,1) = str2double(resp_split(1)); 

             resp_value(count,2) = str2double(resp_split(3)); 

             count = count + 1;
         end
    end

    % here we down scale puls sinces its 200 Hz and upscale resps since its 50 Hz
    % get them both to a sample rate of 100Hz
    % then here we can start and since we know the rate of recording and we
    % know the TR we can get rid of the first 11 seconds of samples which
    % is TR*volumes removed which is 5
    timeRemoval = volRemoved*TR;
    
%     puls_value_down = puls_value(1:2:end, :); %puls_value(timeRemoval*200:2:end, :);
%     resp_value_up = interp(resp_value(:, 2), 2); %
   
    puls_value_down = puls_value(timeRemoval*200:2:end, :);
    resp_value_up = interp(resp_value(timeRemoval*50:end, 2), 2); %
    sizeI = min(length(puls_value_down), length(resp_value_up));

    physio_data = horzcat(puls_value_down(1:sizeI, 1:2), resp_value_up(1:sizeI , :));
    
    % this adds trgiger every 2.2 s
    physio_data(1:220:end, 4) = 1;
    
    % save the file
    writematrix(physio_data(:, 2:end), fullfile(direc, physio_folders(folder).name, 'physio.txt'), 'Delimiter', 'space')


    % please what does 1 tic mean. WHy do I have so many gaps between
    % ptrig = diff(puls_value(find(puls_value(:,3)), 1));
    % rtrig = diff(puls_value(find(resp_value(:,3)), 1));

    % puls is every 2 tics
    % resp is every 8 tics

    % what does 1 tic mean? 
    % 2200 ms for TR

    % so the scan is 398 seconds in total from when the trigger starts to end. The
    % scan last 360 seconds so we ignore the first 40 seconds of data and then 
    % take the last 360 seconds worth of samples 
    % we know puls is 200 Hz samples while resp is 50 Hz


    % [idx,loc] = ismember(sliceTimings(1, 3), puls_value);
    % out = loc(idx);

    % find(sliceTimings(1, 3) == puls_value)

    % fid = fopen('D:\SBSN\SBSN_S_001_raw\MRI\flywheel\pirondinilab\SBSN\SBSN_S_001\ses-01\func-bold_task-sbsn_run-02__1.5mmaniso_PhysioLog\1.3.12.2.1107.5.2.43.166114.2022032411283075813970568.0.0.0.dicom\smooth_popp_card.txt','rt');
    % checking2 = textscan(fid, '%[^\n]','HeaderLines', 3);
    % fclose(fid);

    % lTotalScanTimeSec = 398;
    % 0.0025 samples/second
    % 400 Hz
    % 200 Hz for puls
    % 50 Hz for resp
    % final_values = horzcat(upscaled_resp2, check4(:,1));

    % writematrix(final_values, 'D:\SBSN\SBSN_S_001_raw\MRI\flywheel\pirondinilab\SBSN\SBSN_S_001\ses-01\func-bold_task-sbsn_run-02__1.5mmaniso_PhysioLog\1.3.12.2.1107.5.2.43.166114.2022032411283075813970568.0.0.0.dicom\sbsn.txt', 'Delimiter', 'space')
end