% File watcher and replacer
% Continuously monitors a folder for a trigger file, then replaces it

% Define paths
watchFolder = 'D:\SBSN\Data\Spine\SBSN_S_008\func\func4\physio_denoised.feat';           % Folder to watch
triggerFile = 'mask.nii.gz';                       % File to watch for
sourceFolder = 'D:\SBSN\Data\Spine\SBSN_S_008\func\func3\physio_denoised.feat';         % Folder to get replacement file
replacementFile = 'mask.nii.gz';               % File to copy in

% Full paths
triggerPath = fullfile(watchFolder, triggerFile);
sourcePath  = fullfile(sourceFolder, replacementFile);
destPath    = fullfile(watchFolder, replacementFile);

fprintf('Watching folder: %s\n', watchFolder);

while true
    if exist(triggerPath, 'file') == 2
        fprintf('Trigger file detected');
        
        % Delete the trigger file
        delete(triggerPath);
        fprintf('Deleted trigger file.\n');

        % Copy the replacement file
        copyfile(sourcePath, destPath);

        fprintf('Replaced with file: %s\n\n', replacementFile);

        exit
    end
    
    pause(0.04); % Polling interval (reduce CPU load)
end

% addpath('/Users/pirondinilab/spinalcordtoolbox/cbiNifti');
addpath('D:\NHP_code\cbiNifti')

gunzip('D:\SBSN\Data\Spine\SBSN_S_008\func\func4\fmri_spine_moco.nii.gz');
[spineLevels1, ~] = cbiReadNifti('D:\SBSN\Data\Spine\SBSN_S_008\func\func4\fmri_spine_moco.nii');
% watchFolder = 'D:\SBSN\Data\Spine\SBSN_S_008\func\func4\physio_denoised.feat';           % Folder to watch

gunzip('D:\SBSN\Data\Spine\SBSN_S_008\func\func3\fmri_spine_moco.nii.gz');
[spineLevels2, ~] = cbiReadNifti('D:\SBSN\Data\Spine\SBSN_S_008\func\func3\fmri_spine_moco.nii');