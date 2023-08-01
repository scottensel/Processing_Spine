addpath('D:\NHP_code\cbiNifti')

path = 'D:\SBSN\Data\Spine\SBSN_H_003\func\func1';

% unzip the files from the .gz extention or they dont read into matlab
gunzip(fullfile(path, 'fmri_spine_moco.nii.gz'))
gunzip(fullfile(path, 'fmri_spine_moco_mean_CSF_seg_corr.nii.gz'))


[sliceFileWhole, hdr1] = cbiReadNifti(fullfile(path, 'fmri_spine_moco.nii')); 
[sliceFile, hdr2] = cbiReadNifti(fullfile(path, 'fmri_spine_moco_mean_CSF_seg_corr.nii'));
size(sliceFileWhole)

hdr2.dim(2) = hdr1.dim(2);
hdr2.dim(3) = hdr1.dim(3);
hdr2.dim(4) = hdr1.dim(4);
sliceFile2 = sliceFile(:,1:133,1:53);
cbiWriteNifti(fullfile(path, 'fmri_spine_moco_mean_CSF_seg_corr.nii'), sliceFile2, hdr2);  

gzip(fullfile(path, 'fmri_spine_moco_mean_CSF_seg_corr.nii'))