clear all

%%% SPINE

% addpath('/Users/pirondinilab/spinalcordtoolbox/cbiNifti');
addpath('D:\NHP_code\cbiNifti')

% varibales to set up before
subName = {'SBSN_H_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_007','SBSN_H_008','SBSN_H_010',...
    'SBSN_H_101', 'SBSN_H_102', 'SBSN_H_103', 'SBSN_H_104',...
    'SBSN_S_001', 'SBSN_S_002', 'SBSN_S_003', 'SBSN_S_004', 'SBSN_S_055', 'SBSN_S_066', 'SBSN_S_077'}; 

%% THINGS TO ADD
for i = 1:length(subName)
    disp(subName{i})
    gunzip(['D:\SBSN\Data\Spine\' subName{i} '\anat\t2_roots.nii.gz']);
    [roots, ~] = cbiReadNifti(['D:\SBSN\Data\Spine\' subName{i} '\anat\t2_roots.nii']);

    gunzip(['D:\SBSN\Data\Spine\' subName{i} '\anat\t2_seg.nii.gz']);
    [seg, param] = cbiReadNifti(['D:\SBSN\Data\Spine\' subName{i} '\anat\t2_seg.nii']);

    gunzip(['D:\SBSN\Data\Spine\' subName{i} '\anat\t2_seg_dilate.nii.gz']);
    [segD, ~] = cbiReadNifti(['D:\SBSN\Data\Spine\' subName{i} '\anat\t2_seg_dilate.nii']);

    % multiplys the levels by segment that so each x,y that is overlapping with a
    % label is now a level
    disp(size(roots))
    if size(roots, 3) == 64
        newLevel = max(max(roots,[],1),[],3).*segD;
        newLevel(newLevel < 1) = 0;
        
        newLevel(seg~=1) = 0;
    
        seg(newLevel < 1) = 0;
        
        seg = seg.*newLevel;
    
        cbiWriteNifti(['D:\SBSN\Data\Spine\' subName{i} '\anat\t2_spinal_levels.nii'], seg, param);  
        gzip(['D:\SBSN\Data\Spine\' subName{i} '\anat\t2_spinal_levels.nii'])
    else
        newLevel = max(max(roots,[],1),[],2).*segD;
        newLevel(newLevel < 1) = 0;
        
        newLevel(seg~=1) = 0;
    
        seg(newLevel < 1) = 0;
        
        seg = seg.*newLevel;
    
        cbiWriteNifti(['D:\SBSN\Data\Spine\' subName{i} '\anat\t2_spinal_levels.nii'], seg, param);  
        gzip(['D:\SBSN\Data\Spine\' subName{i} '\anat\t2_spinal_levels.nii'])
    end
end
