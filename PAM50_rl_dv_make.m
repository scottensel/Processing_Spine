gunzip('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_cord.nii.gz');
[spine, ~] = cbiReadNifti('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_cord.nii');

% PAM50_cervical_cord_all.nii.gz

gunzip('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii.gz');
[spineLevels, param] = cbiReadNifti('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii');
% spineLevels(spineLevels < 5) = 0;
spineLevels(spineLevels > 8) = 0;
% spineLevels(spineLevels > 0) = spineLevels(spineLevels > 0) - 4;

spineLevels(spineLevels>=1) = 1;

spine(1:71, :, :) = 2;
spine(71:end, :, :) = 1;

spine(spineLevels == 0) = 0;

cbiWriteNifti('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_all_rl.nii', spine, param);  
gzip('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_all_rl.nii')

%%
gunzip('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_cord.nii.gz');
[spine, ~] = cbiReadNifti('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_cord.nii');

gunzip('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii.gz');
[spineLevels, param] = cbiReadNifti('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii');
% spineLevels(spineLevels < 5) = 0;
spineLevels(spineLevels > 8) = 0;
% spineLevels(spineLevels > 0) = spineLevels(spineLevels > 0) - 4;

spineLevels(spineLevels >= 1) = 1;

spine(:, 1:71, :) = 2;
spine(:, 72:end, :) = 1;

spine(spineLevels == 0) = 0;

cbiWriteNifti('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_all_dv.nii', spine, param);  
gzip('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_all_dv.nii')

% %%
% gunzip('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii.gz');
% [spineLevels, param] = cbiReadNifti('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels.nii');
% % spineLevels(spineLevels < 5) = 0;
% spineLevels(spineLevels > 8) = 0;
% cbiWriteNifti('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels_5_8.nii', spineLevels, param);  
% gzip('D:\SMA\MRI_data_upper_limb\Spine\template\PAM50_spinal_levels_5_8.nii')
