subject = 'SBSN_H_010';

direc_start = 'D:\SBSN\Data\Spine\';


mkdir(fullfile(direc_start, subject, 'anat'))
mkdir(fullfile(direc_start, subject, 'func'))
mkdir(fullfile(direc_start, subject, 'physio'))
mkdir(fullfile(direc_start, subject, 'task'))
for i = 1:6
    
    mkdir(fullfile(direc_start, subject, 'func', ['func' num2str(i)]))
    mkdir(fullfile(direc_start, subject, 'physio', ['physio' num2str(i)]))
    mkdir(fullfile(direc_start, subject, 'task', ['task' num2str(i)]))

end
    