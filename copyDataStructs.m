startPath = "P:/projects/human/SCS_SMA/DATA/";
endPathCopy = 'D:\SMA\DATA';

subject = ["SMA01", "SMA02","SMA03"];
dates = ["20220919","20220926","20221003","20221011";
    "20230227","20230307","20230313","20230320";
    "20230522","20230529","20230605","20230612"];

endPath = '/DAQ_Data/structs/';

top = dir(startPath);
parfor i = 1:length(subject)

    for j = 1:length(dates)

        copyfile(fullfile(startPath, subject(i), dates(i,j), endPath, '*.pickle'), fullfile(endPathCopy, subject(i), dates(i,j), endPath))

    end
end