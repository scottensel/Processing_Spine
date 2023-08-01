clear all
close all
clc

% load('savedData/allDataPAMcope7.mat')
load('savedData/allDataSmoothPAMcope7.mat')

% load('savedData/allDataRLcope7.mat')
% load('savedData/allDataSmoothRLcope7.mat')
% 
% load('savedData/allDataDVcope7.mat')
% load('savedData/allDataSmoothDVcope7.mat')

allData = allDataSmooth;

for sbj = 1:length(allData)
    for i = 1:5
        for cer = 1:4
            actVx(i,sbj) = sum(cell2mat(allData{sbj, 1}{i, 1}(1,:)));
            zSc{i,sbj}{cer,1} = allData{sbj, 1}{i, 1}{2,cer};
            zScMn{i,sbj}{cer,1} = mean(zSc{i,sbj}{cer,1});
            runs(i,1) = str2num(allData{1, 1}{i, 2});
        end

    end

end
clear sbj
%actVx = actVx(:,[1 2 4 5 6 7]);

%avRuns = [runs actVx];
totVox = 34947;
% totVox = 35133;

%% 
% %sub = 5;
% twoR = zScMn([1 6 10 13 15],:);
% threeR = zScMn([2 7 11 14],:);
% fourR = zScMn([3 8 12],:);
% fiveR = zScMn([4 9],:);
% sixR = zScMn(5,:);
% 
twoR = zScMn(1,:);
threeR = zScMn(2,:);
fourR = zScMn(3,:);
fiveR = zScMn(4,:);
sixR = zScMn(5,:);
% clear i
% runsCom = {};
% for i = 1:5
%     runsCom = zScMn(i,:);
% end

%%
combs = fourR;
% combs = sixR;
for cer = 1
    for sub = 1:7
        %cerAll{cer,sub} = [combs{1, sub}{cer, 1}];
        c4(sub,cer) = [combs{cer, sub}{1, 1}];  %; combs{c, 2}{1, 1}; combs{c, 3}{1, 1}; combs{c, 4}{1, 1}; combs{c, 5}{1, 1}; combs{c, 6}{1, 1}; combs{c, 7}{1, 1}];
        c5(sub,cer) = [combs{cer, sub}{2, 1}];
        c6(sub,cer) = [combs{cer, sub}{3, 1}];
        c7(sub,cer) = [combs{cer, sub}{4, 1}];
        
    end
end


%%
% close all
figure
c4mn = mean(c4, 2,'omitnan');
c5mn = mean(c5, 2,'omitnan');
c6mn = mean(c6, 2,'omitnan');
c7mn = mean(c7, 2,'omitnan');
%c7 = mean(cell2mat(cerAll(4,:)),'omitnan');

%sp = [cell2mat(cerAll(4,:)); cell2mat(cerAll(3,:)); cell2mat(cerAll(2,:)); cell2mat(cerAll(2,:));];
%grp = [ones(size(cell2mat(cerAll(4,:)))); 2.*ones(size(cell2mat(cerAll(3,:)))); 3.*ones(size(cell2mat(cerAll(2,:)))); 4.*ones(size(cell2mat(cerAll(1,:))))];
sp = [c7;c6;c5;c4];
grp = [ones(size(c7)); 2.*ones(size(c6)); 3.*ones(size(c5)); 4.*ones((size(c4)))];
col = linspace(1,10,length(sp));
allMeans = [mean(c7,'omitnan'),mean(c6,'omitnan'),mean(c5,'omitnan'),mean(c4,'omitnan')];
%allMeans = [mean(cell2mat(cerAll(4,:)),'omitnan'),mean(cell2mat(cerAll(3,:)),'omitnan'),mean(cell2mat(cerAll(2,:)),'omitnan'),mean(cell2mat(cerAll(1,:)),'omitnan')];
% light blue 4R [0.3010 0.7450 0.9330] %dark blue 6R [0 0.4470 0.7410]
barh(allMeans,'FaceColor',[0.3010 0.7450 0.9330]);hold on  %4 runs
%barh(allMeans,'FaceColor',[0 0.4500 0.750]);hold on  % All runs
sub=scatter(sp, grp,20,'k','o','filled'); 
xlim([0 8])
%s=scatter(grp,sp,[],'r','.'); 

sub.YJitter = 'rand'; sub.YJitterWidth = 0.5;
yticks(1:length(1:4)); yticklabels({'C7','C6','C5','C4'})
ylabel('Cervical level')
xlabel('z Score');
title(sprintf('Average z-Scores 4 Runs combined'));
%savefig('SBSN_figures/spLvl_zScore_4R_all_sub.fig')


