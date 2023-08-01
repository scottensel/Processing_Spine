clear all
% close all
% clc
% 
% load('savedData/allDataPAMcope1.mat')
% load('savedData/allDataSmoothPAMcope7.mat')

% load('savedData/allDataRLcope7.mat')
load('savedData/allDataSmoothRLcope7.mat')
% 
% load('savedData/allDataDVcope7.mat')
% load('savedData/allDataSmoothDVcope7.mat')
allData = allDataSmooth;

totVox = 35133;
% totDV = [18594 16539];
%RL
totDV = [18514 16619];

for sbj = 1:length(allData)
    for i = 1:5
        for cer = 1:2
            actVx(i,sbj) = sum(cell2mat(allData{sbj, 1}{i, 1}(1,:)));
            zSc{i,sbj}{cer,1} = allData{sbj, 1}{i, 1}{2,cer};
            acVxCer{i,sbj}{cer,1} = length(zSc{i, sbj}{cer, 1});
            acVxCer{i,sbj}{cer,2} = totDV(cer); 
            acVxCer{i,sbj}{cer,3} = (length(zSc{i, sbj}{cer, 1})/totDV(cer)*100);
            runs(i,1) = str2num(allData{1, 1}{i, 2});
        end
    end
end

clear sbj


%% Using only series starting with 1
twoR = actVx(1,:);
threeR = actVx(2,:);
fourR = actVx(3,:);
fiveR = actVx(4,:);
sixR = actVx(5,:);

percR = zeros(5,7);
percR(1,:) = (twoR/totVox)*100;    %(mean((twoR/totVox)*100))';
percR(2,:) = (threeR/totVox)*100;                           %(mean((threeR/totVox)*100))';
percR(3,:) = (fourR/totVox)*100;                           %(mean((fourR/totVox)*100))';
percR(4,:) = (fiveR/totVox)*100;                           %(mean((fiveR/totVox)*100))';
percR(5,:) = (sixR/totVox)*100;                           %(mean((sixR/totVox)*100))';
                 
allMeans= mean(percR,2)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plotting code - BE CAREFUL

% 
% close all
% figure
% clear i
% for i  = 1: size(percR,1)
%     SEM(i) = std(percR(i,:))/sqrt(length(percR(i,:)));
% end
% 
% for i = 1:length(allMeans)
%     errorbar(i*ones(size(allMeans(i))), allMeans(i), SEM(i),'LineWidth',1.5,'Color','0.5 0.5 0.5'); hold on
%     %scatter(i*ones(size(allMeans(i))), allMeans(i),75,'red','square','filled');
% end
% 
% hold on
% bar(allMeans,'FaceColor',[.7 .7 .7]);
% xticks(1:length(allMeans))
% xticklabels({'1-2','1-2-3','1-2-3-4','1-2-3-4-5','All'})
% ylabel('Percentge of Active voxels'); xlabel('Number of Runs');
% title('Percentage of active voxels per combination of runs averaged across all subjects')
% % savefig('SBSN_figures/actVox_all_sub.fig')

%% spinal segment active voxels DON'T HARD CODE
% close all
avRuns = [runs actVx];
a = [];
for sbj = 1:7
    for i = 3
        % avg cer activations 
        a(:, sbj) = cell2mat(acVxCer{i, sbj}(:, end));
        for c = 1:2
            cerAvg(c,:)= mean(a(c,:));
        end

    end
end


% %
% figure
% % close all
% h = bar(cerAvg,1,'FaceColor', [.85 .85 .85]); %4 runs
% %h = bar(flip(cerAvg,1),'FaceColor',[.35 .35 .35]); %all runs
% title (sprintf('All subjects - average : %d runs combined',i+1))
% xticks(1:2)
% ylim([1 16])
% xticklabels({'Dorsal','Ventral'})
% % xticklabels({'Right','Left'})
% ylabel('Percentge of Active voxels/side'); xlabel('Vertebral level');

figure
% close all
h = bar(cerAvg,1,'FaceColor', [.85 .85 .85]); %4 runs
%h = bar(flip(cerAvg,1),'FaceColor',[.35 .35 .35]); %all runs
title (sprintf('All subjects - average : %d runs combined',i+1))
xticks(1:2)
ylim([1 16])
xticklabels({'Left','Right'})
ylabel('Percentge of Active voxels/side'); xlabel('Vertebral level');
