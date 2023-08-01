% clear all
% close all
% clc
% 
% load('savedData/allDataPAMcope1.mat')
load('savedData/allDataSmoothPAMcope1.mat')


% load('savedData/allDataRLcope7.mat')
% load('savedData/allDataSmoothRLcope7.mat')
% 
% load('savedData/allDataDVcope7.mat')
% load('savedData/allDataSmoothDVcope7.mat')

allData = allDataSmooth;

totVox = 34947;
totCer = [10750 9438 7870 6889];
% totVox = 35133;
% totCer = [18594 16539];

for sbj = 1:length(allData)
    for i = 1:5
        for cer = 1:4 % 1:2
%         for cer = 1:length(totCer)
            actVx(i,sbj) = sum(cell2mat(allData{sbj, 1}{i, 1}(1,:)));
            zSc{i,sbj}{cer,1} = allData{sbj, 1}{i, 1}{2,cer};
            acVxCer{i,sbj}{cer,1} = length(zSc{i, sbj}{cer, 1});
            acVxCer{i,sbj}{cer,2} = totCer(cer); 
            acVxCer{i,sbj}{cer,3} = (length(zSc{i, sbj}{cer, 1})/totCer(cer)*100);
            runs(i,1) = str2num(allData{1, 1}{i, 2});
        end

    end

end

clear sbj


%%
avRuns = [runs actVx];

% % %% OLD combs
% % %sub = 7;
% 
% twoR = actVx([1 10 14 17 19],:);
% threeR = actVx([2 11 15 18],:);
% fourR = actVx([3 12 16],:);
% fiveR = actVx([4 13],:);
% sixR = actVx(5,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% close all
figure;
clear i
for i  = 1:size(percR,1)
    SEM(i) = std(percR(i,:))/sqrt(length(percR(i,:)));
end

for i = 1:length(allMeans)
    errorbar(i*ones(size(allMeans(i))), allMeans(i), SEM(i),'LineWidth',1.5,'Color','0.5 0.5 0.5'); hold on
    %scatter(i*ones(size(allMeans(i))), allMeans(i),75,'red','square','filled');
end

hold on
bar(allMeans,'FaceColor',[.7 .7 .7]);
xticks(1:length(allMeans))
xticklabels({'1-2','1-2-3','1-2-3-4','1-2-3-4-5','All'})
ylabel('Percentge of Active voxels'); xlabel('Number of Runs');
% title('Percentage of active voxels per combination of runs averaged across all subjects')
% % savefig('SBSN_figures/actVox_all_sub.fig')

%% spinal segment active voxels DON'T HARD CODE
% close all
avRuns = [runs actVx];
% totVox = 34947;
% totCer = [10750 9438 7870 6889];

a = [];
for sbj = 1:7
    for i = 3
        %figure
        %bar(cell2mat(acVxCer{i,sbj}(:,1)));
        %title (sprintf('SBSN-H-00%d : %d runs combined',sbj,i+1)
        %saveas(h,sprintf('SBSN-H-00%d_%druns.png',sbj,i));
        %saveas(h,sprintf('SBSN-H-00%d_%druns.fig',sbj,i));
        % avg cer activations 
        a(:,sbj) = cell2mat(acVxCer{i,sbj}(:,end));
%         for c = 1:2
        for c = 1:length(totCer)
            %cerAvg4R(c,:)= mean(a(c,:));
            cerAvg(c,:)= mean(a(c,:));
%             h = bar(cerAvg);
%             title (sprintf('All subjects - average : %d runs combined',i+1))
        end

    end
end


%%
% close all
figure;
% h = barh(flip([cerAvg4R,cerAvg6R],1));
% legend('4R','6R','Location','eastoutside')
%h = barh(flip(cerAvg,1),'FaceColor',[.85 .85 .85]); %4 runs
h = barh(flip(cerAvg,1),'FaceColor',[.35 .35 .35]); %all runs
title (sprintf('All subjects - average : %d runs combined', i+1))
yticks(1:4)
xlim([1 24])
yticklabels({'C7','C6','C5','C4'})
xlabel('Percentge of Active voxels/segment'); ylabel('Vertebral level');
%saveas(h,sprintf('AllSub_SegAvg%druns.png',i+1));
%saveas(h,sprintf('AllSub_SegAvg%druns.fig',i+1));
%title('Percentage of active voxels per vertebral level - 4R,')


