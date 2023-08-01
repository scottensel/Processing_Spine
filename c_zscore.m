clear all
close all
clc

% load('savedData/allDataPAMcope1.mat')
load('savedData/allDataSmoothPAMcope7.mat')


% load('savedData/allDataRLcope7.mat')
load('savedData/allDataSmoothRLcope7.mat')
% 
% load('savedData/allDataDVcope7.mat')
% load('savedData/allDataSmoothDVcope7.mat')
allData = allDataSmooth;

for sbj = 1:length(allData)
    for i = 1:5
        for cer = 1:2
            actVx(i,sbj) = sum(cell2mat(allData{sbj, 1}{i, 1}(1,:)));
            zSc{i,sbj}{cer,1} = allData{sbj, 1}{i, 1}{2,cer};
            zScMn{i,sbj}{cer,1} = mean(zSc{i,sbj}{cer,1});
            runs(i,1) = str2num(allData{1, 1}{i, 2});
        end
    end
end

clear sbj

totVox = 35133;
totDV = [18594 16539];

%%

zThresh = 1.5;
for sbj = 1:length(allData)
    for i = 1:5
        avgZs(i,sbj) = mean(allData{sbj, 1}{i, 3}(find(allData{sbj, 1}{i, 3}>zThresh)));
        allavgZsc{i,sbj} = allData{sbj, 1}{i, 3}(find(allData{sbj, 1}{i, 3}>zThresh));
        
    end
end


%%
figure
allMeans = mean(avgZs,2);
clear i
for i  = 1: size(avgZs,1)
    SEM(i) = std(avgZs(i,:))/sqrt(length(avgZs(i,:)));
end
for i = 1:length(allMeans)
    errorbar(i*ones(size(allMeans(i))), allMeans(i), SEM(i),'LineWidth',1.5,'Color','0.5 0.5 0.5'); hold on
    %scatter(i*ones(size(allMeans(i))), allMeans(i),75,'red','square','filled');
end

bar(allMeans,'FaceColor',[.7 .7 .7]);
xticks(1:length(allMeans))
xticklabels({'1-2','1-2-3','1-2-3-4','1-2-3-4-5','All'}) 
ylabel('Average z Score'); xlabel('Number of Runs');
title('Combination of runs staring with 1: All subjects');


%% commented from here
figure
% close all
trii=diff(avgZs);
l= 4; % young
green = [0.4660 0.6740 0.1880];
light_green = [0.5000 0.950 0.750];
colors_p = [linspace(green(1),light_green(1),l)', linspace(green(2),light_green(2),l)', linspace(green(3),light_green(3),l)'];

for i =[1 2 3 4]  %1:4
    plot(trii(:,i),'o','MarkerSize',8,'MarkerEdgeColor',colors_p(i,:),'MarkerFaceColor',colors_p(i,:)); hold on
    plot(trii(:,i),'Color',colors_p(i,:))
end


l2= 3; % old
black = [0 0 0];
light_black = [0.65 0.65 0.65];
colors_p2 = [linspace(black(1),light_black(1),l2)', linspace(black(2),light_black(2),l2)', linspace(black(3),light_black(3),l2)'];

for j = 1:3
    plot(trii(:,j+4),'o','MarkerSize',8,'MarkerEdgeColor',colors_p2(j,:),'MarkerFaceColor',colors_p2(j,:)); hold on
    plot(trii(:,j+4),'Color',colors_p2(j,:))
end

xlim([0.5 4.5]);
%ylim([-1.75 3.75]);

xticks(1:4);
xticklabels({'2-3','3-4','4-5','5-6'})

legend('Y-1','','Y-2','','Y-3','','Y-4','','O-7','','O-8','','O-10','Location','eastoutside') 
% legend('Y-1','','Y-2','','Y-4','','O-7','','O-8','','O-10') %expt 3
%legend('Y-1','','Y-2','','Y-4') %young
%legend('O-7','','O-8','','O-10') %old 

%legend('Sub 1','','Sub 2','','Sub 4','','Sub 7','','Sub 8','','Sub 10') 
xlabel('Runs involved'); ylabel('Diff in percentage of mean z scores');
title('Diff in mean z scores in incremental runs')

