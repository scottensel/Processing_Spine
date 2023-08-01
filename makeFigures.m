%% template matching
clear all

% %%% SPINE
addpath('D:\NHP_code\cbiNifti')
% 
load('allData.mat')
% varibales to set up before
% subName = {'SBSN\_H\_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_007','SBSN_H_008','SBSN_H_010'};
% subName = {'SBSN\_H\_001','SBSN_H_002','SBSN_H_003','SBSN_H_004','SBSN_H_007','SBSN_H_010'};
subName = {'Y1','Y2','Y4','O1','O2'};

% subName = {'SBSN_H_010','SBSN_H_008'};
zScore = 1.5;

vertLevels = [10750, 9438, 7870, 6889];
totalVoxels = sum(vertLevels);
%% now can loop through allData and make the plots
arrayData = {};
for i = 1:length(allData)
    arrayData{i,1} = {};
    
    counter = ones(5,1);
    for j = 1:length(allData{i,1})
        
        % get runs in run
        runNumber = split(allData{i,1}{j,2}, '');
        runNumber = runNumber(2:end-1);
        runNum = str2num(cell2mat(runNumber));
        
        idx = length(runNum) - 1;
        
        % only take consecutive runs
        if mean(diff(runNum)) == 1
        
            subZscores = {allData{i,1}{j,1}{2,:}};
            subZscores = cell2mat(subZscores');
            mZscores = mean(subZscores);
            
            arrayData{i,1}{idx, counter(idx, 1)} = mZscores;
            
            
            subActive = (mean(cell2mat({allData{i,1}{j,1}{1,:}}))/totalVoxels)*100;
            
            arrayData{i,3}{idx, counter(idx, 1)} = subActive;
            
            
            arrayData{i,5}{idx, counter(idx, 1)} = allData{i,1}{j,3};
            
            counter(idx, 1) = counter(idx, 1) + 1;
        end
    end
    
    shit = [];
    shit2 = [];
    for k = 1:5
        
        shit(k) = mean([arrayData{i,1}{k,:}]);
        shit2(k) = mean([arrayData{i,3}{k,:}]);
        

    end
    
    alldice = zeros(1,4);
    idx = [5, 4, 3 , 2];
    % doing all the dice coeff
    for j = 1:4
        shit3 = [];
        for k = 1:idx(j)
            for kk = 1:idx(j)

                shit3(k, kk) = dice(arrayData{i,5}{j,k}>zScore, arrayData{i,5}{j,kk}>zScore);

            end
        end
        shit3(logical(tril(ones(size(shit3))))) = nan;
        totalDice = mean(mean(shit3 , 'omitnan'),'omitnan');
        alldice(1,j) = totalDice;
    end
    
    
    arrayData{i, 2} = shit;
    arrayData{i, 4} = shit2;
    arrayData{i, 6} = alldice;
end

for i = 1:length(arrayData)
    plotZscore(i,:) = [arrayData{i,2}];
    plotVoxels(i,:) = [arrayData{i,4}];
    plotDice(i,:) = [arrayData{i,6}];
end


%calculate dice coeff


plotZscore(6,:) = [];
plotZscore(3, :) = [];

figure;
h = plot(plotZscore', '--o');
xlim([0.5 5.5])
set(h(1), 'color', '#045275');
set(h(2), 'color', '#089099');
set(h(3), 'color', '#7CCBA2');
set(h(4), 'color', '#FEB24C');
set(h(5), 'color', '#F0746E');
% set(h(6), 'color', '#7C1D6F');
legend(subName,'location','best')
make_pretty
title('Z-Score')
%     plot((1:length(durationsAvg{i,2}(B)))+offset(i), 2.4*durationsAvg{i,2}(B), '-', 'linewidth', 2, 'Color', [colorList{i}, 0.25]);
%     hold on
%     plot((1:length(durationsAvg{i,2}(B)))+offset(i), 2.4*durationsAvg{i,2}(B), 'o', 'linewidth', 15, 'Color', colorList{i});


figure;
b = bar(mean(plotZscore)');
xlim([0.5 5.5])
ylim([2.5 3.5])
% make_pretty
hold on
ngroups = 5; 
nbars = 1;
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x', mean(plotZscore)', (std(plotZscore)/sqrt(length(plotZscore))'), 'k', 'linestyle', 'none');
% set(gca,'xticklabel', conditions)
title('Z-score')
make_pretty


% active voxels
plotVoxels(6,:) = [];
plotVoxels(3, :) = [];
% make this avg or all subjects instead of each
% with std error on them
figure;
b = bar(mean(plotVoxels)');
xlim([0.5 5.5])
% make_pretty
hold on
ngroups = 5; 
nbars = 1;
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x', mean(plotVoxels)', (std(plotVoxels)/sqrt(length(plotVoxels))'), 'k', 'linestyle', 'none');
% set(gca,'xticklabel', conditions)
title('Mean % Active Voxels')
make_pretty


plotDice(6,:) = [];
plotDice(3, :) = [];

figure;
h = plot(plotDice', '--o');
xlim([0.5 4.5])
set(h(1), 'color', '#045275');
set(h(2), 'color', '#089099');
set(h(3), 'color', '#7CCBA2');
set(h(4), 'color', '#FEB24C');

% set(h(5), 'color', '#F0746E');
% set(h(6), 'color', '#7C1D6F');
legend(subName,'location','best')
make_pretty
title('dice')

figure;
b = bar(mean(plotDice)');
xlim([0.5 4.5])
% ylim([2.5 3.5])
% make_pretty
hold on
ngroups = 4; 
nbars = 1;
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x', mean(plotDice)', (std(plotDice)/sqrt(length(plotDice))'), 'k', 'linestyle', 'none');
% set(gca,'xticklabel', conditions)
title('dice')
make_pretty

%% now can loop through allData and make the plots
vertData = {};
allActive = [];
for i = 1:length(allData)
    vertData{i,1} = {};
    
    counter = 1;
    for j = 1:length(allData{i,1})
        
        % get runs in run
        runNumber = split(allData{i,1}{j,2}, '');
        runNumber = runNumber(2:end-1);
        runNum = str2num(cell2mat(runNumber));
        
        idx = length(runNum) - 1;
        
        % only take consecutive runs
        if length(runNum) == 4 && mean(diff(runNum)) == 1
        
            subZscores = {allData{i,1}{j,1}{2,:}};
            
            vertData{i,1}{1, counter} = subZscores;
            
            subActive = cell2mat({allData{i,1}{j,1}{1,:}});
            
            allActive(size(allActive,1)+1, :) = subActive;
            
%             counter(idx, 1) = counter(idx, 1) + 1;
            counter = counter + 1;
        end
        
    end
    
end

counter = 1;
for i = 1:length(vertData)
   
   for j = 1:length(vertData{i})
       
       
       scores = [];
       for k = 1:4
           
           scores(k) = mean(vertData{i}{j}{k});
           
           
       end
       
       scores2(counter, :) = scores;
       counter = counter + 1;
   end
    
end

counter = 1;
for i = 1:3:length(scores2)
    
    actScores(counter,:) = mean(allActive(i:i+2, :),  'omitnan');
    
    scores3(counter,:) = mean(scores2(i:i+2, :),  'omitnan');
    counter = counter + 1;

end

% 
% scores3(6, :) = []; 
% scores3(3, :) = [];
stutter = ones(size(scores3));
for i = 1:4
 stutter(:,i) = (stutter(:,i)*i) + randi([-30,30], size(scores3,1), 1)./100;
end


figure;
b = bar(mean(scores3));
hold on
h = plot(stutter', scores3', '.k');
view([90 -90])
set ( gca, 'xdir', 'reverse' )
xlim([0.5 4.5])
ylim([2 4])
xticks([1 2 3 4])
xticklabels({'C4','C5','C6','C7'})
% set(h(1), 'color', '#045275');
% set(h(2), 'color', '#089099');
% set(h(3), 'color', '#7CCBA2');
% set(h(4), 'color', '#FEB24C');
% set(h(5), 'color', '#F0746E');
% set(h(6), 'color', '#DC3977');
% set(h(6), 'color', '#7C1D6F');
% legend(subName,'location','best')
make_pretty
% 
% % make this avg or all subjects instead of each

% % with std error on them
% actScores(6, :) = [];
% actScores(3, :) = [];

actScores2 = actScores./repmat(vertLevels, length(actScores), 1)*100;

figure;
b = bar(mean(actScores2));
hold on
h = plot(stutter', actScores2', '.k');
view([90 -90])
set ( gca, 'xdir', 'reverse' )
xlim([0.5 4.5])
xticks([1 2 3 4])
% set(h(1), 'color', '#045275');
% set(h(2), 'color', '#089099');
% set(h(3), 'color', '#7CCBA2');
% set(h(4), 'color', '#FEB24C');
% set(h(5), 'color', '#F0746E');
% set(h(6), 'color', '#DC3977');
% set(h(6), 'color', '#7C1D6F');
xticklabels({'C4','C5','C6','C7'})
% legend(subName,'location','best')
make_pretty

% 
% scores2(2, :) = []; 
% 
% figure;
% h = plot(scores2', '--o');
% view([90 -90])
% set ( gca, 'xdir', 'reverse' )
% xlim([0.5 4.5])
% xticks([1 2 3 4])
% xticklabels({'C4','C5','C6','C7'})
% set(h(1), 'color', '#045275');
% set(h(2), 'color', '#089099');
% set(h(3), 'color', '#7CCBA2');
% set(h(4), 'color', '#FEB24C');
% set(h(5), 'color', '#F0746E');
% % set(h(6), 'color', '#DC3977');
% set(h(6), 'color', '#7C1D6F');
% legend(subName,'location','best')
% make_pretty
% % 
% % % make this avg or all subjects instead of each
% 
% % % with std error on them
% allActive(6, :) = [];
% 
% figure;
% h = plot(allActive', '--o');
% view([90 -90])
% set ( gca, 'xdir', 'reverse' )
% xlim([0.5 4.5])
% xticks([1 2 3 4])
% set(h(1), 'color', '#045275');
% set(h(2), 'color', '#089099');
% set(h(3), 'color', '#7CCBA2');
% set(h(4), 'color', '#FEB24C');
% set(h(5), 'color', '#F0746E');
% % set(h(6), 'color', '#DC3977');
% set(h(6), 'color', '#7C1D6F');
% xticklabels({'C4','C5','C6','C7'})
% legend(subName,'location','best')
% make_pretty