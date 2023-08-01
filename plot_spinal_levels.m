% y = rand(1,10); % your mean vector;
% x = 1:numel(y);
% std_dev = 1;
% 
% curve1 = y + std_dev;
% curve2 = y - std_dev;
% x2 = [x, fliplr(x)];
% 
% inBetween = [curve1, fliplr(curve2)];
% fill(x2, inBetween, 'g','FaceAlpha',0.5);
% hold on;
% plot(x, y, 'r', 'LineWidth', 2);
% 
clear all

% load('savedData/allDataSLcope1.mat')
load('savedData/allDataSmoothSLcope7.mat')
allData = allDataSmooth;

for sbj = 1:length(allData)
    allVoxels(sbj, 1) = allData{sbj, 1}{3, 1};
end

for i = 1:length(allVoxels)
    cAll(:, i) = allVoxels{i, 1};

end

% cAM = mean(cAll(735:870,:), 2);
% cAS = std(cAll(735:870,:), 0, 2)/sqrt(7);
cAll = smoothdata(cAll, 'gaussian', 20);

cAM = mean(cAll, 2);
cAS = std(cAll, 0, 2)/sqrt(7);

% figure
% plot(c5)
%             717:888
figure
% plot(c5M)
% hold on
% plot(c6M)
% plot(c7M)
% plot(c8M)

% add smoothing to the lines to create it better
curve1 = cAM + cAS;
curve2 = cAM - cAS;
h=area([curve2,curve1]);
set(h(1),'FaceColor','w','EdgeColor','none', 'FaceAlpha',0.2)
hold on
plot(cAM)
ylim([0 165])
xlim([726 880])


