% function Problem4_12(Data)
clc;
out = pmdDataSetup(Data);
firingRates = out(:,2:end);
numN = size(firingRates, 2);
numTrials = size(firingRates,1);
%%
[spike, pos, velo] = PosVeloTime(Data);

%%
X = [];
y = []; 
xTest = [];
yTest = [];
plotXpos = [];
plotYpos = [];
for i=1:numTrials
    temp = [pos.x{i}(1:end-1) pos.y{i}(1:end-1) velo.x{i}(1:end-1) velo.y{i}(1:end-1)];
    if (max(max(isnan(temp)))==1)
        continue;
    end
    if i <= floor(numTrials*0.8)
%             X = spike{i}(bin,:);
%             y = [pos.x{i}(1:end-1) pos.y{i}(1:end-1) velo.x{i}(1:end-1) velo.y{i}(1:end-1)];
%             result =(X' *X)^(-1)*(X'*y);
            result = regress([pos.x{i}(1:end-1) pos.y{i}(1:end-1) velo.x{i}(1:end-1) velo.y{i}(1:end-1)], spike{i});

%             predictedTraj = (yTest*(result)')/50;
        

%         X = [X;spike{i}];
%         y= [y;temp];
    else
%         xTest = [xTest;spike{i}];
%         yTest = [yTest;temp];
%         plotXpos = [plotXpos;pos.x{i}(:)];
%         plotYpos = [plotYpos;pos.y{i}(:)];
    end
    
end


%%
result = regress(y, X);
result =(X' *X)^(-1)*(X'*y);
predictedTraj = (yTest*(result)')/50;
plot(plotXpos, plotYpos)
hold on
plot(predictedTraj(:,3), predictedTraj(:,4), 'Color', 'r')
hold off