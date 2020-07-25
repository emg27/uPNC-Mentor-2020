function predictedValues =VelolinearRegression(data, spike,velo,wantComputerToDie)
firingRates = data(:,2:end);
numN = size(firingRates, 2);
numTrials = size(firingRates,1);
X = [];
y = []; 
xTest = [];
yTest = [];
plotXVelo = [];
plotYVelo = [];
for i=1:numTrials
    temp = [velo.x{i}(1:end-1) velo.y{i}(1:end-1)];
    if (max(max(isnan(temp)))==1)
        continue;
    end
    if i <= floor(numTrials*0.8)
        X = [X;spike{i}];
        y= [y;temp];
    else
        xTest = [xTest;spike{i}];
        yTest = [yTest;temp];
        plotXVelo = [plotXVelo;velo.x{i}(1:end-1)];
        plotYVelo = [plotYVelo;velo.y{i}(1:end-1)];
    end
    
end


%%
result =(X' *X)^(-1)*(X'*y);
predictedTraj = (xTest*(result));


%%
if wantComputerToDie
    plot(plotXVelo, plotYVelo)
    hold on
    plot(predictedTraj(:,1), predictedTraj(:,2), 'Color', 'r')
    hold off
end
%%
errorX = mean((abs((plotXVelo-predictedTraj(:,1)))/plotXVelo));
errorY = mean((abs((plotYVelo-predictedTraj(:,2)))/plotYVelo));
predictedValues.Errorperformance = mean([errorX errorY]);

