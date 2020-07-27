function predictedValues =VelolinearRegression(spikes,velo,wantComputerToDie, independent)
numTrials = size(spikes,2);
X = [];
y = []; 
xTest = [];
yTest = [];
plotXVelo = [];
plotYVelo = [];
testingCoeff = 0.80;
testingsize = floor(length(spikes)*testingCoeff);
while rem(testingsize,8)~=0
    testingsize = testingsize +1;
end
rng('default')
rng(independent)
sortedIndex = randperm(numTrials);
for i=1:numTrials
    temp = [velo.x{sortedIndex(i)} velo.y{sortedIndex(i)}];
    if (max(max(isnan(temp)))==1)
        continue;
    end
    if i <= testingsize
        X = [X;spikes{sortedIndex(i)}];
        y= [y;temp];
    else
        xTest = [xTest;spikes{sortedIndex(i)}];
        yTest = [yTest;temp];
        plotXVelo = [plotXVelo;velo.x{sortedIndex(i)}];
        plotYVelo = [plotYVelo;velo.y{sortedIndex(i)}];
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
errorX = mean((abs((plotXVelo'-predictedTraj(:,1)')/plotXVelo')))*100;
errorY = mean((abs((plotYVelo'-predictedTraj(:,2)')/plotYVelo')))*100;
predictedValues.Errorperformance = mean([errorX errorY]);

