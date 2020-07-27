function predictedValues =linearRegression(spikes, posData, veloData,wantComputerToDie, independent)
numTrials = size(spikes,2);
X = [];
y = []; 
xTest = [];
yTest = [];
plotXpos = [];
plotYpos = [];
testingCoeff = 0.80;
testingsize = floor(length(spikes)*testingCoeff);
while rem(testingsize,8)~=0
    testingsize = testingsize +1;
end
rng('default')
rng(independent)
sortedIndex = randperm(numTrials);

for i=1:numTrials
    temp = [posData.x{sortedIndex(i)} posData.y{sortedIndex(i)} veloData.x{sortedIndex(i)} veloData.y{sortedIndex(i)}];
    if (max(max(isnan(temp)))==1||max(max(isnan(spikes{sortedIndex(i)}))))
        continue;
    end
    if i <= testingsize
        X = [X;spikes{sortedIndex(i)}];
        y= [y;temp];
    else
        xTest = [xTest;spikes{sortedIndex(i)}];
        yTest = [yTest;temp];
        plotXpos = [plotXpos;posData.x{sortedIndex(i)}];
        plotYpos = [plotYpos;posData.y{sortedIndex(i)}];
    end
    
end


%%
result =(X' *X)^(-1)*(X'*y);
predictedTraj = (xTest*(result));


%%
if wantComputerToDie
    plot(plotXpos, plotYpos)
    hold on
    plot(predictedTraj(:,1), predictedTraj(:,2), 'Color', 'r')
    hold off
end
%%
differencePosX = zeros(length(numTrials));
differencePosY = zeros(length(numTrials));
for trial= 1:length(numTrials)
    if trial ==1
        spikeIndex = size(spikes{i},1);
        differencePosX(trial) = mean(abs(plotXpos(spikeIndex) - predictedTraj(spikeIndex,1)));
        differencePosY(trial) = mean(abs(plotYpos(spikeIndex) - predictedTraj(spikeIndex,2)));
    else
        spikeIndex = spikeIndex +length(spikes{i});
        differencePosX(trial) = mean(abs(plotXpos(spikeIndex) - predictedTraj(spikeIndex,1)));
        differencePosY(trial) = mean(abs(plotYpos(spikeIndex) - predictedTraj(spikeIndex,2)));
    end
end 
predictedValues.Distanceperformance = mean([differencePosX; differencePosY]);
errorX = mean((abs((plotXpos'-predictedTraj(:,1)')/plotXpos')))*100;
errorY = mean((abs((plotYpos'-predictedTraj(:,2)')/plotYpos')))*100;
predictedValues.Errorperformance = mean([errorX errorY]);