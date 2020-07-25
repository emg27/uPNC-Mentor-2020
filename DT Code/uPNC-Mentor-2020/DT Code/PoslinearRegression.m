function predictedValues =PoslinearRegression(data, spike, pos,wantComputerToDie)
firingRates = data(:,2:end);
numN = size(firingRates, 2);
numTrials = size(firingRates,1);
X = [];
y = []; 
xTest = [];
yTest = [];
plotXpos = [];
plotYpos = [];
for i=1:numTrials
    temp = [pos.x{i}(1:end-1) pos.y{i}(1:end-1)];
    if (max(max(isnan(temp)))==1)
        continue;
    end
    if i <= floor(numTrials*0.8)
        X = [X;spike{i}];
        y= [y;temp];
    else
        xTest = [xTest;spike{i}];
        yTest = [yTest;temp];
        plotXpos = [plotXpos;pos.x{i}(1:end-1)];
        plotYpos = [plotYpos;pos.y{i}(1:end-1)];
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
        spikeIndex = length(spike{i});
        differencePosX(trial) = mean(abs(plotXpos(spikeIndex) - predictedTraj(spikeIndex,1)));
        differencePosY(trial) = mean(abs(plotYpos(spikeIndex) - predictedTraj(spikeIndex,2)));
    else
        spikeIndex = spikeIndex +length(spike{i});
        differencePosX(trial) = mean(abs(plotXpos(spikeIndex) - predictedTraj(spikeIndex,1)));
        differencePosY(trial) = mean(abs(plotYpos(spikeIndex) - predictedTraj(spikeIndex,2)));
    end
end 
predictedValues.Distanceperformance = mean([differencePosX; differencePosY]);
%%
errorX = mean((abs((plotXPos-predictedTraj(:,1)))/plotXPos));
errorY = mean((abs((plotYPos-predictedTraj(:,2)))/plotYPos));
predictedValues.Errorperformance = mean([errorX errorY]);