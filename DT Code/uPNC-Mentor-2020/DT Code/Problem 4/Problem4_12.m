% function problem4_12(data)
clc; close all;
out = pmdDataSetup(Data);
firingRates = out(:,2:end);
numN = size(firingRates, 2);
numTrials = size(firingRates,1);
%%
[spike, pos, velo] = PosVeloTime(Data);


posX = [];
posY = [];
veloX = [];
veloY = []; 
toggle = false;
for i=1:1222
        if i== 1
            X = [spike{i}];
            y= [pos.x{i} pos.y{i} velo.x{i} velo.y{i}];
        else
            X = [X;spike{i}];
            posX = [posX;pos.x{i}];
            posY = [posY;pos.y{i}];
            veloX = [veloX;velo.x{i}];
            veloY = [veloY;velo.y{i}]; 

        end
end

y = [posX posY veloX veloY];
%%
index = floor(.8*numTrials);
trainX = X(1:index, :);
trainY= y(1:index,:);
testX = X(index+1:end,:);
testY= y(index+1:end,:);
% result = regress(trainY, trainX);
result =(X' *X)\X'*y;
AngularError= abs((testX*result)-testY);
figure;
histogram(AngularError)
title('Histogram of the Absolute Angular Error')

