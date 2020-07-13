function problem4_11(data)
clc; close all;
out = pmdDataSetup(data);
angles = out(:,1);
firingRates = out(:,2:end);
numN = size(firingRates, 2);
numTrials = size(firingRates,1);
x = [ones(numTrials,1), firingRates];
index = floor(.8*numTrials);
trainX = x(1:index, :);
trainAngles = angles(1:index,:);
testX = x(index+1:end,:);
testAngles = angles(index+1:end,:);
% result = (trainX'*trainX)^(-1)*trainX'*trainAngles
result = regress(trainAngles, trainX)
AngularError= abs((testX*result)-testAngles);
figure;
histogram(AngularError)
title('Histogram of the Absolute Angular Error')
