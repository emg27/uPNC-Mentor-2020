% ask for data file
clear; clc; close all;
filename = input('Please input the name of the file here: ', 's');
tic
%load data file with name Data
load(filename)
BC = contains(Data(1).Overview.trialName,'BC','IgnoreCase',true);
if BC
	[spikes,fiftyPos,fiftyVelo] = PosVeloTimeBC(Data);
else
    [spikes,fiftyPos,fiftyVelo] = PosVeloTime(Data);
    data = pmdDataSetup(Data);
end
%%

for i=1:5
neuronLength = length(spikes{1}(2,:));
veloData = fiftyVelo;
posData =fiftyPos;
tempSpikes = spikes;
%%
kalmanErrorPlotting = zeros(neuronLength,2);
kalmanDistancePlotting = zeros(neuronLength,2);

kalPosErrorPlotting = zeros(neuronLength,2);
kalPosDistancePlotting = zeros(neuronLength,2);

kalVeloErrorPlotting = zeros(neuronLength,2);

linearErrorPlotting = zeros(neuronLength,2);
linearDistancePlotting = zeros(neuronLength,2);

linPosErrorPlotting = zeros(neuronLength,2);
linPosDistancePlotting = zeros(neuronLength,2);

linVeloErrorPlotting = zeros(neuronLength,2);
%%
for neuron =neuronLength:-1:226
    if neuron~=neuronLength
        dropppedNeuron = randi([1 size(tempSpikes{neuron},2)],1);
        for trial =1:length(tempSpikes)
            tempSpikes{trial}(:,dropppedNeuron)=[];
        end
    end
    b = randi([1 10000], 1);
    %plot performance vs neurons dropped
    if BC
        outputKalman = kalmanFilterDecoderBC(veloData, posData, tempSpikes, false,b);
        outputKalPos = PosKalmanFilterDecoderBC(posData, tempSpikes,b);
        outputKalVelo = VelokalmanFilterDecoderBC(veloData, tempSpikes,b);
        outputLinear= linearRegression(tempSpikes, posData, veloData,false,b);
        outputLinearPos = PoslinearRegression(tempSpikes, posData,false,b);
        outputLinearVelo = VelolinearRegression(tempSpikes,veloData,false,b);
    else
        outputKalman = kalmanFilterDecoder(data, veloData, posData, tempSpikes, false,b);
        outputKalPos = PosKalmanFilterDecoder(data,posData, tempSpikes,b);
        outputKalVelo = VelokalmanFilterDecoder(data, veloData, tempSpikes,b);
        outputLinear= linearRegression(tempSpikes, posData, veloData,false,b);
        outputLinearPos = PoslinearRegression(tempSpikes, posData,false,b);
        outputLinearVelo = VelolinearRegression(tempSpikes,veloData,false,b);
    end
    %%
    kalmanErrorPlotting(neuron, :) = [neuron outputKalman.Errorperformance];
    kalmanDistancePlotting(neuron, :) = [neuron outputKalman.Distanceperformance];
    %%
    kalPosErrorPlotting(neuron, :) = [neuron outputKalPos.Errorperformance];
    kalPosDistancePlotting(neuron, :) = [neuron outputKalPos.Distanceperformance];
    %%
    kalVeloErrorPlotting(neuron, :) = [neuron outputKalVelo.Errorperformance];
    %%
    linearErrorPlotting(neuron, :) = [neuron outputLinear.Errorperformance];
    linearDistancePlotting(neuron, :) = [neuron outputLinear.Distanceperformance];
    %%
    linPosErrorPlotting(neuron, :) = [neuron outputLinearPos.Errorperformance];
    linPosDistancePlotting(neuron, :) = [neuron outputLinearPos.Distanceperformance];
    %%
    linVeloErrorPlotting(neuron, :) = [neuron outputLinearVelo.Errorperformance];
end
%%

%%
plot(kalmanErrorPlotting(:,1), kalmanErrorPlotting(:,2), 'Color', 'k')
hold on
plot(kalPosErrorPlotting(:,1), kalPosErrorPlotting(:,2), 'Color', 'r')
xlabel('Number of Neurons')
ylabel('Percent Error Performance')
title('The Relationship of Neuron Dropping to Error Performance for the Kalman Filter')
set(gca, 'XDir','reverse')
set(gca, 'YDir','reverse')
grid
hold off
figure;
plot(kalVeloErrorPlotting(:,1), kalVeloErrorPlotting(:,2), 'Color', 'b')
xlabel('Number of Neurons')
ylabel('Percent Error Performance')
title('The Relationship of Neuron Dropping to Error Performance for the Kalman Filter')
set(gca, 'XDir','reverse')
set(gca, 'YDir','reverse')
figure; 
hold on
plot(kalmanDistancePlotting(:,1), kalmanDistancePlotting(:,2), 'Color', 'k')
plot(kalPosDistancePlotting(:,1), kalPosDistancePlotting(:,2), 'Color', 'r')
xlabel('Number of Neurons')
ylabel('Distance Performance')
title('Distance Performance of the Kalman filter as Neurons are Dropped')
set(gca, 'XDir','reverse')
set(gca, 'YDir','reverse')
grid
hold off
%%
figure;
hold on
plot(linearErrorPlotting(:,1), linearErrorPlotting(:,2), 'Color', 'k')
plot(linPosErrorPlotting(:,1), linPosErrorPlotting(:,2), 'Color', 'r')
plot(linVeloErrorPlotting(:,1), linVeloErrorPlotting(:,2), 'Color', 'b')
xlabel('Number of Neurons')
ylabel('Percent Error Performance')
title('Error Performance of the Linear Regression Decoder as Neurons are Dropped')
set(gca, 'XDir','reverse')
set(gca, 'YDir','reverse')
hold off
figure;
hold on
plot(linearDistancePlotting(:,1), linearDistancePlotting(:,2))
plot(linPosDistancePlotting(:,1), linPosDistancePlotting(:,2))
xlabel('Number of Neurons')
ylabel('Distance Performance')
title('The Relationship of Neuron Dropping to Distance Performance using Linear Regression')
set(gca, 'XDir','reverse')
set(gca, 'YDir','reverse')
hold off
end
toc