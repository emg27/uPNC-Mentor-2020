% ask for data file
clear; clc; close all;
filename = input('Please input the name of the file here: ', 's');
tic
%load data file with name Data
load(filename)
switch Data(1).Overview.trialName
    case contains(Data(1).Overview.trialName,'BC','IgnoreCase',true)
        [spikes,fiftyPos,fiftyVelo] = PosVeloTimeBC(Data);
        
    case contains(Data(1).Overview.trialName,'BC','IgnoreCase',true)
end

neuronLength = length(spikes{1}(2,:));
veloData = fiftyVelo;
posData =fiftyPos;
data = pmdDataSetup(Data);
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
for neuron =neuronLength:-1:1
    if neuron~=neuronLength
        dropppedNeuron = randi([1 size(spikes{trial},2)],1);
        for trial =1:length(spikes)
            spikes{trial}(:,dropppedNeuron)=[];
        end
    end
    %plot performance vs neurons dropped
    outputKalman = kalmanFilterDecoder(data, veloData, posData, spikes, false);
    outputKalPos = PosKalmanFilterDecoder(data,posData, spikes);
    outputKalVelo = VelokalmanFilterDecoder(data, veloData, spikes);
    outputLinear= linearRegression(data, spikes, pos, velo,false);
    outputLinearPos = PoslinearRegression(data, spikes, posData,false);
    outputLinearVelo = VelolinearRegression(data, spikes,veloData,false);
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
close all;
plot(kalmanErrorPlotting(:,1), kalmanErrorPlotting(:,2))
xlabel('Number of Neurons')
ylabel('Error Performance of the Kalman filter')
title('The Relationship of Neuron Dropping to Error Performance')
set(gca, 'XDir','reverse')
grid
figure;
plot(kalmanDistancePlotting(:,1), kalmanDistancePlotting(:,2), 'Color', 'r')
xlabel('Number of Neurons')
ylabel('Distance Performance of the Kalman filter')
title('The Relationship of Neuron Dropping to Distance Performance')
set(gca, 'XDir','reverse')
grid
%%
figure;
plot(kalPosErrorPlotting(:,1), kalPosErrorPlotting(:,2))
xlabel('Number of Neurons')
ylabel('Error Performance of the Kalman filter')
title('The Relationship of Neuron Dropping to Error Performance in Position Kalman filter')
set(gca, 'XDir','reverse')
grid
figure;
plot(kalPosDistancePlotting(:,1), kalPosDistancePlotting(:,2), 'Color', 'r')
xlabel('Number of Neurons')
ylabel('Distance Performance of the Position Kalman filter')
title('The Relationship of Neuron Dropping to Distance Performance in Pos')
set(gca, 'XDir','reverse')
grid
%%
figure;
plot(kalVeloErrorPlotting(:,1), kalVeloErrorPlotting(:,2))
xlabel('Number of Neurons')
ylabel('Error Performance of the Velocity Kalman filter')
title('The Relationship of Neuron Dropping to Error Performance in Velocity Kalman')
set(gca, 'XDir','reverse')
grid
%%
% plot(linearErrorPlotting(:,1), linearErrorPlotting(:,2))
% xlabel('Number of Neurons')
% ylabel('Error Performance of the linear regression decoder')
% title('The Relationship of Neuron Dropping to Error Performance')
% figure;
% plot(linearDistancePlotting(:,1), linearDistancePlotting(:,2))
% xlabel('Number of Neurons')
% ylabel('Distance Performance of the linear regression decoder')
% title('The Relationship of Neuron Dropping to Distance Performance')
%%
plot(linPosErrorPlotting(:,1), linPosErrorPlotting(:,2))
xlabel('Number of Neurons')
ylabel('Error Performance of the linear position decoder')
title('The Relationship of Neuron Dropping to Error Performance')
figure;
plot(linPosDistancePlotting(:,1), linPosDistancePlotting(:,2))
xlabel('Number of Neurons')
ylabel('Distance Performance of the linear position decoder')
title('The Relationship of Neuron Dropping to Distance Performance')
%%
plot(linVeloErrorPlotting(:,1), linVeloErrorPlotting(:,2))
xlabel('Number of Neurons')
ylabel('Error Performance of the velocity linear decoder')
title('The Relationship of Neuron Dropping to Error Performance')
toc
