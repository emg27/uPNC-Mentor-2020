clc; close all; clear;
load('Ike20120905_processed_original.mat');
data = pmdDataSetup(Data);
[spikes,posData,veloData] = PosVeloTime(Data);
for neuronNumber = [20 135 201]
    sortedSpikes = problem_1_1(Data, neuronNumber, true);
    outputProblem2 = problem_2(data, neuronNumber, true);
end 
problem_2_2(data)
Problem3_1(data)
Problem3_2(data)
problem4_11(data)
linearRegression(data, spikes, posData, veloData,true)
kalmanFilterDecoder(data, veloData, posData, spikes, true)
