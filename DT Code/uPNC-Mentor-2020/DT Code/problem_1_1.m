function sortedSpikes = problem_1_1(Data, neuronNumber)

clc; close all;
ov = [Data.Overview];
success = [ov.trialStatus] == 1;
data = Data(success);
timeOffset= 500;
% find any trials that do not contain spike data. Exclude these.
noData = false(size(data));
for trl = 1:length(data)
    noData(trl) = isempty(data(trl).TrialData.spikes) || ...
        isnan(data(trl).TrialData.timeMoveOnset);
end
data(noData) = [];

tempSpikes = [];
spikeCounter = [];
sortedSpikes.a = [];
sortedSpikes.b = [];
sortedSpikes.c = [];
sortedSpikes.d = [];
sortedSpikes.e = [];
sortedSpikes.f = [];
sortedSpikes.g = [];
sortedSpikes.h = [];
%%

angles = [0, 45, 90, 135, 180, 225, 270, 315];

for i= 1:length(data)
    Neurons = data(i).TrialData.spikes;
    beforeOnset = data(i).TrialData.timeMoveOnset - timeOffset;
    afterOnset = data(i).TrialData.timeMoveOnset + timeOffset;
    direction = data(i).TrialData.reachAngle;
    temp = Neurons(neuronNumber).timestamps;
        for k = 1:length(temp)
%             spikeCount(i,j) = sum((temp > beforeOnset&temp < afterOnset));
             if temp(k) >= beforeOnset && temp(k) <= afterOnset
                tempSpikes = [tempSpikes, double(temp(k)) - data(i).TrialData.timeMoveOnset];
             end
        end
     switch direction
         case angles(1)
             sortedSpikes.a = [sortedSpikes.a, tempSpikes];
         case angles(2)
             sortedSpikes.b = [sortedSpikes.b, tempSpikes];
         case angles(3)
             sortedSpikes.c= [sortedSpikes.c, tempSpikes];
         case angles(4)
             sortedSpikes.d = [sortedSpikes.d, tempSpikes];
         case angles(5)
             sortedSpikes.e =[sortedSpikes.e, tempSpikes];
         case angles(6)
             sortedSpikes.f =[sortedSpikes.f, tempSpikes];
         case angles(7)
             sortedSpikes.g = [sortedSpikes.g, tempSpikes];
         case angles(8)
             sortedSpikes.h =[sortedSpikes.h, tempSpikes];
         otherwise
     end
     spikeCounter = [spikeCounter, tempSpikes];
     tempSpikes= [];
end
%%
% for i=1:length(angles)
%     indexH = direction(
% end
%%
figure;
subplot(3,3,5)
histogram((spikeCounter), 'binedges', -timeOffset:20:timeOffset)
title(['Peristimulus Time Histogram for Neuron ', num2str(neuronNumber)])
axis([-500 500 0 max(spikeCounter)])
subplot(3,3,6)
histogram((sortedSpikes.a), 'binedges', -timeOffset:20:timeOffset)
title('Angle 0') 
axis([-500 500 0 max(spikeCounter)])
subplot(3,3,3)
histogram((sortedSpikes.b), 'binedges', -timeOffset:20:timeOffset)
title('Angle 45') 
axis([-500 500 0 max(spikeCounter)])
subplot(3,3,2)
histogram((sortedSpikes.c), 'binedges', -timeOffset:20:timeOffset)
title('Angle 90') 
axis([-500 500 0 max(spikeCounter)])
subplot(3,3,1)
histogram((sortedSpikes.d), 'binedges', -timeOffset:20:timeOffset)
title('Angle 135') 
axis([-500 500 0 max(spikeCounter)])
subplot(3,3,4)
histogram((sortedSpikes.e), 'binedges', -timeOffset:20:timeOffset)
title('Angle 180') 
axis([-500 500 0 max(spikeCounter)])
subplot(3,3,7)
histogram((sortedSpikes.f), 'binedges', -timeOffset:20:timeOffset)
title('Angle 225') 
axis([-500 500 0 max(spikeCounter)])
subplot(3,3,8)
histogram((sortedSpikes.g), 'binedges', -timeOffset:20:timeOffset)
title('Angle 270') 
axis([-500 500 0 max(spikeCounter)])
subplot(3,3,9)
histogram((sortedSpikes.h), 'binedges', -timeOffset:20:timeOffset)
title('Angle 315') 
axis([-500 500 0 max(spikeCounter)])





