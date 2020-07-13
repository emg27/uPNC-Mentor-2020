function [spikes,fiftyPos,fiftyVelo] = PosVeloTime(Data)
clc; close all;
%%
% find successful trials. Only include these. 
ov = [Data.Overview];
success = [ov.trialStatus] == 1;
data = Data(success);
timeOffset = 100;
% find any trials that do not contain spike data. Exclude these.
noData = false(size(data));
for trl = 1:length(data)
    noData(trl) = isempty(data(trl).TrialData.spikes) || ...
        isnan(data(trl).TrialData.timeMoveOnset);
end
data(noData) = [];
%%
timers.movementOnset = {length(data)};
timers.movementStart = {length(data)};
timers.movementEnd = {length(data)};
fiftyTimers = timers;
pos.x = {length(data)};
pos.y = {length(data)};
fiftyPos = pos;
velo.x=  {length(data)};
velo.y = {length(data)};
spikes = {length(data)};
fiftyVelo = velo;
time = {length(data)};
refresh = 25/3;

%%
for i=1:length(data)
    time{i} = data(i).TrialData.Marker.rawPositions(:,6);
    timers.movementOnset{i} = round((data(i).TrialData.timeMoveOnset)/50)*50;
    timers.movementStart{i} = round((data(i).TrialData.timeMoveOnset-timeOffset)/50)*50;
    timers.movementEnd{i} = round((data(i).TrialData.timeMoveEnd)/50)*50;
    pos.x{i} = data(i).TrialData.Marker.rawPositions(:,2);
    pos.y{i} = data(i).TrialData.Marker.rawPositions(:,3);
    velo.x{i} =  (pos.x{i}(2:end)-pos.x{i}(1:end-1))/(refresh);
    velo.y{i} =  (pos.y{i}(2:end)-pos.x{i}(1:end-1))/(refresh);


end
%%
for i=1:length(data)
    tempor = timers.movementStart{i}:50:timers.movementEnd{i};
    for j=1:length(tempor)
        if j ==1
             index  = (time{i}>=tempor(j) &time{i}<=tempor(j)+50);
             fiftyPos.x{i} = mean(data(i).TrialData.Marker.rawPositions(index,2));
             fiftyPos.y{i} = mean(data(i).TrialData.Marker.rawPositions(index,3));
             fiftyVelo.x{i} = mean(velo.x{i}(index(1:end-1)));
             fiftyVelo.y{i} = mean(velo.y{i}(index(1:end-1))); 
        else
             index  = (time{i}>=tempor(j) &time{i}<=tempor(j)+50);
             fiftyPos.x{i} = [fiftyPos.x{i};mean(data(i).TrialData.Marker.rawPositions(index,2))];
             fiftyPos.y{i} = [fiftyPos.y{i};mean(data(i).TrialData.Marker.rawPositions(index,3))];
             fiftyVelo.x{i} = [fiftyVelo.x{i};mean(velo.x{i}(index(2:end)))];
             fiftyVelo.y{i} = [fiftyVelo.y{i};mean(velo.y{i}(index(2:end)))];
        end
    end
    spikeCount = zeros(length(tempor), length(data(i).TrialData.spikes));
    for neuron = 1:length(data(i).TrialData.spikes)
        spikeTimes = data(i).TrialData.spikes(neuron).timestamps;
        if isempty(spikeTimes)
            counts = zeros(length(tempor),1);
            spikeCount(:,neuron) = counts;
        else
            counts = histcounts(spikeTimes, 'binedges', [tempor ,tempor(end)]);
            spikeCount(:,neuron) = counts;
        end
    end

    spikes{i} = spikeCount;
end
%%

