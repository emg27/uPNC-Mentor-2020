function [spikes,fiftyPos,fiftyVelo] = Problem4_0_1(Data)
clc;
ov = [Data.Overview];
success = [ov.trialStatus] == 1;
data = Data(success);
timeOffset = 100;
noData = false(size(data));
for trl = 1:length(data)
    noData(trl) = isempty(data(trl).TrialData.spikes) || ...
        isnan(data(trl).TrialData.timeMoveOnset);
end
data(noData) = [];

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
refresh = (25/3)/1000;


for i=1:length(data)
    time{i} = data(i).TrialData.Marker.rawPositions(:,6);
    timers.movementOnset{i} = data(i).TrialData.timeMoveOnset;
    timers.movementStart{i} = data(i).TrialData.timeMoveOnset-timeOffset;
    timers.movementEnd{i} = data(i).TrialData.timeMoveEnd;
    pos.x{i} = data(i).TrialData.Marker.rawPositions(:,2);
    pos.y{i} = data(i).TrialData.Marker.rawPositions(:,3);
    velo.x{i} =  (pos.x{i}(2:end)-pos.x{i}(1:end-1))/(refresh);
    velo.y{i} =  (pos.y{i}(2:end)-pos.y{i}(1:end-1))/(refresh);
end

for i=1:length(data)
    tempor = timers.movementStart{i}:50:timers.movementEnd{i};
    for j=1:length(tempor)
        if j ==1
             index  = (time{i}>=tempor(j) &time{i}<tempor(j)+50);
             fiftyPos.x{i} = mean(data(i).TrialData.Marker.rawPositions(index,2));
             fiftyPos.y{i} = mean(data(i).TrialData.Marker.rawPositions(index,3));
             fiftyVelo.x{i} = mean(velo.x{i}(index));
             fiftyVelo.y{i} = mean(velo.y{i}(index)); 
        else
             index  = (time{i}>=tempor(j) &time{i}<tempor(j)+50);
             fiftyPos.x{i} = [fiftyPos.x{i};mean(data(i).TrialData.Marker.rawPositions(index,2))];
             fiftyPos.y{i} = [fiftyPos.y{i};mean(data(i).TrialData.Marker.rawPositions(index,3))];
             fiftyVelo.x{i} = [fiftyVelo.x{i};mean(velo.x{i}(index))];
             fiftyVelo.y{i} = [fiftyVelo.y{i};mean(velo.y{i}(index))];
        end
    end
    spikeCount = zeros(length(tempor)-1, length(data(i).TrialData.spikes));
    for neuron = 1:length(data(1).TrialData.spikes)
            spikeTimes = data(i).TrialData.spikes(neuron).timestamps;
            if isempty(spikeTimes)
                counts = zeros(length(tempor)-1,1);
                spikeCount(:,neuron) = counts;
            else
                counts = histcounts(spikeTimes, 'binedges', [tempor, tempor(end)+50]);
                spikeCount(:,neuron) = counts(1:end-1);
            end
    end
    spikes{i} = spikeCount;
end
spikeSum = vertcat(spikes{:});
testArray = find(sum(spikeSum(:)))==0;
for i = 1:1222
   spikes{i}(:, 36:37) = []; 
end
