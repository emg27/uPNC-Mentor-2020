 function [spikes,fiftyPos,fiftyVelo] = PosVeloSpikeBC(Data)
clc;
ov = [Data.Overview];
success = [ov.trialStatus] == 1;
data = Data(success);

pos.x = {length(data)};
pos.y = {length(data)};
fiftyPos = pos;
velo.x=  {length(data)};
velo.y = {length(data)};
spikes = {length(data)};
fiftyVelo = velo;
time = {length(data)};
refresh = (25/3)/1000;

%%
for i=1:length(data)  
    time{i} = data(i).TrialData.Decoder.rawDecode(:,2); % get rid of first col,8,9,10
    timers.movementStart{i} = data(i).TrialData.stateTransitions(2,4);
    timers.movementEnd{i} = data(i).TrialData.stateTransitions(2,6);
    pos.x{i} = data(i).TrialData.Decoder.rawDecode(:,3);
    pos.y{i} = data(i).TrialData.Decoder.rawDecode(:,4);
    velo.x{i} =  (pos.x{i}(2:end)-pos.x{i}(1:end-1))/(refresh);
    velo.y{i} =  (pos.y{i}(2:end)-pos.y{i}(1:end-1))/(refresh);
end
%%

for i=1:length(data)
    tempor = timers.movementStart{i}:45:timers.movementEnd{i};
    for j=1:length(tempor)
        if j ==1
             index  = (time{i}>=tempor(j) &time{i}<tempor(j)+45);
             fiftyPos.x{i} = mean(data(i).TrialData.Decoder.rawDecode(index,3));
             fiftyPos.y{i} = mean(data(i).TrialData.Decoder.rawDecode(index,4));
             fiftyVelo.x{i} = mean(velo.x{i}(index));
             fiftyVelo.y{i} = mean(velo.y{i}(index)); 
        else
             index  = (time{i}>=tempor(j) &time{i}<tempor(j)+45);
             fiftyPos.x{i} = [fiftyPos.x{i};mean(data(i).TrialData.Decoder.rawDecode(index,3))];
             fiftyPos.y{i} = [fiftyPos.y{i};mean(data(i).TrialData.Decoder.rawDecode(index,4))];
             fiftyVelo.x{i} = [fiftyVelo.x{i};mean(velo.x{i}(index))];
             fiftyVelo.y{i} = [fiftyVelo.y{i};mean(velo.y{i}(index))];
        end
    end

    % remove ever even col every trial
       spikeCount = zeros(length(tempor)-1, length(data(i).TrialData.spikes));
    for neuron = 1:length(data(1).TrialData.spikes)
            spikeTimes = data(i).TrialData.spikes(neuron).timestamps;
            if isempty(spikeTimes)
                counts = zeros(length(tempor)-1,1);
                spikeCount(:,neuron) = counts;
            else
                counts = histcounts(spikeTimes, 'binedges', [tempor, tempor(end)+45]);
                spikeCount(:,neuron) = counts(1:end-1);
            end
    end
    spikes{i} = spikeCount;
       

end
