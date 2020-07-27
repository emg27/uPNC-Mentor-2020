function [spikes,fiftyPos,fiftyVelo] = PosVeloTimeBC(Data)
clc;
%%
% find successful trials. Only include these. 
ov = [Data.Overview];
success = [ov.trialStatus] == 1;
data = Data(success);
% find any trials that do not contain spike data. Exclude these.
%%
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

%%
for i=1:length(data)
    time{i} = data(i).TrialData.DecodeKinematics.decoderTime;
    timers.movementStart{i} = data(i).TrialData.stateTransitions(2,4);
    timers.movementEnd{i} = data(i).TrialData.stateTransitions(2,6);
    pos.x{i} = data(i).TrialData.DecodeKinematics.position(:,1);
    pos.y{i} =data(i).TrialData.DecodeKinematics.position(:,2);
    velo.x{i} =  1000*data(i).TrialData.DecodeKinematics.velocity(:,1);
    velo.y{i} =  1000*data(i).TrialData.DecodeKinematics.velocity(:,2);


end
%%
for i=1:length(data)
    tempor = timers.movementStart{i}:45:timers.movementEnd{i};
    for j=1:length(tempor)
        if j ==1
             index  = (time{i}>=tempor(j) &time{i}<tempor(j)+45);
             spikes{i} = data(i).TrialData.Decoder.rawSpikeBins(index,1:2:end);
             fiftyPos.x{i} = mean(pos.x{i}(index));
             fiftyPos.y{i} = mean(pos.y{i}(index));
             fiftyVelo.x{i} = mean(velo.x{i}(index));
             fiftyVelo.y{i} = mean(velo.y{i}(index)); 
        else
             index  = (time{i}>=tempor(j) &time{i}<tempor(j)+45);
             spikes{i} = [spikes{i}; data(i).TrialData.Decoder.rawSpikeBins(index,1:2:end)];
             fiftyPos.x{i} = [fiftyPos.x{i};mean(pos.x{i}(index))];
             fiftyPos.y{i} = [fiftyPos.y{i};mean(pos.y{i}(index))];
             fiftyVelo.x{i} = [fiftyVelo.x{i};mean(velo.x{i}(index))];
             fiftyVelo.y{i} = [fiftyVelo.y{i};mean(velo.y{i}(index))];
        end
    end
end
%%
checking = 0;
for i=1:length(data)
    checking = sum(spikes{i},1)+checking;
end
index = checking==0;
for i=1:length(data)
    spikes{i}(:,index) = [];
end