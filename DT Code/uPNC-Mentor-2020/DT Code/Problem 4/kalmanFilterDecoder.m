% function decodedData = kalmanFilterDecoder(data)
% M  = parameters %velo pos etc
% D = numNeurons;
% time = []
%  II and V are sample mean and variance from z1
% Z_t = []; %state variable = velocity
% 
% X_t = []; %observation ??spike count vector?
% 
% %% STATE MODEL
% AmeanZt = [A*Zt-1]; 
% Qgausscovar = []; %gaussian noise
% 
% model(AmeanZt,Qgausscovar)
% 
% z_t = [A11,A12;A21,A22]*z_t1 +Qgausscovar
% xt = xt_1*A11 + A12*z_t1
% %% OBSERVATION MODEL
% A11 =1;
% A12 = delta(t)
% post = A11*post-1 + A12*velt-1
% 
% model(z1)*Prod(model(z = model(z-1))*(prod(model(x-1 = model(z-1)))))
[spikes,fiftyPos,fiftyVelo] = PosVeloTime(Data);
data = pmdDataSetup(Data)
totaldata = Data;
% veloX = fiftyVelo.x{1}(1);
% veloY = fiftyVelo.y{1}(1);
% posX = fiftyPos.x{1}(1);
% posY = fiftyPos.y{1}(1);
% z_t1 = [veloX posX posY veloY]';
% T = length(fiftyVelo.x{1});
% veloX2 = fiftyVelo.x{1}(2);
% veloY2 = fiftyVelo.y{1}(2);
% posX2 = fiftyPos.x{1}(2);
% posY2 = fiftyPos.y{1}(2);
% z_t = [veloX2 posX2 posY2 veloY2]';
% x_t = spikes{1}(2,:);
D = length(x_t);
M = 4;
%%
angles = [0 45 90 135 180 225 270 315];
reachAngles = data(:,1);
count = 0;
trainingData.veloX = {808};
trainingData.veloY = {808};
trainingData.posX = {808};
trainingData.posY = {808};
trainingData.spikes= {808};

sortedIndex = [];
for i =1:length(angles)
    index = reachAngles==angles(i);
    for j=1:length(index)
        if index(j)==1&&count<101
            sortedIndex = [sortedIndex j];
            count = count +1;
        elseif index(j)==1&&count>=101
            count = 0;
            break;
        end
    end
end


%%
for i = 1:length(sortedIndex)
    trainingData.veloX{i} = fiftyVelo.x{sortedIndex(i)}(:);
    fiftyVelo.x{sortedIndex(i)}(:) = [];
    trainingData.veloY{i} = fiftyVelo.y{sortedIndex(i)}(:);
    fiftyVelo.y{sortedIndex(i)}(:) = [];
    trainingData.posX{i} = fiftyPos.x{sortedIndex(i)}(:);
    fiftyPos.x{sortedIndex(i)}(:) = [];
    trainingData.posY{i} = fiftyPos.y{sortedIndex(i)}(:);
    fiftyPos.y{sortedIndex(i)}(:) = [];
    trainingData.spikes{i} = spikes{sortedIndex(i)}(:);
    spikes{i}(:) = [];
end
%%
fnV = fieldnames(fiftyVelo);
% fnP = fieldnames(fiftyPos)
tfV = cellfun(@(c) isempty(fiftyVelo.(c)), fnV);
S2 = rmfield(fiftyVelo, fnV(tfV))
%%
Asum= zeros(M,M);
Qsum= zeros(M,M);
Csum = zeros(D,M);
Rsum = zeros(D,D);
%% Prediction 
for trial =1:1222
    
    trialveloX = fiftyVelo.x{trial};
    trialposX = fiftyPos.x{trial};
    trialveloY = fiftyVelo.y{trial};
    trialposY = fiftyPos.y{trial};
    T = length(trialveloX);
        for bin =2:T
            veloX = fiftyVelo.x{trial}(bin-1);
            veloY = fiftyVelo.y{trial}(bin-1);
            posX = fiftyPos.x{trial}(bin-1);
            posY = fiftyPos.y{trial}(bin-1);
            z_t1 = [veloX posX posY veloY]';
            veloX2 = fiftyVelo.x{trial}(bin);
            veloY2 = fiftyVelo.y{trial}(bin);
            posX2 = fiftyPos.x{trial}(bin);
            posY2 = fiftyPos.y{trial}(bin);
            x_t = spikes{trial}(bin,:);
            z_t = [veloX2 posX2 posY2 veloY2]';
            A= (z_t*z_t1'T)*(z_t1*z_t1'T)^(-1);
            C = (x_t'*z_t'T)*(z_t*z_t'T)^-1;
            Asum = A+ Asum;
            Csum = C+Csum;
        end
end

for trial =1:1222
    trialveloX = fiftyVelo.x{trial};
    trialposX = fiftyPos.x{trial};
    trialveloY = fiftyVelo.y{trial};
    trialposY = fiftyPos.y{trial};
    T = length(trialveloX);
        for bin =1:T
            veloX = fiftyVelo.x{trial}(bin);
            veloY = fiftyVelo.y{trial}(bin);
            posX = fiftyPos.x{trial}(bin);
            posY = fiftyPos.y{trial}(bin);
            x_t = spikes{trial}(bin,:);
            z_t = [veloX2 posX2 posY2 veloY2]';
            Q= (1/(T-1)) * (z_t-A*z_t1)*(z_t-A*z_t1)'^(-1);
            R = (1/T)* (x_t' - C*z_t)*(x_t'-C*z_t)';
            Qsum = Q + Qsum;
            Rsum = R + Rsum;
        end
end

% for trial =1:1222
%     for bin = 1:T
%         z_t
%     end
% end
% A= (z_t*z_t1'.^T)*(z_t1*z_t1'.^T)^(-1);
% Q= (1/(T-1)) * (z_t-A*z_t1)*(z_t-A*z_t1)'.^(-1);

% R = (1/T)* (x_t' - C*z_t)*(x_t'-C*z_t)'.^T;
% %% testing phase