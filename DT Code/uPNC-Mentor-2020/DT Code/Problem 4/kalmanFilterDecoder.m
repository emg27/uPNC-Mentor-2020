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
%%
testSpikes = spikes;
data = pmdDataSetup(Data);
totaldata = Data;
%%
D = length(spikes{1}(2,:));
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
    trainingData.veloX{i} = fiftyVelo.x{sortedIndex(i)};
    fiftyVelo.x{sortedIndex(i)}(:) = [];
    trainingData.veloY{i} = fiftyVelo.y{sortedIndex(i)};
    fiftyVelo.y{sortedIndex(i)}(:) = [];
    trainingData.posX{i} = fiftyPos.x{sortedIndex(i)};
    fiftyPos.x{sortedIndex(i)}(:) = [];
    trainingData.posY{i} = fiftyPos.y{sortedIndex(i)};
    fiftyPos.y{sortedIndex(i)}(:) = [];
    trainingData.spikes{i} = spikes{sortedIndex(i)};
    spikes{sortedIndex(i)} = [];
end
%%
Asum= zeros(M,M);
Qsum= zeros(M,M);
Csum = zeros(D,M);
Rsum = zeros(D,D);
buffer = false;
%Prediction 
for trial =1:length(trainingData.veloX)
    
    trialveloX = trainingData.veloX{trial};
    trialposX =  trainingData.posX{trial};
    trialveloY = trainingData.veloY{trial};
    trialposY = trainingData.posY{trial};
    trialSpikes = trainingData.spikes{trial};  
    T = length(trialveloX);
        for bin =2:T
            if ~buffer
                if (isnan(trialveloX(bin))||isnan(trialveloY(bin))||isnan(trialposX(bin))||isnan(trialposY(bin))||isnan(trialSpikes(bin)))
                    buffer = true;
                    continue;
                elseif (isnan(trialveloX(bin-1))||isnan(trialveloY(bin-1))||isnan(trialposX(bin-1))||isnan(trialposY(bin-1))||isnan(trialSpikes(bin-1)))
                    buffer = true;
                    continue;      
                else

                    veloX = trialveloX(bin-1);
                    veloY = trialveloY(bin-1);
                    posX = trialposX(bin-1);
                    posY = trialposY(bin-1);
                    z_t1 = [veloX posX posY veloY]';
                    veloX2 = trialveloX(bin);
                    veloY2 = trialveloY(bin);
                    posX2 = trialposX(bin);
                    posY2 = trialposY(bin);
                    x_t = trialSpikes(bin,:);
                    z_t = [veloX2 posX2 posY2 veloY2]';
                    A= (z_t*z_t1')*(z_t1*z_t1')^(-1);
                    C = (x_t'*z_t')*(z_t*z_t')^-1;
                    Asum = Asum+A;
                    Csum = C+Csum;
                end
            else
                buffer = false;
            end
        end
end
%%
buffer = false;
for trial =1:length(trainingData.veloX)
    trialveloX = trainingData.veloX{trial};
    trialposX =  trainingData.posX{trial};
    trialveloY = trainingData.veloY{trial};
    trialposY = trainingData.posY{trial};
    trialSpikes = trainingData.spikes{trial};
    T = length(trialveloX);
        for bin =1:T
            if ~buffer
                if (isnan(trialveloX(bin))||isnan(trialveloY(bin))||isnan(trialposX(bin))||isnan(trialposY(bin))||isnan(trialSpikes(bin)))
                    buffer = true;
                    continue;
                else
                    veloX = trialveloX(bin);
                    veloY = trialveloY(bin);
                    posX = trialposX(bin);
                    posY = trialposY(bin);
                    x_t = trialSpikes(bin,:);
                    z_t = [veloX posX posY veloY]';
                    Q= (1/(T-1)) * (z_t-Asum*z_t1)*(z_t-Asum*z_t1)';
                    R = (1/T)* (x_t' - Csum*z_t)*(x_t'-Csum*z_t)';
                    Qsum = Q + Qsum;
                    Rsum = R + Rsum;
                end
            else
                buffer = false;
            end
        end
end

%% testing phase

