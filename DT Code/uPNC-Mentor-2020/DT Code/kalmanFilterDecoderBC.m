function predictedValues = kalmanFilterDecoderBC(data,veloData, posData, spikes, wantComputerToDie, independent)
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
%%
D = length(spikes{1}(2,:));
M = 4;

%%
testingCoeff = 0.80;
testingsize = floor(length(veloData.x)*testingCoeff);
while rem(testingsize,8)~=0
    testingsize = testingsize +1;
end

trainingData.veloX = {testingsize};
trainingData.veloY = {testingsize};
trainingData.posX = {testingsize};
trainingData.posY = {testingsize};
trainingData.spikes= {testingsize};
rng('default')
rng(independent)
varyingTrials = randperm(length(veloData.x));
sortedIndex = varyingTrials(1:testingsize);
for i = 1:length(sortedIndex)
    trainingData.veloX{i} = veloData.x{sortedIndex(i)};
    veloData.x{sortedIndex(i)}= [];
    veloData.x{sortedIndex(i)} = 0;
    trainingData.veloY{i} = veloData.y{sortedIndex(i)};
    veloData.y{sortedIndex(i)} = [];
    veloData.y{sortedIndex(i)} = 0;
    trainingData.posX{i} = posData.x{sortedIndex(i)};
    posData.x{sortedIndex(i)} = [];
    posData.x{sortedIndex(i)} = 0;
    trainingData.posY{i} = posData.y{sortedIndex(i)};
    posData.y{sortedIndex(i)} = [];
    posData.y{sortedIndex(i)} = 0;
    trainingData.spikes{i} = spikes{sortedIndex(i)};
    spikes{sortedIndex(i)} = [];
    spikes{sortedIndex(i)} = 0;
end
%%
testingPosX = posData.x(cellfun(@(j) ~isequal(j,0), posData.x));
testingPosY = posData.y(cellfun(@(j) ~isequal(j,0), posData.y));
testingVeloX = veloData.x(cellfun(@(j) ~isequal(j,0), veloData.x));
testingVeloY = veloData.y(cellfun(@(j) ~isequal(j,0), veloData.y));
testingSpikes = spikes(cellfun(@(j) ~isequal(j,0), spikes));


%%

QSum= zeros(M,M);
QSum1= zeros(M,M);
firstTermQ= 0;
firstTermR= 0;
Rsum = zeros(D,D);
%Prediction 
for trial =1:length(trainingData.spikes)
    
    trialveloX = trainingData.veloX{trial};
    trialposX =  trainingData.posX{trial};
    trialveloY = trainingData.veloY{trial};
    trialposY = trainingData.posY{trial};
    trialSpikes = trainingData.spikes{trial};  
    T = size(trialSpikes,1);
        for bin =1:T
                if (isnan(trialveloX(bin))||isnan(trialveloY(bin))||isnan(trialposX(bin))||isnan(trialposY(bin))||isnan(trialSpikes(bin)))
                    continue;
                end
                if bin~=1&&(isnan(trialveloX(bin-1))||isnan(trialveloY(bin-1))||isnan(trialposX(bin-1))||isnan(trialposY(bin-1))||isnan(trialSpikes(bin-1)))
                    continue; 
                end
                if bin ==1
                    veloX2 = trialveloX(bin);
                    veloY2 = trialveloY(bin);
                    posX2 = trialposX(bin);
                    posY2 = trialposY(bin);
                    x_t = trialSpikes(bin,:);
                    z_t = [posX2 posY2 veloX2 veloY2]';
                    Ctop = (x_t'*z_t');
                    Cbottom = (z_t*z_t');
                elseif bin ==2
                    veloX = trialveloX(bin-1);
                    veloY = trialveloY(bin-1);
                    posX = trialposX(bin-1);
                    posY = trialposY(bin-1);
                    z_t1 = [posX posY veloX veloY]';
                    veloX2 = trialveloX(bin);
                    veloY2 = trialveloY(bin);
                    posX2 = trialposX(bin);
                    posY2 = trialposY(bin);
                    x_t = trialSpikes(bin,:);
                    z_t = [posX2 posY2 veloX2 veloY2]';
                    Atop =(z_t*z_t1');
                    Abottom =(z_t1*z_t1');
                    Ctop = Ctop+ (x_t'*z_t');
                    Cbottom = Cbottom+ (z_t*z_t');
                else
                    veloX = trialveloX(bin-1);
                    veloY = trialveloY(bin-1);
                    posX = trialposX(bin-1);
                    posY = trialposY(bin-1);
                    z_t1 = [posX posY veloX veloY]';
                    veloX2 = trialveloX(bin);
                    veloY2 = trialveloY(bin);
                    posX2 = trialposX(bin);
                    posY2 = trialposY(bin);
                    x_t = trialSpikes(bin,:);
                    z_t = [posX2 posY2 veloX2 veloY2]';
                    Atop = Atop+ (z_t*z_t1');
                    Abottom = Abottom+ (z_t1*z_t1');
                    Ctop = Ctop +(x_t'*z_t');
                    Cbottom = Cbottom+ (z_t*z_t');
                end
        end
end
A = Atop*Abottom^(-1);
C = Ctop*Cbottom^(-1);
for trial =1:length(trainingData.spikes)
    trialveloX = trainingData.veloX{trial};
    trialposX =  trainingData.posX{trial};
    trialveloY = trainingData.veloY{trial};
    trialposY = trainingData.posY{trial};
    trialSpikes = trainingData.spikes{trial};
    T = size(trialSpikes,1);
        for bin =1:T
                if (isnan(trialveloX(bin))||isnan(trialveloY(bin))||isnan(trialposX(bin))||isnan(trialposY(bin))||isnan(trialSpikes(bin)))
                    continue;
                end

                if bin~=1&&(isnan(trialveloX(bin-1))||isnan(trialveloY(bin-1))||isnan(trialposX(bin-1))||isnan(trialposY(bin-1))||isnan(trialSpikes(bin-1)))
                    continue; 
                end
                    if bin==1
                        veloX = trialveloX(bin);
                        veloY = trialveloY(bin);
                        posX = trialposX(bin);
                        posY = trialposY(bin);
                        x_t = trialSpikes(bin,:);
                        z_t = [posX posY veloX veloY]';
                        R = (x_t' - C*z_t)*(x_t'-C*z_t)';
                        Rsum = R + Rsum;
   
                    else
                        veloX = trialveloX(bin-1);
                        veloY = trialveloY(bin-1);
                        posX = trialposX(bin-1);
                        posY = trialposY(bin-1);
                        z_t1 = [posX posY veloX veloY]';
                        veloX2 = trialveloX(bin);
                        veloY2 = trialveloY(bin);
                        posX2 = trialposX(bin);
                        posY2 = trialposY(bin);
                        x_t = trialSpikes(bin,:);
                        z_t = [ posX2 posY2 veloX2 veloY2]';
                        QSum1 = (z_t-A*z_t1)*(z_t-A*z_t1)' +QSum1; %+ QSum1;
                        R =(x_t' - C*z_t)*(x_t'-C*z_t)';
                        Rsum = R + Rsum;
                    end
                    firstTermQ = (T-1)+ firstTermQ;
                    firstTermR = (T)+ firstTermQ;
        end
end

Q= (1/firstTermQ) * QSum1;
R = (1/firstTermR) *Rsum;
%% testing phase
%% Fit step
predictedValues.muPos(length(testingPosX)) = {1};
predictedValues.muVelo(length(testingPosX)) ={1};
predictedValues.sig(length(testingPosX)) = {1};
Z_1 = zeros(length(testingPosX),M);

for trial=1:length(testingPosX)
    Z_1(trial,:) = [testingPosX{trial}(1) testingPosY{trial}(1) testingVeloX{trial}(1) testingVeloY{trial}(1)];
end
Z_1_mean = mean(Z_1)';
Z_1_covar = cov(Z_1);

%%
for trial=1:length(testingPosX)
    T =size(testingSpikes{trial},1);
    Mu_t1 = zeros(size(Z_1_mean));
    Mu_t = Mu_t1;
    Sig_t = zeros(size(Z_1_covar));
    Sig_t1 = Sig_t;
    for bin =1:T
       if bin ==1
           Mu_t1 = Z_1_mean;
           Sig_t1 = Z_1_covar;

       else
           Mu_t1 = A*Mu_t;
           Sig_t1 = A*Sig_t*A'+Q;
       end
       k_t = Sig_t1*C'*(C*Sig_t1*C'+R)^(-1);
       Mu_t= Mu_t1+k_t*(testingSpikes{trial}(bin,:)'-C*Mu_t1);
       Sig_t = Sig_t1-k_t*C*Sig_t1;
       if bin==1
           predictedValues.muPos{trial} = Mu_t(1:2)';
           predictedValues.muVelo{trial} =Mu_t(3:4)';
           predictedValues.sig{trial} = Sig_t;
       else 
           predictedValues.muPos{trial} = [predictedValues.muPos{trial}; Mu_t(1:2)'];
           predictedValues.muVelo{trial} = [predictedValues.muVelo{trial}; Mu_t(3:4)'];
           predictedValues.sig{trial} = [predictedValues.sig{trial}; Sig_t];     
       end
    end 
end
%%
if wantComputerToDie
figure;
    for trial =1:length(predictedValues.muPos)
            tempPosX = zeros(length(predictedValues.muPos{trial}),1);
            tempPosY = zeros(length(predictedValues.muPos{trial}),1);
            for bin = 1:length(predictedValues.muPos{trial})
                if bin ==1
                    tempPosX(1) = predictedValues.muPos{trial}(1,1);
                    tempPosY(1) = predictedValues.muPos{trial}(1,2);
                else
                    tempPosX(bin) = tempPosX(bin-1)+predictedValues.muVelo{trial}(bin,1)*0.045;
                    tempPosY(bin) = tempPosY(bin-1)+predictedValues.muVelo{trial}(bin,2)*0.045;
                end 
            end
            predictedValues.muPos{trial}(:,1) = tempPosX;
            predictedValues.muPos{trial}(:,2) = tempPosY;
            plot(tempPosX, tempPosY, 'Color', 'k')
            xlabel('Position X in mm')
            ylabel('Position Y in mm')
            title('The Estimated Position using the Kalman Filter Plotted over Actual')
            hold on
            plot(testingPosX{trial}, testingPosY{trial}, 'Color', 'g')
    end
end
%% PREFORMANCE CALCULATION

% Postion & Velocity Difference

error.diffXPos(length(testingPosX)) = {1};
error.diffYPos(length(testingPosX)) = {1};
error.diffXVelo(length(testingPosX)) = {1};
error.diffYVelo(length(testingPosX)) = {1};

meanErrorPosX = zeros(length(testingPosX),1);
meanErrorPosY = zeros(length(testingPosX),1);
meanErrorVeloX = zeros(length(testingPosX),1);
meanErrorVeloY = zeros(length(testingPosX),1);

for trial= 1:length(testingPosX)
    
       if (max(max(isnan(predictedValues.muPos{trial}(:,1)))) == 1 || max(max(isnan(predictedValues.muPos{trial}(:,2)))) == 1 || max(max(isnan(testingPosX{trial}(1:end-1)))) == 1 || max(max(isnan(testingPosY{trial}(1:end-1)))) == 1)
           continue;
       end 
       if (max(max(isnan(predictedValues.muVelo{trial}(:,1)))) == 1 || max(max(isnan(predictedValues.muVelo{trial}(:,2)))) == 1 || max(max(isnan(testingVeloX{trial}(1:end-1)))) == 1 || max(max(isnan(testingVeloY{trial}(1:end-1)))) == 1)
           continue;
       end 
       xp = abs(predictedValues.muPos{trial}(:,1) - testingPosX{trial}(1:end))/(testingPosX{trial}(1:end));
       yp = abs(predictedValues.muPos{trial}(:,2) - testingPosY{trial}(1:end))/(testingPosY{trial}(1:end));
       xv = abs(predictedValues.muVelo{trial}(:,1) - testingVeloX{trial}(1:end))/(testingVeloX{trial}(1:end));
       yv = abs(predictedValues.muVelo{trial}(:,2) - testingVeloY{trial}(1:end))/(testingVeloY{trial}(1:end));
       error.diffXPos{trial} = xp((xp ~= 0));
       error.diffYPos{trial} = yp((yp ~= 0));
       error.diffXVelo{trial} = xv((xv ~= 0));
       error.diffYVelo{trial} = yv((yv ~= 0));
       meanErrorPosX(trial) = mean(error.diffXPos{trial},1);
       meanErrorPosY(trial) = mean(error.diffYPos{trial},1);
       meanErrorVeloX(trial) = mean(error.diffXVelo{trial},1);
       meanErrorVeloY(trial) = mean(error.diffYVelo{trial},1);

end
predictedValues.Errorperformance= mean([mean(abs(meanErrorPosX))*100, mean(abs(meanErrorPosY))*100, mean(abs(meanErrorVeloX))*100, mean(abs(meanErrorVeloY))*100]);


%% Target Reach Angle Difference 

% if at end/near end, find the end how far value is from end of postion.
% Taking end of each trial and comparing to. if we know target position

% X and Y pos diff between last value in each bin
differencePosX = zeros(length(testingPosX),1);
differencePosY = zeros(length(testingPosX),1);
for trial= 1:length(testingPosX)
    
         if (max(max(isnan(predictedValues.muPos{trial}(:,1)))) == 1 || max(max(isnan(predictedValues.muPos{trial}(:,2)))) == 1 || max(max(isnan(testingPosX{trial}(1:end-1)))) == 1 || max(max(isnan(testingPosY{trial}(1:end-1)))) == 1)
           continue;
         end 
          
         differencePosX(trial,1) = mean(abs(testingPosX{trial}(end) - predictedValues.muPos{trial}(end,1)));
         
        
         differencePosY(trial,1) = mean(abs(testingPosY{trial}(end) - predictedValues.muPos{trial}(end,2)));

end 
predictedValues.Distanceperformance = mean([differencePosX ;differencePosY]);
