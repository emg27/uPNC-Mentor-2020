function predictedValues = PosKalmanFilterDecoder(data,posData, spikes)
D = length(spikes{1}(2,:));
M = 2;

%%
testingCoeff = 0.90;
testingsize = floor(length(data)*testingCoeff);
while rem(testingsize,8)~=0
    testingsize = testingsize +1;
end
%%
angles = [0 45 90 135 180 225 270 315];
reachAngles = data(:,1);
count = 0;
trainingData.posX = {testingsize};
trainingData.posY = {testingsize};
trainingData.spikes= {testingsize};

sortedIndex = [];
for i =1:length(angles)
    index = reachAngles==angles(i);
    for f=1:length(index)
        if index(f)==1&&count<(testingsize/8)
            sortedIndex = [sortedIndex f];
            count = count +1;
        elseif index(f)==1&&count>=(testingsize/8)
            count = 0;
            break;
        end
    end
end
%%
for i = 1:length(sortedIndex)
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
testingSpikes = spikes(cellfun(@(j) ~isequal(j,0), spikes));


%%

QSum= zeros(M,M);
QSum1= zeros(M,M);
firstTermQ= 0;
firstTermR= 0;
Rsum = zeros(D,D);
%Prediction 
for trial =1:length(trainingData.spikes)
    trialposX =  trainingData.posX{trial};
    trialposY = trainingData.posY{trial};
    trialSpikes = trainingData.spikes{trial};  
    T = size(trialSpikes,1);
        for bin =1:T
                if (isnan(trialposX(bin))||isnan(trialposY(bin))||isnan(trialSpikes(bin)))
                    continue;
                end
                if bin~=1&&(isnan(trialposX(bin-1))||isnan(trialposY(bin-1))||isnan(trialSpikes(bin-1)))
                    continue; 
                end
                if bin ==1
                    posX2 = trialposX(bin);
                    posY2 = trialposY(bin);
                    x_t = trialSpikes(bin,:);
                    z_t = [posX2 posY2]';
                    Ctop = (x_t'*z_t');
                    Cbottom = (z_t*z_t');
                elseif bin ==2
                    posX = trialposX(bin-1);
                    posY = trialposY(bin-1);
                    z_t1 = [posX posY]';
                    posX2 = trialposX(bin);
                    posY2 = trialposY(bin);
                    x_t = trialSpikes(bin,:);
                    z_t = [posX2 posY2]';
                    Atop =(z_t*z_t1');
                    Abottom =(z_t1*z_t1');
                    Ctop = Ctop+ (x_t'*z_t');
                    Cbottom = Cbottom+ (z_t*z_t');
                else
                    posX = trialposX(bin-1);
                    posY = trialposY(bin-1);
                    z_t1 = [posX posY]';
                    posX2 = trialposX(bin);
                    posY2 = trialposY(bin);
                    x_t = trialSpikes(bin,:);
                    z_t = [posX2 posY2]';
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
    trialposX =  trainingData.posX{trial};
    trialposY = trainingData.posY{trial};
    trialSpikes = trainingData.spikes{trial};
    T = size(trialSpikes,1);
        for bin =1:T
                if (isnan(trialposX(bin))||isnan(trialposY(bin))||isnan(trialSpikes(bin)))
                    continue;
                end

                if bin~=1&&(isnan(trialposX(bin-1))||isnan(trialposY(bin-1))||isnan(trialSpikes(bin-1)))
                    continue; 

                end
                    if bin==1
                        posX = trialposX(bin);
                        posY = trialposY(bin);
                        x_t = trialSpikes(bin,:);
                        z_t = [posX posY]';
                        R = (x_t' - C*z_t)*(x_t'-C*z_t)';
                        Rsum = R + Rsum;
   
                    else
                        posX = trialposX(bin-1);
                        posY = trialposY(bin-1);
                        z_t1 = [posX posY]';
                        posX2 = trialposX(bin);
                        posY2 = trialposY(bin);
                        x_t = trialSpikes(bin,:);
                        z_t = [ posX2 posY2]';
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
predictedValues.sig(length(testingPosX)) = {1};
Z_1 = zeros(length(testingPosX),M);

for trial=1:length(testingPosX)
    Z_1(trial,:) = [testingPosX{trial}(1) testingPosY{trial}(1)];
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
           predictedValues.sig{trial} = Sig_t;
       else 
           predictedValues.muPos{trial} = [predictedValues.muPos{trial}; Mu_t(1:2)'];
           predictedValues.sig{trial} = [predictedValues.sig{trial}; Sig_t];     
       end
    end 
end

%% PREFORMANCE CALCULATION

% Postion & Velocity Difference

error.diffXPos(length(testingPosX)) = {1};
error.diffYPos(length(testingPosX)) = {1};
meanErrorPosX = zeros(length(testingPosX),1);
meanErrorPosY = zeros(length(testingPosX),1);

for trial= 1:length(testingPosX)
    
       if (max(max(isnan(predictedValues.muPos{trial}(:,1)))) == 1 || max(max(isnan(predictedValues.muPos{trial}(:,2)))) == 1 || max(max(isnan(testingPosX{trial}(1:end-1)))) == 1 || max(max(isnan(testingPosY{trial}(1:end-1)))) == 1)
           continue;
       end 
       xp = abs(predictedValues.muPos{trial}(:,1) - testingPosX{trial}(1:end-1))/(testingPosX{trial}(1:end-1));
       yp = abs(predictedValues.muPos{trial}(:,2) - testingPosY{trial}(1:end-1))/(testingPosY{trial}(1:end-1));
       error.diffXPos{trial} = xp((xp ~= 0));
       error.diffYPos{trial} = yp((yp ~= 0));
       meanErrorPosX(trial) = mean(error.diffXPos{trial},1);
       meanErrorPosY(trial) = mean(error.diffYPos{trial},1);
end
predictedValues.Errorperformance= mean([mean(abs(meanErrorPosX))*100, mean(abs(meanErrorPosY))*100]);


%% Target Reach Angle Difference 

% if at end/near end, find the end how far value is from end of postion.
% Taking end of each trial and comparing to. if we know target position

% X and Y pos diff between last value in each bin
differencePosX = zeros(length(testingPosX));
differencePosY = zeros(length(testingPosX));
for trial= 1:length(testingPosX)
    
         if (max(max(isnan(predictedValues.muPos{trial}(:,1)))) == 1 || max(max(isnan(predictedValues.muPos{trial}(:,2)))) == 1 || max(max(isnan(testingPosX{trial}(1:end-1)))) == 1 || max(max(isnan(testingPosY{trial}(1:end-1)))) == 1)
           continue;
         end 
          
         differencePosX(trial) = mean(abs(testingPosX{trial}(end) - predictedValues.muPos{trial}(end,1)));
         
        
         differencePosY(trial) = mean(abs(testingPosY{trial}(end) - predictedValues.muPos{trial}(end,2)));

end 
predictedValues.Distanceperformance = mean([differencePosX;differencePosY]);
